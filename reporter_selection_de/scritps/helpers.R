library(tidyverse)
library(plotly)


bh_correction <- function(condition_to_correct, fdr, complete_dataset,
                          fc){
  original_condition_df <- subset(complete_dataset, 
                                  condition==condition_to_correct)
  
  original_condition_df <- original_condition_df %>% 
    mutate(padj = p.adjust(pvalue, method="BH"))
  
  orginal_corrected <- original_condition_df %>% 
    mutate("corrected_de" = case_when(padj <= fdr & log2FoldChange > fc ~ 1,
                                      padj <= fdr & log2FoldChange < (-1*fc) ~ -1,
                                      TRUE ~ 0))
  
  return(orginal_corrected)
}

apply_correction <- function(complete_de_df, 
                             fdr_chosen, fc_chosen){
  
  camp_de <- complete_de_df %>% 
    filter(endsWith(gene_name, "-Campylobacter"))
  
  salm_de <- complete_de_df %>% 
    filter(endsWith(gene_name, "-Salmonella"))
  
  complete_corrected_salmonella <- bind_rows(lapply(unique(salm_de$condition), 
                                                    bh_correction,
                                                    complete_dataset=salm_de,
                                                    fdr=fdr_chosen,
                                                    fc=fc_chosen))
  
  complete_corrected_campylobacter <- bind_rows(lapply(unique(camp_de$condition), 
                                                       bh_correction,
                                                       complete_dataset=camp_de,
                                                       fdr=fdr_chosen,
                                                       fc=fc_chosen))
  
  all_patterns <- bind_rows(complete_corrected_campylobacter, 
                            complete_corrected_salmonella)
  
  return(all_patterns)
}

gather_patterns <- function(de_dataframe, features_df, feature_of_interest,
                            selected_conditions){
  
  # Format the dataframe to get rows of 1s and 0s
  de_patterns.wide <- de_dataframe %>% 
    select(gene_name, corrected_de, condition) %>% 
    filter(condition %in% selected_conditions) %>% 
    pivot_wider(names_from = condition, values_from=corrected_de) %>% 
    column_to_rownames("gene_name")
  
  # Remove genes that have no DE in any condition
  rows_all_zero <- apply(de_patterns.wide, 1, function(x){all(x == 0)})
  joint_de <- de_patterns.wide[!rows_all_zero,]
  
  joint_de <- joint_de %>% 
    rownames_to_column("gene_name") %>% 
    left_join(features_df, by="gene_name") %>% 
    filter(`# feature` %in% feature_of_interest)
  
  return(joint_de)
  
}
  
plot_frequency <- function(joint_de, min_genes, xcoord, ycoord,
                           selected_conditions){
  
  # Count pattern frequency
  pattern_count <- plyr::count(joint_de %>% 
                                 select(-cje_loctag) %>% 
                                 distinct(.keep_all = T), vars = selected_conditions)
  
  freq.df <- pattern_count %>% 
    count(freq) %>% 
    rename(n_genes=freq,
           n_patterns = n) %>% 
    mutate(n_genes_fctr = factor(n_genes, levels = 1:max(n_genes)),
           n_patterns_log = log2(n_patterns),
           selected=if_else(n_genes >= min_genes, "chosen", "not_chosen"),
           ngenes_in_npatterns = n_genes * n_patterns)
  
  # Build the annotation
  n_patterns_selected <- sum(freq.df[freq.df$selected=="chosen", "n_patterns"])
  n_genes_selected <- sum(freq.df[freq.df$selected=="chosen", "ngenes_in_npatterns"])
  annotation <- paste0("There are ", n_patterns_selected, " patterns with at least ", min_genes, " genes\n (", n_genes_selected, " genes represented)")
  
  ggplot(freq.df, aes(x=n_genes_fctr, y=n_patterns_log, fill=selected)) +
    geom_bar(stat="identity", color="black") + 
    geom_text(aes(label=n_patterns), vjust=-0.3, size=3) +
    scale_fill_manual(values=c(alpha("red", 0.5), alpha("grey", 0.5))) +
    annotate("text", x=xcoord, y=ycoord, label=annotation) +
    theme_minimal() +
    theme(legend.position = "none") +
    labs(x="Number of genes with the same pattern", 
         y="Number of Patterns with X genes (log2)",
         title = "Number of patterns followed a number of genes")
  
}

create_pattern_database <- function(joint_de, selected_conditions){
  
  # Count patterns
  pattern_count <- pattern_count <- plyr::count(joint_de %>% 
                                                  select(-cje_loctag) %>% 
                                                  distinct(.keep_all = T), vars = selected_conditions)
  
  pattern.db <- pattern_count %>% 
    rename(n_genes=freq) %>% 
    arrange(desc(n_genes)) %>% 
    mutate(pattern_id = paste0("Pattern ", 1:nrow(pattern_count)),
           patid_ngenes = paste0(pattern_id, " (", n_genes, " genes)"),
           
           pattern_id = factor(pattern_id, levels = pattern_id),
           patid_ngenes = factor(patid_ngenes, levels = patid_ngenes))
  
  pattern.db$pattern_str <- apply(pattern.db[, selected_conditions], 1, paste, 
                                  collapse=',')
  
  return(pattern.db)
}


plot_patterns <- function(pat.db, min_genes, selected_conditions){
  
  selected.patterns <- pat.db %>% 
    filter(n_genes >= min_genes) %>% 
    mutate(Ctrl = 0) %>% 
    gather(c(selected_conditions, "Ctrl"), 
           key = "condition", value="DE") %>% 
    mutate(condition = factor(condition, levels=c("Ctrl", selected_conditions)),
           DE_fctr=factor(DE, levels = c(1,0,-1)))
  
  ggplot(selected.patterns, aes(x=condition, y=DE, group=pattern_id, color=DE_fctr)) +
    geom_line(color=alpha("black", 0.75)) +
    geom_point() +
    scale_color_manual(values=c("#ef8a62", "grey", "#67a9cf")) +
    facet_wrap(~ patid_ngenes) +
    theme_minimal() +
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 90)) +
    labs(x="Condition", y="Direction of Differential Expression",
         title="Differential Expression Patterns")
}

combine_patterns_and_fc <- function(pattern.database.df, pattern_by_gene_df, 
                                    information_df, selected_conditions){
  
  filt_info_df <- information_df %>% 
    filter(condition %in% selected_conditions) %>% 
    select(-locus_tag)
  
  gene_and_patternID <- pattern_by_gene_df %>% 
    mutate(pattern_str = apply(pattern_by_gene_df[,selected_conditions], 1, paste, collapse=',')) %>% 
    left_join(pattern.database.df %>% select(n_genes, pattern_id, patid_ngenes, pattern_str), by="pattern_str")
  
  collected_info <- gene_and_patternID %>% 
    left_join(filt_info_df, by="gene_name") %>% 
    select(-selected_conditions)
  
  return(collected_info)
}


plot_gene_wrap <- function(complete_data, min_genes, color_by_this, selected_conditions){
  
  plot_data <- complete_data %>% 
    filter(n_genes >= min_genes) %>% 
    select(gene_name, regulator_general, pre.selected, clone, patid_ngenes, condition, log2FC.shrink) %>% 
    pivot_wider(names_from = condition, values_from=log2FC.shrink) %>% 
    mutate(Ctrl = 0,
           organism = sapply(gene_name, function(x){tail(strsplit(x, "-")[[1]], n=1)})) %>% 
    gather(c("Ctrl", selected_conditions), key = "condition", value = "log2FC.shrink") %>% 
    mutate(condition = factor(condition, levels = c("Ctrl", selected_conditions)))
  
  if(color_by_this == "organism"){
    camp.data <- subset(plot_data, organism=="Campylobacter")
    salm.data <- subset(plot_data, organism=="Salmonella")
    
    campy_highlight <- ggplot() +
      geom_line(data=salm.data, aes(x=condition, y=log2FC.shrink, group=gene_name),
                color=alpha("grey", 0.5)) +
      geom_line(data=camp.data, aes(x=condition, y=log2FC.shrink, group=gene_name),
                color=alpha("steelblue", 0.5)) +
      facet_wrap(~ patid_ngenes, scales="free") +
      theme_minimal() +
      theme(legend.position = "none",
            axis.text.x = element_text(angle = 90)) +
      labs(x="Condition", y="log Fold Change",
           title="Fold Change Patterns",
           subtitle = "Salmonella in Grey\nCampylobacter in blue")
    
    return(campy_highlight)
    
  }else if(color_by_this == "Salmonella pre-selection"){
    salm.data <- subset(plot_data, organism=="Salmonella")
    salm.pre <- subset(salm.data, pre.selected == "pre.selected")
    
    salmonella.pre.wrap <- ggplot() +
      geom_line(data=salm.data, aes(x=condition, y=log2FC.shrink, group=gene_name),
                color=alpha("grey", 0.5)) +
      geom_line(data=salm.pre, aes(x=condition, y=log2FC.shrink, group=gene_name),
                color=alpha("red", 0.5)) +
      facet_wrap(~ patid_ngenes, scales="free") +
      theme_minimal() +
      theme(legend.position = "none",
            axis.text.x = element_text(angle = 90)) +
      labs(x="Condition", y="log Fold Change",
           title="Salmonella pre-selection",
           subtitle = "Pre-selected in red")
    
    return(salmonella.pre.wrap)
    
  }else if(color_by_this == "Campylobacter pre-selection"){
    camp.data <- subset(plot_data, organism=="Campylobacter")
    camp.pre <- subset(camp.data, pre.selected == "pre.selected")
    camp.clo <- subset(camp.pre, clone =="yes")
    
    campylobacter.pre.wrap <- ggplot() +
      geom_line(data=camp.data, aes(x=condition, y=log2FC.shrink, group=gene_name),
                color=alpha("grey", 0.5)) +
      geom_line(data=camp.pre, aes(x=condition, y=log2FC.shrink, group=gene_name),
                color=alpha("steelblue", 1)) +
      geom_line(data=camp.clo, aes(x=condition, y=log2FC.shrink, group=gene_name),
                color=alpha("red", 1)) +
      facet_wrap(~ patid_ngenes, scales="free") +
      theme_minimal() +
      theme(legend.position = "none",
            axis.text.x = element_text(angle = 90)) +
      labs(x="Condition", y="log Fold Change",
           title="Campylobacter pre-selection",
           subtitle = "Pre-selected in blue\nCloned in red")
    
    return(campylobacter.pre.wrap)
    
  }else if(color_by_this == "Transcription Factors"){
    camp.tfs <- subset(plot_data, organism=="Campylobacter" & regulator_general=="TF")
    salm.tfs <- subset(plot_data, organism=="Salmonella" & regulator_general=="TF")
    
    tf.wrap <- ggplot() + 
      geom_line(data=plot_data, aes(x=condition, y=log2FC.shrink, group=gene_name),
                color=alpha("grey", 0.5)) +
      geom_line(data=salm.tfs, aes(x=condition, y=log2FC.shrink, group=gene_name),
                color=alpha("#d95f02", 1)) +
      geom_line(data=camp.tfs, aes(x=condition, y=log2FC.shrink, group=gene_name),
                color=alpha("#7570b3", 1)) +
      facet_wrap(~ patid_ngenes, scales="free") +
      theme_minimal() +
      theme(legend.position = "none",
            axis.text.x = element_text(angle = 90)) +
      labs(x="Condition", y="log Fold Change",
           title="Transcription Factors",
           subtitle = "Campylobacter TFs in purple\nSalmonella TFs in red")
    
    return(tf.wrap)
      
  }
}

focused_profile <- function(which_organism, which_cluster, which_gene, 
                            dataframe_for_plot, selected_conditions){
  
  dataframe_for_plot <- dataframe_for_plot %>% 
    filter(pattern_id == which_cluster) %>% 
    select(gene_name, condition, log2FC.shrink) %>% 
    pivot_wider(names_from = condition, values_from=log2FC.shrink) %>% 
    mutate(Ctrl = 0,
           organism = sapply(gene_name, function(x){tail(strsplit(x, "-")[[1]], n=1)})) %>% 
    gather(c("Ctrl", selected_conditions), key = "condition", value = "log2FC.shrink") %>% 
    mutate(condition = factor(condition, levels = c("Ctrl", selected_conditions)))
  

  if(which_organism == "Both organisms"){
    f.df <- dataframe_for_plot
    
  }else if(which_organism == "Only Salmonella"){
    f.df <- dataframe_for_plot %>% 
      filter(organism == "Salmonella")
    
  }else if(which_organism == "Only Campylobacter"){
    f.df <- dataframe_for_plot %>% 
      filter(organism == "Campylobacter")
  }
  
  f.plot <- ggplot(f.df, aes(x=condition, y=log2FC.shrink, group=gene_name)) +
    geom_line(color=alpha("darkgrey", 0.5))
  
  if(which_gene != "None"){
    f.plot <- f.plot + 
      geom_line(data=f.df %>% filter(gene_name %in% which_gene),
                aes(x=condition, y=log2FC.shrink, group=gene_name, color=gene_name))
  }
  ggplotly(f.plot)
}


focused_profile_de <- function(pat.db, which_pattern, selected_conditions){
  
  pat.lot <- pat.db %>% 
    filter(pattern_id == which_pattern) %>% 
    mutate(Ctrl = 0) %>% 
    gather(c("Ctrl",selected_conditions), 
           key = "condition", value="DE") %>% 
    mutate(condition = factor(condition, levels=c("Ctrl",selected_conditions)),
           DE_fctr=factor(DE, levels = c(1,0,-1)))
  
  p <- ggplot(pat.lot, aes(x=condition, y=DE, group=1, color=DE_fctr)) +
    geom_line(color=alpha("black", 0.75)) +
    geom_point() +
    scale_color_manual(values=c("#ef8a62", "grey", "#67a9cf")) +
    theme_minimal() +
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 90)) +
    labs(x="Condition", y="Direction of Differential Expression")
  
  return(p)
}

pattern_table <- function(complete_dataframe, which_pattern){
  complete_dataframe %>% 
  select(gene_name, locus_tag, cje_loctag, symbol, name, 
         regulator_general, pre.selected, clone, bbh, pattern_id) %>% 
    filter(pattern_id == which_pattern) %>% 
    distinct()
}

all_genes_table <- function(complete_dataframe){
  complete_dataframe %>% 
    select(gene_name, locus_tag, cje_loctag, symbol, name, 
           regulator_general, pre.selected, clone, bbh, pattern_id, n_genes, pattern_str) %>% 
    distinct()
}
  









