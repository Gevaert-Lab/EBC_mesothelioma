
dep_volcano <- function ( label, data ,  level , FC_thr  ){
  cmp = gsub('Group', '',label) 
  all_res <-  rowData(data[[level]])[[label]]
  all_res <- all_res %>% rownames_to_column(var = 'Uniprot_id' )

  #perc_field <- rowData(pe[[level]]) %>% colnames() %>%  stringr::str_subset('perc') 
  temp <- as.data.frame(rowData(data[[level]])) %>% rownames_to_column('Uniprot_id') %>%  
    dplyr::select(Uniprot_id,Protein.Names,  Genes, Protein.Names ) 

  all_res <-  all_res %>% left_join( temp, by=join_by(Uniprot_id)) 
  
  log_info(paste0(cmp,'#by MSqrob: ', dim(all_res)[1]))
  all_res <- all_res[ ! is.na(all_res$adjPval),]
  log_info(paste0(cmp,'# by MSqrob after p-adj Null filt.: ', dim(all_res)[1]))
  
  #all_res$Protein.names <- rowData(pe[["proteinRS"]])[['Protein.names']]
  all_res$differential_expressed <- "NO"
  all_res$differential_expressed[all_res$logFC >= FC_thr & all_res$adjPval <= 0.05 ] <- "UP"
  all_res$differential_expressed[all_res$logFC < FC_thr & all_res$adjPval <= 0.05] <- "DOWN"
  
  #browser()
  #sum(all_res$adjPval<0.05, na.rm = T)/(nrow(all_res)) *0.05
  #  0.05
  p1 <- ggplot(data = all_res , aes(x = logFC, y = -log10(pval) ,col=differential_expressed , 
                                      text = sprintf("Protein_name: %s <br> Gene %s" , all_res$Protein.Names, all_res$Genes )   )  )  +
      geom_point() +
      theme_minimal() +
      #geom_text_repel() +
      geom_vline(xintercept = c(- FC_thr , FC_thr),col="grey") +
      geom_hline(yintercept = -log10(0.05),col="grey") +
      scale_color_manual(values=c("DOWN"="blue","NO"="black", "UP"="red")) +
      ggtitle( cmp )
    
    #log_info({c('Uniprot_id',  "Protein.Names" , "Genes", "adjPval","pval","logFC", "differential_expressed",perc_field)})
    #log_info(head(all_res))
    DEall <- all_res[!is.na(all_res$adjPval) , c('Uniprot_id', 'Protein.Names'
                                                 ,'Genes', "adjPval","pval","logFC", "differential_expressed")]
  
  ## volcano annotate with gene name 
  str(all_res)
  all_res_file <- all_res  %>% mutate( label_DE =  ifelse(abs(logFC) > 0.65, Genes, NA)
                                          )
  p_toFile <- ggplot(data = all_res_file , aes(x = logFC, y = -log10(pval) ,col=differential_expressed , 
                                               label= label_DE  )  )  +
    geom_point() +
    geom_text_repel(max.overlaps = Inf, size = 3) +
    geom_vline(xintercept = c(- FC_thr, FC_thr),col="grey") +
    geom_hline(yintercept = -log10(0.05),col="grey") +
    scale_color_manual(values=c("DOWN"="blue","NO"="black", "UP"="red"))+
    ggtitle(cmp) +
    theme(legend.position = "none")
  
  
  
  return ( list( toptable =DEall , volcano = p1, volcano2file =p_toFile ) )
}

prot_boxplot_global <- function(id_uniprot, id_gene_name,  pe , layer ){
  
  ab <- as.data.frame(assay(pe[[layer]])[id_uniprot,])
  base <- as.data.frame(assay(pe[['proteinB']])[id_uniprot,])
  base <- base %>%
    mutate(across(everything(),   ~ case_when(
      !is.na(.) ~ "not imputed",
      is.na(.)  ~ "imputed"
    ) )) 
  
  base <-base %>% rownames_to_column('Protein') %>% 
    gather(Sample, Imputed ,- c(Protein ))
  
  ab <- ab %>% rownames_to_column('Protein') %>% 
    gather(Sample, Intensity ,- c(Protein )) %>% 
    left_join( as.data.frame(colData(pe)) %>% rownames_to_column('Sample'), join_by(Sample) )
  
  ab <- ab %>% left_join( base,  join_by(Sample, Protein) )
  
  ab <- ab %>% mutate(Gene = id_gene_name[match(Protein, id_uniprot)])
  
  
  
  dep_plot_crip <- ab %>%
    ggplot( aes(x= Group , y= Intensity )) +
    geom_boxplot( outliers = FALSE) +
    geom_jitter(aes(alpha = Imputed, color = Imputed), size = 1) +
    scale_alpha_manual(values = c("not imputed" = 0.8, "imputed" = 0.7)) +
    scale_color_manual(values = c("not imputed" = "black", "imputed" = "#EFB118FF")) +
    ylab('Summarized Protein Intensity') +
    xlab('')+
    ggtitle('') + 
    facet_wrap( . ~ Gene  , nrow = 3, ncol = 3  ) +
    guides(
    color = guide_legend(title = NULL),
  alpha = guide_legend(title = NULL)
    )
  
  
  #dep_plot_crip
  
  return (dep_plot_crip)
}

prot_boxplot <- function ( ttop, g2search , comp_vec , pe , layer , title_ ){
   print(g2search)
   id <- ttop %>% filter(Genes %in% g2search) %>% pull(Uniprot_id )
   
   sample <- as.data.frame(colData(pe)) %>% filter(Group == comp_vec[1] | Group == comp_vec[2] ) %>% rownames()
  
   ab <- as.data.frame(assay(pe[[layer]])[id,sample])
   
   print(head(ab))
   ann <-  ttop  %>%   filter ( Uniprot_id  %in% id ) %>% dplyr::select(Uniprot_id ,Genes ,logFC,adjPval) %>%  
     mutate( Protein_info= str_wrap( paste(Genes,'\n','LogFC:',round(logFC,2),'adjPval',format(adjPval,scientific= TRUE) ,sep=' '), width =15) )  %>% 
     dplyr::select(- Genes, , -adjPval)
   
   ab <-  ab %>% rownames_to_column('Protein') %>% 
     left_join(ann, join_by(Protein  == Uniprot_id)) %>% 
     gather(Sample, Intensity ,- c(Protein ,Protein_info,logFC)) %>% 
     left_join( as.data.frame(colData(pe)) %>% rownames_to_column('Sample'), join_by(Sample) )
   
   print(ab %>% head())
     
   dep_plot_crip <- ab %>%
     ggplot( aes(x= Group , y= Intensity, fill= Group )) +
     geom_boxplot() +
     scale_fill_manual( values= c( "#6d0286", "#35B779FF") ) +
     geom_jitter(color="black", size=1, alpha=0.8)  +
     ylab('Summarized Protein Intensity ') +
     ggtitle(title_) +
     facet_wrap( . ~ reorder(Protein_info, logFC), nrow = 2  )
   return(dep_plot_crip)
  }
