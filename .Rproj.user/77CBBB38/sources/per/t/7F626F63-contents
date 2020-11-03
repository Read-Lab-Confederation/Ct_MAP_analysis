read_MAP <- function(path_mapfile,mfn) {
  mapdf <- read.table(path_mapfile, row.names = NULL, sep = "\t")
  colnames(mapdf) <- c("Position","Refnt","Coverage","E","E.1","Aa","Gg","Cc","Tt","Ref_cov","Map_value")
  mapdf$Position <- as.integer(mapdf$Position)
  # find all places where the numbers add up. The "regular" column means that the ref coverage + all minor bases == total coverage
  mapdf <- mapdf %>% 
    add_column(strain = mfn) %>%
    select(-E,-E.1) %>%
    mutate(regular = (Coverage - Ref_cov - Aa - Gg - Cc - Tt)) %>%
    mutate(minor_counts = Aa + Gg + Cc + Tt) %>%
    mutate(SNP_percent = round(Ref_cov/(Ref_cov + Aa + Gg + Cc + Tt)*100,digits=1))
  return(mapdf)
}

plot_MAP_histo <- function(mapDF,map_file_out,mfn) {
  plotDF <- mapDF %>%
    filter(Coverage > 10) %>%
    filter(SNP_percent < 99) %>%
    filter(regular > -1 & regular < 1)
    ggplot(plotDF,aes(SNP_percent)) + 
    geom_histogram(binwidth=1, fill = "red") +
    ggtitle(mfn) +
    theme_bw()
  plot_fn <- paste0(mfn,"_histo.pdf")
  ggsave(filename= plot_fn, path=map_file_out,device = pdf)
  
}

plot_MAP_coverage <- function(mapDF,map_file_out,mfn) {
  # remove 90% of positions for speed
  plotDF <- mapDF %>%
    filter(Position%%10 == 0) 
  ggplot(plotDF,aes(Position, Coverage)) + 
    geom_line(aes(y=rollmean(Coverage, 1000, fill = NA))) +
    geom_smooth() +
    ggtitle(mfn) +
    theme_bw()
  plot_fn <- paste0(mfn,"_covplot.pdf")
  ggsave(filename= plot_fn, path=map_file_out,device = pdf)
}

write_potential_minors <- function(mapDF,map_file_out,mfn) {
  potential_true_minors <- mapDF %>%
    filter(regular > -1 & regular < 1) %>%
    filter(Coverage > 10) %>%
    filter(SNP_percent >= 10 & SNP_percent < 90) %>%
    select(strain,Position,Refnt,Coverage,Ref_cov,Aa,Gg,Cc,Tt)
  map_file_out_path <- paste0(map_file_out,mfn,"_potential_minors.tsv")
  write_tsv(potential_true_minors,map_file_out_path)
}

write_potential_SNPs <- function(mapDF,map_file_out,mfn) {
  potential_true_minors <- mapDF %>%
    filter(regular > -1 & regular < 1) %>%
    filter(Coverage > 10) %>%
    filter(SNP_percent < 10) %>%
    select(strain,Position,Refnt,Coverage,Ref_cov,Aa,Gg,Cc,Tt)
  map_file_out_path <- paste0(map_file_out,mfn,"_potential_SNPs.tsv")
  write_tsv(potential_true_minors,map_file_out_path)
}

MAP_stats_summary <- function(mapDF,map_file_out,mfn) {
  total_nt <- mapDF %>% 
    tally(name="total_nt") 
  reg0 <- mapDF %>% 
    filter(regular == 0) %>%
    tally(name="reg0")
  map_value1 <- mapDF %>% 
    filter(Map_value == 1) %>%
    tally(name="map_value1")
  reg_more_than_1 <- mapDF %>%
    filter(regular < -1 | regular > 1) %>%
    tally(name="reg_more_than_1")
  map_value_more_than1 <- mapDF %>% 
    filter(Map_value > 1) %>%
    tally(name="map_value_more_than1")
  cov10 <- mapDF %>% 
    filter(Coverage > 10) %>%
    tally(name="cov10")
  pot_true_SNPs <- mapDF %>%
    filter(regular > -1 & regular < 1) %>%
    filter(Coverage > 10) %>%
    filter(SNP_percent < 10) %>%
    tally(name="pot_true_SNP")
  pot_true_minors <- mapDF %>%
    filter(regular > -1 & regular < 1) %>%
    filter(Coverage > 10) %>%
    filter(SNP_percent >= 10 & SNP_percent < 90) %>%
    tally(name="pot_true_minors")
  summary_df <- data_frame(total_nt,reg0,map_value1,reg_more_than_1,map_value_more_than1,cov10,pot_true_minors,pot_true_SNPs) %>%
    add_column(strain = mfn,.before = "total_nt")
  map_file_out_path <- paste0(map_file_out,mfn,"_stats.tsv")
  write_tsv(summary_df,map_file_out_path)
  return(summary_df)
  
}