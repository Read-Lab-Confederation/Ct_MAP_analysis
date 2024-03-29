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

write_potential_minors <- function(mapDF,map_file_out,mfn,cov_cutoff=10) {
  potential_true_minors <- mapDF %>%
    filter(regular > -1 & regular < 1) %>%
    filter(Coverage > cov_cutoff) %>%
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

recomb_coords <- function(fiji_df,samp_name){
  tmp_df <- fiji_df %>% 
    filter(c1 == samp_name | c2 == samp_name | c3== samp_name) %>% 
    select(Start,End)
  if (nrow(tmp_df)==0) {
    return(vector())
  }
  else {
    result <- vector()
    for (i in 1:nrow(tmp_df)) {
      result <- c(result,seq(tmp_df$Start[i],tmp_df$End[i]))
    }
    return(result)
  }
}

SNP_analysis_pipeline <- function(subject_id, rare_cutoff = 3,cov_cutoff = 10) {
  print(paste("SUBJECT_ID is",subject_id))
  
  library(tidyverse)
  source("MAP_functions.R")
  library(assertthat)
  #load Pos_count
  #rare_cutoff <- 3
  Pos_count <- read_delim("Pos_count.tsv", show_col_types = FALSE)
  assert_that(nrow(Pos_count) > 0)
  rare_iSNPs <- Pos_count %>% filter(n<=rare_cutoff) %>%
    .$Position
  common_iSNPs <- Pos_count %>% filter(n>rare_cutoff) %>%
    .$Position
  num_rare = length(rare_iSNPs)
  num_common = length(common_iSNPs)
  print(paste("Rare iSNP cutoff = ",rare_cutoff))
  print(paste("num rare iSNPs = ",num_rare))
  print(paste("num common iSNPs = ",num_common))
  #using the "iSNP" designation that I later changed to "SNV"
  
  #load recomb_coords - recombination blocks cutoffs = 1
  fastgear_fiji_recombs <- read.delim("~/GitHub/Ct_MAP_analysis/fastgear_fiji_recombs") %>%
    select(Start,End,...1,...2,...3)
  print("Recombination block cutoff = 3")
  assert_that(nrow(fastgear_fiji_recombs) > 0)
  colnames(fastgear_fiji_recombs) <- c("Start","End","c1","c2","c3")
  R_name = paste0(subject_id,"R")
  V_name = paste0(subject_id,"V")
  C_name = paste0(subject_id,"C")
  C_recomb <- recomb_coords(fastgear_fiji_recombs,C_name)
  R_recomb <- recomb_coords(fastgear_fiji_recombs,R_name)
  V_recomb <- recomb_coords(fastgear_fiji_recombs,V_name)
  C_recomb_bases <- length(C_recomb)
  R_recomb_bases <- length(R_recomb)
  V_recomb_bases <- length(V_recomb)
  print(paste("Num bases in C recomb blocks= ",C_recomb_bases))
  print(paste("Num bases in R recomb = blocks",R_recomb_bases))
  print(paste("Num bases in V recomb = blocks",V_recomb_bases))
  
  # load SNPs from MAP analysis (NOTE cutoff cov=10)
  #note default cov_cutoff = 10
  print(paste("Coverage cutoff = ",cov_cutoff))
  R_MAP <- read_MAP(paste0("~/Documents/2021-Fiji-Ct_paper/MAP_results/MAP_files/",subject_id,"R_MAP.txt"),"R") %>%
    #filter(regular > -1 & regular < 1) %>%
    filter(Coverage > cov_cutoff) %>%
    select(Position,SNP_percent,regular)
  assert_that(nrow(R_MAP) > 0)
  R_bases <- nrow(R_MAP)
  print(paste("Numbers of R bases above cutoff is ",R_bases))
  C_MAP <- read_MAP(paste0("~/Documents/2021-Fiji-Ct_paper/MAP_results/MAP_files/",subject_id,"C_MAP.txt"),"C") %>%
    #filter(regular > -1 & regular < 1) %>%
    filter(Coverage > cov_cutoff) %>%
    select(Position,SNP_percent,regular)
  assert_that(nrow(C_MAP) > 0)
  C_bases <- nrow(C_MAP)
  print(paste("Numbers of C bases above cutoff is ",C_bases))
  V_MAP <- read_MAP(paste0("~/Documents/2021-Fiji-Ct_paper/MAP_results/MAP_files/",subject_id,"V_MAP.txt"),"V") %>%
    #filter(regular > -1 & regular < 1) %>%
    filter(Coverage > cov_cutoff) %>%
    select(Position,SNP_percent,regular)
  assert_that(nrow(V_MAP) > 0)
  V_bases <- nrow(V_MAP)
  print(paste("Numbers of V bases above cutoff is ", V_bases))
  
# Join all 3 together and parse out the rare sites (ie < 90% in any of the 3)  
  Total_iSNPs <- inner_join(R_MAP,V_MAP, by = "Position") %>%
    inner_join(.,C_MAP, by = "Position") %>%
    filter(!(SNP_percent.x > 90 & SNP_percent.y > 90 & SNP_percent > 90)) %>%
    mutate(SNP = case_when(
      Position %in% rare_iSNPs ~ "rare_iSNP",
      Position %in% common_iSNPs ~ "common_iSNP",
      TRUE ~ "other_SNP"
    )) %>%
    mutate(fastGEAR = ifelse(Position %in% c(R_recomb,C_recomb,V_recomb),"fastGEAR_block","other")) %>%
    add_column(subject_id = subject_id)
  colnames(Total_iSNPs) <- c("Position","R_percent","R_regular","V_percent","V_regular","C_percent","C_regular","SNP","fastGEAR_block", "subject_id")
  outfile2 <- paste0("./iSNPs_by_subject/iSNPs_subject_",subject_id,".tsv")
  write_tsv(Total_iSNPs,outfile2)
  assert_that(nrow(Total_iSNPs) > 0)
  print(paste("Number of iSNP positions = ",nrow(Total_iSNPs)))
  #at this stage could get the same results as first pass by filtering for regular = 0
  
  # R only fixed SNPs
  RonlyF <- Total_iSNPs %>%
    filter(R_percent < 10 & V_percent >= 90 & C_percent >= 90) %>%
    nrow()
  RonlyFrare <- Total_iSNPs %>%
    filter(R_percent < 10 & V_percent >= 90 & C_percent >= 90) %>%
    filter(SNP == "other_SNP" | SNP == "rare_iSNP") %>%
    nrow()
  RonlyFrecomb <- Total_iSNPs %>%
    filter(R_percent < 10 & V_percent >= 90 & C_percent >= 90) %>%
    filter(fastGEAR_block == "fastGEAR_block") %>%
    nrow()
  print(paste("R only fixed SNPs (RonlyF) = ",RonlyF))
  print(paste("R only rare fixed SNPs (RonlyFrare) = ",RonlyFrare))
  print(paste("R only fixed SNPs in recombinant regions (RonlyFrecomb) = ",RonlyFrecomb))
  
  #R only iSNPs
  RonlyRI <- Total_iSNPs %>%
    filter(R_percent >= 10 & R_percent < 90 & V_percent >= 90 & C_percent >= 90) %>%
    nrow()
  RonlyRIrare <- Total_iSNPs %>%
    filter(R_percent >= 10 & R_percent < 90 & V_percent >= 90 & C_percent >= 90) %>%
    filter(SNP == "other_SNP" | SNP == "rare_iSNP") %>%
    nrow()
  RonlyRIrecomb <- Total_iSNPs %>%
    filter(R_percent >= 10 & R_percent < 90 & V_percent >= 90 & C_percent >= 90) %>%
    filter(fastGEAR_block == "fastGEAR_block") %>%
    nrow()
  print(paste("R only iSNPs (RonlyRI) = ",RonlyRI))
  print(paste("R only rare iSNPs (RonlyRIrare) = ",RonlyRIrare))
  print(paste("R only iSNPs in recombinant regions (RonlyRIrecomb) = ",RonlyRIrecomb))
  
  #C only fixed
  ConlyF <-  Total_iSNPs %>%
    filter(C_percent < 10 & V_percent >= 90 & R_percent >= 90) %>%
    nrow()
  ConlyFrare <- Total_iSNPs %>%
    filter(C_percent < 10 & V_percent >= 90 & R_percent >= 90) %>%
    filter(SNP == "other_SNP" | SNP == "rare_iSNP") %>%
    nrow()
  ConlyFrecomb <- Total_iSNPs %>%
    filter(C_percent < 10 & V_percent >= 90 & R_percent >= 90) %>%
    filter(fastGEAR_block == "fastGEAR_block") %>%
    nrow()
  print(paste("C only fixed SNPs (ConlyF) = ",ConlyF))
  print(paste("C only rare fixed SNPs (ConlyFrare) = ",ConlyFrare))
  print(paste("C only fixed SNPs in recombinant regions (ConlyFrecomb) = ",ConlyFrecomb))
  
  #C only iSNPs
  ConlyRI <- Total_iSNPs %>%
    filter(C_percent >= 10 & C_percent < 90 & V_percent >= 90 & R_percent >= 90) %>%
    nrow()
  ConlyRIrare <- Total_iSNPs %>%
    filter(C_percent >= 10 & C_percent < 90 & V_percent >= 90 & R_percent >= 90) %>%
    filter(SNP == "other_SNP" | SNP == "rare_iSNP") %>%
    nrow()
  ConlyRIrecomb <- Total_iSNPs %>%
    filter(C_percent >= 10 & C_percent < 90 & V_percent >= 90 & R_percent >= 90) %>%
    filter(fastGEAR_block == "fastGEAR_block") %>%
    nrow()
  print(paste("C only iSNPs (ConlyRI) = ",ConlyRI))
  print(paste("C only rare iSNPs (ConlyRIrare) = ",ConlyRIrare))
  print(paste("C only iSNPs in recombinant regions (ConlyRIrecomb) = ",ConlyRIrecomb))
  
  #V only fixed
  VonlyF <- Total_iSNPs %>%
    filter(V_percent < 10 & C_percent >= 90 & R_percent >= 90) %>%
    nrow()
  VonlyFrare <- Total_iSNPs %>%
    filter(V_percent < 10 & C_percent >= 90 & R_percent >= 90) %>%
    filter(SNP == "other_SNP" | SNP == "rare_iSNP") %>%
    nrow()
  VonlyFrecomb <- Total_iSNPs %>%
    filter(V_percent < 10 & C_percent >= 90 & R_percent >= 90) %>%
    filter(fastGEAR_block == "fastGEAR_block") %>%
    nrow()
  print(paste("V only fixed SNPs (VonlyF) = ",VonlyF))
  print(paste("V only rare fixed SNPs (VonlyFrare) = ",VonlyFrare))
  print(paste("V only fixed SNPs in recombinant regions (VonlyFrecombb) = ",VonlyFrecomb))
  
  #V only iSNPs
  VonlyRI <- Total_iSNPs %>%
    filter(V_percent >= 10 & V_percent < 90 & C_percent >= 90 & R_percent >= 90) %>%
    nrow()
  VonlyRIrare <- Total_iSNPs %>%
    filter(V_percent >= 10 & V_percent < 90 & C_percent >= 90 & R_percent >= 90) %>%
    filter(SNP == "other_SNP" | SNP == "rare_iSNP") %>%
    nrow()
  VonlyRIrecomb <- Total_iSNPs %>%
    filter(V_percent >= 10 & V_percent < 90 & C_percent >= 90 & R_percent >= 90) %>%
    filter(fastGEAR_block == "fastGEAR_block") %>%
    nrow()
  print(paste("V only iSNPs (VonlyRI) = ",VonlyRI))
  print(paste("V only rare iSNPs (VonlyRIrare) = ",VonlyRIrare))
  print(paste("V only iSNPs in recombinant regions (VonlyRIrecomb) = ",VonlyRIrecomb))
  
  #RC only fixed
  RConlyF <- Total_iSNPs %>%
    filter(V_percent >= 90 & C_percent < 10 & R_percent < 10) %>%
    nrow()
  RConlyFrare <- Total_iSNPs %>%
    filter(V_percent >= 90 & C_percent < 10 & R_percent < 10) %>%
    filter(SNP == "other_SNP" | SNP == "rare_iSNP") %>%
    nrow()
  RConlyFrecomb <- Total_iSNPs %>%
    filter(V_percent >= 90 & C_percent < 10 & R_percent < 10) %>%
    filter(fastGEAR_block == "fastGEAR_block") %>%
    nrow()
  print(paste("RC only fixed SNPs (RConlyF) = ",RConlyF))
  print(paste("RC only rare fixed SNPs (RConlyFrare) = ",RConlyFrare))
  print(paste("RC only fixed SNPs in recombinant regions (RConlyFrecomb) = ",RConlyFrecomb))
  
  #RC only iSNPs
  RConlyRI <- Total_iSNPs %>%
    filter(V_percent >= 90 & C_percent >= 10 & C_percent < 90 & R_percent >= 10 & R_percent < 90) %>%
    nrow()
  RConlyRIrare <- Total_iSNPs %>%
    filter(V_percent >= 90 & C_percent >= 10 & C_percent < 90 & R_percent >= 10 & R_percent < 90) %>%
    filter(SNP == "other_SNP" | SNP == "rare_iSNP") %>%
    nrow()
  RConlyRIrecomb <- Total_iSNPs %>%
    filter(V_percent >= 90 & C_percent >= 10 & C_percent < 90 & R_percent >= 10 & R_percent < 90) %>%
    filter(fastGEAR_block == "fastGEAR_block") %>%
    nrow()
  print(paste("RC only iSNPs (RConlyRI) = ",RConlyRI))
  print(paste("RC only rare iSNPs (RConlyRIrare) = ",RConlyRIrare))
  print(paste("RC only iSNPs in recombinant regions (RConlyRIrecomb) = ",RConlyRIrecomb))
  
  #RV only fixed
  RVonlyF <- Total_iSNPs %>%
    filter(C_percent >= 90 & V_percent < 10 & R_percent < 10) %>%
    nrow()
  RVonlyFrare <- Total_iSNPs %>%
    filter(C_percent >= 90 & V_percent < 10 & R_percent < 10) %>%
    filter(SNP == "other_SNP" | SNP == "rare_iSNP") %>%
    nrow()
  RVonlyFrecomb <- Total_iSNPs %>%
    filter(C_percent >= 90 & V_percent < 10 & R_percent < 10) %>%
    filter(fastGEAR_block == "fastGEAR_block") %>%
    nrow()
  print(paste("RV only fixed SNPs (RVonlyF) = ",RVonlyF))
  print(paste("RV only rare fixed SNPs (RVonlyFrare) = ",RVonlyFrare))
  print(paste("RV only fixed SNPs in recombinant regions (RVonlyFrecomb) = ",RVonlyFrecomb))
  
  #RV only iSNPs
  RVonlyRI <- Total_iSNPs %>%
    filter(C_percent >= 90 & V_percent >= 10 & V_percent < 90 & R_percent >= 10 & R_percent < 90) %>%
    nrow()
  RVonlyRIrare <- Total_iSNPs %>%
    filter(C_percent >= 90 & V_percent >= 10 & V_percent < 90 & R_percent >= 10 & R_percent < 90) %>%
    filter(SNP == "other_SNP" | SNP == "rare_iSNP") %>%
    nrow()
  RVonlyRIrecomb <- Total_iSNPs %>%
    filter(C_percent >= 90 & V_percent >= 10 & V_percent < 90 & R_percent >= 10 & R_percent < 90) %>%
    filter(fastGEAR_block == "fastGEAR_block") %>%
    nrow()
  print(paste("RV only iSNPs (RVonlyRI) = ",RVonlyRI))
  print(paste("RV only rare iSNPs (RVonlyRIrare) = ",RVonlyRIrare))
  print(paste("RV only iSNPs in recombinant regions (RVonlyRIrecomb) = ",RVonlyRIrecomb))
  
  #VC only fixed
  VConlyF <- Total_iSNPs %>%
    filter(R_percent >= 90 & V_percent < 10 & C_percent < 10) %>%
    nrow()
  VConlyFrare <- Total_iSNPs %>%
    filter(R_percent >= 90 & V_percent < 10 & C_percent < 10) %>%
    filter(SNP == "other_SNP" | SNP == "rare_iSNP") %>%
    nrow()
  VConlyFrecomb <- Total_iSNPs %>%
    filter(R_percent >= 90 & V_percent < 10 & C_percent < 10) %>%
    filter(fastGEAR_block == "fastGEAR_block") %>%
    nrow()
  print(paste("VC only fixed SNPs (VConlyF) = ",VConlyF))
  print(paste("VC only rare fixed SNPs (VConlyFrare) = ",VConlyFrare))
  print(paste("VC only fixed SNPs in recombinant regions (VConlyFrecomb) = ",VConlyFrecomb))
  
  #VC only iSNPs
  VConlyRI <- Total_iSNPs %>%
    filter(R_percent >= 90 & V_percent >= 10 & V_percent < 90 & C_percent >= 10 & C_percent < 90) %>%
    nrow()
  VConlyRIrare <- Total_iSNPs %>%
    filter(R_percent >= 90 & V_percent >= 10 & V_percent < 90 & C_percent >= 10 & C_percent < 90) %>%
    filter(SNP == "other_SNP" | SNP == "rare_iSNP") %>%
    nrow()
  VConlyRIrecomb <- Total_iSNPs %>%
    filter(R_percent >= 90 & V_percent >= 10 & V_percent < 90 & C_percent >= 10 & C_percent < 90) %>%
    filter(fastGEAR_block == "fastGEAR_block") %>%
    nrow()
  print(paste("VC only iSNPs (VConlyRI) = ",VConlyRI))
  print(paste("VC only rare iSNPs (VConlyRIrare) = ",VConlyRIrare))
  print(paste("VC only iSNPs in recombinant regions (VConlyRIrecomb) = ",VConlyRIrecomb))
  
  #all_fixed
  AllF <- Total_iSNPs %>%
    filter(R_percent <10 & V_percent < 10 & C_percent < 10) %>%
    nrow()
  AllFrare <- Total_iSNPs %>%
    filter(R_percent <10 & V_percent < 10 & C_percent < 10) %>%
    filter(SNP == "other_SNP" | SNP == "rare_iSNP") %>%
    nrow()
  AllFrecomb <- Total_iSNPs %>%
    filter(R_percent <10 & V_percent < 10 & C_percent < 10) %>%
    filter(fastGEAR_block == "fastGEAR_block") %>%
    nrow()
  print(paste("All fixed SNPs (AllF) = ",AllF))
  print(paste("All rare fixed SNPs (AllFrare) = ",AllFrare))
  print(paste("All fixed SNPs in recombinant regions (AllFrecomb) = ",AllFrecomb))
  
  #all iSNPs
  AllRI <- Total_iSNPs %>%
    filter(R_percent >= 10 & R_percent < 90 & V_percent >= 10 & V_percent < 90 & C_percent >= 10 & C_percent < 90) %>%
    nrow()
  AllRIrare <- Total_iSNPs %>%
    filter(R_percent >= 10 & R_percent < 90 & V_percent >= 10 & V_percent < 90 & C_percent >= 10 & C_percent < 90) %>%
    filter(SNP == "other_SNP" | SNP == "rare_iSNP") %>%
    nrow()
  AllRIrecomb <- Total_iSNPs %>%
    filter(R_percent >= 10 & R_percent < 90 & V_percent >= 10 & V_percent < 90 & C_percent >= 10 & C_percent < 90) %>%
    filter(fastGEAR_block == "fastGEAR_block") %>%
    nrow()
  print(paste("All iSNPs (AllRI) = ",AllRI))
  print(paste("All rare  iSNPs (AllRIrare) = ",AllRIrare))
  print(paste("All iSNPs in recombinant regions (AllRIrecomb) = ",AllRIrecomb))
  
  #fixed in one compartment, intermediate in others
  
  RonlyF_VIrare <- Total_iSNPs %>%
    filter(R_percent < 10 & V_percent >= 10 & V_percent < 90 & C_percent >= 90) %>%
    filter(SNP == "other_SNP" | SNP == "rare_iSNP") %>%
    nrow()
  
  RonlyF_CIrare <- Total_iSNPs %>%
    filter(R_percent < 10 & C_percent >= 10 & C_percent < 90 & V_percent >= 90) %>%
    filter(SNP == "other_SNP" | SNP == "rare_iSNP") %>%
    nrow()
  
  RonlyF_CandVIrare <- Total_iSNPs %>%
    filter(R_percent < 10 & C_percent >= 10 & C_percent < 90 & V_percent >= 10 & V_percent < 90) %>%
    filter(SNP == "other_SNP" | SNP == "rare_iSNP") %>%
    nrow()
  
  ConlyF_VIrare <- Total_iSNPs %>%
    filter(C_percent < 10 & V_percent >= 10 & V_percent < 90 & R_percent >= 90) %>%
    filter(SNP == "other_SNP" | SNP == "rare_iSNP") %>%
    nrow()
  
  ConlyF_RIrare <- Total_iSNPs %>%
    filter(C_percent < 10 & R_percent >= 10 & R_percent < 90 & V_percent >= 90) %>%
    filter(SNP == "other_SNP" | SNP == "rare_iSNP") %>%
    nrow()
  
  ConlyF_RandVIrare <- Total_iSNPs %>%
    filter(C_percent < 10 & R_percent >= 10 & R_percent < 90 & V_percent >= 10 & V_percent < 90) %>%
    filter(SNP == "other_SNP" | SNP == "rare_iSNP") %>%
    nrow()
  
  VonlyF_RIrare <- Total_iSNPs %>%
    filter(V_percent < 10 & R_percent >= 10 & R_percent < 90 & C_percent >= 90) %>%
    filter(SNP == "other_SNP" | SNP == "rare_iSNP") %>%
    nrow()
  
  VonlyF_CIrare <- Total_iSNPs %>%
    filter(V_percent < 10 & C_percent >= 10 & C_percent < 90 & R_percent >= 90) %>%
    filter(SNP == "other_SNP" | SNP == "rare_iSNP") %>%
    nrow()
  
  VonlyF_CandRIrare <- Total_iSNPs %>%
    filter(V_percent < 10 & R_percent >= 10 & R_percent < 90 & C_percent >= 10 & C_percent < 90) %>%
    filter(SNP == "other_SNP" | SNP == "rare_iSNP") %>%
    nrow()
  
  # fixed in two sites, intermediate in the other
  
  RCF_VIrare <- Total_iSNPs %>%
    filter(R_percent < 10 & V_percent >= 10 & V_percent < 90 & C_percent < 10) %>%
    filter(SNP == "other_SNP" | SNP == "rare_iSNP") %>%
    nrow()
  
  RVF_CIrare <- Total_iSNPs %>%
    filter(R_percent < 10 & C_percent >= 10 & C_percent < 90 & V_percent < 10) %>%
    filter(SNP == "other_SNP" | SNP == "rare_iSNP") %>%
    nrow()
  
  CVF_RIrare <- Total_iSNPs %>%
    filter(V_percent < 10 & R_percent >= 10 & R_percent < 90 & C_percent < 10) %>%
    filter(SNP == "other_SNP" | SNP == "rare_iSNP") %>%
    nrow()
  
  
  print(paste("R is fixed, V only is intermediate, and rare (RonlyF_VIrare) = ",RonlyF_VIrare))
  print(paste("R is fixed, C only is intermediate, and rare (RonlyF_CIrare) = ",RonlyF_CIrare))
  print(paste("R is fixed, V and C are intermediate, and rare (RonlyF_CandVIrare) = ",RonlyF_CandVIrare))
  print(paste("C is fixed, V only is intermediate, and rare (ConlyF_VIrare) = ",ConlyF_VIrare))
  print(paste("C is fixed, R only is intermediate, and rare (ConlyF_RIrare) = ",ConlyF_RIrare))
  print(paste("C is fixed, R and V are intermediate, and rare (ConlyF_RandVIrare) = ",ConlyF_RandVIrare))
  print(paste("V is fixed, C only is intermediate, and rare (VonlyF_CIrare) = ",VonlyF_CIrare))
  print(paste("V is fixed, R only is intermediate, and rare (VonlyF_RIrare) = ",VonlyF_RIrare))
  print(paste("V is fixed, R and C are intermediate, and rare (VonlyF_CandRIrare) = ",VonlyF_CandRIrare))
  
  print(paste("R and C are fixed, V is intermediate, and rare (RCF_VIrare) = ",RCF_VIrare))
  print(paste("R and V are fixed, C is intermediate, and rare (RVF_CIrare) = ",RVF_CIrare))
  print(paste("V and C are fixed, R is intermediate, and rare (CVF_RIrare) = ",CVF_RIrare))
  
  #summarise and write SNP table
  result_table <- data.frame(subject_name = subject_id, R_bases, C_bases, V_bases,R_recomb_bases,C_recomb_bases,V_recomb_bases,AllF,AllFrare,AllFrecomb,AllRI,AllRIrare,AllRIrecomb,RonlyF,RonlyFrare,RonlyFrecomb,RonlyRI,RonlyRIrare,RonlyRIrecomb,ConlyF,ConlyFrare,ConlyFrecomb,ConlyRI,ConlyRIrare,ConlyRIrecomb,VonlyF,VonlyFrare,VonlyFrecomb,VonlyRI,VonlyRIrare,VonlyRIrecomb,RConlyF,RConlyFrare,RConlyFrecomb,RConlyRI,RConlyRIrare,RConlyRIrecomb,RVonlyF,RVonlyFrare,RVonlyFrecomb,RVonlyRI,RVonlyRIrare,RVonlyRIrecomb,VConlyF,VConlyFrare,VConlyFrecomb,VConlyRI,VConlyRIrare,VConlyRIrecomb)
  result_table <- data.frame(result_table,RonlyF_VIrare,RonlyF_CIrare,RonlyF_CandVIrare,ConlyF_VIrare,ConlyF_RIrare,ConlyF_RandVIrare,VonlyF_CIrare,VonlyF_RIrare,VonlyF_CandRIrare)
  result_table <- data.frame(result_table,RCF_VIrare,RVF_CIrare,CVF_RIrare)
  
  outfile <- paste0("./iSNP_analysis_table/iSNPs_",subject_id,".tsv")
  write_tsv(result_table,outfile)
  
} #end SNP_analysis_pipeline function block


SNP_analysis_for_pub <- function(subject_id) {
  cov_cutoff <- 10
  print(paste("SUBJECT_ID is",subject_id))
  print(paste("Coverage cutoff =",cov_cutoff))
 
  ## replace with the table created from the above function 
  library(tidyverse)
  source("MAP_functions.R")
  library(assertthat)
  
  in_file_name <- paste0("./iSNPs_by_subject/iSNPs_subject_",subject_id,".tsv")
  Total_iSNPs <- read.delim(in_file_name)
  assert_that(nrow(Total_iSNPs) > 0)
  print(paste("Number of iSNP positions = ",nrow(Total_iSNPs)))
  
  # #Tryptich 1
  # library(cowplot)
  # RCplot <- ggplot(filter(Total_iSNPs,!(R_percent > 90 & C_percent > 90)), aes(x=R_percent, y=C_percent, color = SNP)) + 
  #   geom_point() +
  #   xlab("Ref allele % rectal") +
  #   ylab("Ref allele % cervical") +
  #   xlim(0,100) + 
  #   scale_color_manual(values = c("gainsboro","orange","red","blue")) +
  #   theme_bw() 
  # RVplot <- ggplot(filter(Total_iSNPs,!(R_percent > 90 & V_percent > 90)), aes(x=R_percent, y=V_percent, color = SNP)) + 
  #   geom_point() +
  #   xlab("Ref allele % rectal") +
  #   ylab("Ref allele % vaginal") +
  #   xlim(0,100) +
  #   scale_color_manual(values = c("gainsboro","orange","red","blue")) +
  #   theme_bw() 
  # VCplot <- ggplot(filter(Total_iSNPs,!(V_percent > 90 & C_percent > 90)), aes(x=V_percent, y=C_percent, color = SNP)) + 
  #   geom_point() +
  #   xlab("Ref allele % vaginal") +
  #   ylab("Ref allele % cervical") +
  #   xlim(0,100) +
  #   scale_color_manual(values = c("gainsboro","orange","red","blue")) +
  #   theme_bw() 
  # try1_title <- ggdraw() +
  #   theme(plot.margin = margin(0, 0, 0, 7))
  # tryptich1 <- plot_grid(
  #   RCplot, RVplot, VCplot,
  #   labels = c("e","f","g")
  # )
  
  
  #canonical SNP analysis
  T1T2_denovo_snps <- read_lines("T1T2_snps_denovo.list")
  
  R_MAP <- read_MAP(paste0("~/Documents/2021-Fiji-Ct_paper/MAP_results/MAP_files/",subject_id,"R_MAP.txt"),"R") %>%
    #filter(regular > -1 & regular < 1) %>%
    filter(Coverage > cov_cutoff) %>%
    select(Position,SNP_percent)
  assert_that(nrow(R_MAP) > 0)
  R_bases <- nrow(R_MAP)
  print(paste("Numbers of R bases above cutoff is ",R_bases))
  C_MAP <- read_MAP(paste0("~/Documents/2021-Fiji-Ct_paper/MAP_results/MAP_files/",subject_id,"C_MAP.txt"),"C") %>%
    #filter(regular > -1 & regular < 1) %>%
    filter(Coverage > cov_cutoff) %>%
    select(Position,SNP_percent)
  assert_that(nrow(C_MAP) > 0)
  C_bases <- nrow(C_MAP)
  print(paste("Numbers of C bases above cutoff is ",C_bases))
  V_MAP <- read_MAP(paste0("~/Documents/2021-Fiji-Ct_paper/MAP_results/MAP_files/",subject_id,"V_MAP.txt"),"V") %>%
    #filter(regular > -1 & regular < 1) %>%
    filter(Coverage > cov_cutoff) %>%
    select(Position,SNP_percent)
  assert_that(nrow(V_MAP) > 0)
  
  
  T1T2_denovo_df <- full_join(R_MAP,V_MAP, by = "Position") %>%
    full_join(.,C_MAP, by = "Position") %>%
    filter(Position %in% T1T2_denovo_snps)
  colnames(T1T2_denovo_df) <- c("Position","Rectum","Vagina","Endocervix")
  
  T1T2_denovo_df_long <- gather(T1T2_denovo_df,"Body_site","SNP_percent",-Position)
  #write_tsv(T1T2_denovo_df_long,file = "T1T2_denovo_df_long.tsv")
 
  library(RColorBrewer)
  cust_pal <- brewer.pal(3,"Accent")
  snps_ref1 <- ggplot(filter(T1T2_denovo_df_long,Body_site == "Endocervix"), aes(x= Position, y= SNP_percent, color = Body_site)) +
    geom_point(alpha = 0.7, size = 0.5) + 
    scale_color_manual(values = cust_pal[1]) +
    xlim(1,1042504)+
    ylab("% ref") +
    theme_bw() +
    theme(axis.text = element_text(size=9))
  
  snps_ref2 <- ggplot(filter(T1T2_denovo_df_long,Body_site == "Rectum"), aes(x= Position, y= SNP_percent, color = Body_site)) +
    geom_point(alpha = 0.7, size = 0.5) + 
    scale_color_manual(values = cust_pal[2]) +
    xlim(1,1042504)+
    labs(y = "% ref") +
    theme_bw() +
    theme(axis.text = element_text(size=9))
  
  snps_ref3 <- ggplot(filter(T1T2_denovo_df_long,Body_site == "Vagina"), aes(x= Position, y= SNP_percent, color = Body_site)) +
    geom_point(alpha = 0.7, size = 0.5) + 
    scale_color_manual(values = cust_pal[3]) +
    xlim(1,1042504)+
    labs(y = "% ref") +
    theme_bw()+
    theme(axis.text = element_text(size=9))
  
  snps_ref4 <- ggplot(T1T2_denovo_df_long, aes(y=SNP_percent,x=Body_site,color = Body_site)) +
    geom_boxplot() + 
    scale_color_brewer(palette = "Accent") +
    scale_y_continuous(name="% ref", limits=c(0, 100)) +
    theme_bw() +
    theme(legend.position="none") 
  
  CSNP_plot2 <- plot_grid(snps_ref1,snps_ref2,snps_ref3, ncol = 1, align = 'v',labels = c('a', 'b','c'),label_y = 0.95)
  
  CSNP_plot3 <- plot_grid(snps_ref4, ncol = 2, nrow=2, align = 'v',labels = c('d'),label_y = 0.98,label_x = 0.015) 
  
  #combined_plot <- plot_grid(CSNP_plot2,CSNP_plot3,tryptich1,align = 'v')
  plot_filename = paste0(subject_id,"_combined.jpg")
  #ggsave2(filename=plot_filename,device = "jpeg",dpi = 300,units = "in", width = 8, height = 3)
  
  
  # return all the plots (except tryptich)
  return(list(CSNP_plot2,CSNP_plot3))
  print("")
} #end SNP_analysis_for_pub function block


#### ORIGINAL SNP pipeline figures

# #Tryptich 1
# library(cowplot)
# RCplot <- ggplot(filter(Total_iSNPs,!(R_percent > 90 & C_percent > 90)), aes(x=R_percent, y=C_percent, color = SNP)) + 
#   geom_point() +
#   xlab("Ref allele % rectal") +
#   ylab("Ref allele % cervical") +
#   xlim(0,100) + 
#   scale_color_manual(values = c("gainsboro","orange","red","blue")) +
#   theme_bw() 
# RVplot <- ggplot(filter(Total_iSNPs,!(R_percent > 90 & V_percent > 90)), aes(x=R_percent, y=V_percent, color = SNP)) + 
#   geom_point() +
#   xlab("Ref allele % rectal") +
#   ylab("Ref allele % vaginal") +
#   xlim(0,100) +
#   scale_color_manual(values = c("gainsboro","orange","red","blue")) +
#   theme_bw() 
# VCplot <- ggplot(filter(Total_iSNPs,!(V_percent > 90 & C_percent > 90)), aes(x=V_percent, y=C_percent, color = SNP)) + 
#   geom_point() +
#   xlab("Ref allele % vaginal") +
#   ylab("Ref allele % cervical") +
#   xlim(0,100) +
#   scale_color_manual(values = c("gainsboro","orange","red","blue")) +
#   theme_bw() 
# try1_title <- ggdraw() +
#   draw_label(paste(subject_id), 
#              fontface = 'bold',
#              x = 0,
#              hjust = 0) +
#   theme(plot.margin = margin(0, 0, 0, 7))
# tryptich1 <- plot_grid(
#   try1_title,RCplot, RVplot, VCplot,
#   labels = c("","a","b","c")
# )
# 
# #Tryptich 2
# RCplot2 <- ggplot(filter(Total_iSNPs,fastGEAR_block == "fastGEAR_block"), aes(x=R_percent, y=C_percent, color = fastGEAR_block)) + 
#   geom_point() +
#   xlab("Ref allele % rectal") +
#   ylab("Ref allele % cervical") +
#   xlim(0,100) + 
#   ylim(0,100) +
#   scale_color_manual(values = c("green","gainsboro")) +
#   theme_bw() 
# RVplot2 <- ggplot(filter(Total_iSNPs,fastGEAR_block == "fastGEAR_block"),aes(x=R_percent, y=V_percent, color = fastGEAR_block)) + 
#   geom_point() +
#   xlab("Ref allele % rectal") +
#   ylab("Ref allele % vaginall") +
#   xlim(0,100) + 
#   ylim(0,100) +
#   scale_color_manual(values = c("green","gainsboro")) +
#   theme_bw() 
# VCplot2 <- ggplot(filter(Total_iSNPs,fastGEAR_block == "fastGEAR_block"), aes(x=V_percent, y=C_percent, color = fastGEAR_block)) + 
#   geom_point() +
#   xlab("Ref allele % vaginal") +
#   ylab("Ref allele % cervical") +
#   xlim(0,100) + 
#   ylim(0,100) +
#   scale_color_manual(values = c("green","gainsboro")) +
#   theme_bw()
# try2_title <- ggdraw() +
#   draw_label(paste("Recombination: subject",subject_id), 
#              fontface = 'bold',
#              x = 0,
#              hjust = 0) +
#   theme(plot.margin = margin(0, 0, 0, 7))
# tryptich2 <- plot_grid(
#   try2_title,RCplot2, RVplot2, VCplot2,
#   labels = c("","d","e","f")
# )
# #plots of each pair in Tryptich 1
# RCplot_title <- RCplot +
#   labs(title = paste("Minor SNPs in R(x) versus C(y). Subject:",subject_id))
# RVplot_title <- RVplot +
#   labs(title = paste("Minor SNPs in R(x) versus V(y). Subject:",subject_id))
# VCplot_title <- VCplot +
#   labs(title = paste("Minor SNPs in V(x) versus C(y). Subject:",subject_id))
# 
# #parallel coords plot
# library(GGally)
# parcoord_title <- paste("Subject",subject_id,", number iSNPs =",nrow(Total_iSNPs))
# parplot <- ggparcoord(data = Total_iSNPs,
#                       columns = 2:4,
#                       scale = "globalminmax",
#                       alphaLines = 0.05,
#                       title = parcoord_title)
# 
# #canonical SNP analysis
# T1T2_denovo_snps <- read_lines("T1T2_snps_denovo.list")
# 
# T1T2_denovo_df <- full_join(R_MAP,V_MAP, by = "Position") %>%
#   full_join(.,C_MAP, by = "Position") %>%
#   filter(Position %in% T1T2_denovo_snps)
# colnames(T1T2_denovo_df) <- c("Position","Rectum","Vagina","Endocervix")
# 
# T1T2_denovo_df_long <- gather(T1T2_denovo_df,"Body_site","SNP_percent",-Position)
# 
# CSNP_plot1 <- ggparcoord(data = T1T2_denovo_df,
#                          columns = 2:4,
#                          scale = "globalminmax",
#                          alphaLines = 0.2,
#                          title = paste0("cDenovo_SNPs for" ,subject_id))
# 
# library(RColorBrewer)
# cust_pal <- brewer.pal(3,"Accent")
# snps_ref1 <- ggplot(filter(T1T2_denovo_df_long,Body_site == "Endocervix"), aes(x= Position, y= SNP_percent, color = Body_site)) +
#   geom_point(alpha = 0.7, size = 0.5) + 
#   scale_color_manual(values = cust_pal[1]) +
#   xlim(1,1042504)+
#   ylab("% reference") +
#   theme_bw()
# 
# snps_ref2 <- ggplot(filter(T1T2_denovo_df_long,Body_site == "Rectum"), aes(x= Position, y= SNP_percent, color = Body_site)) +
#   geom_point(alpha = 0.7, size = 0.5) + 
#   scale_color_manual(values = cust_pal[2]) +
#   xlim(1,1042504)+
#   ylab("% reference") +
#   theme_bw()
# 
# snps_ref3 <- ggplot(filter(T1T2_denovo_df_long,Body_site == "Vagina"), aes(x= Position, y= SNP_percent, color = Body_site)) +
#   geom_point(alpha = 0.7, size = 0.5) + 
#   scale_color_manual(values = cust_pal[3]) +
#   xlim(1,1042504)+
#   ylab("% reference") +
#   theme_bw()
# 
# snps_ref4 <- ggplot(T1T2_denovo_df_long, aes(x=SNP_percent,color = Body_site)) +
#   geom_boxplot() + 
#   scale_color_brewer(palette = "Accent") +
#   xlab("Percent_cDenovo_reference") +
#   theme_bw() +
#   theme(axis.text.y = element_blank())
# 
# CSNP_plot2 <- plot_grid(snps_ref1,snps_ref2,snps_ref3, ncol = 1, align = 'v',labels = c('d', 'e','f'),vjust = 0.1)
# 
# CSNP_plot3 <- plot_grid(snps_ref4, ncol = 1, align = 'v',labels = c('g'),vjust = 0.2,hjust = 0.2) 
# 
# print(paste("Number cDenovo Ref SNPs for C = ",nrow(filter(T1T2_denovo_df,Endocervix > 90))))
# print(paste("Number cDenovo Non-Ref SNPs for C = ",nrow(filter(T1T2_denovo_df,Endocervix < 10))))
# print(paste("Number cDenovo iSNPs for C = ",nrow(filter(T1T2_denovo_df,Endocervix >= 10 & Endocervix <= 90))))
# 
# print(paste("Number cDenovo Ref SNPs for R = ",nrow(filter(T1T2_denovo_df,Rectum > 90))))
# print(paste("Number cDenovo Non-Ref SNPs for R = ",nrow(filter(T1T2_denovo_df,Rectum < 10))))
# print(paste("Number cDenovo iSNPs for R = ",nrow(filter(T1T2_denovo_df,Rectum >= 10 & Rectum <= 90))))
# 
# print(paste("Number cDenovo Ref SNPs for V = ",nrow(filter(T1T2_denovo_df,Vagina > 90))))
# print(paste("Number cDenovo Non-Ref SNPs for V = ",nrow(filter(T1T2_denovo_df,Vagina < 10))))
# print(paste("Number cDenovo iSNPs for V = ",nrow(filter(T1T2_denovo_df,Vagina >= 10 & Vagina <= 90))))
# 
# 
# # return all the plots
# return(list(tryptich1, tryptich2,RCplot_title,RVplot_title,VCplot_title,parplot,CSNP_plot1,CSNP_plot2,CSNP_plot3))
# print("")
# 
# 
# 
