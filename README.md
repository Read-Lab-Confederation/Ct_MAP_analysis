# Ct MAP analysis

Code for analysis of minor allele distributions in Chlamydia trachomatis pull-down genome sequencing project

The original input files were created from mpileup output of Ct genomes mapped against the serovar D genome

```MAP_file_pipeline``` script was used to create plots 
Inputs - mpileup files
Outputs ```MAP_files```.
One of the columns in the MAP_file report is called "regular", which reports the discrepancy between the total coverage and the sum of the reads mapped to the referecne+ other allels.  This is essentially a report of possibly mismatched bases.

```analyzing_MAP_pipeline_results.Rmd```  
Inputs ```MAP_files``` and reports stats about iSNPs(SNVs).
Outputs the ```Pos_count.tsv``` file for the numbers of times each iSNP occurs in the 87 original high_quality samples.

```snippy_analysis.Rmd``` identifies the canonical (cSNPs) and 
outputs lists of positions that discriminate T1 and T2 clades including ```T1T2_denovo_snps```

```fastGEAR_ancrestral_recombination``` and ```fastgear_recent_recombination``` pipelines analyse the output of fastGEAR analysis and make plots and gather stats. Als produce the ```fastgear_fiji_recombs``` file of recombinant positions.

```SNP_pipelines``` analyzes SNPs for each strain
INputs are ``Pos_count.tsv```, ```fastgear_fiji_recombs```,```MAP_files```, and ```T1T2_denovo_snps```, and ```2021_08_24-Grand_Summary - MAP_analysis_summary.tsv```. 
Outputs are ```iSNPs_by_subject``` and ```iSNP_analysis_table```., and various plots

```SNP_analysis_figures``` is similar to SNP_pipelines but updated figures for publications.

```common_iSNPs.Rmd```: Makes some plots.
inputs are ```/iSNPs_by_subject/```.

```counting_suppl_table1``` scripts to analyze plasmids and ompA types in Supplementary Table 1.

