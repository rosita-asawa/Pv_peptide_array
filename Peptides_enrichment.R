# Script 2: Enrichment analyses
# coordinating custom scripts for Plos Pathogens manuscript, published October 2026
  # Proteome-wide characterization of Plasmodium vivax antigens using a high-density peptide array
  # https://journals.plos.org/plospathogens/article?id=10.1371/journal.ppat.1014621

library(tidyverse)
library(data.table)
library(plotly)
library(ggtext)
library(ggpubr)
library(seqinr)
library(rtracklayer)
library(Biostrings)

install.packages("")

# ---------------- Symptomatic Pv patient enrichment
# ---------------------- Numbers: disorder enrichment within immunoreactive segments -----------
# look for enrichment of disordered segment (>25% of segment is disordered) vs immunoreactivity
# removed pseudogenes (n = 152) as they cannot be ran on PONDR
# upload percent disordered for each 100aa protein segment in genome from PONDR
data <- fread("/Path/pondr_segments_output.txt")

    # rename protein.name column
    names(data)[names(data) == "Accession"] <- "merged"
    
    # adjust data ID to match ID in data
    # for weird genes
    others <- data %>% filter(str_detect(merged, 'PVP01_API|PVP01_MIT'))  
    
        # select first 14 characters of merged
        others <- others %>% mutate(protein.name = substr(merged, 1, 14))
        
    # for main data
    # remove weird genes from data now that they're pulled
    data <- data[!grepl("PVP01_API|PVP01_MIT", merged)] 
        
        # take first 13 characters from ID
        data <- data %>% mutate(protein.name = substr(merged, 1, 13))
      
        # merge others back into data
        data <- rbind(data, others)
        
        # select columns
        data <- data[ , c('protein.name', 'merged', 'Annotation',
                          'Sequence Length', 'Percent Disordered', 'Average Score')]
      
# add immunoreactivity        
# upload list of most reactive segment per gene
# code for list found in "Peptides_immunogenicity.R"
immunoreac <- fread('/Path/symptomatic_pv_most_reactive_seg_per_protein.txt')

    # select columns
    immunoreac <- immunoreac[,-2]
    
    # change name in data to match immunoreac
    data$protein.name <- paste0(data$protein.name, ".1-p1")

    # merge immunoreac and data
    all <- left_join(immunoreac, data, by = "protein.name")
    
    # remove empty rows (found in immunoreac but no disorder info)
    all <- na.omit(all)
    
        # 45389 segments available with BOTH immunoreactivity and disorder info

          # select immunoreactive segments
          reac <- all %>% filter(total > 6)
        
              # 4546 reactive segments
          
                  # set cut off at >25% is disordered
                  # count disordered segments in reactive
                  reac.disord.25 <- reac %>% filter(`Percent Disordered` > 24.9999999)
                  
                            #  2937 segments disordered in reactive (>25% threshold)
                      
                        # count ordered segments in reactive
                        reac.ord.25 <- reac %>% filter(`Percent Disordered` < 24.9999999)
                      
                              # 1609 segments ordered in reactive (>25% threshold)
                              # total of 4546 segments counted reactive genes - CORRECT
                          
                  # set cut off at >50% is disordered
                  # count disordered segments in reactive
                  reac.disord.50 <- reac %>% filter(`Percent Disordered` > 49.9999999)
                  
                          #  1695 segments disordered in reactive (>50% threshold)
                      
                  # count ordered segments in reactive
                  reac.ord.50 <- reac %>% filter(`Percent Disordered` < 49.9999999)
                
                      # 2851 segments ordered in reactive (>50% threshold)
                      # total of 4546 segments counted in reactive segments - CORRECT
                       
          # select NON-reactive segments
          nonreac <- all %>% filter(total < 7)
        
              # 40843 NON-reactive segments
          
                  # set cut off at >25% is disordered
                  # count disordered segments in NON-reactive
                  nonreac.disord.25 <- nonreac %>% filter(`Percent Disordered` > 24.9999999)
                  
                            #  20968 segments disordered in NON-reactive (>25% threshold)
                      
                        # count ordered segments in NON-reactive
                        nonreac.ord.25 <- nonreac %>% filter(`Percent Disordered` < 24.9999999)
                      
                              # 19875 segments ordered in NON-reactive (>25% threshold)
                              # total of 40843 segments counted in NON-reactive segments - CORRECT
                          
                  # set cut off at >50% is disordered
                  # count disordered segments in NON-reactive
                  nonreac.disord.50 <- nonreac %>% filter(`Percent Disordered` > 49.9999999)
                  
                          #  10720 segments disordered in NON-reactive (>50% threshold)
                      
                  # count ordered segments in NON-reactive
                  nonreac.ord.50 <- nonreac %>% filter(`Percent Disordered` < 49.9999999)
                
                      # 30123 segments ordered in NON-reactive (>50% threshold)
                      # total of 40843 segments counted in NON-reactive segments - CORRECT

                          # extra note: looks like 21 new genes were added to P01 genome
                          # since the array was produced (missing in reac file, found in disorder file)

# make chi square table for >25% is disordered
    # reactive segments: 4546
          # ordered: 1609
          # disordered: 2937
    # NON-reactive segments: 40843
          # ordered: 19875
          # disordered: 20968                
                          
# chi square test - disorder (>25%) vs immunoreactivity
    # create object containing values for chi square test
    # tests: proportion of disorder within reactive and nonreactive segments
    # put in this order: ordered-reac, ordered-nonreac, disorder-reac, disorder-nonreac
    chi <- matrix(c(
      1609, 19875,                     # ordered: reactive, non-reactive
      2937, 20968                      # disordered: reactive, non-reactive
      ), nrow = 2, byrow = TRUE)
    dimnames(chi) <- list(c("Ordered", "Disordered"), c("Immunoreactive", "Nonimmunoreactive"))
    
    # chi square test for disorder vs immunoreactivity
    chi <- chisq.test(chi, correct = TRUE)
          # Pearson's Chi-squared test with Yates' continuity correction
          # data:  chi
          # X-squared = 288.35, df = 1, p-value < 2.2e-16
          # IMPORTANT: 2.2e-16 is a bunk value from being blown out in R
          # pull real p value
          chi$p.value
                
            # >25% disorder p value = 1.138934e-64
          
                # for >25% is disordered    
                # question: are disordered proteins enriched within immunoreactive proteins? 
                    # reactive group
                          # ordered = 1695/4546 = 35.4%
                          # disordered = 2937/4546 = 64.6%
                    # non-reactive group
                          # ordered = 19875/40843= 48.7%
                          # disordered = 20968/40843 = 51.3%
    
                    # proportion higher of disordered proteins in immunoreactive group
                          # 64.6% / 51.3% = 1.6
                
                          # interpretation: immunoreactive proteins contain 1.6-fold higher proportion 
                          # of disordered proteins compared to non-immunoreactive proteins
          

# make chi square table for >50% is disordered
    # reactive segments: 4546
          # ordered: 2851
          # disordered: 1695
    # NON-reactive segments: 40843
          # ordered: 30123
          # disordered: 10720      
    
# chi square test - disorder (>50%) vs immunoreactivity
    # create object containing values for chi square test
    # tests: proportion of disorder within reactive and nonreactive segments
    # put in this order: ordered-reac, ordered-nonreac, disorder-reac, disorder-nonreac
    chi <- matrix(c(
      2851, 30123,                     # ordered: reactive, non-reactive
      1695, 10720                      # disordered: reactive, non-reactive
      ), nrow = 2, byrow = TRUE)
    dimnames(chi) <- list(c("Ordered", "Disordered"), c("Immunoreactive", "Nonimmunoreactive"))
    
    # chi square test for disorder vs immunoreactivity
    chi <- chisq.test(chi, correct = TRUE)
          # Pearson's Chi-squared test with Yates' continuity correction
          # data:  chi
          # X-squared = 250.29, df = 1, p-value < 2.2e-16
          # IMPORTANT: 2.2e-16 is a bunk value from being blown out in R
          # pull real p value
          chi$p.value
                
            # <50% disorder p value = 2.240037e-56    
                
                # for >50% is disordered    
                # question: are disordered proteins enriched within immunoreactive proteins? 
                    # reactive group
                          # ordered = 2851/4546 = 62.7%
                          # disordered = 1695/4546 = 37.3%
                    # non-reactive group
                          # ordered = 30123/40843 = 73.6%
                          # disordered = 10720/40843 = 26.2%
    
                    # proportion higher of disordered proteins in immunoreactive group
                          # 37.3% / 26.2% = 1.4
                
                          # interpretation: immunoreactive proteins contain 1.4-fold higher proportion 
                          # of disordered proteins compared to non-immunoreactive proteins

# ---------------------- Numbers: repeated peptide enrichment within immunoreactive segments ---------     
# Question: are repeat segments enriched within immunoreactive segments?
# determine repeat/unique status for SEGMENTS (pv only)
# repeat segment = >50% repeat peptides

# create segments in data
# load ORDERED dataset (code to create ordered dataset found in "Peptides_immunogenicity.R")
data <- fread("/Path/peparray_ordered.txt")

    # select columns
    data <- data[ , c('protein.name', 'order', 'PROBE_SEQUENCE')]
    
    # IMPORTANT: filter for pv only
    data <- data %>% filter(grepl('PVP01', protein.name))
    
    # make df
    data$order <- as.numeric(data$order)
    
    # IMPORTANT: add in segment info
    data <- data %>%
      group_by(protein.name) %>%
      mutate(id = as.integer(order/100))
    data <- data %>%
      unite('merged', protein.name, id, remove = FALSE)
      
# determine how many times a peptide is found in dataset
count <- data %>%
  group_by(PROBE_SEQUENCE) %>%
  summarise(Count = n())

    # add repeat count to data
    data <- full_join(data, count, by = "PROBE_SEQUENCE")

    # recode count as a binary
    data$rep.pep <- ifelse(data$Count == 1, 'unique', 'repeat')  
    
    # select columns
    data <- data[, c("merged", "order", "PROBE_SEQUENCE", "rep.pep")]
    
    # make df
    data <- as.data.frame(data)
    
# add immunoreactivity data for UM patients - PVP01 only, time 1, 0 shared included
  hits <- fread('/Users/rosita.asawa/Library/CloudStorage/OneDrive-UniversityofMarylandSchoolofMedicine/Graduate School/Peptide array/R code/postchunk_all_PvProbe_vs_PvSamp_7consec_0shared_included.txt')
    
      # filter for segments with 7+ patients reactive in time 1
        hits <- hits %>% filter(total>6)

# determine repeat status of REACTIVE segments
# merge "hits" with "data" by left_join to take only reactive segments
  reac <- left_join(hits, data, by = "merged")
    
  # IMPORTANT: remove duplicate segments
  reac <- distinct(reac, merged, .keep_all= TRUE)
  
  # table it
    table(reac$rep.pep)
         
        # all reactive segment segments (315 pv proteins only)
            # repeat unique 
            #  8    307
            
# determine repeat status of NONREACTIVE segments
  # remove rows in "data" that are found in "hits"
  nonreac <- anti_join(data, hits, by = "merged")
    
  # IMPORTANT: remove duplicate segments
  nonreac <- distinct(nonreac, merged, .keep_all= TRUE)
    
  # table it
    table(nonreac$rep.pep)
    
          # all reactive segment segments (315 segments)
              # repeat unique 
              #  8    307
      
          # all NON-reactive segment segments (44907 segments)
              # repeat  unique 
              # 1342  43565

# table for chi square
    # reactive segments: 315
        # repeat segments: 8
        # unique segments:307
    # NON-reactive segments: 44907
        # repeat segments: 1342
        # unique segments:43565    
    
# chi square test  
  # create object containing values for chi square test
  chi <- matrix(c(
    8, 1342,                        # repeat: reac, nonreac
    307,  43565                     # unique: reac, nonreac
    ), nrow = 2)
  dimnames(chi) <- list(c("Repeat", "Unique"), c("Immunoreactive", "Nonimmunoreactive"))
  
  # chi square test
  # numbers are large so turn off yates correction (correct = FALSE)
  chisq.test(chi, correct = FALSE)
  
        # Pearson's Chi-squared test with Yates' continuity correction
        # data:  chi
        # X-squared = 0.21747, df = 1, p-value = 0.641
            
        # Interpretation: repeated peptides are NOT significantly enriched within
        # immunoreactive segments (p = 0.641)
              
# ---------------------- Numbers: NUP enrichment within immunoreactive proteins  ------------
# enrichment analysis
# question: among immunoreactive proteins, are NUPs enriched in symptomatic patients? 
# get following numbers
    # reactive genes
        # Nup = ?
        # NON-Nup = ?
    # NON-reactive genes
        # Nup = ?
        # NON-Nup = ?

# upload immunoreactivity list of all segments
immunoreac <- fread('/symptomatic_pv_immunogenicity_all_seg.txt')

        # has 6670 unique genes - CORRECT

# adjust files for merge
  # add protein.name back in - numbers indicate character range to grab
  # for other assessed proteins
  others <- immunoreac %>% filter(str_detect(protein.name, 'PVP01_API|PVP01_MIT'))  
  
      # select first 14 characters of Gene
      others <- others %>% mutate(protein.name = substr(protein.name, 1, 14))
      
    # for main list
    # remove weird genes from list now that they're pulled
    main <- immunoreac %>% filter(!str_detect(protein.name, "PVP01_API|PVP01_MIT"))
    
        # take first 13 characters from ID
        main <- main %>% mutate(protein.name = substr(protein.name, 1, 13))
          
        # merge others back into list
        all <- rbind(main, others)

# add gene name in
# upload pv genome
genome <- fread('/Users/rosita.asawa/Library/CloudStorage/OneDrive-UniversityofMarylandSchoolofMedicine/Graduate School/Peptide array/R code/pv_genome.txt')

    # select columns 
    genome <- genome[,c(1:2)]
    
    # rename ID column
    colnames(genome)[colnames(genome) == 'ID'] <- 'protein.name'
    
    # remove duplicates (where from?)
    genome <- distinct(genome, protein.name, .keep_all= TRUE)
    
    # add gene decription to all
    final <- left_join(all, genome, by = "protein.name")
    
# separate reactive and non-reactive
# filter for reactive
reac <- final %>% filter(total > 6)

        # 283 reactive proteins

              # select Nups in reactive
              reac.nup <- reac %>% filter(grepl('nucleoporin', Description))
              
                    # 4 "nucleoporin" in description for reactive
              
              # select NON-Nups in reactive
              reac.non.nup <- reac %>% filter(!grepl('nucleoporin', Description))
          
                    # 279 without "nucleoporin" in description for reactive
              
# filter for NON-reactive
nonreac <- final %>% filter(total < 7)

        # 6387 non-reactive proteins
    
              # select Nups in reactive
              nonreac.nup <- nonreac %>% filter(grepl('nucleoporin', Description))
              
                  # 6 "nucleoporin" in description for reactive
              
              # select NON-Nups in reactive
              nonreac.non.nup <- nonreac %>% filter(!grepl('nucleoporin', Description))    
                
                  # 6381 without 'nucleoporin' in description for NON-reactive

# set up chi square table
    # reactive genes (n = 283)
        # Nup = 4
        # Non-Nup = 279
    # NON-reactive genes (n = 6387)
        # Nup = 6
        # Non-Nup = 6381    
          
# chi square test
    # create object containing values for chi square test
    # tests: are NUPs enriched within immunoreactive proteins of symptomatic patients?
    chi <- matrix(c(
      4, 6,                           # nups: reac, nonreac
      279, 6381),                      # NON-nup: reac, nonreac
      nrow = 2, byrow = TRUE)
    dimnames(chi) <- list(c("Nup", "Non-Nup"), c("Immunoreactive", "Not immunoreactive"))
    
    # chi square test for 7+ shared segments
    # correct = TRUE gives Yates (conservative)
    # correct = FALSE gives Pearson's chisq (less conservative)
    chisq.test(chi, correct = TRUE)
          # Pearson's Chi-squared test with Yates' continuity correction
          # data:  chi
          # X-squared = 23.319, df = 1, p-value = 1.372e-06
              
              # find proportions
                # reactive genes (n=283)
                    # Nup = 4/283 = 1.4%
                    # NON-Nup = 279/283 = 98.6%
                # NON-reactive genes (n=6387)
                    # Nup = 6/6387 = 0.094%
                    # NON-Nup = 6381/6387 = 99.9%
    
                          # p = 1.37 x 10e-6
                          # fold difference = 1.4%/0.01% = 15x
               
    # Interpretation: NUPs are significantly enriched in immunoreactive proteins compared to non-NUPs (p = 1.37 x 10e-6)

# ---------------------- Numbers: asparagine repeat enrichment within immunoreactive segments ----------
# upload kmer analysis
data <- fread('/Path/kmer_analysis.txt')
summary <- fread('/Path/kmer_summary.txt')

# rename columns
names(data)[names(data) == "In negative chunks"] <- "neg_abun"
names(data)[names(data) == "In positive chunks"] <- "pos_abun"

# make numeric
data$neg_abun <- as.numeric(data$neg_abun)
data$pos_abun <- as.numeric(data$pos_abun)

# add total seg
data$total.seg <- data$neg_abun + data$pos_abun

# find percent pos out of each motif
data$perc.pos <- data$pos_abun / data$total.seg

# find total pos/neg segments
sum(data$pos_abun) 
sum(data$neg_abun) 

        # 419 total positive segments
        # 41240 total negative segments

# table for fishers - compare motifs vs NON-NAANA motifs in pos/neg segments
  
                          # NAANA       # NON-NAANA
      # counts in pos seg    22       (419 - 22) =  397                 
      # counts in neg seg    52       (41240 - 52) = 41188                                              
 
                          # ANAAN      # NON-ANAAN
      # counts in pos seg    34       (419 - 34) = 385                          
      # counts in neg seg    82       (41240 - 82) = 41158                                               

                          # AANAA     # NON-AANAA
      # counts in pos seg    21       (419 - 21) = 398
      # counts in neg seg    76       (41240 - 76) = 41164                                  

# make NAANA matrix
fish <- matrix(c(22, 397, 52, 41188), nrow = 2)
colnames(fish) <- c("Positive segments", "Negative segments")
rownames(fish) <- c("NAANA", "NON_NAANA")

      # Fisher's exact test
      # alternative = greater specifies one-sided test (test is pos assoc, OR >1)
      # "greater" tests for over-representation (enrichment)
      fisher.test(fish, alternative = "greater")
      test <- fisher.test(fish, alternative = "greater")
      test$p.value
      
      # Fisher's exact for NAANA (p = 1.48e-26, OR = 43.8)
              # p-value 1.475992e-26
              # alternative hypothesis: true odds ratio is greater than 1
              # 95 percent confidence interval:
              #  27.4191     Inf
              # sample estimates:
              # odds ratio 
              #  43.84832 

# make ANAAN matrix
fish <- matrix(c(34, 385, 82, 41158), nrow = 2)
colnames(fish) <- c("Positive segments", "Negative segments")
rownames(fish) <- c("ANAAN", "NON_ANAAN")

      # Fisher's exact test
      # alternative = greater specifies one-sided test (test is pos assoc, OR >1)
      # "greater" tests for over-representation (enrichment)
      fisher.test(fish, alternative = "greater")
      test <- fisher.test(fish, alternative = "greater")
      test$p.value
      
      # Fisher's exact for ANAAN (p = 3.59e-40, OR = 44.3)
            # p-value = 3.590702e-40
            # alternative hypothesis: true odds ratio is greater than 1
            # 95 percent confidence interval:
            # 30.47307      Inf
            # sample estimates:
            # odds ratio 
            #  44.29232
  
# make AANAA matrix
fish <- matrix(c(21, 398, 76, 41164), nrow = 2)
colnames(fish) <- c("Positive segments", "Negative segments")
rownames(fish) <- c("AANAA", "NON_AANAA")

      # Fisher's exact test
      # alternative = greater specifies one-sided test (test is pos assoc, OR >1)
      # "greater" tests for over-representation (enrichment)
      fisher.test(fish, alternative = "greater")
      test <- fisher.test(fish, alternative = "greater")
      test$p.value
      
      # Fisher's exact for AANAA (p = 3.39e-22, OR = 28.6)  
              # p-value = 3.390714e-22
              # alternative hypothesis: true odds ratio is greater than 1
              # 95 percent confidence interval:
              #  18.07139      Inf
              # sample estimates:
              # odds ratio 
              #   28.55134 

# ---------------------- Fig S5: correlation plot of kmers by antigenicity ------------------
# Kmer analysis ran on perl to determine enrichment of repetitive motifs
# upload data
data <- fread('/Path/kmer_analysis.txt')
summary <- fread('/Path/kmer_summary.txt')

# rename columns
names(data)[names(data) == "In negative chunks"] <- "neg_abun"
names(data)[names(data) == "In positive chunks"] <- "pos_abun"

# make numeric
data$neg_abun <- as.numeric(data$neg_abun)
data$pos_abun <- as.numeric(data$pos_abun)

# plot with ggtext
ggplot(data, aes(x=neg_abun, y=pos_abun)) +
  geom_point(color = "black") + 
  theme_classic() +
  labs(title = "Kmer analysis of reptitive motifs",
     x = "Abundance in non-antigenic segments", 
     y = "Abundance in antigenic segments") +
  theme(plot.caption = element_markdown(lineheight = 1.2),
     plot.title = element_markdown(size=15, face="bold", hjust=0.5, color="#666666"),
     axis.title.x = element_markdown(face="bold", size = 12),
     axis.title.y = element_markdown(face="bold", size = 12),
     axis.text.x = element_markdown(size = 11),
     axis.text.y = element_markdown(size = 11)) +
  scale_x_continuous(limits=c(0, 2500), oob = scales::oob_keep, breaks = scales::pretty_breaks(n = 5)) +
  scale_y_continuous(limits=c(0, 40), oob = scales::oob_keep)

# ---------------- Asymptomatic Pv patient enrichment ---------
# ---------------------- Numbers: PIR enrichment within immunoreactive proteins ---------------
# enrichment analysis
# question: among immunoreactive proteins, are PIRs enriched in asymptomatic patients? 
# get following numbers
    # reactive genes in symptomatic
        # PIRs = ?
        # NON-PIRS = ?
    
    # reactive genes in asymptomatic
        # PIRs = ?
        # NON-PIRS = ?
    
# SYMPTOMATIC numbers
# upload immunoreactivity list for SYMPTOMATIC patients - most reactive segment per protein
immunoreac.symp <- fread('/symptomatic_pv_most_reactive_seg_per_protein.txt.txt')

        # has 6670 unique genes - CORRECT

# adjust files for merge
  # add protein.name back in - numbers indicate character range to grab
    # for weird genes
    others <- immunoreac.symp %>% filter(str_detect(protein.name, 'PVP01_API|PVP01_MIT'))  
    
        # select first 14 characters of Gene
        others <- others %>% mutate(protein.name = substr(protein.name, 1, 14))
        
    # for main list
    # remove weird genes from list now that they're pulled
    main <- immunoreac.symp %>% filter(!str_detect(protein.name, "PVP01_API|PVP01_MIT"))
    
        # take first 13 characters from ID
        main <- main %>% mutate(protein.name = substr(protein.name, 1, 13))
          
        # merge others back into list
        all <- rbind(main, others)

# add gene name in
# upload pv genome
genome <- fread('/Path/pv_genome.txt')

    # select columns 
    genome <- genome[,c(1:2)]
    
    # rename ID column
    colnames(genome)[colnames(genome) == 'ID'] <- 'protein.name'
    
    # remove duplicates (where from?)
    genome <- distinct(genome, protein.name, .keep_all= TRUE)
    
    # add gene decription to all
    immunoreac.symp <- left_join(all, genome, by = "protein.name")
    
# IMPROTANT: select immunoreactive genes only
immunoreac.symp <- immunoreac.symp %>% filter(total > 6)

        # 282 reactive proteins in UM
    
# separate PIRs and non-PIRs into different dfs      
# filter for PIR genes
pir <- immunoreac.symp %>% filter(grepl('PIR', Description))
    
        # 20 PIR proteins in immunoreactive UM proteins

    # filter for NON-PIR genes
    non.pir <- immunoreac.symp %>% filter(!grepl('PIR', Description))
    
        # 263 NON-PIR proteins in immunoreactive UM proteins

# ASYMPTOMATIC numbers
# upload immunoreactivity list for SYMPTOMATIC patients - most reactive segment per protein
immunoreac.asymp <- fread('/asymptomatic_pv_most_reactive_seg_per_protein.txt.txt')

        # has 6670 unique genes - CORRECT

# add gene name in
# upload pv genome
genome <- fread('/Path/pv_genome.txt')

    # select columns 
    genome <- genome[,c(1:2)]
    
    # rename ID column
    colnames(genome)[colnames(genome) == 'ID'] <- 'protein.name'
    
    # remove duplicates (where from?)
    genome <- distinct(genome, protein.name, .keep_all= TRUE)
    
    # add gene description to all
    immunoreac.asymp <- left_join(immunoreac.asymp, genome, by = "protein.name")

# IMPORTANT: select for reactive proteins only
immunoreac.asymp <- immunoreac.asymp %>% filter(total == 5)    

        # 414 proteins reactive in AM
    
# separate PIRs and non-PIRs into different dfs      
# filter for PIR genes
pir <- immunoreac.asymp %>% filter(grepl('PIR', Description))
    
        # 65 PIR proteins are reactive in AM

    # filter for NON-PIR genes
    non.pir <- immunoreac.asymp %>% filter(!grepl('PIR', Description))
    
        # 349 NON-PIR proteins in immunoreactive AM proteins

# table for chi square 
    # reactive genes in symptomatic (n = 283)
        # PIRs = 20
        # NON-PIRS = 263
    # reactive genes in asymptomatic (n = 414)
        # PIRs = 65
        # NON-PIRS = 349
    
# chi square test
    # create object containing values for chi square test
    # tests: are asymptomatic or symptomatic individuals more likely to have PIRs that are reactive?
    chi <- matrix(c(
      20, 65,                       # reactive pirs: symp, asymp
      263, 349                      # reactive non-pirs: symp, asymp
      ), nrow = 2, byrow = TRUE)
    dimnames(chi) <- list(c("Immunoreactive PIRs", "Immunoreactive non-PIRs"), c("Symptomatic", "Asymptomatic"))
    
    # chi square test for 7+ shared segments
    chisq.test(chi, correct = TRUE)
          # Pearson's Chi-squared test with Yates' continuity correction
          # data:  chi
          # X-squared = 10.908, df = 1, p-value = 0.0009574
          
              # find proportions
                  # reactive genes in symptomatic
                      # PIRs = 20/283 = 7.1%
                      # NON-PIRS = 263/283 =9 2.9%
                  # reactive genes in asymptomatic
                      # PIRs = 65/414 = 15.7%
                      # NON-PIRS = 349/414 = 84.3%

              # fold enrichment of immunoreactive PIRs in asympto
                  # 15.7%/7.1% = 2.2 fold difference
    
                            # Interpretation: immunoreactive PIRs are significantly more immunoreactive in 
                            # asymptomatic patients compared to symptomatic (p =0.001)
  