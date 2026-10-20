# Script 1: Immunogenicity analysis
# coordinating custom scripts for Plos Pathogens manuscript, published October 2026
  # Proteome-wide characterization of Plasmodium vivax antigens using a high-density peptide array
  # https://journals.plos.org/plospathogens/article?id=10.1371/journal.ppat.1014621

library(tidyverse)
library(data.table)
library(plyr)
library(ggplot2)
library(ggtext)
library(ggpubr)

install.packages("")

# ----------- Order peptide array dataset by amino acid position ------------------
# load PROCESSED data faster using data.table package
# appears to be processed dataset given by Nimble with proprietary normalization
# probably normalized across patient arrays (one array/patient) to account for plate effect 
data <- fread("/Path/Processed_aggregate_data.txt")

# change headers for a PROCESSED data (normalized is changed already)
setnames(data, c("PROBE_DESIGN_ID",	"CONTAINER",	"DESIGN_NOTE",	"SELECTION_CRITERIA",	"SEQ_ID",	"PROBE_SEQUENCE",	"MISMATCH",	"MATCH_INDEX",	"FEATURE_ID",	"ROW_NUM",	"COL_NUM",	"PROBE_CLASS",	"PROBE_ID", "POSITION",	"DESIGN_ID",	"X",	"Y",	"CYCLES",	"DMD",	
                 "Pv_UMc_5.2",	"Pv_UMc_5.1",	"Pv_UMr_5.1",	"Pv_UMc_2.2",	"Pv_UMr_1.1",	"Pv_UMr_3.1",	"Pv_UMc_4.1",	"Pv_AM_1",	"Pv_AM_2",	"Pv_AM_3",	"Pv_UMr_3.2",	"Pv_AM_4",	"Pv_UMr_4.2",	"Pv_UMc_3.1",	"Pv_UMr_4.1",	"Pv_UMc_3.2",	"Pv_UMr_5.2",	"Pv_UMc_1.2",	"Pv_UMr_2.1",	"Pv_UMr_1.2",	"NA_ctrl_1",	"NA_ctrl_2",	"NA_ctrl_3",	"NA_ctrl_4",	"NA_ctrl_5",	"NA_ctrl_6",	"Pv_AM_5",	"Pv_UMr_2.2",	"NA_ctrl_7",	"Pv_UMc_4.2",	"Pv_monkey_t2",	"Pv_monkey_t1",	"NA_ctrl_8",	"NA_ctrl_9",	"NA_ctrl_10",	"Pv_UMc_1.1",	"Pv_UMc_2.1",	"Pf_pool",	"Pv_monkey_t3"))

    # Nomenclature key
      # NA_ctrl = north american malaria-naive control
      # Pv_UMc = uncomplicated malaria, cleared of parasite at time (symptomatic pv patients)
      # Pv_UMr = uncomplicated malaria, relapsed with parasite at time (symptomatic pv patients)
      # Pv_AM = asymptomatic malaria (asymptomatic pv patients)

# upload probe duplication file
dup <- fread("/Users/rosita.asawa/Library/CloudStorage/OneDrive-UniversityofMarylandSchoolofMedicine/Graduate School/Peptide array/Peparray raw data/pvivax_peptide_array_duplicate_info.txt")

# merge dup and data by sequence: "all=TRUE" contains extra peptides with no data (remove later)
combo <- merge(data, dup, by = "PROBE_SEQUENCE", all = TRUE)
       
# remove blanks from combo (?)
combo <- combo[-which(is.na(combo$PROBE_DESIGN_ID))]

# deduplicate combo into new row for each protein
dedup <- combo %>% separate_rows("PROTEIN", sep = ",")

# create protein.name column
dedup <- dedup %>% separate(PROTEIN, into = c("protein.name", "order"),sep = "_(?!.*_)")

# remove blanks from dedup due to missing keys
dedup <- dedup[-which(is.na(dedup$protein.name)), ]

# convert order character to number
dedup$order <- as.numeric(dedup$order)

# order by protein name and by amino acid position (order) 
ordered <- dedup %>%
  group_by(protein.name) %>%
  arrange(order, .by_group = TRUE, numeric = TRUE)

# check for NAs
sum(is.na(ordered$protein.name))
sum(is.na(ordered$Pv_UMc_1.1))

# export to excel
write.table(ordered,"/Path/peparray_ordered.txt", sep = "\t", row.names = TRUE, col.names = NA)

# ------------------- Create list of repeat peptides across PVP01 proteins -----------------
# load ORDERED dataset
data <- fread("/Path/peparray_ordered.txt")

# select UM initial time point and NA only for ordered dataset
data <- data[ , c('PROBE_SEQUENCE', 'protein.name', 'order')]

# IMPORTANT: remove duplicates of pf genes ran twice (16mer/15 and 16mer/13)
# keep 16mer/15 (PF3D7_1133400.1) and remove 16mer/13 (PF3D7_1133400.1-p1)
data <- data %>% filter(!grepl('^PF3D7_1133400.1-p1$', protein.name))        # PfAMA1
data <- data %>% filter(!grepl('^PF3D7_0930300.1-p1$', protein.name))        # PfMSP1
data <- data %>% filter(!grepl('^PF3D7_1035400.1-p1$', protein.name))        # PfMSP3
data <- data %>% filter(!grepl('^PF3D7_0620400.1-p1$', protein.name))        # PfAMSP10
data <- data %>% filter(!grepl('^PF3D7_1364100.1-p1$', protein.name))        # Pf 6-cysteine protein P92
data <- data %>% filter(!grepl('^PF3D7_0612700.1-p1$', protein.name))        # Pf 6-cysteine protein P12
data <- data %>% filter(!grepl('^PF3D7_0104200.1-p1$', protein.name))        # Pf StAR-related lipid transfer protein

# count number of times a peptide is found in dataset
count <- data %>%
  select(PROBE_SEQUENCE, protein.name, order) %>%
  group_by(PROBE_SEQUENCE) %>%
  summarise(pep.count = n())

# recode count as a binary
count$peptide <- ifelse(count$pep.count == 1, 'unique', 'repeat')

    # select columns
    count <- count[c('PROBE_SEQUENCE', 'peptide')]

    # sort dfs to see if it worked
    count <- count %>% arrange(PROBE_SEQUENCE)
    data.sort <- data %>% arrange(PROBE_SEQUENCE)

# merge 'data' with 'count' 
all <- merge(data, count, by.x = "PROBE_SEQUENCE", by.y = "PROBE_SEQUENCE", 
             all.x = TRUE, all.y = FALSE)
 
# check single gene (PfAMA1 16mer/15)
gene <- all %>% filter(grepl('PF3D7_1133400.1', protein.name))
   
    # yes it worked!

# select columns
all <- all[ , c('protein.name', 'order', 'PROBE_SEQUENCE', 'peptide')]

# export list of repeat peptides with useful information
write.table(all,"/Path/repeat_peptide_list.txt", sep = "\t", row.names = TRUE, col.names = NA)

# ----------- Immunogenicity set up - Symptomatic Pv patients, Pv probes, time 1 -------------------------
# upload data
data <- fread("/Path/peparray_ordered.txt")

# filter for Pv only
data <- data %>% filter(grepl('PVP01', protein.name))

# create 500 MFI column as NA_11 as extra cutoff
data$NA_ctrl_11 <- 500

# select UM initial time point and NA only for ordered dataset
data <- data[ , c('protein.name', 'order',
                  'NA_ctrl_1','NA_ctrl_2','NA_ctrl_3','NA_ctrl_4','NA_ctrl_5','NA_ctrl_6','NA_ctrl_7','NA_ctrl_8','NA_ctrl_9','NA_ctrl_10', 'NA_ctrl_11',
                  'Pv_UMc_1.1',	'Pv_UMc_2.1',	'Pv_UMc_3.1',	'Pv_UMc_4.1',	'Pv_UMc_5.1',	
                  'Pv_UMr_1.1',	'Pv_UMr_2.1',	'Pv_UMr_3.1',	'Pv_UMr_4.1',	'Pv_UMr_5.1')]

# create max control and max UM and 500 MFI column
data$NA_max <- pmax(data$NA_ctrl_1, data$NA_ctrl_2, data$NA_ctrl_3, data$NA_ctrl_4, data$NA_ctrl_5, data$NA_ctrl_6, data$NA_ctrl_7, data$NA_ctrl_8, data$NA_ctrl_9, data$NA_ctrl_10, data$NA_ctrl_11)

# create column to indicate if sample is greater than max NA and less than 100 mfi in NA 
# ordered dataset
data$UMc_1.1 <- with(data, ifelse(data$Pv_UMc_1.1 > data$NA_max, 1, 0))
data$UMc_2.1 <- with(data, ifelse(data$Pv_UMc_2.1 > data$NA_max, 1, 0))
data$UMc_3.1 <- with(data, ifelse(data$Pv_UMc_3.1 > data$NA_max, 1, 0))
data$UMc_4.1 <- with(data, ifelse(data$Pv_UMc_4.1 > data$NA_max, 1, 0))
data$UMc_5.1 <- with(data, ifelse(data$Pv_UMc_5.1 > data$NA_max, 1, 0))
data$UMr_1.1 <- with(data, ifelse(data$Pv_UMr_1.1 > data$NA_max, 1, 0))
data$UMr_2.1 <- with(data, ifelse(data$Pv_UMr_2.1 > data$NA_max, 1, 0))
data$UMr_3.1 <- with(data, ifelse(data$Pv_UMr_3.1 > data$NA_max, 1, 0))
data$UMr_4.1 <- with(data, ifelse(data$Pv_UMr_4.1 > data$NA_max, 1, 0))
data$UMr_5.1 <- with(data, ifelse(data$Pv_UMr_5.1 > data$NA_max, 1, 0))

# select new columns only
data <- data[ , c('protein.name', 'order',
                  'UMc_1.1', 'UMc_2.1', 'UMc_3.1', 'UMc_4.1', 'UMc_5.1',
                  'UMr_1.1', 'UMr_2.1', 'UMr_3.1', 'UMr_4.1', 'UMr_5.1')]

# make numeric
data <- as.data.frame(data)

# determine consecutively reactive epitopes longer than >6 peptides
UMc_1.1 <- data %>%
  select(UMc_1.1, order, protein.name) %>%
  mutate(grp = data.table::rleid(UMc_1.1)) %>%
  filter(UMc_1.1 ==1) %>%
  group_by(protein.name, grp) %>%
  mutate(consec = row_number())
UMc_1.1 <- UMc_1.1 %>%
  select(UMc_1.1, order, protein.name, grp, consec) %>%
  filter(consec > 6) %>%
  summarise_all(last)           

UMc_2.1 <- data %>%
  select(UMc_2.1, order, protein.name) %>%
  mutate(grp = data.table::rleid(UMc_2.1)) %>%
  filter(UMc_2.1 ==1) %>%
  group_by(protein.name, grp) %>%
  mutate(consec = row_number())
UMc_2.1 <- UMc_2.1 %>%
  select(UMc_2.1, order, protein.name, grp, consec) %>%
  group_by(grp) %>% 
  filter(consec > 6) %>%
  summarise_all(last)  

UMc_3.1 <- data %>%
  select(UMc_3.1, order, protein.name) %>%
  mutate(grp = data.table::rleid(UMc_3.1)) %>%
  filter(UMc_3.1 ==1) %>%
  group_by(protein.name, grp) %>%
  mutate(consec = row_number())
UMc_3.1 <- UMc_3.1 %>%
  select(UMc_3.1, order, protein.name, grp, consec) %>%
  group_by(grp) %>% 
  filter(consec > 6) %>%
  summarise_all(last)  

UMc_4.1 <- data %>%
  select(UMc_4.1, order, protein.name) %>%
  mutate(grp = data.table::rleid(UMc_4.1)) %>%
  filter(UMc_4.1 ==1) %>%
  group_by(protein.name, grp) %>%
  mutate(consec = row_number())
UMc_4.1 <- UMc_4.1 %>%
  select(UMc_4.1, order, protein.name, grp, consec) %>%
  group_by(grp) %>% 
  filter(consec > 6) %>%
  summarise_all(last)  

UMc_5.1 <- data %>%
  select(UMc_5.1, order, protein.name) %>%
  mutate(grp = data.table::rleid(UMc_5.1)) %>%
  filter(UMc_5.1 ==1) %>%
  group_by(protein.name, grp) %>%
  mutate(consec = row_number())
UMc_5.1 <- UMc_5.1 %>%
  select(UMc_5.1, order, protein.name, grp, consec) %>%
  group_by(grp) %>% 
  filter(consec > 6) %>%
  summarise_all(last)  

UMr_1.1 <- data %>%
  select(UMr_1.1, order, protein.name) %>%
  mutate(grp = data.table::rleid(UMr_1.1)) %>%
  filter(UMr_1.1 ==1) %>%
  group_by(protein.name, grp) %>%
  mutate(consec = row_number())
UMr_1.1 <- UMr_1.1 %>%
  select(UMr_1.1, order, protein.name, grp, consec) %>%
  group_by(grp) %>% 
  filter(consec > 6) %>%
  summarise_all(last)  

UMr_2.1 <- data %>%
  select(UMr_2.1, order, protein.name) %>%
  mutate(grp = data.table::rleid(UMr_2.1)) %>%
  filter(UMr_2.1 ==1) %>%
  group_by(protein.name, grp) %>%
  mutate(consec = row_number())
UMr_2.1 <- UMr_2.1 %>%
  select(UMr_2.1, order, protein.name, grp, consec) %>%
  group_by(grp) %>% 
  filter(consec > 6) %>%
  summarise_all(last)

UMr_3.1 <- data %>%
  select(UMr_3.1, order, protein.name) %>%
  mutate(grp = data.table::rleid(UMr_3.1)) %>%
  filter(UMr_3.1 ==1) %>%
  group_by(protein.name, grp) %>%
  mutate(consec = row_number())
UMr_3.1 <- UMr_3.1 %>%
  select(UMr_3.1, order, protein.name, grp, consec) %>%
  group_by(grp) %>% 
  filter(consec > 6) %>%
  summarise_all(last) 

UMr_4.1 <- data %>%
  select(UMr_4.1, order, protein.name) %>%
  mutate(grp = data.table::rleid(UMr_4.1)) %>%
  filter(UMr_4.1 ==1) %>%
  group_by(protein.name, grp) %>%
  mutate(consec = row_number())
UMr_4.1 <- UMr_4.1 %>%
  select(UMr_4.1, order, protein.name, grp, consec) %>%
  group_by(grp) %>% 
  filter(consec > 6) %>%
  summarise_all(last)  

UMr_5.1 <- data %>%
  select(UMr_5.1, order, protein.name) %>%
  mutate(grp = data.table::rleid(UMr_5.1)) %>%
  filter(UMr_5.1 ==1) %>%
  group_by(protein.name, grp) %>%
  mutate(consec = row_number())
UMr_5.1 <- UMr_5.1 %>%
  select(UMr_5.1, order, protein.name, grp, consec) %>%
  group_by(grp) %>% 
  filter(consec > 6) %>%
  summarise_all(last)

# chunk proteins into 100aa segments and merge "protein.name" and "id"
UMc_1.1 <- UMc_1.1 %>%
  group_by(protein.name) %>%
  mutate(id = as.integer(order/100))
UMc_1.1 <- UMc_1.1 %>%
  unite('merged', protein.name, id, remove = FALSE)

UMc_2.1 <- UMc_2.1 %>%
  group_by(protein.name) %>%
  mutate(id = as.integer(order/100))
UMc_2.1 <- UMc_2.1 %>%
  unite('merged', protein.name, id, remove = FALSE)

UMc_3.1 <- UMc_3.1 %>%
  group_by(protein.name) %>%
  mutate(id = as.integer(order/100))
UMc_3.1 <- UMc_3.1 %>%
  unite('merged', protein.name, id, remove = FALSE)

UMc_4.1 <- UMc_4.1 %>%
  group_by(protein.name) %>%
  mutate(id = as.integer(order/100))
UMc_4.1 <- UMc_4.1 %>%
  unite('merged', protein.name, id, remove = FALSE)

UMc_5.1 <- UMc_5.1 %>%
  group_by(protein.name) %>%
  mutate(id = as.integer(order/100))
UMc_5.1 <- UMc_5.1 %>%
  unite('merged', protein.name, id, remove = FALSE)

UMr_1.1 <- UMr_1.1 %>%
  group_by(protein.name) %>%
  mutate(id = as.integer(order/100))
UMr_1.1 <- UMr_1.1 %>%
  unite('merged', protein.name, id, remove = FALSE)

UMr_2.1 <- UMr_2.1 %>%
  group_by(protein.name) %>%
  mutate(id = as.integer(order/100))
UMr_2.1 <- UMr_2.1 %>%
  unite('merged', protein.name, id, remove = FALSE)

UMr_3.1 <- UMr_3.1 %>%
  group_by(protein.name) %>%
  mutate(id = as.integer(order/100))
UMr_3.1 <- UMr_3.1 %>%
  unite('merged', protein.name, id, remove = FALSE)

UMr_4.1 <- UMr_4.1 %>%
  group_by(protein.name) %>%
  mutate(id = as.integer(order/100))
UMr_4.1 <- UMr_4.1 %>%
  unite('merged', protein.name, id, remove = FALSE)

UMr_5.1 <- UMr_5.1 %>%
  group_by(protein.name) %>%
  mutate(id = as.integer(order/100))
UMr_5.1 <- UMr_5.1 %>%
  unite('merged', protein.name, id, remove = FALSE)

# ------------------- Table S3: Immunogenicity of symptomatic Pv patients across P01 genome ------------------------------
# create file reporting immunogenicity of all segments
# keep relevant columns for SEGMENTS and PROTEINS
xUMc_1 <- UMc_1.1[ , c('merged', 'protein.name', 'UMc_1.1')]
xUMc_2 <- UMc_2.1[ , c('merged', 'protein.name', 'UMc_2.1')]
xUMc_3 <- UMc_3.1[ , c('merged', 'protein.name', 'UMc_3.1')]
xUMc_4 <- UMc_4.1[ , c('merged', 'protein.name', 'UMc_4.1')]
xUMc_5 <- UMc_5.1[ , c('merged', 'protein.name', 'UMc_5.1')]
xUMr_1 <- UMr_1.1[ , c('merged', 'protein.name', 'UMr_1.1')]
xUMr_2 <- UMr_2.1[ , c('merged', 'protein.name', 'UMr_2.1')]
xUMr_3 <- UMr_3.1[ , c('merged', 'protein.name', 'UMr_3.1')]
xUMr_4 <- UMr_4.1[ , c('merged', 'protein.name', 'UMr_4.1')]
xUMr_5 <- UMr_5.1[ , c('merged', 'protein.name', 'UMr_5.1')]

    # remove duplicates from segments
    xUMc_1 <- distinct(xUMc_1, merged, .keep_all= TRUE)
    xUMc_2 <- distinct(xUMc_2, merged, .keep_all= TRUE)
    xUMc_3 <- distinct(xUMc_3, merged, .keep_all= TRUE)
    xUMc_4 <- distinct(xUMc_4, merged, .keep_all= TRUE)
    xUMc_5 <- distinct(xUMc_5, merged, .keep_all= TRUE)
    xUMr_1 <- distinct(xUMr_1, merged, .keep_all= TRUE)
    xUMr_2 <- distinct(xUMr_2, merged, .keep_all= TRUE)
    xUMr_3 <- distinct(xUMr_3, merged, .keep_all= TRUE)
    xUMr_4 <- distinct(xUMr_4, merged, .keep_all= TRUE)
    xUMr_5 <- distinct(xUMr_5, merged, .keep_all= TRUE)
    
    # keep relevant columns for SEGMENTS
    xUMc_1 <- xUMc_1[ , c('merged', 'UMc_1.1')]
    xUMc_2 <- xUMc_2[ , c('merged', 'UMc_2.1')]
    xUMc_3 <- xUMc_3[ , c('merged', 'UMc_3.1')]
    xUMc_4 <- xUMc_4[ , c('merged', 'UMc_4.1')]
    xUMc_5 <- xUMc_5[ , c('merged', 'UMc_5.1')]
    xUMr_1 <- xUMr_1[ , c('merged', 'UMr_1.1')]
    xUMr_2 <- xUMr_2[ , c('merged', 'UMr_2.1')]
    xUMr_3 <- xUMr_3[ , c('merged', 'UMr_3.1')]
    xUMr_4 <- xUMr_4[ , c('merged', 'UMr_4.1')]
    xUMr_5 <- xUMr_5[ , c('merged', 'UMr_5.1')]
    
    # join dfs and convert NAs to 0
    all <- full_join(xUMc_1, xUMc_2, by = "merged")
    all <- full_join(all, xUMc_3, by = "merged")
    all <- full_join(all, xUMc_4, by = "merged")
    all <- full_join(all, xUMc_5, by = "merged")
    all <- full_join(all, xUMr_1, by = "merged")
    all <- full_join(all, xUMr_2, by = "merged")
    all <- full_join(all, xUMr_3, by = "merged")
    all <- full_join(all, xUMr_4, by = "merged")
    all <- full_join(all, xUMr_5, by = "merged")
    all[is.na(all)] <- 0
    
    # filter for Pv only
    all <- all %>% filter(grepl('PVP01', merged))
    
    # add total column
    all$total <- rowSums(all[, c("UMc_1.1", "UMc_2.1", "UMc_3.1", "UMc_4.1", "UMc_5.1", 
                                 "UMr_1.1", "UMr_2.1", "UMr_3.1", "UMr_4.1", "UMr_5.1")])
    
    # check to make sure the numbers are the same as previous code
    table(all$total)
          
          #    1     2     3     4     5     6     7     8     9    10 
          # 12941  9401  5584  2794  1285   527   221    75    16     3  
    
    # select columns
    all <- all[, c("merged", "total")]

# create file of all segments
# load ORDERED dataset
data <- fread("/Path/peparray_ordered.txt")

    # select UM initial time point and NA only for ordered dataset
    data <- data[ , c('protein.name', 'order')]
    
    # IMPORTANT: filter for pv only
    data <- data %>% filter(grepl('PVP01', protein.name))
    
    # IMPORTANT: add in segment info
    data <- data %>%
      group_by(protein.name) %>%
      mutate(id = as.integer(order/100)) %>%
      unite('merged', protein.name, id, remove = FALSE)
    
    # remove duplicates (one segment each)
    data <- distinct(data, merged, .keep_all= TRUE)
    
    # select columns
    data <- data[, c("protein.name", "merged")]

# add reactivity data to segment data
both <- left_join(data, all, by = "merged")
    
    # turn NAs to 0
    both[is.na(both)] <- 0
    
    # check gene number
    test <- distinct(both, merged, .keep_all= TRUE)
    test <- distinct(both, protein.name, .keep_all= TRUE)
          
          # total of 45222 segments
          # total of 6670 genes

# export all segments 
write.table(both, '/Path/symptomatic_pv_immunogenicity_all_seg.txt', sep = "\t", row.names = FALSE, col.names = TRUE)
        
# ------------------- Numbers: Immunoreactive peptides per sample ------------
# upload "ordered" dataset in which proteins are ordered by amino acid position
data <- fread("/Path/peparray_ordered.txt")

# filter for Pv only
data <- data %>% filter(grepl('PVP01', protein.name))

# create 500 MFI column as NA_11 as extra cutoff
data$NA_ctrl_11 <- 500

# select UM initial time point and NA only for ordered dataset
data <- data[ , c('protein.name', 'order',
                  'NA_ctrl_1','NA_ctrl_2','NA_ctrl_3','NA_ctrl_4','NA_ctrl_5','NA_ctrl_6','NA_ctrl_7','NA_ctrl_8','NA_ctrl_9','NA_ctrl_10', 'NA_ctrl_11',
                  'Pv_UMc_1.1',	'Pv_UMc_2.1',	'Pv_UMc_3.1',	'Pv_UMc_4.1',	'Pv_UMc_5.1',	
                  'Pv_UMr_1.1',	'Pv_UMr_2.1',	'Pv_UMr_3.1',	'Pv_UMr_4.1',	'Pv_UMr_5.1'), with = FALSE]

# create max control column
data$NA.max <- pmax(data$NA_ctrl_1, data$NA_ctrl_2, data$NA_ctrl_3, data$NA_ctrl_4, data$NA_ctrl_5, data$NA_ctrl_6, data$NA_ctrl_7, data$NA_ctrl_8, data$NA_ctrl_9, data$NA_ctrl_10, data$NA_ctrl_11)

# create column to indicate if sample is greater than max NA and less than 500 mfi in NA 
# ordered dataset
data$UMc_1.1 <- with(data, ifelse(data$Pv_UMc_1.1 > data$NA.max, 1, 0))
data$UMc_2.1 <- with(data, ifelse(data$Pv_UMc_2.1 > data$NA.max, 1, 0))
data$UMc_3.1 <- with(data, ifelse(data$Pv_UMc_3.1 > data$NA.max, 1, 0))
data$UMc_4.1 <- with(data, ifelse(data$Pv_UMc_4.1 > data$NA.max, 1, 0))
data$UMc_5.1 <- with(data, ifelse(data$Pv_UMc_5.1 > data$NA.max, 1, 0))
data$UMr_1.1 <- with(data, ifelse(data$Pv_UMr_1.1 > data$NA.max, 1, 0))
data$UMr_2.1 <- with(data, ifelse(data$Pv_UMr_2.1 > data$NA.max, 1, 0))
data$UMr_3.1 <- with(data, ifelse(data$Pv_UMr_3.1 > data$NA.max, 1, 0))
data$UMr_4.1 <- with(data, ifelse(data$Pv_UMr_4.1 > data$NA.max, 1, 0))
data$UMr_5.1 <- with(data, ifelse(data$Pv_UMr_5.1 > data$NA.max, 1, 0))

# select new columns only
data <- data[ , c('protein.name', 'order',
                  'UMc_1.1',	'UMc_2.1',	'UMc_3.1',	'UMc_4.1',	'UMc_5.1',	
                  'UMr_1.1',	'UMr_2.1',	'UMr_3.1',	'UMr_4.1',	'UMr_5.1'), with = FALSE]

# determine immunogenicity - how many immunoreactive peptides are there per sample
table(data$UMc_1.1)
    # 0       1 
    # 4050209  118833 
table(data$UMc_2.1)
    # 0       1 
    # 3961393  207649 
table(data$UMc_3.1)
    # 0       1 
    # 3927529  241513 
table(data$UMc_4.1)
    # 0       1 
    # 4052153  116889 
table(data$UMc_5.1)
    # 0       1 
    # 3988775  180267 
table(data$UMr_1.1)
    # 0       1 
    # 4072902   96140 
table(data$UMr_2.1)
    # 0       1 
    # 4069712   99330 
table(data$UMr_3.1)
    # 0       1 
    # 3772387  396655 
table(data$UMr_4.1)
    # 0       1 
    # 3982001  187041
table(data$UMr_5.1)
    # 0       1 
    # 3763078  405964 

# average immunoreactive Pv peptides per sample
# pivot longer
data.long <- data %>% pivot_longer(cols = -c(protein.name, order), names_to = 'Samples', values_to = 'Reac.pep')

# find average reactive peptides
ave <- data.long %>%
          select(Samples, Reac.pep) %>%
          group_by(Samples) %>%
          summarise(mean(sum(Reac.pep)))

# total immunoreactive Pv peptides = 2,050,281 peptides
# range of immunoreactive Pv peptides = 96,140 - 405,964 peptides

# ------------------- Numbers: Immunoreactive regions per sample ------------
# start with immunogenicity set up - UM, PVP01 filtered
# keep relevant columns for SEGMENTS
xUMc_1 <- UMc_1.1[ , c('merged', 'UMc_1.1', 'order', 'consec')]
xUMc_2 <- UMc_2.1[ , c('merged', 'UMc_2.1', 'order', 'consec')]
xUMc_3 <- UMc_3.1[ , c('merged', 'UMc_3.1', 'order', 'consec')]
xUMc_4 <- UMc_4.1[ , c('merged', 'UMc_4.1', 'order', 'consec')]
xUMc_5 <- UMc_5.1[ , c('merged', 'UMc_5.1', 'order', 'consec')]
xUMr_1 <- UMr_1.1[ , c('merged', 'UMr_1.1', 'order', 'consec')]
xUMr_2 <- UMr_2.1[ , c('merged', 'UMr_2.1', 'order', 'consec')]
xUMr_3 <- UMr_3.1[ , c('merged', 'UMr_3.1', 'order', 'consec')]
xUMr_4 <- UMr_4.1[ , c('merged', 'UMr_4.1', 'order', 'consec')]
xUMr_5 <- UMr_5.1[ , c('merged', 'UMr_5.1', 'order', 'consec')]

# combine order and consec into one column
xUMc_1$id.order.consec <- do.call(paste, c(xUMc_1[c("merged", "order", "consec")], sep = "-"))
xUMc_2$id.order.consec <- do.call(paste, c(xUMc_2[c("merged", "order", "consec")], sep = "-"))
xUMc_3$id.order.consec <- do.call(paste, c(xUMc_3[c("merged", "order", "consec")], sep = "-"))
xUMc_4$id.order.consec <- do.call(paste, c(xUMc_4[c("merged", "order", "consec")], sep = "-"))
xUMc_5$id.order.consec <- do.call(paste, c(xUMc_5[c("merged", "order", "consec")], sep = "-"))
xUMr_1$id.order.consec <- do.call(paste, c(xUMr_1[c("merged", "order", "consec")], sep = "-"))
xUMr_2$id.order.consec <- do.call(paste, c(xUMr_2[c("merged", "order", "consec")], sep = "-"))
xUMr_3$id.order.consec <- do.call(paste, c(xUMr_3[c("merged", "order", "consec")], sep = "-"))
xUMr_4$id.order.consec <- do.call(paste, c(xUMr_4[c("merged", "order", "consec")], sep = "-"))
xUMr_5$id.order.consec <- do.call(paste, c(xUMr_5[c("merged", "order", "consec")], sep = "-"))

# keep relevant columns for immunoreactive REGIONS (id.order.consec)
xUMc_1 <- xUMc_1[ , c('UMc_1.1', 'id.order.consec')]
xUMc_2 <- xUMc_2[ , c('UMc_2.1', 'id.order.consec')]
xUMc_3 <- xUMc_3[ , c('UMc_3.1', 'id.order.consec')]
xUMc_4 <- xUMc_4[ , c('UMc_4.1', 'id.order.consec')]
xUMc_5 <- xUMc_5[ , c('UMc_5.1', 'id.order.consec')]
xUMr_1 <- xUMr_1[ , c('UMr_1.1', 'id.order.consec')]
xUMr_2 <- xUMr_2[ , c('UMr_2.1', 'id.order.consec')]
xUMr_3 <- xUMr_3[ , c('UMr_3.1', 'id.order.consec')]
xUMr_4 <- xUMr_4[ , c('UMr_4.1', 'id.order.consec')]
xUMr_5 <- xUMr_5[ , c('UMr_5.1', 'id.order.consec')]

# keep duplicates from merged to identify ALL reactive regions

# join dfs in long format so that each epitope is a new line
# pivot longer
xUMc_1 <- xUMc_1 %>% pivot_longer(cols = -c(id.order.consec), names_to = 'Samples', values_to = 'Reactive')
xUMc_2 <- xUMc_2 %>% pivot_longer(cols = -c(id.order.consec), names_to = 'Samples', values_to = 'Reactive')
xUMc_3 <- xUMc_3 %>% pivot_longer(cols = -c(id.order.consec), names_to = 'Samples', values_to = 'Reactive')
xUMc_4 <- xUMc_4 %>% pivot_longer(cols = -c(id.order.consec), names_to = 'Samples', values_to = 'Reactive')
xUMc_5 <- xUMc_5 %>% pivot_longer(cols = -c(id.order.consec), names_to = 'Samples', values_to = 'Reactive')
xUMr_1 <- xUMr_1 %>% pivot_longer(cols = -c(id.order.consec), names_to = 'Samples', values_to = 'Reactive')
xUMr_2 <- xUMr_2 %>% pivot_longer(cols = -c(id.order.consec), names_to = 'Samples', values_to = 'Reactive')
xUMr_3 <- xUMr_3 %>% pivot_longer(cols = -c(id.order.consec), names_to = 'Samples', values_to = 'Reactive')
xUMr_4 <- xUMr_4 %>% pivot_longer(cols = -c(id.order.consec), names_to = 'Samples', values_to = 'Reactive')
xUMr_5 <- xUMr_5 %>% pivot_longer(cols = -c(id.order.consec), names_to = 'Samples', values_to = 'Reactive')

# merge all dfs
merge <- rbind(xUMc_1, xUMc_2, xUMc_3, xUMc_4, xUMc_5, xUMr_1, xUMr_2, xUMr_3, xUMr_4, xUMr_5) 

# total observations in "merged" is the total reactive regions = 81,362 total reactive regions

        # split "merged" back into protein.name
        # split at the SECOND hyphen - see regex
        # find explanation of regex at https://stackoverflow.com/questions/41028681/split-string-last-delimiter
        merge <- merge %>% 
          separate(id.order.consec, into = c("protein.name", "merged"), sep = "_(?!.*_)", remove = FALSE)
        
        # determine total number of proteins that are reactive in at least one sample
        merge <- distinct(merge, protein.name, .keep_all= TRUE)  
        
                # 6363 total proteins with at least one epitope 

# determine total reactive region count AND range of region
# total observations in "UMx_x.x" is the total reactive regions for that sample
# max and min value within "consec" column of "UMx_x.x" is the range of reactive regions
    # UMc_1.1 - 4280, 7-142
    # UMc_2.1 - 8617, 7-145
    # UMc_3.1 - 9788, 7-145
    # UMc_4.1 - 4231, 7-48
    # UMc_5.1 - 7641, 7-86
    # UMr_1.1 - 3332, 7-46
    # UMr_2.1 - 3700, 7-63
    # UMr_3.1 - 16197, 7-121
    # UMr_4.1 - 7685, 7-69
    # UMr_5.1 - 15891, 7-89
        
    # total reactive regions: 81,362 
    # range of reactive regions: 3,332-16,197

# ------------------- Numbers: Immunoreactive regions per gene ------------
# start with immunogenicity set up
# keep relevant columns for PROTEINS
xUMc_1 <- UMc_1.1[ , c('protein.name', 'UMc_1.1', 'order', 'consec')]
xUMc_2 <- UMc_2.1[ , c('protein.name', 'UMc_2.1', 'order', 'consec')]
xUMc_3 <- UMc_3.1[ , c('protein.name', 'UMc_3.1', 'order', 'consec')]
xUMc_4 <- UMc_4.1[ , c('protein.name', 'UMc_4.1', 'order', 'consec')]
xUMc_5 <- UMc_5.1[ , c('protein.name', 'UMc_5.1', 'order', 'consec')]
xUMr_1 <- UMr_1.1[ , c('protein.name', 'UMr_1.1', 'order', 'consec')]
xUMr_2 <- UMr_2.1[ , c('protein.name', 'UMr_2.1', 'order', 'consec')]
xUMr_3 <- UMr_3.1[ , c('protein.name', 'UMr_3.1', 'order', 'consec')]
xUMr_4 <- UMr_4.1[ , c('protein.name', 'UMr_4.1', 'order', 'consec')]
xUMr_5 <- UMr_5.1[ , c('protein.name', 'UMr_5.1', 'order', 'consec')]

# count reactive regions per gene
xUMc_1 <- xUMc_1 %>% group_by(protein.name) %>% count(UMc_1.1, sort = TRUE)
xUMc_2 <- xUMc_2 %>% group_by(protein.name) %>% count(UMc_2.1, sort = TRUE)
xUMc_3 <- xUMc_3 %>% group_by(protein.name) %>% count(UMc_3.1, sort = TRUE)
xUMc_4 <- xUMc_4 %>% group_by(protein.name) %>% count(UMc_4.1, sort = TRUE)
xUMc_5 <- xUMc_5 %>% group_by(protein.name) %>% count(UMc_5.1, sort = TRUE)
xUMr_1 <- xUMr_1 %>% group_by(protein.name) %>% count(UMr_1.1, sort = TRUE)
xUMr_2 <- xUMr_2 %>% group_by(protein.name) %>% count(UMr_2.1, sort = TRUE)
xUMr_3 <- xUMr_3 %>% group_by(protein.name) %>% count(UMr_3.1, sort = TRUE)
xUMr_4 <- xUMr_4 %>% group_by(protein.name) %>% count(UMr_4.1, sort = TRUE)
xUMr_5 <- xUMr_5 %>% group_by(protein.name) %>% count(UMr_5.1, sort = TRUE)

# rename column: UMx_x = number of epitopes per gene
xUMc_1 <- xUMc_1 %>% mutate(c1 = n) %>% select(-UMc_1.1, -n)
xUMc_2 <- xUMc_2 %>% mutate(c2 = n) %>% select(-UMc_2.1, -n)
xUMc_3 <- xUMc_3 %>% mutate(c3 = n) %>% select(-UMc_3.1, -n)
xUMc_4 <- xUMc_4 %>% mutate(c4 = n) %>% select(-UMc_4.1, -n)
xUMc_5 <- xUMc_5 %>% mutate(c5 = n) %>% select(-UMc_5.1, -n)
xUMr_1 <- xUMr_1 %>% mutate(r1 = n) %>% select(-UMr_1.1, -n)
xUMr_2 <- xUMr_2 %>% mutate(r2 = n) %>% select(-UMr_2.1, -n)
xUMr_3 <- xUMr_3 %>% mutate(r3 = n) %>% select(-UMr_3.1, -n)
xUMr_4 <- xUMr_4 %>% mutate(r4 = n) %>% select(-UMr_4.1, -n)
xUMr_5 <- xUMr_5 %>% mutate(r5 = n) %>% select(-UMr_5.1, -n)

# List of your data frames
list <- list(xUMc_1, xUMc_2, xUMc_3, xUMc_4, xUMc_5, xUMr_1, xUMr_2, xUMr_3, xUMr_4, xUMr_5)

# merge them all by protein
merged.list <- reduce(list, full_join, by = "protein.name")

# add average column at end (ave epitope per gene, across samples)
merged.list$ave <- rowMeans(merged.list[, c("c1", "c2", "c3", "c4", "c5", "r1", "r2", "r3", "r4", "r5")])

    # turn NAs to 0s
    merged.list[is.na(merged.list)] <- 0

    # average reactive regions across patients across genes
    mean(merged.list$ave)
    
        # average reactive region per gene, across patients = 0.388 regions      
    
# add average column across columns (average epitope per sample, across genes)
    # create list of column averages
    colMeans(merged.list[,c("c1", "c2", "c3", "c4", "c5", "r1", "r2", "r3", "r4", "r5")], na.rm = TRUE)
    
        # average reactive region per sample, across all genes
        #        c1        c2        c3        c4        c5        r1        r2        r3        r4        r5 
        #      0.6726387 1.3542354 1.5382681 0.6649379 1.2008487 0.5236524 0.5814867 2.5454974 1.2077636 2.4974069 

        # average reactive regions across samples
        mean(0.6726387, 1.3542354, 1.5382681, 0.6649379, 1.2008487, 0.5236524, 0.5814867, 2.5454974, 1.2077636, 2.4974069 )
    
        # average reactive region per gene, across samples = 0.673 regions
    
# ------------------- Numbers: Immunoreactive segments/proteins shared across patients ------------------------------
# start with immunogenicity set up
# keep relevant columns for PROTEINS
xUMc_1 <- UMc_1.1[ , c('merged', 'UMc_1.1')]
xUMc_2 <- UMc_2.1[ , c('merged', 'UMc_2.1')]
xUMc_3 <- UMc_3.1[ , c('merged', 'UMc_3.1')]
xUMc_4 <- UMc_4.1[ , c('merged', 'UMc_4.1')]
xUMc_5 <- UMc_5.1[ , c('merged', 'UMc_5.1')]
xUMr_1 <- UMr_1.1[ , c('merged', 'UMr_1.1')]
xUMr_2 <- UMr_2.1[ , c('merged', 'UMr_2.1')]
xUMr_3 <- UMr_3.1[ , c('merged', 'UMr_3.1')]
xUMr_4 <- UMr_4.1[ , c('merged', 'UMr_4.1')]
xUMr_5 <- UMr_5.1[ , c('merged', 'UMr_5.1')]

# remove duplicates from merged so there are no duplicate segments in each sample
xUMc_1 <- distinct(xUMc_1, merged, .keep_all= TRUE)
xUMc_2 <- distinct(xUMc_2, merged, .keep_all= TRUE)
xUMc_3 <- distinct(xUMc_3, merged, .keep_all= TRUE)
xUMc_4 <- distinct(xUMc_4, merged, .keep_all= TRUE)
xUMc_5 <- distinct(xUMc_5, merged, .keep_all= TRUE)
xUMr_1 <- distinct(xUMr_1, merged, .keep_all= TRUE)
xUMr_2 <- distinct(xUMr_2, merged, .keep_all= TRUE)
xUMr_3 <- distinct(xUMr_3, merged, .keep_all= TRUE)
xUMr_4 <- distinct(xUMr_4, merged, .keep_all= TRUE)
xUMr_5 <- distinct(xUMr_5, merged, .keep_all= TRUE)

# join dfs and convert NAs to 0
all <- full_join(xUMc_1, xUMc_2, by = "merged")
all <- full_join(all, xUMc_3, by = "merged")
all <- full_join(all, xUMc_4, by = "merged")
all <- full_join(all, xUMc_5, by = "merged")
all <- full_join(all, xUMr_1, by = "merged")
all <- full_join(all, xUMr_2, by = "merged")
all <- full_join(all, xUMr_3, by = "merged")
all <- full_join(all, xUMr_4, by = "merged")
all <- full_join(all, xUMr_5, by = "merged")
all[is.na(all)] <- 0

# add total column
all$total <- rowSums(all[, c("UMc_1.1", "UMc_2.1", "UMc_3.1", "UMc_4.1", "UMc_5.1", 
                             "UMr_1.1", "UMr_2.1", "UMr_3.1", "UMr_4.1", "UMr_5.1")])

# how many are recognized by sample - reactive segments, PVP01 only
table(all$total)

        #        1     2     3     4     5     6     7     8     9    10 
        #     12941  9401  5584  2794  1285   527   221    75    16     3  

# split "merged" back into protein.name
# split at the SECOND hyphen - see regex
# find explanation of regex at https://stackoverflow.com/questions/41028681/split-string-last-delimiter
all <- all %>% 
  separate(merged, into = c("protein.name", "segment"), sep = "_(?!.*_)", remove = FALSE)

# select columns
all <- all[ , c('protein.name', 'merged', 'total')]

    # remove duplicates from proteins, keeping most reactive segment per protein
    all <- all %>% group_by(protein.name) %>% slice_max(total, n=1, with_ties = FALSE)
    
    # table most reactive segment per gene
    table(all$total)
    
        # table of genes, keeping only most reactive segment per gene
            #   1    2    3    4    5    6    7    8    9   10 
            # 860 1388 1516 1186  748  382  192   72   16    3 
    
        # total of 6363 genes with at least one epitope across samples
        
# add 0 shared back in
# upload the ORDERED PROCESSED dataset
data <- fread("/Path/peparray_ordered.txt")

    # select gene column
    data <- data[,'protein.name']
    
    # remove duplicates    
    data <- distinct(data, protein.name, .keep_all= TRUE)
        
    # filter for all pv genes ran on array
    pv <- data %>% filter(grepl('PVP01', protein.name))
    
        # should be 6670 genes on array - CORRECT
    
    # find genes in pv (all genes) and not in all (immunoreactive genes)
    extra <- anti_join(pv, all)
    
        # 307 genes that are on peparray with no reactive segments (0 shared)

# add 0 shared genes back to list of all
merged <- rbind(all, extra)
    
# turn NAs to 0
merged <- merged %>% mutate(across(where(is.numeric), ~replace_na(.x, 0)))
    
    # should be 6670 genes again - CORRECT 

# filter for significantly reactive segments (7+ samples reactive within same segment, see FDR for calculation)
merged <- all %>% filter(total > 6)
  
    # check for duplicate proteins - no duplicates 
    merged[duplicated(merged$protein.name),]
              
    # 315 reactive segments total
    # 283 reactive proteins total (most reactive segment per gene)
    
# export most reactive segments per protein 
write.table(merged, '/Path/symptomatic_pv_most_reactive_seg_per_protein.txt', sep = "\t", row.names = FALSE, col.names = TRUE)

# ------------------- Numbers: FDR for significance --------
# determine significance of shared samples
# find probability of a segment being immunoreactive

# run immunogenicity set up first
# total pv segments = 45,222

# find average number of immunoreactive segments across sample
# keep relevant columns for PROTEINS
xUMc_1 <- UMc_1.1[ , c('merged', 'UMc_1.1')]
xUMc_2 <- UMc_2.1[ , c('merged', 'UMc_2.1')]
xUMc_3 <- UMc_3.1[ , c('merged', 'UMc_3.1')]
xUMc_4 <- UMc_4.1[ , c('merged', 'UMc_4.1')]
xUMc_5 <- UMc_5.1[ , c('merged', 'UMc_5.1')]
xUMr_1 <- UMr_1.1[ , c('merged', 'UMr_1.1')]
xUMr_2 <- UMr_2.1[ , c('merged', 'UMr_2.1')]
xUMr_3 <- UMr_3.1[ , c('merged', 'UMr_3.1')]
xUMr_4 <- UMr_4.1[ , c('merged', 'UMr_4.1')]
xUMr_5 <- UMr_5.1[ , c('merged', 'UMr_5.1')]

# remove duplicates from merged so there are no duplicate segments in each sample
xUMc_1 <- distinct(xUMc_1, merged, .keep_all= TRUE)
xUMc_2 <- distinct(xUMc_2, merged, .keep_all= TRUE)
xUMc_3 <- distinct(xUMc_3, merged, .keep_all= TRUE)
xUMc_4 <- distinct(xUMc_4, merged, .keep_all= TRUE)
xUMc_5 <- distinct(xUMc_5, merged, .keep_all= TRUE)
xUMr_1 <- distinct(xUMr_1, merged, .keep_all= TRUE)
xUMr_2 <- distinct(xUMr_2, merged, .keep_all= TRUE)
xUMr_3 <- distinct(xUMr_3, merged, .keep_all= TRUE)
xUMr_4 <- distinct(xUMr_4, merged, .keep_all= TRUE)
xUMr_5 <- distinct(xUMr_5, merged, .keep_all= TRUE)

# number of observations in "UMx_x.x" is number of reactive segments

# find average reactive segment
num <- c(4005, 7760, 8441, 4020, 6814,
     3197, 3500, 13627, 6962, 13243)
mean(num)

      # mean reactive segments = 7156.9

# find probability of being immunoreactive (mean reac segment/total segment)
7156.9/45222

      # probability of being immunoreactive = 0.1582615

# create df for trial runs
trial <- tibble(
  sucesses = c(1,2,3,4,5,6,7,8,9,10),
  fails = c(9,8,7,6,5,4,3,2,1,0)
)

# add probability
trial$prob <- dbinom(x = trial$sucesses, size = 10, prob = 0.1582615)

# add expected immunoreactive segments
trial$expected <- trial$prob * 45222

# add observed immunoreactive segments
# join dfs and convert NAs to 0
all <- full_join(xUMc_1, xUMc_2, by = "merged")
all <- full_join(all, xUMc_3, by = "merged")
all <- full_join(all, xUMc_4, by = "merged")
all <- full_join(all, xUMc_5, by = "merged")
all <- full_join(all, xUMr_1, by = "merged")
all <- full_join(all, xUMr_2, by = "merged")
all <- full_join(all, xUMr_3, by = "merged")
all <- full_join(all, xUMr_4, by = "merged")
all <- full_join(all, xUMr_5, by = "merged")
all[is.na(all)] <- 0

# add total column
all$total <- rowSums(all[, c("UMc_1.1", "UMc_2.1", "UMc_3.1", "UMc_4.1", "UMc_5.1", 
                             "UMr_1.1", "UMr_2.1", "UMr_3.1", "UMr_4.1", "UMr_5.1")])

# how many are recognized by sample - reactive segments, PVP01 only
table(all$total)

    #        1     2     3     4     5     6     7     8     9    10 
    #     12941  9401  5584  2794  1285   527   221    75    16     3  

# add observed to trial df
trial$observed <- c(12941, 9401, 5584, 2794, 1285, 527, 221, 75, 16, 3)

# add FDR (expected/observed)
trial$FDR <- trial$expected/trial$observed

    # FDR cutoff is <0.05 = segments with 7+ patients reactive are significant

# ------------------- Fig S3: histogram of reactive region length ------------
# start with immunogenicity set up
# keep relevant columns for SEGMENTS
xUMc_1.1 <- UMc_1.1[ , c('consec')]
xUMc_2.1 <- UMc_2.1[ , c('consec')]
xUMc_3.1 <- UMc_3.1[ , c('consec')]
xUMc_4.1 <- UMc_4.1[ , c('consec')]
xUMc_5.1 <- UMc_5.1[ , c('consec')]
xUMr_1.1 <- UMr_1.1[ , c('consec')]
xUMr_2.1 <- UMr_2.1[ , c('consec')]
xUMr_3.1 <- UMr_3.1[ , c('consec')]
xUMr_4.1 <- UMr_4.1[ , c('consec')]
xUMr_5.1 <- UMr_5.1[ , c('consec')]

# join dfs and convert NAs to 0
epi <- rbind(xUMc_1.1, xUMc_2.1, xUMc_3.1, xUMc_4.1, xUMc_5.1,
             xUMr_1.1, xUMr_2.1, xUMr_3.1, xUMr_4.1, xUMr_5.1)

# histogram of epitope lengths across all epitopes of all UM samples
gghistogram(epi, x = "consec", fill = "black", bins =100, size = 1.2) +
      labs(title = "Number of consecutively positive peptides <br> across symptomatic *P. vivax* patients",
           x = "No. of consecutive seoreactive peptides", y = "Count") +
      theme(plot.caption = element_markdown(lineheight = 1.5),
           plot.title = element_markdown(size=15, face="bold", hjust=0.5, lineheight = 1.1, color="#666666"),
           axis.title.x = element_markdown(face="bold"),
           axis.title.y = element_markdown(face="bold")) +
      scale_x_continuous(breaks = scales::pretty_breaks(n = 10), limits = c(7, 56)) +
      scale_y_continuous(limits = c(0, 15000))

# ------------------- Fig S3: IQR of histogram of reactive region length ------------
# determine IQR (difference between the 75th and 25th percentiles) for epitope length
# run immunogenicity set up first
# keep relevant columns for SEGMENTS
xUMc_1.1 <- UMc_1.1[ , c('consec')]
xUMc_2.1 <- UMc_2.1[ , c('consec')]
xUMc_3.1 <- UMc_3.1[ , c('consec')]
xUMc_4.1 <- UMc_4.1[ , c('consec')]
xUMc_5.1 <- UMc_5.1[ , c('consec')]
xUMr_1.1 <- UMr_1.1[ , c('consec')]
xUMr_2.1 <- UMr_2.1[ , c('consec')]
xUMr_3.1 <- UMr_3.1[ , c('consec')]
xUMr_4.1 <- UMr_4.1[ , c('consec')]
xUMr_5.1 <- UMr_5.1[ , c('consec')]

# join dfs and convert NAs to 0
epi <- rbind(xUMc_1.1, xUMc_2.1, xUMc_3.1, xUMc_4.1, xUMc_5.1,
             xUMr_1.1, xUMr_2.1, xUMr_3.1, xUMr_4.1, xUMr_5.1)

# find IQR
# pull upper and lower quantile
q1 <- quantile(epi$consec, 0.25)
q3 <- quantile(epi$consec, 0.75)

# Calculate the IQR manually
q3 - q1

      # IQR = 3
      # Interpretation: middle 50% of data has spread of 3 units 

# ----------- Immunogenicity set up - Asymptomatic Pv patients, Pv probes -------------------------
# upload ORDERED dataset
data <- fread("/Path/peparray_ordered.txt")

# select AM initial time point and NA only for ordered dataset
data <- data[ , c('protein.name', 'order',
                  'NA_ctrl_1','NA_ctrl_2','NA_ctrl_3','NA_ctrl_4','NA_ctrl_5','NA_ctrl_6','NA_ctrl_7','NA_ctrl_8','NA_ctrl_9','NA_ctrl_10',
                  'Pv_AM_1', 'Pv_AM_2', 'Pv_AM_3', 'Pv_AM_4', 'Pv_AM_5'), with = FALSE]

# option to filter for pv
data <- data %>% filter(grepl("PVP01", protein.name))

# create 500 MFI column as NA_11 as extra cutoff
data$NA_ctrl_11 <- 500

# create max control and max UM and 500 MFI control column
data$NA.max <- pmax(data$NA_ctrl_1, data$NA_ctrl_2, data$NA_ctrl_3, data$NA_ctrl_4, data$NA_ctrl_5, data$NA_ctrl_6, data$NA_ctrl_7, data$NA_ctrl_8, data$NA_ctrl_9, data$NA_ctrl_10, data$NA_ctrl_11)

# create column to indicate if sample is greater than max NA and less than 100 mfi in NA 
# ordered dataset
data$AM_1 <- with(data, ifelse(data$Pv_AM_1 > data$NA.max, 1, 0))
data$AM_2 <- with(data, ifelse(data$Pv_AM_2 > data$NA.max, 1, 0))
data$AM_3 <- with(data, ifelse(data$Pv_AM_3 > data$NA.max, 1, 0))
data$AM_4 <- with(data, ifelse(data$Pv_AM_4 > data$NA.max, 1, 0))
data$AM_5 <- with(data, ifelse(data$Pv_AM_5 > data$NA.max, 1, 0))

# select new columns only
data <- data[ , c('protein.name', 'order',
                  'AM_1', 'AM_2', 'AM_3', 'AM_4', 'AM_5'), with = FALSE]

# make df
data <- as.data.frame(data)

# determine consecutively reactive epitopes longer than >6 peptides
AM_1 <- data %>%
  select(AM_1, order, protein.name) %>%
  mutate(grp = data.table::rleid(AM_1)) %>%
  filter(AM_1 ==1) %>%
  group_by(protein.name, grp) %>%
  mutate(consec = row_number())
AM_1 <- AM_1 %>%
  select(AM_1, order, protein.name, grp, consec) %>%
  filter(consec > 6) %>%
  summarise_all(last)  

AM_2 <- data %>%
  select(AM_2, order, protein.name) %>%
  mutate(grp = data.table::rleid(AM_2)) %>%
  filter(AM_2 ==1) %>%
  group_by(protein.name, grp) %>%
  mutate(consec = row_number())
AM_2 <- AM_2 %>%
  select(AM_2, order, protein.name, grp, consec) %>%
  filter(consec > 6) %>%
  summarise_all(last)  

AM_3 <- data %>%
  select(AM_3, order, protein.name) %>%
  mutate(grp = data.table::rleid(AM_3)) %>%
  filter(AM_3 ==1) %>%
  group_by(protein.name, grp) %>%
  mutate(consec = row_number())
AM_3 <- AM_3 %>%
  select(AM_3, order, protein.name, grp, consec) %>%
  filter(consec > 6) %>%
  summarise_all(last)  

AM_4 <- data %>%
  select(AM_4, order, protein.name) %>%
  mutate(grp = data.table::rleid(AM_4)) %>%
  filter(AM_4 ==1) %>%
  group_by(protein.name, grp) %>%
  mutate(consec = row_number())
AM_4 <- AM_4 %>%
  select(AM_4, order, protein.name, grp, consec) %>%
  filter(consec > 6) %>%
  summarise_all(last)  

AM_5 <- data %>%
  select(AM_5, order, protein.name) %>%
  mutate(grp = data.table::rleid(AM_5)) %>%
  filter(AM_5 ==1) %>%
  group_by(protein.name, grp) %>%
  mutate(consec = row_number())
AM_5 <- AM_5 %>%
  select(AM_5, order, protein.name, grp, consec) %>%
  filter(consec > 6) %>%
  summarise_all(last)  

# chunk proteins into 100aa segments and merge "protein.name" and "id"
AM_1 <- AM_1 %>%
  group_by(protein.name) %>%
  mutate(id = as.integer(order/100))
AM_1 <- AM_1 %>%
  unite('merged', protein.name, id, remove = FALSE)

AM_2 <- AM_2 %>%
  group_by(protein.name) %>%
  mutate(id = as.integer(order/100))
AM_2 <- AM_2 %>%
  unite('merged', protein.name, id, remove = FALSE)

AM_3 <- AM_3 %>%
  group_by(protein.name) %>%
  mutate(id = as.integer(order/100))
AM_3 <- AM_3 %>%
  unite('merged', protein.name, id, remove = FALSE)

AM_4 <- AM_4 %>%
  group_by(protein.name) %>%
  mutate(id = as.integer(order/100))
AM_4 <- AM_4 %>%
  unite('merged', protein.name, id, remove = FALSE)

AM_5 <- AM_5 %>%
  group_by(protein.name) %>%
  mutate(id = as.integer(order/100))
AM_5 <- AM_5 %>%
  unite('merged', protein.name, id, remove = FALSE)

# ------------------- Table S4: immunogenicity of asymptomatic Pv patients across P01 genome -------------
# start from immunogenicity set up - AM, PVP01 filtered
# keep relevant columns for SEGMENTS
xAM_1 <- AM_1[ , c('merged', 'AM_1')]
xAM_2 <- AM_2[ , c('merged', 'AM_2')]
xAM_3 <- AM_3[ , c('merged', 'AM_3')]
xAM_4 <- AM_4[ , c('merged', 'AM_4')]
xAM_5 <- AM_5[ , c('merged', 'AM_5')]

# remove duplicates from merged so there are no duplicate segments per sample
xAM_1 <- distinct(xAM_1, merged, .keep_all= TRUE)
xAM_2 <- distinct(xAM_2, merged, .keep_all= TRUE)
xAM_3 <- distinct(xAM_3, merged, .keep_all= TRUE)
xAM_4 <- distinct(xAM_4, merged, .keep_all= TRUE)
xAM_5 <- distinct(xAM_5, merged, .keep_all= TRUE)

# join dfs and convert NAs to 0
all <- full_join(xAM_1, xAM_2, by = "merged")
all <- full_join(all, xAM_3, by = "merged")
all <- full_join(all, xAM_4, by = "merged")
all <- full_join(all, xAM_5, by = "merged")
all[is.na(all)] <- 0

# add total column
all$total <- rowSums(all[, c("AM_1", "AM_2", "AM_3", "AM_4", "AM_5")])

# how many are recognized by sample - reactive segments, PVP01 only
table(all$total)
            
      # reactive pv segments = 29,407
            #    1     2     3     4     5 
            #  14439  9109  4027  1347   483 

# add in protein.name
  # for weird genes
  others <- all %>% filter(str_detect(merged, 'PVP01_API|PVP01_MIT'))  
  
        # select first 14 characters of protein.name (numbers indicate character range to grab)
        others <- others %>% mutate(protein.name = substr(merged, 1, 14))
        
  # for main list
  # remove weird genes from list now that they're pulled
  main <- all %>% filter(!str_detect(merged, "PVP01_API|PVP01_MIT"))
  
        # take first 13 characters from ID
        main <- main %>% mutate(protein.name = substr(merged, 1, 13))
          
        # merge others back into list
        all <- rbind(main, others)

# add in 0 shared
# reupload data
data <- fread("/Path/peparray_ordered.txt")

    # filter for Pv only
    data <- data %>% filter(grepl('PVP01', protein.name))
    
    # select UM initial time point and NA only for ordered dataset
    data <- data[ , c('protein.name', 'order')]
     
    # split into 100aa segments
    data <- data %>%
      group_by(protein.name) %>%
      mutate(id = as.integer(order/100))
    data <- data %>%
      unite('merged', protein.name, id, remove = FALSE)

    # remove duplicate segments    
    data <- distinct(data, merged, .keep_all= TRUE)

    # adjust array for merge
    # for weird genes
    others <- data %>% filter(str_detect(protein.name, 'PVP01_API|PVP01_MIT'))  
    
          # select first 14 characters of protein.name (numbers indicate character range to grab)
          others <- others %>% mutate(protein.name = substr(protein.name, 1, 14))
              
    # for main list
    # remove weird genes from list now that they're pulled
    main <- data %>% filter(!str_detect(protein.name, "PVP01_API|PVP01_MIT"))
    
      # take first 13 characters from ID
      main <- main %>% mutate(protein.name = substr(protein.name, 1, 13))
            
      # merge others back into list
      data <- rbind(main, others)
  
    # find genes in all (reactive genes) that are not in array (all genes)
    extra <- anti_join(data, all, by = "merged")
    
           # 15,815 segments that are on peparray with no reactive segments (0 shared)

        # add 0 shared genes back to list of all reactive segments
        final <- rbind.fill(all, extra)

        # turn NAs to 0
        final$total[is.na(final$total)] <- 0
        
            # 45222 total pv segments
        
    # filter for significantly reactive segments (7+ samples reactive within same segment, see FDR for calculation)
    final <- all %>% filter(total > 6)
    
# export most reactive segment per protein for ASYMPTOMATIC patients
write.table(merged, '/Path/asymptomatic_pv_most_reactive_seg_per_protein.txt', sep = "\t", row.names = FALSE, col.names = TRUE)

# ------------------- Numbers: AM vs random UM segment reactivity rates -----------------------------------
# first run symptomatic and asymptomatic immunogenicity set up, Pv probes
# remove duplicate reactive segments
UMc_1.1 <- distinct(UMc_1.1, merged, .keep_all= TRUE)
UMc_2.1 <- distinct(UMc_2.1, merged, .keep_all= TRUE)
UMc_3.1 <- distinct(UMc_3.1, merged, .keep_all= TRUE)
UMc_4.1 <- distinct(UMc_4.1, merged, .keep_all= TRUE)
UMc_5.1 <- distinct(UMc_5.1, merged, .keep_all= TRUE)
UMr_1.1 <- distinct(UMr_1.1, merged, .keep_all= TRUE)
UMr_2.1 <- distinct(UMr_2.1, merged, .keep_all= TRUE)
UMr_3.1 <- distinct(UMr_3.1, merged, .keep_all= TRUE)
UMr_4.1 <- distinct(UMr_4.1, merged, .keep_all= TRUE)
UMr_5.1 <- distinct(UMr_5.1, merged, .keep_all= TRUE)
AM_1 <- distinct(AM_1, merged, .keep_all= TRUE)
AM_2 <- distinct(AM_2, merged, .keep_all= TRUE)
AM_3 <- distinct(AM_3, merged, .keep_all= TRUE)
AM_4 <- distinct(AM_4, merged, .keep_all= TRUE)
AM_5 <- distinct(AM_5, merged, .keep_all= TRUE)

# select columns
c1 <- UMc_1.1[ , c('merged', 'UMc_1.1')]
c2 <- UMc_2.1[ , c('merged', 'UMc_2.1')]
c3 <- UMc_3.1[ , c('merged', 'UMc_3.1')]
c4 <- UMc_4.1[ , c('merged', 'UMc_4.1')]
c5 <- UMc_5.1[ , c('merged', 'UMc_5.1')]
r1 <- UMr_1.1[ , c('merged', 'UMr_1.1')]
r2 <- UMr_2.1[ , c('merged', 'UMr_2.1')]
r3 <- UMr_3.1[ , c('merged', 'UMr_3.1')]
r4 <- UMr_4.1[ , c('merged', 'UMr_4.1')]
r5 <- UMr_5.1[ , c('merged', 'UMr_5.1')]
a1 <- AM_1[ , c('merged', 'AM_1')]
a2 <- AM_2[ , c('merged', 'AM_2')]
a3 <- AM_3[ , c('merged', 'AM_3')]
a4 <- AM_4[ , c('merged', 'AM_4')]
a5 <- AM_5[ , c('merged', 'AM_5')]

# join dfs and convert NAs to 0
um <- full_join(c1, c2, by = "merged")
um <- full_join(um, c3, by = "merged")
um <- full_join(um, c4, by = "merged")
um <- full_join(um, c5, by = "merged")
um <- full_join(um, r1, by = "merged")
um <- full_join(um, r2, by = "merged")
um <- full_join(um, r3, by = "merged")
um <- full_join(um, r4, by = "merged")
um <- full_join(um, r5, by = "merged")
um[is.na(um)] <- 0

am <- full_join(a1, a2, by = "merged")
am <- full_join(am, a3, by = "merged")
am <- full_join(am, a4, by = "merged")
am <- full_join(am, a5, by = "merged")
am[is.na(am)] <- 0

# randomly select 5 um columns for comparison to am

    # run 5 times to get random columns for sympto
    # set seed for reproducibility (pseudo-random generation, set base)
    set.seed(123)
    
        # list of random columns
        um.cols <- sample(2:ncol(um), 5, replace = FALSE)
    
        # filter for selected random columns
        um.new <- um %>%
          select(all_of(um.cols))
    
        # add total column
        um.new$total.um <- rowSums(um.new[, c(1:5)])
        
        # remove columns with 0 reactive
        um.new <- um.new %>% filter(total.um >0)

# add total to am
am$total.am <- rowSums(am[, c(2:6)])

  # COMPARE SEGMENT REACTIVITY RATES
      # number of reactive segments for AM samples
      table(am$total.am)
            #        1     2     3     4     5 
            #      14439  9109  4027  1347   483 

      # reactivity rates for 5 random UM samples (set 1)
      table(um.new$total.um)
            #     1     2     3     4     5 
            #  15802  8704  3167   789    82 
      
      # reactivity rates for 5 random UM samples (set 2)
      table(um.new$total.um)
            #     1     2     3     4     5 
            #  14625  5714  1454   272    32 
            
      # reactivity rates for 5 random UM samples (set 3)
      table(um.new$total.um)
            #     1     2     3     4     5 
            #  15809  8403  3211   788   146 
      
      # reactivity rates for 5 random UM samples (set 4)
      table(um.new$total.um)
            #     1     2     3     4     5 
            #   15221  7730  2799   838   160 
            
      # reactivity rates for 5 random UM samples (set 5)
      table(um.new$total.um)
            #     1     2     3     4     5 
            #   12570  4302  1242   269    29 
            
# find average 5 shared in random UM sets
mean(82, 32, 146, 160, 29)

    # average average 5 shared of random sets = 82

# find fold diff
483/82

    # asymptomatic have a 5.89-fold increase in immunoreactive segments compared to symptomatic

# ----------- Immunogenicity plots: Symptomatic Pv patients, time 1, repeat peptide marker --------------------------------------
# run immunogenicity set up for symptomatic patients first
# prepare repeat peptide information
# upload repeat peptide info
rep <- fread('/Path/repeat_peptide_list.txt')

# select peptide for gene of interest
gene <- rep %>% filter(grepl('PVP01_0000170', protein.name))
  
  # select repeat peptides only
  gene <- gene %>% filter(grepl('repeat', peptide))
  
  # add peptide information as a binary
  gene$peptide <- ifelse(gene$peptide == 'unique', 0, 1)
  
  # select columns
  gene <- gene[ , c('protein.name', 'order', 'peptide')]
  
# start from immunogenicity info with UMc_1, UMc_2, ect 
# pull out reactivity info from "UMc_" file from above
c1 <- UMc_1.1 %>% filter(grepl('PVP01_0000170', protein.name))
c2 <- UMc_2.1 %>% filter(grepl('PVP01_0000170', protein.name))
c3 <- UMc_3.1 %>% filter(grepl('PVP01_0000170', protein.name))
c4 <- UMc_4.1 %>% filter(grepl('PVP01_0000170', protein.name))
c5 <- UMc_5.1 %>% filter(grepl('PVP01_0000170', protein.name))
r1 <- UMr_1.1 %>% filter(grepl('PVP01_0000170', protein.name))
r2 <- UMr_2.1 %>% filter(grepl('PVP01_0000170', protein.name))
r3 <- UMr_3.1 %>% filter(grepl('PVP01_0000170', protein.name))
r4 <- UMr_4.1 %>% filter(grepl('PVP01_0000170', protein.name))
r5 <- UMr_5.1 %>% filter(grepl('PVP01_0000170', protein.name))
pep <- gene %>% filter(grepl('PVP01_0000170', protein.name))

  # create row id for each consecutive peptide
  c1 <-c1 %>% uncount(consec, .id = "id")
  c2 <-c2 %>% uncount(consec, .id = "id")
  c3 <-c3 %>% uncount(consec, .id = "id")
  c4 <-c4 %>% uncount(consec, .id = "id")
  c5 <-c5 %>% uncount(consec, .id = "id")
  r1 <-r1 %>% uncount(consec, .id = "id")
  r2 <-r2 %>% uncount(consec, .id = "id")
  r3 <-r3 %>% uncount(consec, .id = "id")
  r4 <-r4 %>% uncount(consec, .id = "id")
  r5 <-r5 %>% uncount(consec, .id = "id")
  
  # create new order to indicate which amino acid
  c1$order <- (c1$order - c1$id)
  c2$order <- (c2$order - c2$id)
  c3$order <- (c3$order - c3$id)
  c4$order <- (c4$order - c4$id)
  c5$order <- (c5$order - c5$id)
  r1$order <- (r1$order - r1$id)
  r2$order <- (r2$order - r2$id)
  r3$order <- (r3$order - r3$id)
  r4$order <- (r4$order - r4$id)
  r5$order <- (r5$order - r5$id)
  
  # select columns
  c1 <- c1[ , c('protein.name', 'order', 'grp', 'UMc_1.1')]
  c2 <- c2[ , c('protein.name', 'order', 'grp', 'UMc_2.1')]
  c3 <- c3[ , c('protein.name', 'order', 'grp', 'UMc_3.1')]
  c4 <- c4[ , c('protein.name', 'order', 'grp', 'UMc_4.1')]
  c5 <- c5[ , c('protein.name', 'order', 'grp', 'UMc_5.1')]
  r1 <- r1[ , c('protein.name', 'order', 'grp', 'UMr_1.1')]
  r2 <- r2[ , c('protein.name', 'order', 'grp', 'UMr_2.1')]
  r3 <- r3[ , c('protein.name', 'order', 'grp', 'UMr_3.1')]
  r4 <- r4[ , c('protein.name', 'order', 'grp', 'UMr_4.1')]
  r5 <- r5[ , c('protein.name', 'order', 'grp', 'UMr_5.1')]
  pep <- gene[ , c('protein.name', 'order','peptide')]
  
  # pivot longer
  c1 <- c1 %>% pivot_longer(cols = -c(protein.name, grp, order), names_to = 'Samples', values_to = 'MFI')
  c2 <- c2 %>% pivot_longer(cols = -c(protein.name, grp, order), names_to = 'Samples', values_to = 'MFI')
  c3 <- c3 %>% pivot_longer(cols = -c(protein.name, grp, order), names_to = 'Samples', values_to = 'MFI')
  c4 <- c4 %>% pivot_longer(cols = -c(protein.name, grp, order), names_to = 'Samples', values_to = 'MFI')
  c5 <- c5 %>% pivot_longer(cols = -c(protein.name, grp, order), names_to = 'Samples', values_to = 'MFI')
  r1 <- r1 %>% pivot_longer(cols = -c(protein.name, grp, order), names_to = 'Samples', values_to = 'MFI')
  r2 <- r2 %>% pivot_longer(cols = -c(protein.name, grp, order), names_to = 'Samples', values_to = 'MFI')
  r3 <- r3 %>% pivot_longer(cols = -c(protein.name, grp, order), names_to = 'Samples', values_to = 'MFI')
  r4 <- r4 %>% pivot_longer(cols = -c(protein.name, grp, order), names_to = 'Samples', values_to = 'MFI')
  r5 <- r5 %>% pivot_longer(cols = -c(protein.name, grp, order), names_to = 'Samples', values_to = 'MFI')
  pep <- pep %>% pivot_longer(cols = -c(protein.name, order), names_to = 'Samples', values_to = 'MFI')
  
  # merge all dfs
  merge <- rbind(c1, c2, c3, c4, c5, r1, r2, r3, r4, r5,pep) 
  
  # add "MFI" values to samples 
  merge$MFI <- ifelse(merge$Samples == 'peptide', -250,
                      ifelse(merge$Samples == 'UMc_1.1', -500,
                             ifelse(merge$Samples == 'UMc_2.1', -700,
                                    ifelse(merge$Samples == 'UMc_3.1', -900,
                                           ifelse(merge$Samples == 'UMc_4.1', -1100,
                                                  ifelse(merge$Samples == 'UMc_5.1', -1300,
                                                         ifelse(merge$Samples == 'UMr_1.1', -1500,
                                                                ifelse(merge$Samples == 'UMr_2.1', -1700,
                                                                       ifelse(merge$Samples == 'UMr_3.1', -1900,
                                                                              ifelse(merge$Samples == 'UMr_4.1', -2100,
                                                                                     ifelse(merge$Samples == 'UMr_5.1', -2300, NA)))))))))))
  
  # create gene plots with the location of reactivity below the plot
  data <- fread("/Users/rosita.asawa/Library/CloudStorage/OneDrive-UniversityofMarylandSchoolofMedicine/Graduate School/Peptide array/For R (old)/peparray_ordered.txt")
  
  # select UM initial time point and NA only for sliding/ordered dataset, SMOOTHED dataset
  data <- data[ , c('protein.name', 'order',
                    'NA_ctrl_1','NA_ctrl_2','NA_ctrl_3','NA_ctrl_4','NA_ctrl_5','NA_ctrl_6','NA_ctrl_7','NA_ctrl_8','NA_ctrl_9','NA_ctrl_10',
                    'Pv_UMc_1.1',	'Pv_UMc_2.1',	'Pv_UMc_3.1',	'Pv_UMc_4.1',	'Pv_UMc_5.1',	
                    'Pv_UMr_1.1',	'Pv_UMr_2.1',	'Pv_UMr_3.1',	'Pv_UMr_4.1',	'Pv_UMr_5.1'), with = FALSE]
  
  # add in columns with mfi data
  p1 <- data %>% filter(grepl('PVP01_0000170', protein.name))
  
  # pivot longer 
  p1 <- p1 %>% pivot_longer(cols = -c(protein.name, order), names_to = 'Samples', values_to = 'MFI')
  
        # check protein length
        #p1.dedup <- distinct(p1, order, .keep_all= TRUE)
  
  # merge dfs
  p1 <- p1[ , c('protein.name', 'order', 'Samples', 'MFI')]
  merge <- merge[ , c('protein.name', 'order', 'Samples', 'MFI')]
  together <- rbind(p1, merge)
  
  # pivot wider
  together <- together %>% pivot_wider(names_from = 'Samples', values_from = 'MFI')
  
  # plot it out
  # must grey out samples with 0 observations
ggplot(together, aes(x = order)) +
    geom_line(aes(y = as.numeric(NA_ctrl_1)), color = "blue") +
    geom_line(aes(y = as.numeric(NA_ctrl_2)), color = "blue") +
    geom_line(aes(y = as.numeric(NA_ctrl_3)), color = "blue") +
    geom_line(aes(y = as.numeric(NA_ctrl_4)), color = "blue") +
    geom_line(aes(y = as.numeric(NA_ctrl_5)), color = "blue") +
    geom_line(aes(y = as.numeric(NA_ctrl_6)), color = "blue") +
    geom_line(aes(y = as.numeric(NA_ctrl_7)), color = "blue") +
    geom_line(aes(y = as.numeric(NA_ctrl_8)), color = "blue") +
    geom_line(aes(y = as.numeric(NA_ctrl_9)), color = "blue") +
    geom_line(aes(y = as.numeric(NA_ctrl_10)), color = "blue") +
    geom_line(aes(y = as.numeric(Pv_UMc_1.1)), color = "darkgoldenrod1") +
    geom_line(aes(y = as.numeric(Pv_UMc_2.1)), color = "orange3") +
    geom_line(aes(y = as.numeric(Pv_UMc_3.1)), color = "orange4") +
    geom_line(aes(y = as.numeric(Pv_UMc_4.1)), color = "darkorange2") +
    geom_line(aes(y = as.numeric(Pv_UMc_5.1)), color = "darkorange3") +
    geom_line(aes(y = as.numeric(Pv_UMr_1.1)), color = "chocolate") +
    geom_line(aes(y = as.numeric(Pv_UMr_2.1)), color = "chocolate1") +
    geom_line(aes(y = as.numeric(Pv_UMr_3.1)), color = "chocolate3") +
    geom_line(aes(y = as.numeric(Pv_UMr_4.1)), color = "orange") +
    geom_line(aes(y = as.numeric(Pv_UMr_5.1)), color = "darkgoldenrod") +
    geom_line(aes(y = as.numeric(peptide)), color = "black", size =1) +
    geom_line(aes(y = as.numeric(UMc_1.1)), color = "darkgoldenrod1", size =1) +
    geom_line(aes(y = as.numeric(UMc_2.1)), color = "orange4", size =1) +
    geom_line(aes(y = as.numeric(UMc_3.1)), color = "darkorange2", size =1) +
    geom_line(aes(y = as.numeric(UMc_4.1)), color = "darkorange3", size =1) +
    geom_line(aes(y = as.numeric(UMc_5.1)), color = "orange3", size =1) +
    geom_line(aes(y = as.numeric(UMr_1.1)), color = "chocolate", size =1) +
    geom_line(aes(y = as.numeric(UMr_2.1)), color = "chocolate1", size =1) +
    geom_line(aes(y = as.numeric(UMr_3.1)), color = "chocolate3", size =1) +
    geom_line(aes(y = as.numeric(UMr_4.1)), color = "orange", size =1) +
    geom_line(aes(y = as.numeric(UMr_5.1)), color = "darkgoldenrod", size =1) +
    xlab("Amino Acid Position") +
    ylab("Fluorescence") + 
    ggtitle("Symptomatic<br>PROTEIN OF INTEREST") +
    theme_classic(base_size = 14) +
    theme(plot.title = element_text(size=15, face="bold", hjust=0.5, color="#666666")) +
    theme(plot.caption = element_markdown(lineheight = 1.5),
        plot.title = element_markdown(size=15, face="bold", hjust=0.5, lineheight = 1.1, color="#666666"),
        axis.title.x = element_markdown(face="bold"), 
        axis.title.y = element_markdown(face="bold")) +
    scale_x_continuous(breaks = scales::pretty_breaks(n = 10)) +
    scale_y_continuous(limits = c(-3000, 30000))

# ----------- Immunogenicity plots: Asymptomatic Pv patients, time 1, repeat peptide marker  --------------------------------------
# run immunogenicity set up for asymptomatic patients first
# prepare repeat peptide information
# upload repeat peptide info
rep <- fread('/Path/repeat_peptide_list.txt')

# select peptides for gene of interest
gene <- rep %>% filter(grepl('PVP01_0000170', protein.name))

# select repeat peptides only
gene <- gene %>% filter(grepl('repeat', peptide))

# add peptide information as a binary
gene$peptide <- ifelse(gene$peptide == 'unique', 0, 1)

# select columns
gene <- gene[ , c('protein.name', 'order', 'peptide')]

# start from immunogenicity info with UMc_1, UMc_2, ect 
# pull out reactivity info from "UMc_" file from above
a1 <- AM_1 %>% filter(grepl('PVP01_0000170', protein.name))
a2 <- AM_2 %>% filter(grepl('PVP01_0000170', protein.name))
a3 <- AM_3 %>% filter(grepl('PVP01_0000170', protein.name))
a4 <- AM_4 %>% filter(grepl('PVP01_0000170', protein.name))
a5 <- AM_5 %>% filter(grepl('PVP01_0000170', protein.name))
pep <- gene %>% filter(grepl('PVP01_0000170', protein.name))

# create row id for each consecutive peptide
a1 <-a1 %>% uncount(consec, .id = "id")
a2 <-a2 %>% uncount(consec, .id = "id")
a3 <-a3 %>% uncount(consec, .id = "id")
a4 <-a4 %>% uncount(consec, .id = "id")
a5 <-a5 %>% uncount(consec, .id = "id")

# create new order to indicate which amino acid
a1$order <- (a1$order - a1$id)
a2$order <- (a2$order - a2$id)
a3$order <- (a3$order - a3$id)
a4$order <- (a4$order - a4$id)
a5$order <- (a5$order - a5$id)

# select columns
a1 <- a1[ , c('protein.name', 'order', 'AM_1')]
a2 <- a2[ , c('protein.name', 'order', 'AM_2')]
a3 <- a3[ , c('protein.name', 'order', 'AM_3')]
a4 <- a4[ , c('protein.name', 'order', 'AM_4')]
a5 <- a5[ , c('protein.name', 'order', 'AM_5')]
pep <- gene[ , c('protein.name', 'order', 'peptide')]

# pivot longer
a1 <- a1 %>% pivot_longer(cols = -c(protein.name, order), names_to = 'Samples', values_to = 'MFI')
a2 <- a2 %>% pivot_longer(cols = -c(protein.name, order), names_to = 'Samples', values_to = 'MFI')
a3 <- a3 %>% pivot_longer(cols = -c(protein.name, order), names_to = 'Samples', values_to = 'MFI')
a4 <- a4 %>% pivot_longer(cols = -c(protein.name, order), names_to = 'Samples', values_to = 'MFI')
a5 <- a5 %>% pivot_longer(cols = -c(protein.name, order), names_to = 'Samples', values_to = 'MFI')
pep <- pep %>% pivot_longer(cols = -c(protein.name, order), names_to = 'Samples', values_to = 'MFI')

# merge all dfs
merge <- rbind(a1, a2, a3, a4, a5, pep) 

# add "MFI" values to samples 
merge$MFI <- ifelse(merge$Samples == 'peptide', -300,
                    ifelse(merge$Samples == 'AM_1', -500,
                           ifelse(merge$Samples == 'AM_2', -700,
                                  ifelse(merge$Samples == 'AM_3', -900,
                                         ifelse(merge$Samples == 'AM_4', -1100,
                                                ifelse(merge$Samples == 'AM_5', -1300, 'NA'))))))

# create gene plots with the location of reactivity below the plot
data <- fread("/Users/rosita.asawa/Library/CloudStorage/OneDrive-UniversityofMarylandSchoolofMedicine/Graduate School/Peptide array/For R (old)/peparray_ordered.txt")

# select UM initial time point and NA only for sliding/ordered dataset, SMOOTHED dataset
data <- data[ , c('protein.name', 'order',
                  'NA_ctrl_1','NA_ctrl_2','NA_ctrl_3','NA_ctrl_4','NA_ctrl_5','NA_ctrl_6','NA_ctrl_7','NA_ctrl_8','NA_ctrl_9','NA_ctrl_10',
                  'Pv_AM_1', 'Pv_AM_2', 'Pv_AM_3', 'Pv_AM_4', 'Pv_AM_5'), with = FALSE]

# add in columns with mfi data
p1 <- data %>% filter(grepl('PVP01_0000170', protein.name))

# pivot longer 
p1 <- p1 %>% pivot_longer(cols = -c(protein.name, order), names_to = 'Samples', values_to = 'MFI')

# merge dfs
p1 <- p1[ , c('protein.name', 'order', 'Samples', 'MFI')]
merge <- merge[ , c('protein.name', 'order', 'Samples', 'MFI')]
all <- rbind(p1, merge)

# pivot wider
all <- all %>% pivot_wider(names_from = 'Samples', values_from = 'MFI')

# plot it out
# must grey out samples with 0 observations
ggplot(all, aes(x = order)) +
    geom_line(aes(y = as.numeric(NA_ctrl_1)), color = "blue") +
    geom_line(aes(y = as.numeric(NA_ctrl_2)), color = "blue") +
    geom_line(aes(y = as.numeric(NA_ctrl_3)), color = "blue") +
    geom_line(aes(y = as.numeric(NA_ctrl_4)), color = "blue") +
    geom_line(aes(y = as.numeric(NA_ctrl_5)), color = "blue") +
    geom_line(aes(y = as.numeric(NA_ctrl_6)), color = "blue") +
    geom_line(aes(y = as.numeric(NA_ctrl_7)), color = "blue") +
    geom_line(aes(y = as.numeric(NA_ctrl_8)), color = "blue") +
    geom_line(aes(y = as.numeric(NA_ctrl_9)), color = "blue") +
    geom_line(aes(y = as.numeric(NA_ctrl_10)), color = "blue") +
    geom_line(aes(y = as.numeric(Pv_AM_1)), color = "green") +
    geom_line(aes(y = as.numeric(Pv_AM_2)), color = "lightgreen") +
    geom_line(aes(y = as.numeric(Pv_AM_3)), color = "green3") +
    geom_line(aes(y = as.numeric(Pv_AM_4)), color = "green4") +
    geom_line(aes(y = as.numeric(Pv_AM_5)), color = "darkgreen") +
    geom_line(aes(y = as.numeric(peptide)), color = "black", size =1) +
    geom_line(aes(y = as.numeric(AM_1)), color = "green", size =1) +
    geom_line(aes(y = as.numeric(AM_2)), color = "lightgreen", size =1) +
    geom_line(aes(y = as.numeric(AM_3)), color = "green3", size =1) +
    geom_line(aes(y = as.numeric(AM_4)), color = "green4", size =1) +
    geom_line(aes(y = as.numeric(AM_5)), color = "darkgreen", size =1) +
    xlab("Amino Acid Position") +
    ylab("Fluorescence") + 
    ggtitle("Asymptomatic<br>PROTEIN OF INTEREST") +
    theme_classic(base_size = 14) +
    theme(plot.title = element_text(size=15, face="bold", hjust=0.5, color="#666666")) +
    theme(plot.caption = element_markdown(lineheight = 1.5),
        plot.title = element_markdown(size=15, face="bold", hjust=0.5, lineheight = 1.1, color="#666666"),
        axis.title.x = element_markdown(face="bold"), 
        axis.title.y = element_markdown(face="bold")) +
    scale_x_continuous(breaks = scales::pretty_breaks(n = 10)) +
    scale_y_continuous(limits = c(-3000, 30000))
