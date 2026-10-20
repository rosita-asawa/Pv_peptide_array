# Script 3: Bulk RNA-seq expression analyses
# coordinating custom scripts for Plos Pathogens manuscript, published October 2026
  # Proteome-wide characterization of Plasmodium vivax antigens using a high-density peptide array
  # https://journals.plos.org/plospathogens/article?id=10.1371/journal.ppat.1014621

library(tidyverse)
library(data.table)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(reshape2)
library(ggthemes)
library(DT)
library(binsreg)
library(ggtext)

install.packages("")

# --------------- Symptomatic pv patient expression
# --------------------- Numbers: normalize to CPM and rank expression ---------
# bulk RNA-seq for cambodian symptomatic pv patients
# normalize to CPM
  # upload raw counts of cambodian UM patients
  # 6861 genes with RAW expression
  # 6670 genes with peparray reactivity
  data <- fread('/Path/Pv_CountTable_RAW.txt')
    
    # check normalization - add total to sample - NOT normalized  
    colSums(as.data.frame(data[, 7:16]))
    
    # select columns
    data <- data[,c(1, 7:16)]
    
    # transform 
    data.t <- as.data.frame(t(data))
    
    # make first row header
    colnames(data.t) <- data.t[1,]
    
    # remove first column
    data.t <- data.t[-1,]
    
    # make numeric
    data.t <- data.t %>%
      mutate_if(is.character, as.numeric) %>%
      mutate_if(is.factor, as.numeric)
        
    # add total
    data.t$sum <- rowSums(data.t[,])
    
    # formula for normalization
    # sum sample COLUMNS and normalize to CPM (counts per million)
    norm <- data.t %>%
      mutate(across(starts_with("PVP01"), ~ (.x / sum) * 1000000)) %>%
      mutate(across(starts_with("PVP01"), ~ replace_na(.x, 0))) %>%
      select(-sum)
    
    # check normalization - add total - all normalized to 1000
    rowSums(norm)
    
    # transform back
    norm.final <- as.data.frame(t(norm))
    
    # check normalization - CORRECT
    colSums(norm.final)   
    
    # make rownames first column 
    norm.final$Gene <- rownames(norm.final)

# add in immunoreactivity
# import list of immunoreactive segments and how many are shared
# code to create list found in "Peptides_immunogenicity.R"
list <- fread('/Path/symptomatic_pv_most_reactive_seg_per_protein.txt')

    # rename protein.name column
    names(list)[names(list) == "protein.name"] <- "Gene"
    names(list)[names(list) == "total"] <- "shared"
    
    # adjust list ID to match ID in data
    # select weird genes
    others <- list %>% filter(str_detect(Gene, 'PVP01_API|PVP01_MIT'))  
    
        # select first 14 characters of Gene
        others <- others %>% mutate(Gene = substr(Gene, 1, 14))
        
    # for main list
    # remove weird genes from list now that they're pulled
    list <- list[!grepl("PVP01_API|PVP01_MIT", Gene)] 
        
      # take first 13 characters from ID
      list <- list %>% mutate(Gene = substr(Gene, 1, 13))
    
    # merge others back into list
    list <- rbind(list, others)
    
    # select columns
    list <- list[ , c('Gene', 'merged', 'shared')]
    
    # remove duplicate reactive genes - leaving most reactive segment
    list <- list %>% group_by(Gene) %>% slice_max(shared, with_ties = FALSE)
    
      # find duplicate rows - no duplicates 
      list[duplicated(list$Gene),]

            # n = 6670 genes on array
    
# add shared to data.ave
shared <- left_join(norm.final, list, by = "Gene")
    
    # IMPORTANT: DO NOT removes any rows with NA
    # keep genes with expression data but NOT reactivity data

          # n = 6861 PVP01 (MIT/API included) genes with expression
          # n = 6670 PVP01 (MIT/API included) genes with expression data AND reactivity data
          # all genes on array have expression data

    # check normalization
    colSums(shared[,1:10])

        # normalization 1M = CORRECT
    
# determine expression rank
shared$rank <- rowSums(shared[,c(1:10)])

        # can now order by rank expression in object
   
# --------------------- Fig 2: binned expression vs immunoreactivity ---------
# bulk RNA-seq of symptomatic pv patients
# run expression normalization (CPM) + immunoreactivity first (found in "Peptides_immunogenicity.R")

        # 6861 genes with RAW expression
        # 6670 genes with peparray reactivity 
        # file contains 6670 genes with expression (for most) and reactivity (for all)

    # check normalization - CORRECT
      colSums(shared[,1:10])

    # turn NAs to 0s
    shared[is.na(shared)] <- 0
          
#adjust data for bins - MEAN EXPRESSION
    # bin by ranked expression DECENDING order
    # largest expression = smallest
    binned <- shared %>% mutate(bin = ntile(dplyr::desc(rank), 100))
    
    # add count of reactive per bin
    # add count of total genes per bin
    # add average of total for each bin
    # returns one row per bin
    binned <- binned %>%
      select(Gene, shared, rank, bin) %>%
      group_by(bin) %>%
      dplyr::summarise(
        bin.gene.reac = sum(shared > 6),        # count rows reactive per bin
        bin.gene.nonreac = sum(shared < 7),     # count rows nonreactive per bin
        bin.gene.sum = n(),                     # count rows (genes) per bin
        bin.express.ave = mean(rank))           # average expression per bin
    
    # make bins numeric
    binned <- binned %>% 
      mutate_at(vars(bin, bin.gene.reac, bin.gene.nonreac, bin.gene.sum, bin.express.ave), as.numeric)
   
# BINNED MEAN EXPRESSION
    # add trend line of mean expression across bins
    # make object with mean expression
    express.mean <- binned %>% mutate(mean(bin.express.ave))
    
    # make scale factor
    scale_factor_mean <- max(binned$bin.gene.reac, na.rm = TRUE) /
                         max(binned$bin.express.ave, na.rm = TRUE)
      
    # plot with markdown - BINNED + MEAN EXPRESSION
    # average CPM expression across patients, bin 100, plot number of genes reactive per bin
    ggplot(binned, aes(x = as.numeric(bin), y = bin.gene.reac)) +
      geom_bar(stat="identity", fill = "#666666") +
      geom_line(aes(y = bin.express.ave * scale_factor_mean), color = "black", size = 1.2) +
      theme_classic() +
      labs(title = "Shared immunoreactivity of *P. vivax* proteins, <br> ranked by gene expression",
          x = "*P. vivax* genes, ranked by binned expression (descending)", 
          y = "No. immunoreactive proteins") +
      theme(plot.caption = element_markdown(lineheight = 1.2),
          plot.title = element_markdown(size=15, face="bold", hjust=0.5, color="#666666"),
          axis.title.x = element_markdown(face="bold", size = 12),
          axis.title.y = element_markdown(face="bold", size = 12),
          axis.title.y.right = element_text(margin = margin(l = 10)),
          axis.text.x = element_markdown(size = 12),
          axis.text.y = element_markdown(size = 12)) +
      coord_cartesian(xlim = c(-1, 100)) +
      coord_cartesian(ylim = c(0, 11.5)) +
      scale_y_continuous(name = "No. immunoreactive proteins", sec.axis = sec_axis(~ . / scale_factor_mean, name = "Mean expression value (CPM)"))

# --------------------- Numbers: logistic regression of mean expression vs immunoreactivity -------------- 
# IMPORTANT: DO NOT BIN FOR LOGISTIC REGRESSION
# start from "shared" object creation in normalization section
shared
    
    # set binary
    shared$reac <- ifelse(shared$shared > 6, 1, 0)

            # check data
            table(shared$reac)
    
            # 283 genes are reactive - CORRECT
            # NON-reac (0)    reac (1) 
            #     6387          283 
    
    # select columns
    shared <- shared[, c('Gene', 'rank', 'reac')]
    
    # make numeric
    shared$rank <- as.numeric(shared$rank)

        # run logistic regression on LOGGED DATA (if needed?)
        # log transform
        shared$rank.log <- log2(shared$rank + 1)
    
        # select columns for LOGGED CPM 
        shared.log <- shared[, c("Gene", "rank.log", "reac")]

        # find logit for LOGGED CPM
        model.log <- glm(reac ~ rank.log, data = shared.log, family = binomial)
        sum <- summary(model.log)
        
                # Coefficients:
                #            Estimate Std. Error z value Pr(>|z|)    
                # (Intercept) -5.11599    0.27276 -18.756  < 2e-16 ***
                # total.log    0.21610    0.02702   7.997 1.28e-15 ***
              
                    # extract real p value
                    coef(summary(model.log))[,'Pr(>|z|)']
                      
                              # p = 1.276434e-15 

                    # find OR
                    exp(coef(model.log))
            
                          # (Intercept)   total.log 
                          # 0.006000064   1.241232604 
                  
                          # Interpretation: increased odds of highly expressed
                          # genes being immunoreactive (OR = 1.24) 
                        
                    # find range
                    summary(shared.log$total.log)
          
                            #   Min.   1st Qu.  Median   Mean   3rd Qu.  Max. 
                            #  0.000   7.183    9.252   8.570   10.433  17.132 
        
                            # Interpretation
                            # Range = 17.132 − 0 ≈ 17
                                  # means most highly expressed gene is ~17 log2 units higher than the lowest expressed gene
                                  # in log2 scale: +1 = double expression, +17 = 17 doublings  
                                  # so highest gene expressed 2^17 ≈ 131,072 more than lowest
                                  # each +1 unit increase in log2(expression) multiplies the odds by 1.24
                                  # so OR = 1.24^17 ≈ 32.1
                                      # Interpretation: the highest expressed genes have ~32× higher odds 
                                      # of being immunoreactive than the lowest expressed genes (OR= 1.24, p = 1.28e-15)
                                                  # Each doubling of expression = +24% odds
                                                  # Double expression 17 times → the effect compounds
                                                  # Result ≈ 32× higher odds
                    
                # Interpretation: expression correlates with immunoreactivity (p value = 1.28e-15, OR = 32)
                          