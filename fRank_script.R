rm(list = ls())
lapply(c("dplyr", "tibble"), library, character.only=T)

# Set your working directory
setwd("~/Desktop/Stuff/AMeyer_Lab/Experiment_results/5seq/TSSPredIremOutput/")

# Read your TSSpreditor file output.
BigData <- read.csv("ChlData_Cleaned.csv", header = T)
BigData <- BigData %>% 
  add_column(TSSUTR ="", Locus ="", termUTR ="", .after = "SuperPos")
top <- BigData %>% filter(SuperStrand == "+")
comp <- BigData %>% filter(SuperStrand == "-")

# Read your annotation file
genes <- read.delim("NC_003028.v3.17.ncrna.genes", header = F)
genes <- genes[complete.cases(genes),]
genes.top <- genes %>% filter(V4 == "+")
colnames(genes.top) <- c("genome", "from", "to", "strand", "name", "old.name", "new.name", "WP", "rfam", "other" )
genes.comp <- genes %>% filter(V4 == "-")
colnames(genes.comp) <- c("genome", "to", "from", "strand", "name", "old.name", "new.name", "WP", "rfam", "other" )

# List of columns to convert to numeric
cols_to_numeric <- c("from", "to")
# Loop through the selected columns and convert them to numeric
genes.comp[cols_to_numeric] <- lapply(genes.comp[cols_to_numeric], as.numeric)
genes.top[cols_to_numeric] <- lapply(genes.top[cols_to_numeric], as.numeric)
comp$SuperPos <- as.numeric(comp$SuperPos)
top$SuperPos <- as.numeric(top$SuperPos)
########################################
# Annotate BigData Top and Comp Strands
########################################

for (i in 1:nrow(top)) {
  for (k in 1:nrow(genes.top)) {
    if(!is.na(genes.top$from[k]) && k < nrow(genes.top)){
    ######       First Annotate 5'UTR     ######
      # Check if there's at least 500 bp between genes
      if (genes.top$from[k + 1] - genes.top$to[k] >= 500) {
        #Then check following conditions are met
        if (
          top$SuperPos[i] > genes.top$to[k] - 50  &&
          top$SuperPos[i] < genes.top$from[k + 1] + 10 &&
          top$SuperPos[i] > genes.top$from[k + 1] - 550
        ) {
          top$TSSUTR[i] <- genes.top$old.name[k + 1]
          break  # Exit the inner loop once a 5' UTR annotation is found
        }
      } else {
        #Annotate everything btw genes as 5' UTR if the distance is <500 bp
        if (
          top$SuperPos[i] > genes.top$to[k] &&
          top$SuperPos[i] < genes.top$from[k + 1] + 10
        ) {
          top$TSSUTR[i] <- genes.top$old.name[k + 1]
          break  # Exit the inner loop once a 5' UTR annotation is found
        }
      }
    }
    # Handle the case for positions before the first gene
    if (k == 1 && top$SuperPos[i] < genes.top$from[k]) {
      top$TSSUTR[i] <- genes.top$old.name[k]
      break
    }
    if (!is.na(genes.top$from[k]) && !is.na(genes.top$to[k]) && (top$TSSUTR[i] == "") && (top$termUTR[i] == "")) {
      ######    Annotate genic regions    ######
      if (
        top$SuperPos[i] >= genes.top$from[k] &&
        top$SuperPos[i] <= genes.top$to[k]
      ) {
        top$Locus[i] <- genes.top$old.name[k]
        break  # Exit the inner loop once a genic annotation is found
      }
    }
  }
}


for (i in 1:nrow(comp)) {
  for (k in 1:nrow(genes.comp)) {
    ######       First Annotate 5'UTR     ######
    # Check if there's at least 500 bp between genes
    if (!is.na(genes.comp$to[k]) && !is.na(genes.comp$from[k + 1])){
      # Check if there's at least 300 bp between genes
      if (genes.comp$to[k + 1] - genes.comp$from[k] >= 500) {
        #Then check following conditions are met
        if (
          comp$SuperPos[i] < genes.comp$to[k+1] + 50  &&
          comp$SuperPos[i] > genes.comp$from[k] - 10 &&
          comp$SuperPos[i] < genes.comp$from[k] + 550
        ) { 
          comp$TSSUTR[i] <- genes.comp$old.name[k]
          break  # Exit the inner loop once a 5' UTR annotation is found
        }
      }  else {
        # #Annotate everything btw genes as 5' UTR if the distance is <500 bp
        if (
          comp$SuperPos[i] > genes.comp$from[k] - 10 &&
          comp$SuperPos[i] < genes.comp$to[k+1]
        ) {
          comp$TSSUTR[i] <- genes.comp$old.name[k]
        }
      }  
    }
    if (!is.na(genes.comp$from[k]) && !is.na(genes.comp$to[k]) && (comp$TSSUTR[i] == "") && (comp$termUTR[i] == "")) {
      ######    Annotate genic regions    ######
      if (
        comp$SuperPos[i] <= genes.comp$from[k] &&
        comp$SuperPos[i] >= genes.comp$to[k]
      ) {
        comp$Locus[i] <- genes.comp$old.name[k]
        break  # Exit the inner loop once a genic annotation is found
      }
    }
  }
}

#Check for antisense in top strand
for(i in 1:nrow(top)){
  for(k in 1:nrow(genes.comp)){
    if(
      top$SuperPos[i] <= genes.comp$from[k] && 
      top$SuperPos[i] >= genes.comp$to[k] 
      ) {
      top$termUTR[i] <- "antisense"
    }
  }
}
#Check for antisense in comp strand
for(i in 1:nrow(comp)){
  for(k in 1:nrow(genes.top)){
    if(
      comp$SuperPos[i] >= genes.top$from[k] && 
      comp$SuperPos[i] <= genes.top$to[k]  
    ) {
      comp$termUTR[i] <- "antisense"
    }
  }
}

# Remove extra columns that came from TSSpredator
BigDatagenes <- rbind(comp, top) %>%
  arrange(SuperPos) %>%
  filter(!Genome == "NDCt0")
rm(comp,top)
BigDatagenes <- BigDatagenes[,-c(15:30,32:33)] 
write.csv(BigDatagenes, file = "ChlData_nofRank.csv", row.names = F, col.names = T)
# Start preparing data for rank analysis -- Remove rows with no gene information 
BigDatagenes <- BigDatagenes %>% 
  filter(!Locus == "" | !TSSUTR == "" | !termUTR == "")

#Paste all locus information into one column for easy data processing
BigDatagenes$Locus <- paste0(BigDatagenes$TSSUTR, BigDatagenes$termUTR, BigDatagenes$Locus) 
BigDatagenes <- BigDatagenes[,-c(2,4)]

#Separate samples based on conditions & remove empty rows
list2env(setNames(lapply(c("NDCt60", "Chl1q_t60", "Chl3q_t60"), function(x)
  BigDatagenes %>% filter(Genome == x)),
  c("BigDataNDC", "BigData1Q", "BigData3Q")), envir = .GlobalEnv)

list2env(set_names(lapply(
  list(BigDataNDC, BigData1Q, BigData3Q),
  function(x) x[complete.cases(x),]),
  c("BigDataNDC", "BigData1Q", "BigData3Q"),), envi = .GlobalEnv)

###################
# Ranking Function
###################
# Get row with max value for a specific primary column (i.e. enrichmentFactor), 
# then max of secondary column (i.e. stepHeight) if there are ties
fRank <- function(df, group_col, primary_col, secondary_col){
  df %>%
    group_by(!!sym(group_col)) %>%
    # First, filter to keep only rows with max value of primary_col
    filter(!!sym(primary_col) == max(!!sym(primary_col))) %>%
    # Then, among these, select the row with max value of secondary_col
    slice(which.max(!!sym(secondary_col))) %>%
    ungroup()
}
################
#  NDC Primary 
################
# Process for Enrichment Factor, then Step Height
BigDataNDCprimary <- BigDataNDC %>%
  fRank("Locus", "enrichmentFactor", "stepHeight") %>%
  arrange(SuperPos)
# Reset row names if needed
rownames(BigDataNDCprimary) <- NULL

###############
# 1Q primary
###############
# Process for Enrichment Factor, then Step Height
BigData1Qprimary <- BigData1Q %>%
  fRank("Locus", "enrichmentFactor", "stepHeight") %>%
  arrange(SuperPos)
# Reset row names if needed
rownames(BigData1Qprimary) <- NULL
###############
# 3Q primary
###############
# Process for Enrichment Factor, then Step Height
BigData3Qprimary <- BigData3Q %>%
  fRank("Locus", "enrichmentFactor", "stepHeight") %>%
  arrange(SuperPos)
# Reset row names if needed
rownames(BigData3Qprimary) <- NULL

###########################
# Annotate BigData primary
###########################
BigDataprimary <- bind_rows(BigDataNDCprimary, BigData1Qprimary, BigData3Qprimary) %>% 
  arrange(SuperPos)
BigDataprimary <- BigDataprimary[,-c(2)]
BigDataprimary <- BigDataprimary %>%
  add_column(TSSUTR ="", Locus ="", termUTR ="", .after = "SuperPos")

top <- BigDataprimary %>% filter(SuperStrand == "+")
comp <- BigDataprimary %>% filter(SuperStrand == "-")

for (i in 1:nrow(top)) {
  for (k in 1:nrow(genes.top)) {
    if(!is.na(genes.top$from[k]) && k < nrow(genes.top)){
      ######       First Annotate 5'UTR     ######
      # Check if there's at least 500 bp between genes
      if (genes.top$from[k + 1] - genes.top$to[k] >= 500) {
        #Then check following conditions are met
        if (
          top$SuperPos[i] > genes.top$to[k] - 50  &&
          top$SuperPos[i] < genes.top$from[k + 1] + 10 &&
          top$SuperPos[i] > genes.top$from[k + 1] - 550
        ) {
          top$TSSUTR[i] <- genes.top$old.name[k + 1]
          break  # Exit the inner loop once a 5' UTR annotation is found
        }
      } else {
        #Annotate everything btw genes as 5' UTR if the distance is <500 bp
        if (
          top$SuperPos[i] > genes.top$to[k] &&
          top$SuperPos[i] < genes.top$from[k + 1] + 10
        ) {
          top$TSSUTR[i] <- genes.top$old.name[k + 1]
          break  # Exit the inner loop once a 5' UTR annotation is found
        }
      }
    }
    # Handle the case for positions before the first gene
    if (k == 1 && top$SuperPos[i] < genes.top$from[k]) {
      top$TSSUTR[i] <- genes.top$old.name[k]
      break
    }
    if (!is.na(genes.top$from[k]) && !is.na(genes.top$to[k]) && (top$TSSUTR[i] == "") && (top$termUTR[i] == "")) {
      ######    Annotate genic regions    ######
      if (
        top$SuperPos[i] >= genes.top$from[k] &&
        top$SuperPos[i] <= genes.top$to[k]
      ) {
        top$Locus[i] <- genes.top$old.name[k]
        break  # Exit the inner loop once a genic annotation is found
      }
    }
  }
}

for (i in 1:nrow(comp)) {
  for (k in 1:nrow(genes.comp)) {
    ######       First Annotate 5'UTR     ######
    # Check if there's at least 500 bp between genes
    if (!is.na(genes.comp$to[k]) && !is.na(genes.comp$from[k + 1])){
      # Check if there's at least 300 bp between genes
      if (genes.comp$to[k + 1] - genes.comp$from[k] >= 500) {
        #Then check following conditions are met
        if (
          comp$SuperPos[i] < genes.comp$to[k+1] + 50  &&
          comp$SuperPos[i] > genes.comp$from[k] - 10 &&
          comp$SuperPos[i] < genes.comp$from[k] + 550
        ) { 
          comp$TSSUTR[i] <- genes.comp$old.name[k]
          break  # Exit the inner loop once a 5' UTR annotation is found
        }
      }  else {
        # #Annotate everything btw genes as 5' UTR if the distance is <500 bp
        if (
          comp$SuperPos[i] > genes.comp$from[k] - 10 &&
          comp$SuperPos[i] < genes.comp$to[k+1]
        ) {
          comp$TSSUTR[i] <- genes.comp$old.name[k]
        }
      }  
    }
    if (!is.na(genes.comp$from[k]) && !is.na(genes.comp$to[k]) && (comp$TSSUTR[i] == "") && (comp$termUTR[i] == "")) {
      ######    Annotate genic regions    ######
      if (
        comp$SuperPos[i] <= genes.comp$from[k] &&
        comp$SuperPos[i] >= genes.comp$to[k]
      ) {
        comp$Locus[i] <- genes.comp$old.name[k]
        break  # Exit the inner loop once a genic annotation is found
      }
    }
  }
}

#Check for antisense in top strand
for(i in 1:nrow(top)){
  for(k in 1:nrow(genes.comp)){
    if(
      top$SuperPos[i] <= genes.comp$from[k] && 
      top$SuperPos[i] >= genes.comp$to[k] 
    ) {
      top$termUTR[i] <- "antisense"
    }
  }
}
#Check for antisense in comp strand
for(i in 1:nrow(comp)){
  for(k in 1:nrow(genes.top)){
    if(
      comp$SuperPos[i] >= genes.top$from[k] && 
      comp$SuperPos[i] <= genes.top$to[k]  
    ) {
      comp$termUTR[i] <- "antisense"
    }
  }
}

BigDataprimary <- rbind(top,comp)
nrow(BigDataprimary %>% filter(Genome == "NDCt60" & TSSUTR != "" & Locus == ""))
nrow(BigDataprimary %>% filter(Genome == "NDCt60" & TSSUTR == "" & Locus != ""))
nrow(BigDataprimary %>% filter(Genome == "Chl1q_t60" & TSSUTR != "" & Locus == ""))
nrow(BigDataprimary %>% filter(Genome == "Chl1q_t60" & TSSUTR == "" & Locus != ""))
nrow(BigDataprimary %>% filter(Genome == "Chl3q_t60" & TSSUTR != "" & Locus == ""))
nrow(BigDataprimary %>% filter(Genome == "Chl3q_t60" & TSSUTR == "" & Locus != ""))
BigDataprimary$rank <- "primary"

                ########################################################
                #########      SECONDARY TSS ANALYSIS         ##########
                ########################################################
                
#Remove Primary positions from BigData  & extra columns that came from TSSpredator 
# BigDataNP: BigDataNotPrimary
BigDataNP <- BigData[!BigData$SuperPos %in% BigDataprimary$SuperPos,]

top <- BigDataNP %>% filter(SuperStrand == "+")
comp <- BigDataNP %>% filter(SuperStrand == "-")
############################
# Annotate BigDataNP
############################
for (i in 1:nrow(top)) {
  for (k in 1:nrow(genes.top)) {
    if(!is.na(genes.top$from[k]) && k < nrow(genes.top)){
      ######       First Annotate 5'UTR     ######
      # Check if there's at least 500 bp between genes
      if (genes.top$from[k + 1] - genes.top$to[k] >= 500) {
        #Then check following conditions are met
        if (
          top$SuperPos[i] > genes.top$to[k] - 50  &&
          top$SuperPos[i] < genes.top$from[k + 1] + 10 &&
          top$SuperPos[i] > genes.top$from[k + 1] - 550
        ) {
          top$TSSUTR[i] <- genes.top$old.name[k + 1]
          break  # Exit the inner loop once a 5' UTR annotation is found
        }
      } else {
        #Annotate everything btw genes as 5' UTR if the distance is <500 bp
        if (
          top$SuperPos[i] > genes.top$to[k] &&
          top$SuperPos[i] < genes.top$from[k + 1] + 10
        ) {
          top$TSSUTR[i] <- genes.top$old.name[k + 1]
          break  # Exit the inner loop once a 5' UTR annotation is found
        }
      }
    }
    # Handle the case for positions before the first gene
    if (k == 1 && top$SuperPos[i] < genes.top$from[k]) {
      top$TSSUTR[i] <- genes.top$old.name[k]
      break
    }
    if (!is.na(genes.top$from[k]) && !is.na(genes.top$to[k]) && (top$TSSUTR[i] == "") && (top$termUTR[i] == "")) {
      ######    Annotate genic regions    ######
      if (
        top$SuperPos[i] >= genes.top$from[k] &&
        top$SuperPos[i] <= genes.top$to[k]
      ) {
        top$Locus[i] <- genes.top$old.name[k]
        break  # Exit the inner loop once a genic annotation is found
      }
    }
  }
}
for (i in 1:nrow(comp)) {
  for (k in 1:nrow(genes.comp)) {
    ######       First Annotate 5'UTR     ######
    # Check if there's at least 500 bp between genes
    if (!is.na(genes.comp$to[k]) && !is.na(genes.comp$from[k + 1])){
      # Check if there's at least 300 bp between genes
      if (genes.comp$to[k + 1] - genes.comp$from[k] >= 500) {
        #Then check following conditions are met
        if (
          comp$SuperPos[i] < genes.comp$to[k+1] + 50  &&
          comp$SuperPos[i] > genes.comp$from[k] - 10 &&
          comp$SuperPos[i] < genes.comp$from[k] + 550
        ) { 
          comp$TSSUTR[i] <- genes.comp$old.name[k]
          break  # Exit the inner loop once a 5' UTR annotation is found
        }
      }  else {
        # #Annotate everything btw genes as 5' UTR if the distance is <500 bp
        if (
          comp$SuperPos[i] > genes.comp$from[k] - 10 &&
          comp$SuperPos[i] < genes.comp$to[k+1]
        ) {
          comp$TSSUTR[i] <- genes.comp$old.name[k]
        }
      }  
    }
    if (!is.na(genes.comp$from[k]) && !is.na(genes.comp$to[k]) && (comp$TSSUTR[i] == "") && (comp$termUTR[i] == "")) {
      ######    Annotate genic regions    ######
      if (
        comp$SuperPos[i] <= genes.comp$from[k] &&
        comp$SuperPos[i] >= genes.comp$to[k]
      ) {
        comp$Locus[i] <- genes.comp$old.name[k]
        break  # Exit the inner loop once a genic annotation is found
      }
    }
  }
}
#Check for antisense in top strand
for(i in 1:nrow(top)){
  for(k in 1:nrow(genes.comp)){
    if(
      top$SuperPos[i] <= genes.comp$from[k] && 
      top$SuperPos[i] >= genes.comp$to[k] 
    ) {
      top$termUTR[i] <- "antisense"
    }
  }
}
#Check for antisense in comp strand
for(i in 1:nrow(comp)){
  for(k in 1:nrow(genes.top)){
    if(
      comp$SuperPos[i] >= genes.top$from[k] && 
      comp$SuperPos[i] <= genes.top$to[k]  
    ) {
      comp$termUTR[i] <- "antisense"
    }
  }
}
# Removing extra columns came from TSSpredator
BigDataNPgenes <- rbind(comp, top) %>%
  arrange(SuperPos) %>%
  filter(!Genome == "NDCt0")
rm(comp,top)
BigDataNPgenes <- BigDataNPgenes[,-c(15:30,32:33)]

# Start preparing data for secondary rank analysis -- Remove rows with no gene information 
BigDataNPgenes <- BigDataNPgenes %>% 
  filter(!Locus == "" | !TSSUTR == "" | !termUTR == "")

# Paste all locus information into one column for easy data processing
BigDataNPgenes$Locus <- paste0(BigDataNPgenes$TSSUTR, BigDataNPgenes$termUTR, BigDataNPgenes$Locus) 
BigDataNPgenes <- BigDataNPgenes[,-c(2,4)]

# Separate samples based on conditions & remove empty rows
list2env(set_names(lapply(c("NDCt60", "Chl1q_t60", "Chl3q_t60"), function(x)
  BigDataNPgenes %>% filter(Genome == x)), 
  c("SecNDC", "Sec1Q", "Sec3Q")), envir = .GlobalEnv )

list2env(set_names(lapply(list(SecNDC, Sec1Q, Sec3Q), function(y) 
  y[complete.cases(y),]), 
  c("SecNDC", "Sec1Q", "Sec3Q"),), envir = .GlobalEnv)

################
#  NDC Secondary 
################
# Process for Enrichment Factor, then Step Height
BigDataNDCSec <- SecNDC %>%
  fRank("Locus", "enrichmentFactor", "stepHeight") %>%
  arrange(SuperPos)
# Reset row names if needed
rownames(BigDataNDCSec) <- NULL

###############
# 1Q secondary
###############
# Process for Enrichment Factor, then Step Height
BigData1Qsec <- Sec1Q %>%
  fRank("Locus", "enrichmentFactor", "stepHeight") %>%
  arrange(SuperPos)
# Reset row names if needed
rownames(BigData1Qsec) <- NULL

###############
# 3Q secondary
###############
# Process for Enrichment Factor, then Step Height
BigData3Qsec <- Sec3Q %>%
  fRank("Locus", "enrichmentFactor", "stepHeight") %>%
  arrange(SuperPos)
# Reset row names if needed
rownames(BigData3Qsec) <- NULL

#############################
# Annotate BigData secondary
#############################
BigDatasecondary <- bind_rows(BigDataNDCSec, BigData1Qsec, BigData3Qsec) %>% 
  arrange(SuperPos)
BigDatasecondary <- BigDatasecondary[,-c(2)]
BigDatasecondary <- BigDatasecondary %>%
  add_column(TSSUTR ="", Locus ="", termUTR ="", .after = "SuperPos")

top <- BigDatasecondary %>% filter(SuperStrand == "+")
comp <- BigDatasecondary %>% filter(SuperStrand == "-")

for (i in 1:nrow(top)) {
  for (k in 1:nrow(genes.top)) {
    if(!is.na(genes.top$from[k]) && k < nrow(genes.top)){
      ######       First Annotate 5'UTR     ######
      # Check if there's at least 500 bp between genes
      if (genes.top$from[k + 1] - genes.top$to[k] >= 500) {
        #Then check following conditions are met
        if (
          top$SuperPos[i] > genes.top$to[k] - 50  &&
          top$SuperPos[i] < genes.top$from[k + 1] + 10 &&
          top$SuperPos[i] > genes.top$from[k + 1] - 550
        ) {
          top$TSSUTR[i] <- genes.top$old.name[k + 1]
          break  # Exit the inner loop once a 5' UTR annotation is found
        }
      } else {
        #Annotate everything btw genes as 5' UTR if the distance is <500 bp
        if (
          top$SuperPos[i] > genes.top$to[k] &&
          top$SuperPos[i] < genes.top$from[k + 1] + 10
        ) {
          top$TSSUTR[i] <- genes.top$old.name[k + 1]
          break  # Exit the inner loop once a 5' UTR annotation is found
        }
      }
    }
    # Handle the case for positions before the first gene
    if (k == 1 && top$SuperPos[i] < genes.top$from[k]) {
      top$TSSUTR[i] <- genes.top$old.name[k]
      break
    }
    if (!is.na(genes.top$from[k]) && !is.na(genes.top$to[k]) && (top$TSSUTR[i] == "") && (top$termUTR[i] == "")) {
      ######    Annotate genic regions    ######
      if (
        top$SuperPos[i] >= genes.top$from[k] &&
        top$SuperPos[i] <= genes.top$to[k]
      ) {
        top$Locus[i] <- genes.top$old.name[k]
        break  # Exit the inner loop once a genic annotation is found
      }
    }
  }
}

for (i in 1:nrow(comp)) {
  for (k in 1:nrow(genes.comp)) {
    ######       First Annotate 5'UTR     ######
    # Check if there's at least 500 bp between genes
    if (!is.na(genes.comp$to[k]) && !is.na(genes.comp$from[k + 1])){
      # Check if there's at least 300 bp between genes
      if (genes.comp$to[k + 1] - genes.comp$from[k] >= 500) {
        #Then check following conditions are met
        if (
          comp$SuperPos[i] < genes.comp$to[k+1] + 50  &&
          comp$SuperPos[i] > genes.comp$from[k] - 10 &&
          comp$SuperPos[i] < genes.comp$from[k] + 550
        ) { 
          comp$TSSUTR[i] <- genes.comp$old.name[k]
          break  # Exit the inner loop once a 5' UTR annotation is found
        }
      }  else {
        # #Annotate everything btw genes as 5' UTR if the distance is <500 bp
        if (
          comp$SuperPos[i] > genes.comp$from[k] - 10 &&
          comp$SuperPos[i] < genes.comp$to[k+1]
        ) {
          comp$TSSUTR[i] <- genes.comp$old.name[k]
        }
      }  
    }
    if (!is.na(genes.comp$from[k]) && !is.na(genes.comp$to[k]) && (comp$TSSUTR[i] == "") && (comp$termUTR[i] == "")) {
      ######    Annotate genic regions    ######
      if (
        comp$SuperPos[i] <= genes.comp$from[k] &&
        comp$SuperPos[i] >= genes.comp$to[k]
      ) {
        comp$Locus[i] <- genes.comp$old.name[k]
        break  # Exit the inner loop once a genic annotation is found
      }
    }
  }
}
#Check for antisense in top strand
for(i in 1:nrow(top)){
  for(k in 1:nrow(genes.comp)){
    if(
      top$SuperPos[i] <= genes.comp$from[k] && 
      top$SuperPos[i] >= genes.comp$to[k] 
    ) {
      top$termUTR[i] <- "antisense"
    }
  }
}
#Check for antisense in comp strand
for(i in 1:nrow(comp)){
  for(k in 1:nrow(genes.top)){
    if(
      comp$SuperPos[i] >= genes.top$from[k] && 
      comp$SuperPos[i] <= genes.top$to[k]  
    ) {
      comp$termUTR[i] <- "antisense"
    }
  }
}
BigDatasecondary <- rbind(top,comp)
nrow(BigDatasecondary %>% filter(Genome == "NDCt60" & TSSUTR != "" & Locus == ""))
nrow(BigDatasecondary %>% filter(Genome == "NDCt60" & TSSUTR == "" & Locus != ""))
nrow(BigDatasecondary %>% filter(Genome == "Chl1q_t60" & TSSUTR != "" & Locus == ""))
nrow(BigDatasecondary %>% filter(Genome == "Chl1q_t60" & TSSUTR == "" & Locus != ""))
nrow(BigDatasecondary %>% filter(Genome == "Chl3q_t60" & TSSUTR != "" & Locus == ""))
nrow(BigDatasecondary %>% filter(Genome == "Chl3q_t60" & TSSUTR == "" & Locus != ""))
BigDatasecondary$rank <- "secondary"

BigDataRanked <- rbind(BigDataprimary, BigDatasecondary) %>% 
  arrange(SuperPos)

write.csv(BigDataRanked, file = "ChlData_fRanked.csv", row.names = F, col.names = T)
