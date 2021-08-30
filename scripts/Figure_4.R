library(plotly)
library(ggtree)
require(cowplot)

#Check to see if plasmid coverage data object exists - if not, run the script which generates it.
if(!exists("plasmid_coverage_table")){
        source("scripts/ST95_all_process.R")
}

`%nin%` <- Negate(`%in%`)

library(highcharter)
library(htmlwidgets)
library(dplyr)

df <- df_small4# %>% select(IncF_RST, ColV, plasmid_map, Revised_Source_Niche)

df$working_name <- rownames(df_small4)

df$plasmid_map <- gsub(".*_70-79.*","Other",df$plasmid_map)

df$plasmid_map <- gsub(".*_60-69.*","Other",df$plasmid_map)

df$plasmid_map <- gsub(".*_50-59.*","Other",df$plasmid_map)

df$plasmid_map <- gsub(".*_80-89.*","",df$plasmid_map)

#df$plasmid_map <- gsub("_80-89","",df$plasmid_map)
df$plasmid_map <- gsub("_90-100","",df$plasmid_map)

df$IncF_RST_all <- df$IncF_RST

pMLST_counts <- df %>% group_by(IncF_RST) %>% summarise(counts = n()) %>% arrange(desc(counts)) %>%  filter(counts > 15) %>% as.data.frame()

df$IncF_RST[df$IncF_RST %nin% pMLST_counts$IncF_RST] <- "Other"

df$Pathogen <- gsub("_pathogen_status_","",df$Pathogen)

df$Pathogen[is.na(df$Pathogen)] <- "Unknown"

df$Revised_Source_Niche[is.na(df$Revised_Source_Niche)] <- "Unknown"

df$Revised_Source_Niche <- gsub("Source_","",df$Revised_Source_Niche)

df$plasmid_map[df$IncF_RST == "F-:A-:B-"] <- "F Null"

df$IncF_RST <- gsub("C4", "F18", df$IncF_RST)

df$IncF_RST <- gsub("A-:", "", df$IncF_RST)

df$plasmid_map[is.na(df$plasmid_map)] <- "Other"
df$plasmid_map[df$plasmid_map == ""] <- "Other"

df$Pathogen[df$Pathogen == "Flora"] <- "Gastrointestinal "

df$Continent[is.na(df$Continent)] <- "Unknown "

#df$plasmid_map[df$plasmid_map == "Other" & df$ColV == 1] <- "Other (ColV Positive)"

#df$plasmid_map[df$plasmid_map == "Other" & df$ColV == 0] <- "Other (ColV Negative)"

#df$plasmid_map[df$IncF_RST == "NA"] <- "Other (ColV Negative)"

df$ColV[df$ColV == 0] <- "ColV Negative"

df$ColV[df$ColV == 1] <- "ColV Positive"

df$MDR <- gsub(" .*","",df$class_res_counts)

df$MDR <- gsub(".*0","",df$MDR)

df$plasmid_map <- gsub("pSF_088_nores","pSF-088*",df$plasmid_map)
df$plasmid_map <- gsub("pEC244_2","pEC244_2*",df$plasmid_map)
df$plasmid_map <- gsub("pACN001_B","pACN001_B*",df$plasmid_map)
df$plasmid_map <- gsub("pUTI89","pUTI89*",df$plasmid_map)
df$plasmid_map <- gsub("pU1_F51_B10","pU1*",df$plasmid_map)
df$plasmid_map <- gsub("pBCE049_1","pBCE049_1*",df$plasmid_map)

df$MDR <- as.numeric(df$MDR)

df$MDR[is.na(df$MDR)] <- 0

df$MDR[df$MDR < 3] <- 0
df$MDR[df$MDR >= 3] <- 1

df$MDR[df$MDR == 0] <- "Not MDR"
df$MDR[df$MDR == 1] <- "MDR"

df <- apply(df, 2, FUN = function(x) paste0(x, " (",as.data.frame(table(x)[x])[,2],"/668)"))

df <- as.data.frame(df)

df$name

#AMR vs ColV vs Plasmid type
AMR_ColV_plas <- df %>% select(class_res_counts, ColV, plasmid_map)

#Source vs ColV vs plasmid type
pMLST_ColV_plasmap <- df %>% select(plasmid_map, Revised_Source_Niche, ColV)

#Source vs ColV vs plasmid type
pMLST_ColV_mdr <- df %>% select(ESBL_ress, MDR, intI1)

plas_count_plasmap_MDR <- df %>% select(IncF_RST, plasmid_map, ColV)

hchart(data_to_sankey(plas_count_plasmap_MDR), "sankey", name = "IncF_RST to ColV and plasmid_map")  %>%
        hc_plotOptions(series = list(dataLabels = list(style = list(fontSize = "20px"))))
