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

df_clean <- apply(X = df, MARGIN = 2, FUN = function(x)gsub(" \\(.*", "", x))

df_clean <- df_clean %>% as.data.frame()

df_clean$class_res_counts <- gsub("classes_res_0([0-9])","\\1", df_clean$class_res_counts)
df_clean$class_res_counts <- as.numeric(df_clean$class_res_counts)

df_clean$plas_counts <- gsub("count_plas_","", df_clean$plas_counts)
df_clean$plas_counts <- as.numeric(df_clean$plas_counts)

df_clean$MDR <- gsub("Not MDR","0", df_clean$MDR)
df_clean$MDR <- gsub("MDR","1", df_clean$MDR)
df_clean$MDR <- as.numeric(df_clean$MDR)

df_clean$plasmid_map_extended <- df_clean$plasmid_map

df_clean$plasmid_map_extended[df_clean$ColV == "ColV Positive"] <- paste( df_clean$plasmid_map_extended[df_clean$ColV == "ColV Positive"], "(ColV Positive)")
df_clean$plasmid_map_extended[df_clean$ColV == "ColV Negative"] <- paste( df_clean$plasmid_map_extended[df_clean$ColV == "ColV Negative"], "(ColV Negative)")
                                                                          

df_clean %>% group_by(ColV) %>% summarise(mean_res = mean(class_res_counts), mean_plas = mean(plas_counts), mean_mdr = mean(MDR)) %>% View()

all_data <- metadata %>% select(-"Revised_Source_Niche", -"Continent", -"Pathogen", -"ColV", -"IncF_RST") %>% left_join(df_clean)

all_data <- Metadata %>% select(working_name, colnames(Metadata)[colnames(Metadata) %nin% colnames(all_data)]) %>% left_join(all_data)

#read in fastbaps
fastbaps_data <- read.csv(file = "fastbaps.csv")

colnames(fastbaps_data) <- c("working_name", "Level_1", "Level_2")

#Add leading zeroes
fastbaps_data$Level_1 <- str_pad(string = fastbaps_data$Level_1, side = "left", pad = "0", width = 2)
fastbaps_data$Level_2 <- str_pad(string = fastbaps_data$Level_2, side = "left", pad = "0", width = 2)

fastbaps_data$Level_1 <- gsub("$","_fastbaps", fastbaps_data$Level_1)
fastbaps_data$Level_2 <- gsub("$","_fastbaps", fastbaps_data$Level_2)


custom_labels <- cbind(c('04_fastbaps','09_fastbaps','07_fastbaps','08_fastbaps','10_fastbaps','05_fastbaps','02_fastbaps','01_fastbaps','03_fastbaps','06_fastbaps'),LETTERS[1:length(unique(fastbaps_data$Level_1))])

custom_labels <- as.data.frame(custom_labels) %>% rename("Level_1" = "V1", "clade_label" = "V2")

fastbaps_data <- left_join(fastbaps_data, custom_labels)

all_data <- left_join(all_data, fastbaps_data)

#Create a function to generate a spectrum of colours for each plasmid type
plasmid_ref_names <- all_data$plasmid_map %>% unique() %>% sort()

plasmid_ref_cols <- c("white", #F Null
                      "white", #Other
                      "red", # pACN
                      #"yellow", # pAPEC 
                      "black", # pBCE
                      "blue", # pEC244
                      "orange", # pSF
                      "purple", # pU1
                      "green" # pUTI89
)

names(plasmid_ref_cols) <- plasmid_ref_names

#### Plot count of AMR classes vs ColV status ####
g <- ggplot(all_data, mapping = aes(x = reorder(ColV, desc(ColV)), y = class_res_counts))


ColV_status_class_res_counts <- g + geom_boxplot() + 
        labs(title="", 
             x="ColV Status",
             y="Count of classes of AMR genes") +
        theme(axis.text.x = element_text(angle=90))


#### Plot count of non incF plasmids vs ColV status ####

g <- ggplot(all_data, mapping = aes(x = reorder(ColV, desc(ColV)), y = plas_counts))


ColV_status_plasmid_counts <- g + geom_boxplot() + 
        labs(title="", 
             x="ColV Status",
             y="Count of plasmids (Non-IncF)") +
        theme(axis.text.x = element_text(angle=90))



#### Plot count of non incF plasmids vs plasmid map ####

g <- ggplot(all_data, mapping = aes(x = reorder(plasmid_map, desc(plas_counts)), y = plas_counts))


plasmid_map_status_plasmid_counts <- g + geom_boxplot() + 
        labs(title="", 
             x="Plasmid",
             y="Count of plasmids (Non-IncF)") +
        theme(axis.text.x = element_text(angle=90))



ggplotly(plasmid_map_status_plasmid_counts)


#### Plot count of AMR classes vs plasmid map ####

g <- ggplot(all_data, mapping = aes(x = reorder(plasmid_map, desc(class_res_counts)), y = class_res_counts))


plasmid_map_status_class_res_counts <- g + geom_boxplot() + 
        labs(title="", 
             x="Plasmid",
             y="Count of classes of AMR genes") +
        theme(axis.text.x = element_text(angle=90))



#### Plot count of AMR classes vs fastbaps ####

g <- ggplot(all_data, mapping = aes(x = reorder(Level_1, desc(class_res_counts)), y = class_res_counts))


plasmid_map_status_class_res_counts <- g + geom_boxplot() + 
        labs(title="", 
             x="Fastbaps group",
             y="Count of classes of AMR genes") +
        theme(axis.text.x = element_text(angle=90))

#### Plot count of AMR classes vs fastbaps ####

g <- ggplot(all_data, mapping = aes(x = reorder(clade_label, desc(class_res_counts)), y = class_res_counts))


test <- g + geom_boxplot() + 
        labs(title="", 
             x="Fastbaps group",
             y="Count of classes of AMR genes")



ggplotly(plasmid_map_status_class_res_counts)

all_data$fimH_type[is.na(all_data$fimH_type)] <- "fimH-"


#### Chi-squared tests ####

#MDR vs ColV Pos/Neg

plasmid_vs_amr <- table(all_data$ColV, all_data$class_res_counts > 2)

colnames(plasmid_vs_amr) <- c("non-MDR", "MDR")
rownames(plasmid_vs_amr) <- c("ColV Neg", "ColV Positive")

chisq.test(plasmid_vs_amr, correct = FALSE)

plasmid_vs_amr




#AMR vs ColV Pos/Neg

plasmid_vs_amr <- table(all_data$ColV, all_data$class_res_counts > 0)

colnames(plasmid_vs_amr) <- c("non-AMR", "AMR")
rownames(plasmid_vs_amr) <- c("ColV Neg", "ColV Positive")

chisq.test(plasmid_vs_amr, correct = FALSE)


plasmid_vs_amr

proportions(plasmid_vs_amr, margin = 1)


#AMR vs pUTI89* Pos/Neg
plasmid_vs_amr <- table(all_data$plasmid_map == "pUTI89*", all_data$class_res_counts > 0)

colnames(plasmid_vs_amr) <- c("non-AMR", "AMR")
rownames(plasmid_vs_amr) <- c("non-pUTI89", "pUTI89*")

chisq.test(plasmid_vs_amr, correct = FALSE)

plasmid_vs_amr

proportions(plasmid_vs_amr, margin = 1)

#MDR vs pUTI89* Pos/Neg
plasmid_vs_amr <- table(all_data$plasmid_map == "pUTI89*", all_data$class_res_counts > 2)

colnames(plasmid_vs_amr) <- c("non-MDR", "MDR")
rownames(plasmid_vs_amr) <- c("non-pUTI89", "pUTI89*")

chisq.test(plasmid_vs_amr, correct = FALSE)

plasmid_vs_amr

proportions(plasmid_vs_amr, margin = 1)


#AMR vs pSF-088* Pos/Neg
plasmid_vs_amr <- table(all_data$plasmid_map == "pSF-088*", all_data$class_res_counts > 0)

colnames(plasmid_vs_amr) <- c("non-AMR", "AMR")
rownames(plasmid_vs_amr) <- c("non-pSF-088*", "pSF-088*")

chisq.test(plasmid_vs_amr, correct = FALSE)

plasmid_vs_amr

proportions(plasmid_vs_amr, margin = 1)


#MDR vs pSF-088* Pos/Neg
plasmid_vs_amr <- table(all_data$plasmid_map == "pSF-088*", all_data$class_res_counts > 2)

colnames(plasmid_vs_amr) <- c("non-MDR", "MDR")
rownames(plasmid_vs_amr) <- c("non-pSF-088*", "pSF-088*")

chisq.test(plasmid_vs_amr, correct = FALSE)

plasmid_vs_amr

proportions(plasmid_vs_amr, margin = 1)


#intI1 vs pSF-088* Pos/Neg
plasmid_vs_amr <- table(all_data$ColV, all_data$ESBL_ress)

colnames(plasmid_vs_amr) <- c("no ESBL_ress", "ESBL_ress")
rownames(plasmid_vs_amr) <- c("non-ColV", "ColV")

chisq.test(plasmid_vs_amr, correct = FALSE)


plasmid_vs_amr

proportions(plasmid_vs_amr, margin = 1)



#intI1 vs pUTI89* Pos/Neg
plasmid_vs_amr <- table(all_data$plasmid_map == "pUTI89*", all_data$intI1)

colnames(plasmid_vs_amr) <- c("no intI1", "intI1")
rownames(plasmid_vs_amr) <- c("non-pUTI89*", "ppUTI89*")

chisq.test(plasmid_vs_amr, correct = FALSE)

plasmid_vs_amr

#intI1 vs AMR

plasmid_vs_amr <- table(all_data$intI1, all_data$class_res_counts > 0)

colnames(plasmid_vs_amr) <- c("non-AMR", "AMR")
rownames(plasmid_vs_amr) <- c("no intI1", "intI1")

chisq.test(plasmid_vs_amr, correct = FALSE)


#intI1 vs MDR

plasmid_vs_amr <- table(all_data$intI1, all_data$class_res_counts > 2)

colnames(plasmid_vs_amr) <- c("non-MDR", "MDR")
rownames(plasmid_vs_amr) <- c("no intI1", "intI1")

chisq.test(plasmid_vs_amr, correct = FALSE)


#intI1 vs MDR

plasmid_vs_amr <- table(all_data$intI1, all_data$ESBL_ress)

colnames(plasmid_vs_amr) <- c("non-ESBL", "ESBL")
rownames(plasmid_vs_amr) <- c("no intI1", "intI1")

chisq.test(plasmid_vs_amr, correct = FALSE)





#plas_counts vs pUTI89* Pos/Neg
plasmid_vs_amr <- table(all_data$plasmid_map == "pUTI89*", all_data$plas_counts > 0)

colnames(plasmid_vs_amr) <- c("no other plas", "Other plas")
rownames(plasmid_vs_amr) <- c("non-pUTI89*", "pUTI89*")

chisq.test(plasmid_vs_amr, correct = FALSE)


#plas_counts vs pSF-088* Pos/Neg
plasmid_vs_amr <- table(all_data$plasmid_map == "pSF-088*", all_data$plas_counts > 0)

colnames(plasmid_vs_amr) <- c("no other plas", "Other plas")
rownames(plasmid_vs_amr) <- c("non-pSF-088*", "pSF-088*")

chisq.test(plasmid_vs_amr, correct = FALSE)


#AMR vs ColV Pos/Neg

plasmid_vs_amr <- table(all_data$ColV, all_data$plas_counts > 0)

colnames(plasmid_vs_amr) <- c("no other plas", "Other plas")
rownames(plasmid_vs_amr) <- c("ColV Neg", "ColV Positive")

chisq.test(plasmid_vs_amr, correct = FALSE)


plasmid_vs_amr

proportions(plasmid_vs_amr, margin = 1)


#intI1 vs plasmid_count

plasmid_vs_amr <- table(all_data$intI1, all_data$plas_counts > 0)

colnames(plasmid_vs_amr) <- c("no other plas", "Other plas")
rownames(plasmid_vs_amr) <- c("no intI1", "intI1")

chisq.test(plasmid_vs_amr, correct = FALSE)





#intI1 vs AMR

all_data$IS26_binary <- all_data$`ISfinder_Feb_2020_IS26:X00011`

all_data$IS26_binary[all_data$IS26_binary > 0] <- 1

plasmid_vs_amr <- table(all_data$IS26_binary, all_data$class_res_counts > 0)

colnames(plasmid_vs_amr) <- c("non-AMR", "AMR")
rownames(plasmid_vs_amr) <- c("no IS26", "IS26")

chisq.test(plasmid_vs_amr, correct = FALSE)


#intI1 vs MDR

plasmid_vs_amr <- table(all_data$IS26_binary, all_data$class_res_counts > 2)

colnames(plasmid_vs_amr) <- c("non-MDR", "MDR")
rownames(plasmid_vs_amr) <- c("no IS26", "IS26")

chisq.test(plasmid_vs_amr, correct = FALSE)


#intI1 vs MDR

plasmid_vs_amr <- table(all_data$IS26_binary, all_data$ESBL_ress)

colnames(plasmid_vs_amr) <- c("non-ESBL", "ESBL")
rownames(plasmid_vs_amr) <- c("no IS26", "IS26")

chisq.test(plasmid_vs_amr, correct = FALSE)











#intI1 vs plasmid_count

for(i in unique(all_data$clade_label)){
        plasmid_vs_amr <- table(all_data$clade_label == i, all_data$class_res_counts > 3)
        
        colnames(plasmid_vs_amr) <- c("not MDR", "MDR")
        rownames(plasmid_vs_amr) <- c("other", i)
        print(paste0("Testing statistical significance of clade ",i, " with MDR"))
        if(chisq.test(plasmid_vs_amr, correct = FALSE)[[3]] < 0.05){
                print(chisq.test(plasmid_vs_amr, correct = FALSE))
                print(proportions(plasmid_vs_amr, margin = 1))
                print(plasmid_vs_amr)
        }else{print("Association not significant")}
        readline("Press enter to see the next association")
}

#### t-tests ####

t.test(all_data$class_res_counts[all_data$ColV == "ColV Positive"], all_data$class_res_counts[all_data$ColV == "ColV Negative"])

res_data <- all_data %>% select(beta_lactam_res, sulfonamide_res, aminoglycoside_res, tetracycline_res, trimethoprim_res, fluoroquinolone_res, ESBL_ress, amphenicol_res, macrolide_lincosamide_res, polymyxin_res, class_res_counts, CIA_resistance, ColV)

res_data$AMR <- res_data$class_res_counts
res_data$MDR <- res_data$class_res_counts

res_data$AMR[res_data$AMR > 0] <- 1 

res_data$MDR[res_data$MDR < 3] <- 0
res_data$MDR[res_data$MDR > 2] <- 1

res_data$ESBL_ress[res_data$ESBL_ress =="ESBL_neg"] <- 0
res_data$ESBL_ress[res_data$ESBL_ress =="ESBL_pos"] <- 1

res_data$ColV[res_data$ColV =="ColV Negative"] <- 0
res_data$ColV[res_data$ColV =="ColV Positive"] <- 1

res_data <- res_data %>% select(-class_res_counts) %>% select(ColV, everything(), CIA_resistance)

res_data <- apply(res_data, MARGIN = 2, function(x) as.numeric(x))

res_data <- as.data.frame(res_data)

for(i in colnames(res_data)[2:ncol(res_data)]){
        ColV_pos_num <- res_data[res_data$ColV == 1, i]
        ColV_neg_num <- res_data[res_data$ColV == 0, i]
        print(paste(sum(ColV_pos_num), sum(ColV_neg_num), i, t.test(ColV_pos_num, ColV_neg_num)[3]))
}

t.test(all_data$class_res_counts[all_data$plasmid_map == "pUTI89*"], all_data$class_res_counts[all_data$ColV != "pUTI89*"])



ColV_neg <- all_data %>% filter(ColV == "ColV Negative")

t.test(ColV_neg$class_res_counts[ColV_neg$plasmid_map == "pUTI89*"], ColV_neg$class_res_counts[ColV_neg$ColV != "pUTI89*"])




notpUTI89 <- all_data %>% filter(plasmid_map != "pUTI89*")

t.test(notpUTI89$class_res_counts[notpUTI89$ColV == "ColV Positive"], notpUTI89$class_res_counts[notpUTI89$ColV == "ColV Negative"])
