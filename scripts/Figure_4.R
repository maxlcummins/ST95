library(cowplot)
library(grid)
library(gridExtra)
library(dplyr)
library(magrittr)
library(readr)
library(ggplot2)

#Read in curated metadata
meta <- read_delim("data/enterobase_diverse/Manual_edit_New_Jan_29_source_New_ColV_Zoetis.txt", 
                   "\t", escape_double = FALSE, trim_ws = TRUE)

#Read in ColV data (partial subset)
geno_meta <- read_delim("data/enterobase_diverse/EnterolordDB_31-3-20.tsv", 
                        "\t")


#Read in ST daata
ST <- read_delim("data/enterobase_diverse/Zoetis_Subset_w-ST.txt", 
                 "\t")

#Read in working_name data
working_names <- read_delim("data/enterobase_diverse/ColV_Zoetis_14-1-20.txt", 
                            "\t")

#Read in pMLST data
pMLST <- read_delim("data/enterobase_diverse/incF_pMLST.txt", 
                    "\t", escape_double = FALSE, trim_ws = TRUE)

pMLST <- pMLST %>% filter(name %in% geno_meta$Assembly_barcode)

#Create supplementary table 2 with the desired columns from our dataframe
Supp_Table_2 <- geno_meta[,420:ncol(geno_meta)]

#Rename the 
Supp_Table_2 <- pMLST %>% rename("Assembly_barcode" = "name") %>% right_join(Supp_Table_2)

#Remove duplicated columns
Supp_Table_2 <- Supp_Table_2 %>% select(-ends_with(".x"), -ends_with(".y"))

#Convert C4s to F18s - these are misidentified by IncF RST due to the competing bitscore values and overlapping ORFs
Supp_Table_2$pMLST <- gsub("^C4:A-:B1$","F18:B1", Supp_Table_2$pMLST)

#Delete A null to clean IncF RSTs
Supp_Table_2$pMLST <- gsub("A-:","", Supp_Table_2$pMLST)

#Unhash below to write the file
#write_delim(Supp_Table_2, "manuscript/Publication/Publication_figures/Supplementary_Table_2.txt",delim = "\t")

#### Generate a list of Plasmids MLSTs with ColV status proportion ####
ColV_status <- geno_meta %>% select(Assembly_barcode, ColV)

#Designate sample names
ColV_status <- ColV_status %>% rename(working_name = Assembly_barcode)

#Combine CoLV status with pMLST data
ColV_status <- left_join(ColV_status, pMLST)

#Isolate columns of interest
ColV_status <- ColV_status %>% select(pMLST, ColV)

#Rename our name column in ST sheet
ST %>% rename("working_name" = Assembly_barcode)


#### Generate a db of pMLST, ST, source, etc ####
STs <- ST %>% select(Uberstrain, ST)

#### Generate a db of pMLST, ST, source, etc ####
pMLST <- pMLST %>% rename("working_name" = "name")

#Generate a table of Uberstrain and Assembly barcodes (names) to bind ST and pMLST data
working_names <- working_names %>% select(Uberstrain, `Assembly barcode`)

#Bind ST and names table
STs <- left_join(STs, working_names, by = "Uberstrain")

#Rename our name columns
STs <- STs %>% rename(working_name = "Assembly barcode")

#Bind our dataframes
DB <- left_join(pMLST, STs, by = "working_name")

#Isolate metadata for the isolates
origins <- meta %>% select(Final_Source, Assembly_barcode, starts_with("Source"), Continent, Country, Collection_Year)

#Clean some designations/typos
origins$Final_Source <- gsub("Compaion", "Companion",origins$Final_Source)
origins$Final_Source <- gsub("Human_Other", "Human_other",origins$Final_Source)

#Rename our name column
origins <- origins %>% rename(working_name = Assembly_barcode)

#Combine our data with our metadata
DB <- left_join(DB, origins)


# Group by Source
Source_sums <- DB %>% group_by(Final_Source) %>% summarise(counts = n())
# Arrange in descending order of source counts
Source_sums <- Source_sums %>% arrange(desc(counts))

# Group by pMLSTs
pMLST_sums <- DB %>% group_by(pMLST) %>% summarise(counts = n())
# Arrange in descending order of pMLST counts
pMLST_sums <- pMLST_sums %>% arrange(desc(counts))

# Group pMLSTs by ST
pMLST_by_ST <- DB %>% group_by(pMLST, ST) %>% summarise(counts = n()) 

# Arrange in descending order of pMLST-ST counts
pMLST_by_ST <- pMLST_by_ST %>% arrange(desc(counts))

# Group pMLSTs by Source and pMLST
pMLST_by_source <- DB %>% group_by(Final_Source, pMLST) %>% summarise(counts = n())

#Create a new column for a simple source with less granular data
pMLST_by_source$simple_source <- pMLST_by_source$Final_Source

#Combine Human ExPEC and Human Other categories
pMLST_by_source$simple_source <- gsub("Human_.*", "Human", pMLST_by_source$simple_source )

#Join our pMLST_by_source tables and pMLST sums table to get total pMLST counts to facilitate generating percentages
pMLST_by_source <- left_join(pMLST_by_source, pMLST_sums, by = "pMLST")

#Join our pMLST_by_source tables and Source sums table to get total Source counts to facilitate generating percentages
pMLST_by_source <- left_join(pMLST_by_source, Source_sums, by = "Final_Source")

#Rename columns appropriately
pMLST_by_source <- pMLST_by_source %>% rename(Count = "counts.x", Total_pMLST = "counts.y", Total_source = "counts")

#Arrange by descending count
pMLST_by_source <- pMLST_by_source %>% arrange(desc(Count))

#Identify top pMLSTs (excluding F null variants)
top_pMLSTs <- pMLST %>% group_by(pMLST) %>% summarise(counts = n()) %>% arrange(desc(counts)) %>% filter(pMLST != "F-:A-:B-") 

#Identify top 10 pMLSTs by source
pMLST_by_source <- pMLST_by_source %>% filter(pMLST != "F-:A-:B-") %>% arrange(pMLST, Final_Source, Count) %>% filter(pMLST %in% top_pMLSTs$pMLST[1:10])

#Generate percentages
pMLST_by_source$percentage <- (pMLST_by_source$Count / pMLST_by_source$Total_pMLST)*100

#Select columns of interest for our figure
df_fig <- pMLST_by_source %>% select(pMLST, Final_Source, percentage)

#Clean pMLSTs as above
df_fig$pMLST <- gsub("^C4:A-:B1$","F18:B1", df_fig$pMLST)
df_fig$pMLST <- gsub("A-:","", df_fig$pMLST)

#Define an order of pMLSTs we later use to pick the order of x axis in our bar graphs
orders <- top_pMLSTs$pMLST[1:10]

#Clean pMLSTs as above
orders <- gsub("^C4:A-:B1$","F18:B1", orders)
orders <- gsub("A-:","", orders)

#Define colours for our sources
source_cols <- c(
        "brown"	,	#	Bovine	#	Source		
        "gold2"	,	#	Canine	#	Source		
        "cyan2"	,	#	Human	#	Source		
        "black"	,	#	Other	#	Source		
        "mediumorchid1", #	Poultry	#	Source
        "DarkSlateGray",
        "DarkMagenta"
)

source_vars <- c(
        "Bovine",
        "Companion_animal",
        "Human_other",
        "Other",
        "Poultry",
        "Human_ExPEC",
        "Porcine")

names(source_cols) <- source_vars

#Filter out our pMLSTs of interest
df_fig <- df_fig %>% filter(pMLST %in% c("F29:B10", "F24:B1", "F2:B1"))
orders <- orders[orders %in% c("F29:B10", "F24:B1", "F2:B1")]


#Figure 4
fig4a <- ggplot(df_fig, aes(fill=Final_Source, y=percentage, x=pMLST)) + 
        geom_bar(position="fill", stat="identity" , color = "black") +
        ylab(label = "Percentage") +
        scale_x_discrete(limits = orders) +
        scale_y_continuous(expand = c(0,0), labels = scales::percent) +
        scale_fill_manual(labels = c("Bovine", "Companion animal", "Human ExPEC", "Human Other", "Porcine", "Poultry"),
                          name = "Source",
                          aesthetics = c("colour", "fill"),
                          values = source_cols,
                          na.value = 'grey'
        ) +
        xlab(label = "IncF RST") +
        theme_classic() + theme(axis.text = element_text(size = 15),
                                axis.title = element_text(size = 15)) +
        theme(plot.margin=unit(c(1,0,1,0.25),"cm"))



#Clean pMLSTs as above
pMLST_by_ST$pMLST <- gsub("^C4:A-:B1$","F18:B1", pMLST_by_ST$pMLST)
pMLST_by_ST$pMLST <- gsub("A-:","", pMLST_by_ST$pMLST)

pMLST_by_ST_top_5 <- pMLST_by_ST %>% filter(pMLST %in% c("F29:B10", "F24:B1", "F2:B1"))

pMLST_sums_ <- pMLST_sums %>% rename('Total_pMLST_count' = counts)

#Clean pMLSTs as above
pMLST_sums_$pMLST <- gsub("^C4:A-:B1$","F18:B1", pMLST_sums_$pMLST)
pMLST_sums_$pMLST <- gsub("A-:","", pMLST_sums_$pMLST)

pMLST_by_ST_top_5 <- left_join(pMLST_by_ST_top_5, pMLST_sums_)

pMLST_by_ST_top_5$percentage <- pMLST_by_ST_top_5$counts/pMLST_by_ST_top_5$Total_pMLST_count

#change ST to other where percentage (as ratio/1) is <0.01
pMLST_by_ST_top_5 <- pMLST_by_ST_top_5 %>% mutate(ST = ifelse(percentage <= 0.02, "Other", as.character(ST)))

#Generate a table for Figure 4
pMLST_by_ST_top_5 <- pMLST_by_ST_top_5 %>% group_by(pMLST, ST) %>% summarise(counts = sum(counts), Total_pMLST_count = sum(Total_pMLST_count), percentage = sum(percentage))
pMLST_by_ST_top_5$pMLST <- gsub("^C4:A-:B1$","F18:B1", pMLST_by_ST_top_5$pMLST)
pMLST_by_ST_top_5$pMLST <- gsub("A-:","", pMLST_by_ST_top_5$pMLST)

#Define X axis order
orders <- gsub("^C4:A-:B1$","F18:B1", orders)

#Figure 4
fig4b <- ggplot(pMLST_by_ST_top_5, aes(fill=gsub("^([0-9])","ST\\1",ST), y=percentage, x=pMLST)) + 
        geom_bar(position="fill", stat="identity", color = "black") +
        scale_x_discrete(limits = orders)+
        scale_x_discrete(limits = orders)+
        geom_text(aes(label = gsub("^([0-9])","ST\\1",ST)), position = position_stack(vjust=0.5), size = 3) +
        ylab(label = "Percentage") +
        scale_y_continuous(expand = c(0,0), labels = scales::percent) +
        scale_fill_manual(
                name = "Sequence Type",
                aesthetics = c("colour", "fill"),
                values = c('#ffffff', '#3cb44b', '#ffe119', '#4363d8', '#f58231', '#911eb4', '#46f0f0', '#f032e6', '#bcf60c', '#fabebe', '#008080', '#e6beff', '#9a6324', '#fffac8', '#800000', '#aaffc3', '#808000', '#ffd8b1', '#4c4c9e', '#e6194b', '#808080', '#000000'),#sample(rainbow(n = 29)),
                na.value = 'grey'
        ) +
        xlab(label = "IncF RST") +
        theme_classic() + theme(axis.text = element_text(size = 15),
                                axis.title = element_text(size = 15)) +
        theme(plot.margin=unit(c(1,0,1,1),"cm"))


#Export legend
legend4a <- get_legend(fig4a)
legend4b <- get_legend(fig4b)
spacer <- ggplot() + theme_minimal()

#nrow should be two but was set to three for easy padding
legends <- grid.arrange(legend4a, spacer, legend4b, ncol=1, nrow=4)

#Combine figure 4a and 4b
fig4a <- fig4a + theme(legend.position = 'none')

fig4b <- fig4b + theme(legend.position = 'none')

plots <- grid.arrange(fig4a, fig4b, ncol=2)

fig4 <- plot_grid(plots, legends, ncol=2, rel_widths = c(4/5, 1/5)) 


#Add annotations        
fig4 + theme(plot.margin=unit(c(3,1,0.75,1.2),"cm")) + 
        annotate(geom = "text", x = 0.015, y = 1, label = "A",
                 parse = TRUE,
                 size = 15) + 
        annotate(geom = "text", x = 0.425, y = 1, label = "B",
                 parse = TRUE,
                 size = 15) +
        annotate(geom = "text", x = 0.125, y = 0.97, label = "n==1078",
                 parse = TRUE,
                 size = 5,
                 hjust = 0) +
        annotate(geom = "text", x = 0.21, y = 0.97, label = "n==611",
                 parse = TRUE,
                 size = 5,
                 hjust = 0) +
        annotate(geom = "text", x = 0.29, y = 0.97, label = "n==491",
                 parse = TRUE,
                 size = 5,
                 hjust = 0) +
        annotate(geom = "text", x = 0.53, y = 0.97, label = "n==1078",
                 parse = TRUE,
                 size = 5,
                 hjust = 0) +
        annotate(geom = "text", x = 0.61, y = 0.97, label = "n==611",
                 parse = TRUE,
                 size = 5,
                 hjust = 0) +
        annotate(geom = "text", x = 0.68, y = 0.97, label = 'n==491',
                 parse = TRUE,
                 size = 5,
                 hjust = 0)

#Generate a dataframe for ColV data
ColV_by_source <- geno_meta %>% group_by(Final_Source, ColV) %>% summarise(counts = n())

#Change ColV binary values to Character based values
ColV_by_source$ColV[ColV_by_source$ColV == 0] <- "no ColV-LP"
ColV_by_source$ColV[ColV_by_source$ColV == "1"] <- "ColV-LP"

#Cleaning column source values
ColV_by_source$Final_Source <- gsub("_", " ", ColV_by_source$Final_Source)
ColV_by_source$Final_Source <- gsub(" animal", " Animal", ColV_by_source$Final_Source)
ColV_by_source$Final_Source <- gsub(" other", " Other", ColV_by_source$Final_Source)

#Generate figure 4c
fig4c <- ggplot(ColV_by_source, aes(fill=as.factor(ColV), y=counts, x=Final_Source)) + 
        geom_bar(position="fill", stat="identity" , color = "black") +
        ylab(label = "Percentage") +
        scale_x_discrete(limits = c("Poultry", "Porcine", "Human ExPEC", "Companion Animal", "Human Other", "Bovine"))+
        scale_y_continuous(expand = c(0,0), labels = scales::percent) +
        scale_fill_manual(labels = c("Yes", "No"),
                          name = "ColV-like plasmid",
                          aesthetics = c("colour", "fill"),
                          values = c("ColV-LP" = rgb(255, 0, 0, max = 255, alpha = 200), "no ColV-LP" = "grey40"),
                          na.value = 'grey'
        ) +
        xlab(label = "Isolate source") +
        theme_classic() + theme(axis.text = element_text(size = 15),
                                axis.title = element_text(size = 15),
                                axis.text.x = element_text(size = 8, angle = 15, vjust = 0.8, hjust = 1)) +
        theme(plot.margin=unit(c(1,1.5,1,1),"cm"))

#Extract the legend
legend4c <- get_legend(fig4c)

#Plot without the legend
fig4c <- fig4c + theme(legend.position = 'none')

#Combine all plots
plots2 <- grid.arrange(fig4a, fig4b, fig4c, ncol=3)

#Combine all legends
legends2 <- grid.arrange(legend4a, legend4b, legend4c, ncol=1, nrow=3)

#Combine plots and legends
fig4 <- plot_grid(plots2, legends2, ncol=2, rel_widths = c(9.5/10, 0.5/10)) 

fig4 <- fig4 + theme(plot.margin=unit(c(1,1.5,0,0.2),"cm"))

#Generate figure 4 
fig4+ #Add annotations for Subfigure labels
        annotate(geom = "text", x = 0.05, y = 1, label = "A",
                 parse = TRUE,
                 size = 12) + 
        annotate(geom = "text", x = 0.385, y = 1, label = "B",
                 parse = TRUE,
                 size = 12) +
        annotate(geom = "text", x = 0.7, y = 1, label = "C",
                 parse = TRUE,
                 size = 12) +
        #Annotations for Figure C
        annotate(geom = "text", x = 0.09, y = 0.97, label = "n==1078",
                 parse = TRUE,
                 size = 5,
                 hjust = 0) +
        annotate(geom = "text", x = 0.17, y = 0.97, label = "n==611",
                 parse = TRUE,
                 size = 5,
                 hjust = 0) +
        annotate(geom = "text", x = 0.25, y = 0.97, label = "n==491",
                 parse = TRUE,
                 size = 5,
                 hjust = 0) +
        #Annotations for Figure B
        annotate(geom = "text", x = 0.425, y = 0.97, label = "n==1078",
                 parse = TRUE,
                 size = 5,
                 hjust = 0) +
        annotate(geom = "text", x = 0.5, y = 0.97, label = "n==611",
                 parse = TRUE,
                 size = 5,
                 hjust = 0) +
        annotate(geom = "text", x = 0.57, y = 0.97, label = 'n==491',
                 parse = TRUE,
                 size = 5,
                 hjust = 0) +
        #Annotations for Figure C
        annotate(geom = "text", x = 0.745, y = 0.97, label = 'n==4254',
                 parse = TRUE,
                 size = 4,
                 angle = 90,
                 hjust = 0.15) +
        annotate(geom = "text", x = 0.7725, y = 0.97, label = 'n==2679',
                parse = TRUE,
                size = 4,
                angle = 90,
                hjust = 0.15) +
        annotate(geom = "text", x = 0.8, y = 0.97, label = 'n==4425',
                 parse = TRUE,
                 size = 4,
                 angle = 90,
                 hjust = 0.15) +
        annotate(geom = "text", x = 0.83, y = 0.97, label = 'n==941',
                 parse = TRUE,
                 size = 4,
                 angle = 90,
                 hjust = 0.15) +
        annotate(geom = "text", x = 0.8575, y = 0.97, label = 'n==14975',
                 parse = TRUE,
                 size = 4,
                 angle = 90,
                 hjust = 0.15) +
        annotate(geom = "text", x = 0.8875, y = 0.97, label = 'n==6902',
                 parse = TRUE,
                 size = 4,
                 angle = 90,
                 hjust = 0.15)+ theme(plot.margin=unit(c(1.25,1.5,0,0.2),"cm"))
        

