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


#### Generate a list of Plasmids MLSTs with ColV status proportion ####
ColV_status <- geno_meta %>% select(Assembly_barcode, ColV)

ColV_status <- ColV_status %>% rename(working_name = Assembly_barcode)

pMLST <- pMLST %>% rename(working_name = "name")

ColV_status <- left_join(ColV_status, pMLST)

ColV_status <- ColV_status %>% select(pMLST, ColV)

ColV_status %>% group_by(pMLST, ColV) %>% summarise(counts = n()) %>% View()

amicolV <- ColV_status %>% group_by(pMLST, ColV) %>% summarise(counts = n()) %>% rename(IncF_RST = pMLST)


#### Generate a db of pMLST, ST, source, etc ####
STs <- ST %>% select(Uberstrain, ST)

working_names <- working_names %>% select(Uberstrain, `Assembly barcode`)

STs <- left_join(STs, working_names, by = "Uberstrain")

STs <- STs %>% rename(working_name = "Assembly barcode")

DB <- left_join(pMLST, STs, by = "working_name")

origins <- meta %>% select(Final_Source, Assembly_barcode, starts_with("Source"), Continent, Country, Collection_Year)

origins$Final_Source <- gsub("Compaion", "Companion",origins$Final_Source)

origins$Final_Source <- gsub("Human_Other", "Human_other",origins$Final_Source)

origins <- origins %>% rename(working_name = Assembly_barcode)

DB <- left_join(DB, origins)


# Group by Source
Source_sums <- DB %>% group_by(Final_Source) %>% summarise(counts = n())

Source_sums <- Source_sums %>% arrange(desc(counts))

# Group by pMLSTs
pMLST_sums <- DB %>% group_by(pMLST) %>% summarise(counts = n())

pMLST_sums <- pMLST_sums %>% arrange(desc(counts))

# Group pMLSTs by ST
pMLST_by_ST <- DB %>% group_by(pMLST, ST) %>% summarise(counts = n()) 

pMLST_by_ST <- pMLST_by_ST %>% arrange(desc(counts))

# Group pMLSTs by Source
pMLST_by_source <- DB %>% group_by(Final_Source, pMLST) %>% summarise(counts = n())

pMLST_by_source$simple_source <- pMLST_by_source$Final_Source

pMLST_by_source$simple_source <- gsub("Human_.*", "Human", pMLST_by_source$simple_source )

pMLST_by_source <- left_join(pMLST_by_source, pMLST_sums, by = "pMLST")

pMLST_by_source <- left_join(pMLST_by_source, Source_sums, by = "Final_Source")

pMLST_by_source <- pMLST_by_source %>% rename(Count = "counts.x", Total_pMLST = "counts.y", Total_source = "counts")

pMLST_by_source <- pMLST_by_source %>% arrange(desc(Count))


top_pMLSTs <- pMLST %>% group_by(pMLST) %>% summarise(counts = n()) %>% arrange(desc(counts)) %>% filter(pMLST != "F-:A-:B-") 

left_join(amicolV) %>% filter()




pMLST_by_source <- pMLST_by_source %>% filter(pMLST != "F-:A-:B-") %>% arrange(pMLST, Final_Source, Count) %>% filter(pMLST %in% top_pMLSTs$pMLST[1:10])

pMLST_by_source$percentage <- (pMLST_by_source$Count / pMLST_by_source$Total_pMLST)*100

df_fig <- pMLST_by_source %>% select(pMLST, Final_Source, percentage)

df_fig

df_fig$pMLST <- gsub("^C4:A-:B1$","F18:B1", df_fig$pMLST)

df_fig$pMLST <- gsub("A-:","", df_fig$pMLST)

orders <- top_pMLSTs$pMLST[1:10]

orders <- gsub("^C4:A-:B1$","F18:B1", orders)

orders <- gsub("A-:","", orders)

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
        theme_classic() + theme(axis.text = element_text(size = 15),
                                axis.title = element_text(size = 15)) +
        theme(plot.margin=unit(c(1,1,1,1.2),"cm"))
#geom_text(aes(label = paste(round(percentage, digits = 1),"%")), position = position_stack(vjust = 0.5), size = 2) #+
#coord_flip()


pMLST_by_ST$pMLST <- gsub("^C4:A-:B1$","F18:B1", pMLST_by_ST$pMLST)

pMLST_by_ST$pMLST <- gsub("A-:","", pMLST_by_ST$pMLST)


#pMLST_by_ST_top_5 <- pMLST_by_ST %>% slice_max(order_by = counts, n = 10) %>% filter(pMLST %in% top_pMLSTs$pMLST[1:10])
pMLST_by_ST_top_5 <- pMLST_by_ST %>% filter(pMLST %in% c("F29:B10", "F24:B1", "F2:B1"))

pMLST_sums_ <- pMLST_sums %>% rename('Total_pMLST_count' = counts)

pMLST_sums_$pMLST <- gsub("^C4:A-:B1$","F18:B1", pMLST_sums_$pMLST)

pMLST_sums_$pMLST <- gsub("A-:","", pMLST_sums_$pMLST)

pMLST_by_ST_top_5 <- left_join(pMLST_by_ST_top_5, pMLST_sums_)

pMLST_by_ST_top_5$percentage <- pMLST_by_ST_top_5$counts/pMLST_by_ST_top_5$Total_pMLST_count

#change ST to other where percentage (as ratio/1) is <0.01
pMLST_by_ST_top_5 <- pMLST_by_ST_top_5 %>% mutate(ST = ifelse(percentage <= 0.02, "Other", as.character(ST)))

pMLST_by_ST_top_5 <- pMLST_by_ST_top_5 %>% group_by(pMLST, ST) %>% summarise(counts = sum(counts), Total_pMLST_count = sum(Total_pMLST_count), percentage = sum(percentage))

pMLST_by_ST_top_5$pMLST <- gsub("^C4:A-:B1$","F18:B1", pMLST_by_ST_top_5$pMLST)

pMLST_by_ST_top_5$pMLST <- gsub("A-:","", pMLST_by_ST_top_5$pMLST)

orders <- gsub("^C4:A-:B1$","F18:B1", orders)

View(pMLST_by_ST_top_5)

set.seed(2)

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
        theme_classic() + theme(axis.text = element_text(size = 15),
                                axis.title = element_text(size = 15)) +
        theme(plot.margin=unit(c(1,1.3,1,1.2),"cm"))


legend4a <- get_legend(fig4a)
legend4b <- get_legend(fig4b)
spacer <- ggplot() + theme_minimal()

#nrow should be two but was set to three for easy padding
legends <- grid.arrange(legend4a, spacer, legend4b, ncol=1, nrow=4)

fig4a <- fig4a + theme(legend.position = 'none')

fig4b <- fig4b + theme(legend.position = 'none')

plots <- grid.arrange(fig4a, fig4b, ncol=2)

fig4 <- plot_grid(plots, legends, ncol=2, rel_widths = c(4/5, 1/5)) 

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





#write_delim(amicolV, "am_i_colv.txt", delim = "\t")

#write_delim(DB, "pMLST_DB.txt", delim = "\t")

#write_delim(pMLST_by_source, "pMLST_by_source.txt", delim = "\t")

#write_delim(pMLST_by_ST, "pMLST_by_ST.txt", delim = "\t")

#write_delim(Source_sums, "Source_sums.txt", delim = "\t")

#write_delim(pMLST_sums, "pMLST_sums.txt", delim = "\t")

DB %>% filter(pMLST == "F2:A-:B1") %>% group_by(ST) %>% summarise(counts = n()) %>% View()











