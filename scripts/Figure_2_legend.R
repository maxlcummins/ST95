library(ggtree)

if(!exists("plasmid_coverage_table")){
        source("scripts/ST95_all_process.R")
}

#Remove all objects from our environment except the ones we need
rm(list=setdiff(ls(), c("plasmid_coverage_table", "metadata", "df_small4")))

#Load in trees
core_tree_path <- "placeholder/output/ST95_all_CGA_snp_sites.tree"
accessory_tree_path <- "placeholder/output/accessory_ST95_all.tree"

#Generate plasmid map bins
plas_names <- plasmid_coverage_table %>% select(starts_with("p"),-plas_counts) %>% colnames()

plas_cuts <- c("0-49", "50-59", "60-69", "70-79", "80-89", "90-100")

plas_cells <- c()

for(i in plas_names){
        for(j in plas_cuts){
                plas_cells <- c(plas_cells, paste(i,j, sep = "_"))
        }
}

#Generate plasmid map colours
plasmid1_col_gen <- colorRampPalette(c("white", "red")) # pACN
plasmid2_col_gen <- colorRampPalette(c("white", "yellow")) # pAPEC 
plasmid3_col_gen <- colorRampPalette(c("white", "black")) # pBCE
plasmid4_col_gen <- colorRampPalette(c("white", "blue")) # pEC244
plasmid5_col_gen <- colorRampPalette(c("white", "orange")) # pSF
plasmid6_col_gen <- colorRampPalette(c("white", "purple")) # pU1
plasmid7_col_gen <- colorRampPalette(c("white", "green")) # pUTI89


plasmid1_cols <- plasmid1_col_gen(6)
plasmid2_cols <- plasmid2_col_gen(6)
plasmid3_cols <- plasmid3_col_gen(6)
plasmid4_cols <- plasmid4_col_gen(6)
plasmid5_cols <- plasmid5_col_gen(6)
plasmid6_cols <- plasmid6_col_gen(6)
plasmid7_cols <- plasmid7_col_gen(6)

#Preview colours
#plot(x = 1:8, y = rep(1, 8), col = class_res_count_cols, pch = 19)

#Bind colours into a list
plas_cols <- c(plasmid1_cols, plasmid2_cols, plasmid3_cols, plasmid4_cols, plasmid5_cols, plasmid6_cols, plasmid7_cols)


#Generate colour list for intI1
intI1_cols <- c(
        "white", # 0 - no gene present
        "black" # 1 - gene present
)

intI1_vars <- c(0,1)

names(intI1_cols) <- intI1_vars

#Generate colour list for ColV carriage
colV_cols <- c(
        "grey40", # ColV_neg - ColV present (Liu criteria)
        "red" # ColV_pos - ColV present (Liu criteria)
)

colv_vars <- c("ColV_neg","ColV_pos")


names(colV_cols) <- colv_vars

#Generate colour list for AMR
amr_col_gen <- colorRampPalette(c("white", "blue", "red"))

amr_count_cols <- amr_col_gen(13)

class_res_count_cols <- amr_col_gen(8)

class_res_count_vars <- c(
        "classes_res_00",
        "classes_res_01",
        "classes_res_02",
        "classes_res_03",
        "classes_res_04",
        "classes_res_05",
        "classes_res_06",
        "classes_res_07"
)

names(class_res_count_cols) <- class_res_count_vars

#Generate colour list for ESBL resistance
ESBL_res_cols <- c("white", # ESBL NEG
                   "black"  # ESBL POS
)

ESBL_res_vars <- c("ESBL_neg", "ESBL_pos")

names(ESBL_res_cols) <- ESBL_res_vars


#Generate colour list for plasmid counts (non-IncF repA)
plas_count_cols <- c("gray1", #   plas_count - res_count_0
                     "gray25", #   plas_count - res_count_1
                     "gray50", #   plas_count - res_count_2
                     "gray75", #   plas_count - res_count_3
                     "gray100" #   plas_count - res_count_4
)

plas_count_vars <- c("count_plas_4",
                     "count_plas_3",
                     "count_plas_2",
                     "count_plas_1",
                     "count_plas_0"
)

names(plas_count_cols) <- plas_count_vars

#Generate colour list for HC200 groups
HC200_cols <- c(        "#60c757"	,	#	1104	#	HC200		
                        "#47037a"	,	#	1108	#	HC200		
                        "#726c00"	,	#	1592	#	HC200		
                        "#ff6ed0"	,	#	4252	#	HC200		
                        "#ff9954"	,	#	44	#	HC200		
                        "#525799"	,	#	55	#	HC200		
                        "#d62678"	,	#	6624	#	HC200		
                        "#d4abff"	,	#	8655	#	HC200
                        "black"
)

HC200_vars <- c(1104,
                1108,
                1592,
                4252,
                44,
                55,
                6624,
                8655,
                "Other"
)

names(HC200_cols) <- HC200_vars


#Generate colour list for Pathogen type
Pathogen_cols <- c("green",  # Pathogen status        # Environmental
                   "black", # Pathogen status         # Flora
                   "grey", # Pathogen status         # Other
                   "pink", # Pathogen status         # Raw Chicken
                   "red", # Pathogen status         # Systemic
                   "yellow") # Pathogen status         # Urine

Pathogen_vars <- c(
        "Environmental",
        "Flora",
        "Other",
        "Raw Chicken",
        "Systemic",
        "Urine"
)

names(Pathogen_cols) <- Pathogen_vars


#Generate colour list for Source type
source_cols <- c(
        "brown"	,	#	Bovine	#	Source		
        "gold2"	,	#	Canine	#	Source		
        "springgreen4"	,#	Environment	#	Source		
        "cyan2"	,	#	Human	#	Source		
        "black"	,	#	Other	#	Source		
        "mediumorchid1" #	Poultry	#	Source
)

source_vars <- c(
        "Source_Bovine",
        "Source_Canine",
        "Source_Environment",
        "Source_Human",
        "Source_Other",
        "Source_Poultry")

names(source_cols) <- source_vars

#Generate colour list for reference plasmids
plas_vars <- c()

#Create all possible combinations of bins and plasmid references so that when we make our legend it includes all options
for(i in sort(plas_names)){
        plas_var <- c(paste0(i, "_0-49"), paste0(i, "_50-69"), paste0(i, "_60-69"), paste0(i, "_70-79"), paste0(i, "_80-89"), paste0(i, "_90-100"))
        plas_vars <- c(plas_vars, plas_var)
}

names(plas_cols) <- plas_vars


var_col <- c(HC200_cols, colV_cols, intI1_cols, plas_cols, class_res_count_cols, ESBL_res_cols, source_cols, Pathogen_cols, plas_count_cols)

#Generate colour list for Continents of origin
continent_cols <- c("#DDCC77", #sand
                    "#EE3377",  #magenta
                    "#009988",   #teal
                    "#4dac26" #dark green
)

continent_names <- c("North America","Europe", "Oceania", "Asia")

names(continent_cols) <- continent_names

#Generate colour list for fastbaps clades
fastbaps_clade_label_cols <- c(
        '#e6194B',
        '#3cb44b',
        #'#ffe119', yellow
        '#4363d8',
        '#f58231',
        '#42d4f4',
        #'#f032e6', magenta
        #'#fabed4', pink
        '#469990',
        #'#dcbeff',
        #'#9A6324',
        #'#fffac8', beige
        '#800000',
        '#aaffc3',
        '#000075',
        'black'
        #'#ffffff', white
)

names(fastbaps_clade_label_cols) <- c('Clade J',
                                      'Clade D',
                                      'Clade E',
                                      'Clade F',
                                      'Clade A',
                                      'Clade B',
                                      'Clade H',
                                      'Clade C',
                                      'Clade G',
                                      'Clade I'
)

var_col <- c(var_col, continent_cols, plas_cols, fastbaps_clade_label_cols)

#Read in trees
core_tree <- read.tree(core_tree_path)
accessory_tree <- read.tree(accessory_tree_path)

#Midpoint root rtrees
core_tree <- midpoint.root(core_tree)

#Clean the IncF RST data
metadata$IncF_RST <- gsub("F-:A-:B-", "", metadata$IncF_RST)

#Plot the tree
tree_core <- ggtree(core_tree,
                    open.angle = 0) %<+%
        metadata +
        geom_tiplab(size = 0.5,
                    align = TRUE,
                    linesize = 0.15,
                    aes(color = as.factor(ColV_pos_neg))
        ) +
        geom_tippoint(size = 0.2, aes(colour=as.factor(HC200_Other))
        )


#Create a small dataframe with the data of interest
fig1_df <- df_small4 %>% select(Revised_Source_Niche, plasmid_map, intI1, class_res_counts, plas_counts, Continent)

fig1_df <- fig1_df %>% select(Revised_Source_Niche)

#Combine all possible cell values into one big column to force all values into our legend
fig1_df$Null_column <- c(rep(names(var_col), 5)[1:668])

#Plot our pseudofigure for legend extraction
fig1a <- gheatmap(
        p = tree_core,
        data = fig1_df,
        colnames_offset_y = -0.1,
        font.size = 1.5,
        hjust = 0,
        colnames = FALSE,
        offset = 0.0020,
        width = 0.25,
        color = rgb(220, 220, 220, max = 255, alpha = 50)
) + 
        theme(legend.position = "bottom",
              legend.box = "vertical",
              legend.key.size = unit(1, "mm"),
              legend.title = element_blank(),
              legend.text = element_text(size = 8),
        ) +
        scale_fill_manual(
                aesthetics = c("colour", "fill"),
                values = var_col,
                na.value = 'grey'
        )

#Extract the legend
require(cowplot)
fig1_theme <- get_legend(fig1a)

plot(fig1_theme)

