library(ggtree)
require(cowplot)

#Create a function for finding objects in list A which are absent from list B
`%nin%` <- Negate(`%in%`)

#Check to see if plasmid coverage data object exists - if not, run the script which generates it.
if(!exists("plasmid_coverage_table")){
        source("scripts/ST95_all_process.R")
}

#Define paths to core tree and accessory tree files
core_tree_path <- "placeholder/output/ST95_all_CGA_snp_sites.tree"
accessory_tree_path <- "placeholder/output/accessory_ST95_all.tree"

#### Define colours ####

#Isolate the names of the plasmids
plas_names <- plasmid_coverage_table %>% select(starts_with("p"),-plas_counts) %>% colnames()

# Provide a list of bins to break the data points into
plas_cuts <- c("0-49", "50-59", "60-69", "70-79", "80-89", "90-100")

#Generate an empty vector
plas_cells <- c()

#For each of the plasmids we have screened for prepend the name of the plasmid with the bin
#eg. plasmid pABC will have elements pABC_0-49, pABC_50-59, pABC_60-69 etc generated
for(i in plas_names){
        for(j in plas_cuts){
                plas_cells <- c(plas_cells, paste(i,j, sep = "_"))
        }
}

#Create a function to generate a spectrum of colours for each plasmid type
plasmid1_col_gen <- colorRampPalette(c("white", "red")) # pACN
plasmid2_col_gen <- colorRampPalette(c("white", "yellow")) # pAPEC 
plasmid3_col_gen <- colorRampPalette(c("white", "black")) # pBCE
plasmid4_col_gen <- colorRampPalette(c("white", "blue")) # pEC244
plasmid5_col_gen <- colorRampPalette(c("white", "orange")) # pSF
plasmid6_col_gen <- colorRampPalette(c("white", "purple")) # pU1
plasmid7_col_gen <- colorRampPalette(c("white", "green")) # pUTI89

#Define the colours using our predefined functions
plasmid1_cols <- plasmid1_col_gen(6)
plasmid2_cols <- plasmid2_col_gen(6)
plasmid3_cols <- plasmid3_col_gen(6)
plasmid4_cols <- plasmid4_col_gen(6)
plasmid5_cols <- plasmid5_col_gen(6)
plasmid6_cols <- plasmid6_col_gen(6)
plasmid7_cols <- plasmid7_col_gen(6)

#If you want to see what your colours look like you can plot a preview like so
#plot(x = 1:8, y = rep(1, 8), col = class_res_count_cols, pch = 19)

#Bind the names of all colours for all plasmid-bin combinations into a vector
plas_cols <- c(plasmid1_cols, plasmid2_cols, plasmid3_cols, plasmid4_cols, plasmid5_cols, plasmid6_cols, plasmid7_cols)

#Generate an empty vector to host variable names for colouring plasmid coverage
plas_vars <- c()
for(i in sort(plas_names)){
        plas_var <- c(paste0(i, "_0-49"), paste0(i, "_50-69"), paste0(i, "_60-69"), paste0(i, "_70-79"), paste0(i, "_80-89"), paste0(i, "_90-100"))
        plas_vars <- c(plas_vars, plas_var)
}

#Assign the plasmid bin combination names to their respective colours
names(plas_cols) <- plas_vars

#Provide potential values to bind to the colours for intI1 carriage (so that 0 = absent = white and 1 = present = black)
intI1_cols <- c(
        "white", # 0 - no gene present
        "black" # 1 - gene present
)
intI1_vars <- c(0,1)
names(intI1_cols) <- intI1_vars

#Provide potential values to bind to the colours for ColV status as per Liu criteria (grey40 = absent, red = present)
colV_cols <- c(
        "grey40", # ColV_neg - ColV present (Liu criteria)
        "red" # ColV_pos - ColV present (Liu criteria)
)
colv_vars <- c("ColV_neg","ColV_pos")
names(colV_cols) <- colv_vars

#Define a colour palette for AMR gene carriage where white indicates no AMR genes, blue indicates mid-range amr gene carriage and red indicates maximal AMR gene carriage
amr_col_gen <- colorRampPalette(c("white", "blue", "red"))
#
amr_count_cols <- amr_col_gen(13)

class_res_count_cols <- amr_col_gen(9)

class_res_count_vars <- c(
        "classes_res_00",
        "classes_res_01",
        "classes_res_02",
        "classes_res_03",
        "classes_res_04",
        "classes_res_05",
        "classes_res_06",
        "classes_res_07",
        "classes_res_08"
)

names(class_res_count_cols) <- class_res_count_vars

#Define colours for ESBL carriage
ESBL_res_cols <- c("white", # ESBL NEG
                   "black"  # ESBL POS
)
ESBL_res_vars <- c("ESBL_neg", "ESBL_pos")
names(ESBL_res_cols) <- ESBL_res_vars

#Define plasmid replicon carriage colours
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

#Define colours for HC200 groupings
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

#Define colours for pathogen status
Pathogen_cols <- c("green",  # Pathogen status        # Environmental
                   "brown", # Pathogen status         # Flora
                   "black", # Pathogen status         # Other
                   "pink", # Pathogen status         # Raw Chicken
                   "red", # Pathogen status         # Systemic
                   "yellow") # Pathogen status         # Urine
Pathogen_vars <- c(
        "Environmental",
        "Gastrointestinal",
        "Other",
        "Raw Chicken",
        "Systemic",
        "Urine"
)
names(Pathogen_cols) <- Pathogen_vars

#Define colours for the sources a strain is from
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

#Bind all our colours together
var_col <- c(HC200_cols, colV_cols, intI1_cols, plas_cols, class_res_count_cols, ESBL_res_cols, source_cols, Pathogen_cols, plas_count_cols)

#Add colours for continent sources
continent_cols <- c("#DDCC77", #sand
                    "#EE3377",  #magenta
                    "#009988",   #teal
                    "#4dac26" #dark green
)

#Define names for the continents
continent_names <- c("North America","Europe", "Oceania", "Asia")



#Bind continent names with their colours
names(continent_cols) <- continent_names

#Combine continent and plasmid colours with our previous colour list
var_col <- c(var_col, continent_cols, plas_cols, Pathogen_cols)


#read in fastbaps
fastbaps_data <- read.csv(file = "fastbaps.csv")

colnames(fastbaps_data) <- c("working_name", "Level_1", "Level_2")

#Add leading zeroes
fastbaps_data$Level_1 <- str_pad(string = fastbaps_data$Level_1, side = "left", pad = "0", width = 2)
fastbaps_data$Level_2 <- str_pad(string = fastbaps_data$Level_2, side = "left", pad = "0", width = 2)

fastbaps_data$Level_1 <- gsub("$","_fastbaps", fastbaps_data$Level_1)
fastbaps_data$Level_2 <- gsub("$","_fastbaps", fastbaps_data$Level_2)

#Define colours for our fastbaps clades
fastbaps_level_1_names <-
        fastbaps_data$Level_1 %>% unique() %>% sort()
fastbaps_level_1_cols <-
        c(
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
                'black',
                #'#ffffff', white
                '#000000'
        )[1:length(fastbaps_level_1_names)]

names(fastbaps_level_1_cols) <- fastbaps_level_1_names

#Bind in our new colours to our colour list
var_col <- c(var_col, fastbaps_level_1_cols)


metadata <- left_join(metadata, fastbaps_data)

#### Tree visualisation ####
#Read in tree files
core_tree <- read.tree(core_tree_path)
accessory_tree <- read.tree(accessory_tree_path)

#Midpoint root core tree
core_tree <- midpoint.root(core_tree)

#Replace NULL IncF_RSTs with a blank space
metadata$IncF_RST <- gsub("F-:A-:B-", "", metadata$IncF_RST)
#Get rid of A NULL notation
metadata$IncF_RST <- gsub(":A-", "", metadata$IncF_RST)
#Replace C4 variant with F18 variant - these sequences overlap but C4 has a slightly higher scoring metric which means it gets selected
metadata$IncF_RST <- gsub("^C4:B$", "F18:B", metadata$IncF_RST)



#Plot the core tree
tree_core <- ggtree(core_tree, layout = 'circular',
                    open.angle = 0) %<+%
        metadata +
        geom_tiplab(size = 0.5,
                    align = TRUE,
                    linesize = 0.15,
                    aes(color = as.factor(ColV_pos_neg))
        ) +
        geom_tippoint(size = 0.2, aes(colour=as.factor(HC200))
        )


#Define our data which we want to generate our rings around the tree with
fig1_df <- df_small4 %>% select(Pathogen, Revised_Source_Niche, plasmid_map, intI1, class_res_counts, plas_counts, Continent)

#Generate our first core tree
fig1a <- gheatmap(
        p = tree_core,
        data = fig1_df,
        colnames_offset_y = -0.1,
        font.size = 1.5,
        hjust = 0,
        colnames = FALSE,
        offset = 0.0015,
        width = 0.3,
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
        ) +
        geom_tiplab2(
                aes(label = IncF_RST),
                align = TRUE,
                linetype = NULL,
                offset = 0.00575,
                size = 0.5
        ) +
        geom_tiplab2(
                aes(label = fimH_type),
                align = TRUE,
                linetype = NULL,
                offset = 0.007,
                size = 0.5
        )+
        geom_tiplab2(
                aes(label = O_type),
                align = TRUE,
                linetype = NULL,
                offset = 0.00775,
                size = 0.5
        ) +
        geom_tiplab2(
                aes(label = H_type),
                align = TRUE,
                linetype = NULL,
                offset = 0.0085,
                size = 0.5
        )


fastbaps_core <- fortify(fig1a$data) %>%
        filter(isTip == TRUE) %>%
        select(label, Level_1, node, y, O_type, H_type, fimH_type) %>%
        arrange(y) %>%
        group_by(Level_1) %>%
        summarise(
                firstocc = first(label),
                lastocc = last(label),
                count = n(),
                top_O = paste0(unique(sort(O_type)),collapse = "/"),
                top_H = paste0(unique(sort(H_type)),collapse = "/"),
                fimH = paste0(unique(sort(fimH_type)),collapse = "/")
        )

#Clade label colours
custom_labels <-
        cbind(
                c(
                        '04_fastbaps',
                        '09_fastbaps',
                        '07_fastbaps',
                        '08_fastbaps',
                        '10_fastbaps',
                        '05_fastbaps',
                        '02_fastbaps',
                        '01_fastbaps',
                        '03_fastbaps',
                        '06_fastbaps'
                ),
                #LETTERS[1:nrow(fastbaps_core)]
                c( 'F',#'04_fastbaps',
                   'G',#'09_fastbaps',
                   'H',#'07_fastbaps',
                   'C',#08_fastbaps,
                   'I',#'10_fastbaps',
                   'A',#'05_fastbaps',
                   'D',#'02_fastbaps',
                   'J',#'01_fastbaps',
                   'E',#'03_fastbaps',
                   'B'#'06_fastbaps'
                )
        )

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
#Level 1s above translate to fastbaps Level 1 groups as below.
#'04_fastbaps',
#'09_fastbaps',
#'07_fastbaps',
#'08_fastbaps',
#'10_fastbaps',
#'05_fastbaps',
#'02_fastbaps',
#'01_fastbaps',
#'03_fastbaps',
#'06_fastbaps'

var_col <- c(var_col, fastbaps_clade_label_cols)
        
custom_labels <- as.data.frame(custom_labels) %>% rename("Level_1" = "V1", "clade_label" = "V2")

fastbaps_core <- left_join(fastbaps_core, custom_labels)

fastbaps_core$color <- fastbaps_level_1_cols

strip_fontsize <- 8
strip_offset <- 0.00885
strip_text_offset <- 0.0037
strip_barsize <- 1
strip_text_hjust <- 0.5

#Plot our core tree
fig1a <- fig1a + geom_strip(
        fastbaps_core$firstocc[1],
        fastbaps_core$lastocc[1],
        barsize = strip_barsize,
        color = fastbaps_core$color[1],
        label = fastbaps_core$clade_label[1],
        fontsize = strip_fontsize,
        offset = strip_offset,
        offset.text = strip_text_offset,
        hjust = strip_text_hjust
) +
        geom_strip(
                fastbaps_core$firstocc[2],
                fastbaps_core$lastocc[2],
                barsize = strip_barsize,
                color = fastbaps_core$color[2],
                label = fastbaps_core$clade_label[2],
                fontsize = strip_fontsize,
                offset = strip_offset,
                offset.text = strip_text_offset,
                hjust = strip_text_hjust
        ) +
        geom_strip(
                fastbaps_core$firstocc[3],
                fastbaps_core$lastocc[3],
                barsize = strip_barsize,
                color = fastbaps_core$color[3],
                label = fastbaps_core$clade_label[3],
                fontsize = strip_fontsize,
                offset = strip_offset,
                offset.text = strip_text_offset,
                hjust = strip_text_hjust
        ) +
        geom_strip(
                fastbaps_core$firstocc[4],
                fastbaps_core$lastocc[4],
                barsize = strip_barsize,
                color = fastbaps_core$color[4],
                label = fastbaps_core$clade_label[4],
                fontsize = strip_fontsize,
                offset = strip_offset,
                offset.text = strip_text_offset,
                hjust = strip_text_hjust
        ) +
        geom_strip(
                fastbaps_core$firstocc[5],
                fastbaps_core$lastocc[5],
                barsize = strip_barsize,
                color = fastbaps_core$color[5],
                label = fastbaps_core$clade_label[5],
                fontsize = strip_fontsize,
                offset = strip_offset,
                offset.text = strip_text_offset,
                hjust = strip_text_hjust
        ) +
        geom_strip(
                fastbaps_core$firstocc[6],
                fastbaps_core$lastocc[6],
                barsize = strip_barsize,
                color = fastbaps_core$color[6],
                label = fastbaps_core$clade_label[6],
                fontsize = strip_fontsize,
                offset = strip_offset,
                offset.text = strip_text_offset,
                hjust = strip_text_hjust
        ) +
        geom_strip(
                fastbaps_core$firstocc[7],
                fastbaps_core$lastocc[7],
                barsize = strip_barsize,
                color = fastbaps_core$color[7],
                label = fastbaps_core$clade_label[7],
                fontsize = strip_fontsize,
                offset = strip_offset,
                offset.text = strip_text_offset,
                hjust = strip_text_hjust
        ) +
        geom_strip(
                fastbaps_core$firstocc[8],
                fastbaps_core$lastocc[8],
                barsize = strip_barsize,
                color = fastbaps_core$color[8],
                label = fastbaps_core$clade_label[8],
                fontsize = strip_fontsize,
                offset = strip_offset,
                offset.text = strip_text_offset,
                hjust = strip_text_hjust
        ) +
        geom_strip(
                fastbaps_core$firstocc[9],
                fastbaps_core$lastocc[9],
                barsize = strip_barsize,
                color = fastbaps_core$color[9],
                label = fastbaps_core$clade_label[9],
                fontsize = strip_fontsize,
                offset = strip_offset,
                offset.text = strip_text_offset,
                hjust = strip_text_hjust
        ) +
        geom_strip(
                fastbaps_core$firstocc[10],
                fastbaps_core$lastocc[10],
                barsize = strip_barsize,
                color =  fastbaps_core$color[10],
                label = fastbaps_core$clade_label[10],
                fontsize = strip_fontsize,
                offset = strip_offset,
                offset.text = strip_text_offset,
                hjust = strip_text_hjust
        )+ theme(legend.position = "none") 
        

metadata <- left_join(metadata, custom_labels)


#Midpoint root our accessory tree
accessory_tree <- midpoint.root(accessory_tree)

#Plot our accessory tree
tree_accessory <- ggtree(accessory_tree, layout = 'circular',
                         open.angle = 0) %<+%
        metadata +
        geom_tiplab(size = 0.5,
                    align = TRUE,
                    linesize = 0.15,
                    aes(color = as.factor(ColV_pos_neg))
        ) +
        geom_tippoint(size = 0.2, aes(colour=as.factor(HC200))
        )


#Plot our accessory tree
fig1b <- gheatmap(
        p = tree_accessory,
        data = fig1_df,
        colnames_offset_y = -0.1,
        font.size = 1.5,
        hjust = 0,
        colnames = FALSE,
        offset = 0.08,
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
        ) +
        geom_tiplab2(
                aes(label = IncF_RST),
                align = TRUE,
                linetype = NULL,
                offset = 0.25,
                size = 0.5
        ) +
        geom_tiplab2(
                aes(label = fimH_type),
                align = TRUE,
                linetype = NULL,
                offset = 0.3,
                size = 0.5
        )+
        geom_tiplab2(
                aes(label = O_type),
                align = TRUE,
                linetype = NULL,
                offset = 0.35,
                size = 0.5
        ) +
        geom_tiplab2(
                aes(label = H_type),
                align = TRUE,
                linetype = NULL,
                offset = 0.375,
                size = 0.5
        ) +
        geom_tiplab2(
                aes(label = gsub("^", "Clade ", clade_label), colour=gsub("^", "Clade ", clade_label)),
                align = TRUE,
                linetype = NULL,
                offset = 0.4,
                size = 0.75
        )        


#Plot our two trees together without their legends
tree_plot <- fig1a + theme(legend.position="none") + fig1b+ theme(legend.position="none") 


Supp_Tab_1 <- metadata

Supp_Tab_1 <- Supp_Tab_1 %>% select(-starts_with(ST95_all.N90.L90.PASS$DATABASE), -ColV, -ColV_pos_neg)

#write_delim(Supp_Tab_1, "manuscript/Publication/Publication_figures/Supplementary_Table_1.txt",delim = "\t")

#write_delim(ST95_all_co_occurence_N90L90, "manuscript/Publication/Publication_figures/Supplementary_Table_4.txt", delim = "\t")
