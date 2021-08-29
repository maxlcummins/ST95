library(ggtree)
require(cowplot)

#Check if the data required for this script is in our environment, or else generate it
if(!exists("metadata")){
        source("scripts/Clade_H2_process.R")
}

#Check that the number of strains is correct and that we aren't importing the full metadata list
if(nrow(metadata) != 51){
        source("scripts/Clade_H2_process.R")
}

#Remove all objects from our environment except the ones we need
rm(list=setdiff(ls(), c("df", "metadata")))

tree_path_2 <- "placeholder/output/overlap_group_snp_sites.treefile"

tree_path_accessory <- "placeholder/output/accessory_overlap_group.tree"

#Define our colours for particular gene types
cols_fig3 <- c("0" = "white",  "1" = "#bebada", "2" = "#80b1d3", "3" = "#fb8072", "4" = "black", "5" = "purple")

#Generate our dataframe we wish to plot next to our tree
df_Fig2 <- df %>% select(-starts_with("colv_zoetis")) %>% select(-starts_with("IS_"))

#Reassign metadata for data manipulation
meta <- metadata

#Replace NULL IncF_RSTs with a blank space
meta$IncF_RST <- gsub("F-:A-:B-", "", meta$IncF_RST)
#Get rid of A NULL notation
meta$IncF_RST <- gsub(":A-", "", meta$IncF_RST)
#Replace C4 variant with F18 variant - these sequences overlap but C4 has a slightly higher scoring metric which means it gets selected
meta$IncF_RST <- gsub("^C4:", "F18:", meta$IncF_RST)


#Read in our tree
tree2 <- read.tree(tree_path_2)

#Extract the column names from the input for scrutiny. I then manually edited some in a text editor, flagging duplicates which were present in multiple databases
# with *** to allow their removal. I also cleaned the names of things with long ugly names and fixed some genes which didn't correctly have gene counts in their names
df_Fig2_colnames <- colnames(df_Fig2) %>% as.data.frame()

#While column names in the next step are manually reassigned I did include a check that the number of columns is correct prior to colname reassignment. This should, in the vast majority of cases,
#throw an error if a different column set is brought in.
if(length(colnames(df_Fig2) == 104)){
        colnames(df_Fig2) <- c("ColV (48/51)",
                               "Col(BS512) (1/51)",
                               "Col(MG828) (17/51)",
                               "Col156 (11/51)",
                               "Col440I (1/51)",
                               "Col8282 (1/51)",
                               "ColpVC (3/51)",
                               "ColRNAI (2/51)",
                               "IncB/O/K/Z (5/51)",
                               "IncFIA (2/51)",
                               "IncFIB (50/51)",
                               "IncFIC (6/51)",
                               "IncFII (8/51)",
                               "IncFII (pUTI89) (1/51)",
                               "IncHI2 (1/51)",
                               "IncHI2A (1/51)",
                               "IncI Gamma (1/51)",
                               "IncI1 Alpha (2/51)",
                               "p0111 (9/51)",
                               "gyrA(D87N) (1/51)",
                               "gyrA(S83L) (3/51)",
                               "aadA2 (1/51)",
                               "aadA5 (1/51)",
                               "APH(3)`-Ib (10/51)",
                               "APH(3)`-Ia (2/51)",
                               "APH(6)-Id (8/51)",
                               "CTX-M-14 (1/51)",
                               "dfrA12 (1/51)",
                               "dfrA14 (1/51)",
                               "dfrA17 (1/51)",
                               "dfrA5 (6/51)",
                               "intI1 (9/51)",
                               "merA (9/51)",
                               "sul1 (6/51)",
                               "sul2 (10/51)",
                               "blaTEM-1 (8/51)",
                               "blaTEM-30 (1/51)",
                               "terA (1/51)",
                               "tet(A) (4/51)",
                               "tet(B) (47/51)",
                               "aslA (51/51)",
                               "astA (2/51)",
                               "cjrA (1/51)",
                               "cjrB (1/51)",
                               "cjrC (1/51)",
                               "cvaA (47/51)",
                               "cvaB (48/51)",
                               "cvaC (45/51)",
                               "cvi (44/51)",
                               "eitA (AY545598) (3/51)",
                               "eitA (KU695535) (1/51)",
                               "etsA (48/51)",
                               "fdeC (51/51)",
                               "fecA (51/51)",
                               "fepA (51/51)",
                               "fepG (51/51)",
                               "fes (51/51)",
                               "fimH (51/51)",
                               "fyuA (51/51)",
                               "gimB marker (27/51)",
                               "gspC (51/51)",
                               "hek (1/51)",
                               "hylF (45/51)",
                               "iha (1/51)",
                               "ireA (51/51)",
                               "iroE (48/51)",
                               "iroN (48/51)",
                               "irp1 (51/51)",
                               "irp2 (51/51)",
                               "iss (51/51)",
                               "iucD (47/51)",
                               "iutA (47/51)",
                               "kpsD (51/51)",
                               "kpsM (51/51)",
                               "*** kpsMT.II._K2.CP000468.1..APEC.O1.3376659.3378102 (47/51)",
                               "kpsT (45/51)",
                               "neuC_CP007275.1 (33/51)",
                               "ompA (51/51)",
                               "ompT (NC_002695) (32/51)",
                               "ompT (HM210637) (48/51)",
                               "papA (51/51)",
                               "papB (51/51)",
                               "***papC (51/51)",
                               "papC (51/51)",
                               "papD (51/51)",
                               "papE (50/51)",
                               "papF (51/51)",
                               "***papG (51/51)",
                               "papGII (51/51)",
                               "papH (51/51)",
                               "papI (51/51)",
                               "papJ (51/51)",
                               "papK (51/51)",
                               "*** senB (2/51)",
                               "senB (2/51)",
                               "shiD (1/51)",
                               "sitA (14/51)",
                               "tia (50/51)",
                               "traT (49/51)",
                               "tsh (1/51)",
                               "usp (51/51)",
                               "vat (47/51)",
                               "yagV/ecpE (51/51)",
                               "ybtA (51/51)")
}else{message("Check the column names are correct as they are hard coded in this script...")}

#Remove genes flagged for removal
df_Fig2 <- df_Fig2 %>% select(-starts_with("***"))

#Plot out tree
p <- ggtree(tree2) %<+%
        meta +
        geom_tiplab(size = 1.75,
                    align = TRUE,
                    linesize = 0.15)

#Plot figure 3
plasmid_tree <- gheatmap(
        p = p,
        data = df_Fig2,
        font.size = 1.5,
        hjust = 0,
        colnames_position = "top",
        colnames = TRUE,
        colnames_angle = 90,
        offset = 0.0175,
        width = 5.5,
        color = "grey",
        low = 'white',
        high = '#fb8072'
)+
        geom_tiplab(
                aes(image = Flag),
                geom = "image",
                size = 0.0175,
                align = TRUE,
                linetype = NULL,
                offset = 0.007
        ) +
        geom_tiplab(
                aes(image = Revised_Source_Niche_img),
                geom = "image",
                size = 0.0135,
                align = TRUE,
                linetype = NULL,
                offset = 0.0095
        ) +
        geom_tiplab(
                aes(image = Pathogen_img),
                geom = "image",
                size = 0.0135,
                align = TRUE,
                linetype = NULL,
                offset = 0.0115
        ) +
        geom_tiplab(
                aes(label = IncF_RST),
                align = TRUE,
                linetype = NULL,
                offset = 0.0135,
                size = 2
                ) +
        geom_tiplab(
                aes(label = HC5),
                align = TRUE,
                linetype = NULL,
                offset = -0.015,
                size = 2
        ) +
        geom_tiplab(
                aes(label = HC2),
                align = TRUE,
                linetype = NULL,
                offset = -0.020,
                size = 2
        ) +
        ggplot2::ylim(NA, 60) +
        scale_fill_manual(
                aesthetics = c("colour", "fill"),
                values = cols_fig3,
                na.value = 'grey'
        )

#Isolate the legend from figure 3
leg1 <- get_legend(plasmid_tree)
#Plot the legend
plot(leg1)


#Plot the plasmid tree (without legend)
plasmid_tree + theme(legend.position="none")


#### Make SNP-matrix tree ####
df <- read.csv(
        "placeholder/output/overlap_group_snp_sites_dists.csv"
        #"placeholder/output/ST95_all_CGA_snp_sites_dists.csv"
)

df2 <- df 
data.csv <- Metadata

#clean row names
rownames(df)<- df$snp.dists.0.6.3
df <- df %>% select(-snp.dists.0.6.3)

##Get order from phylotree
d <- fortify(tree2)
dd <- subset(d, isTip)
colorder <- dd$label[order(dd$y, decreasing = TRUE)]

colordnum <- match(colorder, colnames(df))
rownum <- match(colorder, row.names(df))

#Reset row and column order from SNP matrix to match our tree
df <- df[, colordnum]
df <- df[rownum, ]

#Define SNP matrix colours
snp_names <- c(" g   250+",
               " f   100-250",
               " e   75-100",
               " d   50-75",
               " c   25-50",
               " b   10-25",
               " a   0-10")


snp_cols_fun <- colorRampPalette(c("white", "blue", "red"))
snp_cols <- snp_cols_fun(7)

names(snp_cols) <- snp_names

#Bind our new colours to our old oens
var_col <- c(cols_fig3, snp_cols)

#Define ranges for binning of SNP distance values
mo250 <- df[] > 250
mo100 <- df[] <= 250 & df[] > 100
les100 <- df[] <= 100 & df[] > 75
les75 <- df[] <= 75 & df[] > 50
les50 <- df[] <= 50 & df[] > 25
les25 <- df[] <= 25 & df[] > 10
les10 <- df[] <= 10

#Reassign values to bins
df[mo250] <- " g   250+"
df[mo100] <- " f   100-250"
df[les100] <- " e   75-100"
df[les75] <- " d   50-75"
df[les50] <- " c   25-50"
df[les25] <- " b   10-25"
df[les10] <- " a   0-10"

#Set rownames
rownames(df) <- colnames(df)

#Plot the SNP matrix
snp_matrix <- gheatmap(
        p,
        df,
        width = 10,
        colnames = T,
        color = rgb(220, 220, 220, max = 255, alpha = 25),
        colnames_angle = 90,
        colnames_offset_y = 0,
        colnames_offset_x = 0,
        offset = 0.025,
        colnames_position = "top",
        font.size = 2,
        hjust = 0
) +
        geom_tiplab(
                aes(image = Flag),
                geom = "image",
                size = 0.0175,
                align = TRUE,
                linetype = NULL,
                offset = 0.01
        ) +
        geom_tiplab(
                aes(image = Revised_Source_Niche_img),
                geom = "image",
                size = 0.0135,
                align = TRUE,
                linetype = NULL,
                offset = 0.0135
        ) +
        geom_tiplab(
                aes(image = Pathogen_img),
                geom = "image",
                size = 0.0135,
                align = TRUE,
                linetype = NULL,
                offset = 0.0165
        ) +
        geom_tiplab(
                aes(label = IncF_RST),
                align = TRUE,
                linetype = NULL,
                offset = 0.0195,
                size = 2,
                hjust = 0
        ) +
        geom_tiplab(
                aes(label = HC5),
                align = TRUE,
                linetype = NULL,
                offset = -0.015,
                size = 2
        ) +
        geom_tiplab(
                aes(label = HC2),
                align = TRUE,
                linetype = NULL,
                offset = -0.020,
                size = 2
        ) +
       ggplot2::ylim(NA, 60) +
        scale_fill_manual(
                aesthetics = c("colour", "fill"),
                values = var_col,
                na.value = 'grey'
        )

#require(cowplot)
snp_matrix_theme <- get_legend(snp_matrix)

#Plot below, unhash to visualise and save
#snp_matrix + theme(legend.position="none") 