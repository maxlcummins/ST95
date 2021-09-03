library(ggtree)
library(ape)
library(ggplot2)

if(!exists("plasmid_coverage_table")){
        source("scripts/Figure_2_fastbaps.R")
}

#Read in core genome tree
phylotree <-
        read.tree(
                "analysis/output/ST95_all_CGA_snp_sites.tree"
                #"analysis/output/accessory_ST95_all.tree"
        )

#Midpoint root the tree
phylotree <- midpoint.root(phylotree)

#Read in SNP data
df <-
        read.csv(
                #"analysis/output/overlap_group_snp_sites_dists.csv"
                "analysis/output/ST95_all_CGA_snp_sites_dists.csv"
        )

#Rename Metadata sheet for some reason
data.csv <- Metadata

#clean row names
rownames(df)<- df$snp.dists.0.6.3

#Remove the column with sample names
df <- df %>% select(-snp.dists.0.6.3)

##Extract tip order from tree
d <- fortify(phylotree)
dd <- subset(d, isTip)
colorder <- dd$label[order(dd$y, decreasing = TRUE)]

colordnum <- match(colorder, colnames(df))
rownum <- match(colorder, row.names(df))

#Change order of SNP matrix rows and columns to match the tree
df <- df[, colordnum]
df <- df[rownum, ]

#Set bins
snp_names <- c(" g   250+",
               " f   100-250",
               " e   75-100",
               " d   50-75",
               " c   25-50",
               " b   10-25",
               " a   0-10")

#Define colours
snp_cols_fun <- colorRampPalette(c("white", "blue", "red"))
snp_cols <- snp_cols_fun(7)
        
names(snp_cols) <- snp_names

#Bind our colours to the others we need
var_col <- c(var_col, snp_cols)
        

#Define bins
mo250 <- df[] > 250
mo100 <- df[] <= 250 & df[] > 100
les100 <- df[] <= 100 & df[] > 75
les75 <- df[] <= 75 & df[] > 50
les50 <- df[] <= 50 & df[] > 25
les25 <- df[] <= 25 & df[] > 10
les10 <- df[] <= 10

#Bin the data
df[mo250] <- " g   250+"
df[mo100] <- " f   100-250"
df[les100] <- " e   75-100"
df[les75] <- " d   50-75"
df[les50] <- " c   25-50"
df[les25] <- " b   10-25"
df[les10] <- " a   0-10"

#Paste together the fim, O and H data
metadata$cat_fimH_OH <- paste(metadata$fimH_type, metadata$O_type, metadata$H_type, sep = ":")

#Set the names for our colours (after trimming off a prefix from an earlier step)
names(var_col) <- gsub("Source_","",names(var_col))

#Set names via rownames attribute
fig1_df$working_name <- rownames(fig1_df)

#Select columns with plasmid mapping data
plasmid_map <- fig1_df %>% select(working_name, plasmid_map)

#Bind plasmid mapping data to the metadata/data table
metadata <- left_join(metadata, plasmid_map)

#Plot our tree
ggsim <- ggtree(phylotree, branch.length = "none") %<+% metadata  +
        geom_tiplab(hjust = 0,
                    size = 0.5,
                    offset = 4,
                    aes(color = ColV_pos_neg)
        )  + geom_tippoint(aes(color = plasmid_map), size = 0.05
                           )
rownames(df) <- colnames(df)

#Plot our snp matrix/tree
snp_matrix <- gheatmap(
        ggsim,
        df,
        width = 50,
        colnames = T,
        color = rgb(220, 220, 220, max = 255, alpha = 25),
        colnames_angle = 90,
        colnames_offset_y = 5,
        colnames_offset_x = 0,
        offset = 45,
        colnames_position = "top",
        font.size = 0.5,
        hjust = 0
) +
        geom_tiplab(
                aes(label = HC2),
                align = TRUE,
                linetype = NULL,
                #Accessory
                #offset = 3.3,
                #CGA_snp_sites
                offset = -50,
                size = 0.5
        ) +
        geom_tiplab(
                aes(label = HC5),
                align = TRUE,
                linetype = NULL,
                #Accessory
                #offset = 3.3,
                #CGA_snp_sites
                offset = -75,
                size = 0.5
        ) +
        geom_tiplab(
                aes(label = Pathogen),
                align = TRUE,
                linetype = NULL,
                #Accessory
                #offset = 3.3,
                #CGA_snp_sites
                offset = -100,
                size = 0.5
        ) +
        geom_tiplab(
                aes(label = Revised_Source_Niche, colour = Revised_Source_Niche),
                align = TRUE,
                linetype = NULL,
                #Accessory
                #offset = 3.3,
                #CGA_snp_sites
                offset = -125,
                size = 0.5
        ) +        geom_tiplab(
                aes(label = cat_fimH_OH),
                align = TRUE,
                linetype = NULL,
                #Accessory
                #offset = 3.3,
                #CGA_snp_sites
                offset = -176,
                size = 0.5
        )+
        ylim(NA, 700) +
        scale_fill_manual(
                aesthetics = c("colour", "fill"),
                values = var_col,
                na.value = 'grey'
        )

#Extract the legend
require(cowplot)

legend_snp_matrix <- gheatmap(
        ggsim,
        df,
        width = 50,
        colnames = T,
        color = rgb(220, 220, 220, max = 255, alpha = 25),
        colnames_angle = 90,
        colnames_offset_y = 5,
        colnames_offset_x = 0,
        offset = 45,
        colnames_position = "top",
        font.size = 0.5,
        hjust = 0
)+
        scale_fill_manual(
                aesthetics = c("colour", "fill"),
                values = var_col,
                na.value = 'grey'
        )
snp_matrix_theme <- get_legend(legend_snp_matrix)

#Plot the legend
plot(snp_matrix_theme)

#Define variables for SNP matrix
strip_fontsize <- 5
strip_offset <- -200
strip_text_offset <- -25
strip_barsize <- 1
strip_text_hjust <- 0.5

#Plot snp Matrix
snp_matrix + geom_strip(
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


#Plot the legend
plot(snp_matrix_theme)

