library(ggtree)

if(!exists("metadata")){
        source("scripts/Clade_H2_process.R")
}

if(nrow(metadata) != 51){
        quit()
}

#Remove all objects from our environment except the ones we need
rm(list=setdiff(ls(), c("df", "metadata")))

tree_path_2 <- "placeholder/output/overlap_group_snp_sites.treefile"
#tree_path_2 <- "placeholder/output/accessory_overlap_group.tree"

cols_fig3 <- c("0" = "white",  "1" = "#bebada", "2" = "#80b1d3", "3" = "#fb8072", "4" = "black", "5" = "purple")

df_Fig2 <- df %>% select(-starts_with("colv_zoetis")) %>% select(-starts_with("IS_"))

df_plas <- df %>% select(-starts_with("colv")) %>% select(starts_with("Col"), starts_with("Inc"), starts_with("p0111"),starts_with("rep"))

df_IS <- df %>% select(starts_with("IS_"))

df_noplas <- df %>% select(-starts_with("IS_")) %>% select(-starts_with("colv")) %>% select(-starts_with("Col"), -starts_with("Inc"), -starts_with("p0111"),-starts_with("rep"))

tree_path_2

tree2 <- read.tree(tree_path_2)

df_Fig2_colnames <- colnames(df_Fig2) %>% as.data.frame()

if(length(colnames(df_Fig2) == 104)){
        colnames(df_Fig2) <- c("ColV (48/51)",
                               "Col(BS512) (1/51)",
                               "Col(MG828) (17/51)",
                               "Col156 (11/51)",
                               "Col440I (1/51)",
                               "Col8282 (1/51)",
                               "ColpVC (3/51)",
                               "ColRNAI (2/51)",
                               "IncB/O/K/Z",
                               "IncFIA (2/51)",
                               "IncFIB (50/51)",
                               "IncFIC (6/51)",
                               "IncFII (8/51)",
                               "IncFII (pUTI89) (1/51)",
                               "IncHI2 (1/51)",
                               "IncHI2A (1/51)",
                               "IncI Gamma  (1/51)",
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

df_Fig2 <- df_Fig2 %>% select(-starts_with("***"))

p <- ggtree(tree2) %<+%
        metadata +
        geom_tiplab(size = 1.75,
                    align = TRUE,
                    linesize = 0.15)

listy <- read_csv("placeholder/delims/pSF_088_nores_plasmid_coverage.csv")

rownames(listy) <- listy$X1

listy2 <- listy %>% select(-X1)

rownames(listy2) <- listy$X1

plasmid_tree <- gheatmap(
        p = p,
        data = df_Fig2,
        #colnames_offset_y = -0.1,
        font.size = 1.5,
        hjust = 0,
        colnames_position = "top",
        colnames = TRUE,
        colnames_angle = 90,
        #Accessory
        #offset = 0.65,
        #CGA_snps
        offset = 0.0175,
        #Accessory
        #width = 3.5,
        #CGA_snps
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
                #Accessory
                #offset = 0.35
                #CGA_snp_sites
                offset = 0.007
        ) +
        geom_tiplab(
                aes(image = Revised_Source_Niche_img),
                geom = "image",
                size = 0.0135,
                align = TRUE,
                linetype = NULL,
                #Accessory
                #offset = 0.45
                #CGA_snp_sites
                offset = 0.0095
        ) +
        geom_tiplab(
                aes(image = Pathogen_img),
                geom = "image",
                size = 0.0135,
                align = TRUE,
                linetype = NULL,
                #Accessory
                #offset = 0.55
                #CGA_snp_sites
                offset = 0.0115
        ) +
        #        geom_tiplab(
        #                aes(label = Revised_Collection_Year),
        #                align = TRUE,
        #                linetype = NULL,
        #                offset = 19,
        #                size = 0.5
        #        ) +
        #geom_tiplab(
        #        aes(label = ColV_percent_hit),
        #        align = TRUE,
        #        linetype = NULL,
#        offset = 0.018,
#        size = 0.5
#) +
geom_tiplab(
        aes(label = IncF_RST),
        align = TRUE,
        linetype = NULL,
        #Accessory
        #offset = 3.3,
        #CGA_snp_sites
        offset = 0.0135,
        size = 2
) +
        geom_tiplab(
                aes(label = HC5),
                align = TRUE,
                linetype = NULL,
                #Accessory
                #offset = 3.3,
                #CGA_snp_sites
                offset = -0.015,
                size = 2
        ) +
        geom_tiplab(
                aes(label = HC2),
                align = TRUE,
                linetype = NULL,
                #Accessory
                #offset = 3.3,
                #CGA_snp_sites
                offset = -0.020,
                size = 2
        ) +
#        theme(legend.position = "none") #+
ggplot2::ylim(NA, 60) +
        scale_fill_manual(
                aesthetics = c("colour", "fill"),
                values = cols_fig3,
                na.value = 'grey'
        )# +
#geom_image(x = 104.5, y = 48, image = "../Linear_Map.png", size = 0.630) 

require(cowplot)
leg1 <- get_legend(plasmid_tree)

plot(leg1)

plasmid_tree + theme(legend.position="none")
