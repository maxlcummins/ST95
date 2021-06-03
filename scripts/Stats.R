
if(!exists("plasmid_coverage_table")){
        source("scripts/Publication/New_plas_Script_1_Combine_and_postprocess.R")
}

ColV_vs_AMR <- Metadata  %>% select(class_res_counts, `ColV (328/668)`, working_name, plas_counts, ESBL_ress, `intI1 (98/668)`) %>% rename(ColV = `ColV (328/668)`)

plasmids <- df_small2 %>% select(starts_with("p"), intI1)

plasmids$working_name <- rownames(df_small4)

ColV_vs_AMR <- left_join(ColV_vs_AMR, plasmids, by = "working_name")

IncF_RSTs <- pMLST %>% select(name, IncF_RST) %>% rename("working_name" = "name")

ColV_vs_AMR <- left_join(ColV_vs_AMR, IncF_RSTs, by = "working_name")

ColV_vs_AMR$pU1_F51_B10[ColV_vs_AMR$pU1_F51_B10 < 80] <- 0
ColV_vs_AMR$pU1_F51_B10[ColV_vs_AMR$pU1_F51_B10 > 1] <- 1

ColV_vs_AMR$pUTI89[ColV_vs_AMR$pUTI89 < 80] <- 0
ColV_vs_AMR$pUTI89[ColV_vs_AMR$pUTI89 > 1] <- 1

ColV_vs_AMR$pACN001_B[ColV_vs_AMR$pACN001_B < 80] <- 0
ColV_vs_AMR$pACN001_B[ColV_vs_AMR$pACN001_B > 1] <- 1

ColV_vs_AMR$pBCE049_1[ColV_vs_AMR$pBCE049_1 < 80] <- 0
ColV_vs_AMR$pBCE049_1[ColV_vs_AMR$pBCE049_1 > 1] <- 1

ColV_vs_AMR$pSF_088_nores[ColV_vs_AMR$pSF_088_nores < 80] <- 0
ColV_vs_AMR$pSF_088_nores[ColV_vs_AMR$pSF_088_nores > 1] <- 1

ColV_vs_AMR$pEC244_2[ColV_vs_AMR$pEC244_2 < 80] <- 0
ColV_vs_AMR$pEC244_2[ColV_vs_AMR$pEC244_2 > 1] <- 1

ColV_vs_AMR$F_plasmid <- ColV_vs_AMR$IncF_RST

ColV_vs_AMR$F_plasmid <- gsub("^[A-Z]-.*[A-Z]-$","0", ColV_vs_AMR$F_plasmid)
ColV_vs_AMR$F_plasmid <- gsub(".*:.*","1", ColV_vs_AMR$F_plasmid)

ColV_vs_AMR <- ColV_vs_AMR %>% rename(plas_counts = "plas_counts.x")

ColV_vs_AMR$class_res_counts <- gsub("classes_res_","",ColV_vs_AMR$class_res_counts)
ColV_vs_AMR$plas_counts <- gsub("count_plas_","",ColV_vs_AMR$plas_counts)

ColV_vs_AMR$ColV <- gsub("ColV_","ColV ",ColV_vs_AMR$ColV)
ColV_vs_AMR$ColV <- gsub("pos","positive",ColV_vs_AMR$ColV)
ColV_vs_AMR$ColV <- gsub("neg","negative",ColV_vs_AMR$ColV)

ColV_vs_AMR$ESBL_ress <- gsub("ESBL_neg","0",ColV_vs_AMR$ESBL_ress)
ColV_vs_AMR$ESBL_ress <- gsub("ESBL_pos","1",ColV_vs_AMR$ESBL_ress)

ColV_vs_AMR$class_res_counts <- as.numeric(ColV_vs_AMR$class_res_counts)
ColV_vs_AMR$plas_counts <- as.numeric(ColV_vs_AMR$plas_counts)

ColV_vs_AMR$plas_counts[ColV_vs_AMR$plas_counts > 1] <- 1

ColV_vs_AMR <- ColV_vs_AMR %>% filter(F_plasmid == 1)


ggplot(data = ColV_vs_AMR,
       aes(x = ColV, y = class_res_counts)) +
        geom_boxplot(fill = "steelblue") +
        labs(y = "Number of classes of AMR determinants",
             x = "ColV Carriage",
             title = "Breadth of antimicrobial resistance versus ColV status")


ColV_vs_AMR$MDR <- ColV_vs_AMR$class_res_counts

ColV_vs_AMR$MDR[ColV_vs_AMR$MDR < 3] <- 0
ColV_vs_AMR$MDR[ColV_vs_AMR$MDR >= 3] <- 1


table_trait_group<- ColV_vs_AMR %>% filter(grepl("[A-Z]", IncF_RST)) %>% group_by(plas_counts, MDR) %>% summarise(counts = n())

plas_counts_mdr<- cbind(c(table_trait_group[1,3]
                          ,table_trait_group[2,3]),
                        c(table_trait_group[3,3],
                          table_trait_group[4,3])
)

plas_counts_mdr <- matrix(unlist(plas_counts_mdr), 2)



#### Check association between plasmid replicons and MDR ####


table_trait_group<- ColV_vs_AMR %>% filter(grepl("[A-Z]", IncF_RST)) %>% group_by(pUTI89, MDR) %>% summarise(counts = n())

pUTI89_mdr<- cbind(c(table_trait_group[1,3]
                     ,table_trait_group[2,3]),
                   c(table_trait_group[3,3],
                     table_trait_group[4,3])
)

pUTI89_mdr <- matrix(unlist(pUTI89_mdr), 2)

table_trait_group<- ColV_vs_AMR %>% filter(grepl("[A-Z]", IncF_RST)) %>% group_by(pSF_088_nores, MDR) %>% summarise(counts = n())

pSF_088_nores_mdr<- cbind(c(table_trait_group[1,3]
                            ,table_trait_group[2,3]),
                          c(table_trait_group[3,3],
                            table_trait_group[4,3])
)

pSF_088_nores_mdr <- matrix(unlist(pSF_088_nores_mdr), 2)

table_trait_group<- ColV_vs_AMR %>% filter(grepl("[A-Z]", IncF_RST)) %>% group_by(pACN001_B, MDR) %>% summarise(counts = n())

pACN001_B_mdr<- cbind(c(table_trait_group[1,3]
                        ,table_trait_group[2,3]),
                      c(table_trait_group[3,3],
                        table_trait_group[4,3])
)

pACN001_B_mdr <- matrix(unlist(pACN001_B_mdr), 2)

table_trait_group<- ColV_vs_AMR %>% filter(grepl("[A-Z]", IncF_RST)) %>% group_by(pBCE049_1, MDR) %>% summarise(counts = n())

pBCE049_mdr<- cbind(c(table_trait_group[1,3]
                      ,table_trait_group[2,3]),
                    c(table_trait_group[3,3],
                      table_trait_group[4,3])
)

pBCE049_mdr <- matrix(unlist(pBCE049_mdr), 2)

table_trait_group<- ColV_vs_AMR %>% filter(grepl("[A-Z]", IncF_RST)) %>% group_by(pU1_F51_B10, MDR) %>% summarise(counts = n())

pU1_F51_B10_mdr <- cbind(c(table_trait_group[1,3]
                           ,table_trait_group[2,3]),
                         c(table_trait_group[3,3],
                           table_trait_group[4,3])
)

pU1_F51_B10_mdr <- matrix(unlist(pU1_F51_B10_mdr), 2)



table_trait_group<- ColV_vs_AMR %>% filter(grepl("[A-Z]", IncF_RST)) %>% group_by(pEC244_2, MDR) %>% summarise(counts = n())

pEC244_2_mdr <- cbind(c(table_trait_group[1,3]
                        ,table_trait_group[2,3]),
                      c(table_trait_group[3,3],
                        table_trait_group[4,3])
)

pEC244_2_mdr <- matrix(unlist(pEC244_2_mdr), 2)




table_trait_group<- ColV_vs_AMR %>% filter(grepl("[A-Z]", IncF_RST)) %>% group_by(intI1, MDR) %>% summarise(counts = n())


intI1_mdr<- cbind(c(table_trait_group[1,3]
                    ,table_trait_group[2,3]),
                  c(table_trait_group[3,3],
                    table_trait_group[4,3])
)

intI1_mdr <- matrix(unlist(intI1_mdr), 2)

table_trait_group<- ColV_vs_AMR %>% filter(grepl("[A-Z]", IncF_RST)) %>% group_by(ColV, MDR) %>% summarise(counts = n())

colv_mdr<- cbind(c(table_trait_group[1,3]
                   ,table_trait_group[2,3]),
                 c(table_trait_group[3,3],
                   table_trait_group[4,3]))

colv_mdr <- matrix(unlist(colv_mdr), 2)


fisher.test(pSF_088_nores_mdr)
fisher.test(pUTI89_mdr)
fisher.test(pU1_F51_B10_mdr)
fisher.test(pACN001_B_mdr)
fisher.test(pBCE049_mdr)
fisher.test(pEC244_2_mdr)

fisher.test(intI1_mdr)
fisher.test(colv_mdr)
fisher.test(plas_counts_mdr)

fisher.test(colv_mdr)


#### Check association between plasmid replicons and ESBL carriage ####

table_trait_group<- ColV_vs_AMR %>% filter(grepl("[A-Z]", IncF_RST)) %>% group_by(plas_counts, ESBL_ress) %>% summarise(counts = n())

plas_counts_ESBL_ress<- cbind(c(table_trait_group[1,3]
                                ,table_trait_group[2,3]),
                              c(table_trait_group[3,3],
                                table_trait_group[4,3])
)

plas_counts_ESBL_ress <- matrix(unlist(plas_counts_ESBL_ress), 2)

table_trait_group<- ColV_vs_AMR %>% filter(grepl("[A-Z]", IncF_RST)) %>% group_by(pUTI89, ESBL_ress) %>% summarise(counts = n())

pUTI89_esbl<- cbind(c(table_trait_group[1,3]
                      ,table_trait_group[2,3]),
                    c(table_trait_group[3,3],
                      table_trait_group[4,3])
)

pUTI89_esbl <- matrix(unlist(pUTI89_esbl), 2)

table_trait_group<- ColV_vs_AMR %>% filter(grepl("[A-Z]", IncF_RST)) %>% group_by(pSF_088_nores, ESBL_ress) %>% summarise(counts = n())

pSF_088_nores_esbl<- cbind(c(table_trait_group[1,3]
                             ,table_trait_group[2,3]),
                           c(table_trait_group[3,3],
                             table_trait_group[4,3])
)

pSF_088_nores_esbl <- matrix(unlist(pSF_088_nores_esbl), 2)

table_trait_group<- ColV_vs_AMR %>% filter(grepl("[A-Z]", IncF_RST)) %>% group_by(pACN001_B, ESBL_ress) %>% summarise(counts = n())

pACN001_B_esbl<- cbind(c(table_trait_group[1,3]
                         ,table_trait_group[2,3]),
                       c(table_trait_group[3,3],
                         table_trait_group[4,3])
)

pACN001_B_esbl <- matrix(unlist(pACN001_B_esbl), 2)

table_trait_group<- ColV_vs_AMR %>% filter(grepl("[A-Z]", IncF_RST)) %>% group_by(pBCE049_1, ESBL_ress) %>% summarise(counts = n())

pBCE049_esbl<- cbind(c(table_trait_group[1,3]
                       ,table_trait_group[2,3]),
                     c(table_trait_group[3,3],
                       table_trait_group[4,3])
)

pBCE049_esbl <- matrix(unlist(pBCE049_esbl), 2)

table_trait_group<- ColV_vs_AMR %>% filter(grepl("[A-Z]", IncF_RST)) %>% group_by(pU1_F51_B10, ESBL_ress) %>% summarise(counts = n())

pU1_F51_B10_esbl <- cbind(c(table_trait_group[1,3]
                            ,table_trait_group[2,3]),
                          c(table_trait_group[3,3],
                            table_trait_group[4,3])
)

pU1_F51_B10_esbl <- matrix(unlist(pU1_F51_B10_esbl), 2)




table_trait_group<- ColV_vs_AMR %>% filter(grepl("[A-Z]", IncF_RST)) %>% group_by(pEC244_2, ESBL_ress) %>% summarise(counts = n())

pEC244_2_esbl <- cbind(c(table_trait_group[1,3]
                         ,table_trait_group[2,3]),
                       c(table_trait_group[3,3],
                         table_trait_group[4,3])
)

pEC244_2_esbl <- matrix(unlist(pEC244_2_esbl), 2)






table_trait_group<- ColV_vs_AMR %>% filter(grepl("[A-Z]", IncF_RST)) %>% group_by(intI1, ESBL_ress) %>% summarise(counts = n())

intI1_esbl<- cbind(c(table_trait_group[1,3]
                     ,table_trait_group[2,3]),
                   c(table_trait_group[3,3],
                     table_trait_group[4,3])
)

intI1_esbl <- matrix(unlist(intI1_esbl), 2)

table_trait_group<- ColV_vs_AMR %>% filter(grepl("[A-Z]", IncF_RST)) %>% group_by(ColV, ESBL_ress) %>% summarise(counts = n())

colv_esbl<- cbind(c(table_trait_group[1,3]
                    ,table_trait_group[2,3]),
                  c(table_trait_group[3,3],
                    table_trait_group[4,3]))

colv_esbl <- matrix(unlist(colv_esbl), 2)


fisher.test(pSF_088_nores_esbl)
fisher.test(pUTI89_esbl)
fisher.test(pU1_F51_B10_esbl)
fisher.test(pACN001_B_esbl)
fisher.test(pBCE049_esbl)
fisher.test(pEC244_2_esbl)

fisher.test(intI1_esbl)
fisher.test(colv_esbl)
fisher.test(plas_counts_ESBL_ress)

fisher.test(colv_esbl)


#### Check association between plasmid replicons and intI1 carriage ####

table_trait_group<- ColV_vs_AMR %>% filter(grepl("[A-Z]", IncF_RST)) %>% group_by(plas_counts, intI1) %>% summarise(counts = n())

plas_counts_intI1<- cbind(c(table_trait_group[1,3]
                            ,table_trait_group[2,3]),
                          c(table_trait_group[3,3],
                            table_trait_group[4,3])
)

plas_counts_intI1 <- matrix(unlist(plas_counts_intI1), 2)

table_trait_group<- ColV_vs_AMR %>% filter(grepl("[A-Z]", IncF_RST)) %>% group_by(pUTI89, intI1) %>% summarise(counts = n())

pUTI89_intI1<- cbind(c(table_trait_group[1,3]
                       ,table_trait_group[2,3]),
                     c(table_trait_group[3,3],
                       table_trait_group[4,3])
)

pUTI89_intI1 <- matrix(unlist(pUTI89_intI1), 2)

table_trait_group<- ColV_vs_AMR %>% filter(grepl("[A-Z]", IncF_RST)) %>% group_by(pSF_088_nores, intI1) %>% summarise(counts = n())

pSF_088_nores_intI1<- cbind(c(table_trait_group[1,3]
                              ,table_trait_group[2,3]),
                            c(table_trait_group[3,3],
                              table_trait_group[4,3])
)

pSF_088_nores_intI1 <- matrix(unlist(pSF_088_nores_intI1), 2)

table_trait_group<- ColV_vs_AMR %>% filter(grepl("[A-Z]", IncF_RST)) %>% group_by(pACN001_B, intI1) %>% summarise(counts = n())

pACN001_B_intI1<- cbind(c(table_trait_group[1,3]
                          ,table_trait_group[2,3]),
                        c(table_trait_group[3,3],
                          table_trait_group[4,3])
)

pACN001_B_intI1 <- matrix(unlist(pACN001_B_intI1), 2)

table_trait_group<- ColV_vs_AMR %>% filter(grepl("[A-Z]", IncF_RST)) %>% group_by(pBCE049_1, intI1) %>% summarise(counts = n())

pBCE049_intI1<- cbind(c(table_trait_group[1,3]
                        ,table_trait_group[2,3]),
                      c(table_trait_group[3,3],
                        table_trait_group[4,3])
)

pBCE049_intI1 <- matrix(unlist(pBCE049_intI1), 2)

table_trait_group<- ColV_vs_AMR %>% filter(grepl("[A-Z]", IncF_RST)) %>% group_by(pU1_F51_B10, intI1) %>% summarise(counts = n())

pU1_F51_B10_intI1 <- cbind(c(table_trait_group[1,3]
                             ,table_trait_group[2,3]),
                           c(table_trait_group[3,3],
                             table_trait_group[4,3])
)

pU1_F51_B10_intI1 <- matrix(unlist(pU1_F51_B10_intI1), 2)

table_trait_group<- ColV_vs_AMR %>% filter(grepl("[A-Z]", IncF_RST)) %>% group_by(ColV, intI1) %>% summarise(counts = n())

colv_intI1<- cbind(c(table_trait_group[1,3]
                     ,table_trait_group[2,3]),
                   c(table_trait_group[3,3],
                     table_trait_group[4,3]))

colv_intI1 <- matrix(unlist(colv_intI1), 2)



table_trait_group<- ColV_vs_AMR %>% filter(grepl("[A-Z]", IncF_RST)) %>% group_by(pEC244_2, intI1) %>% summarise(counts = n())

pEC244_2_intI1 <- cbind(c(table_trait_group[1,3]
                          ,table_trait_group[2,3]),
                        c(table_trait_group[3,3],
                          table_trait_group[4,3])
)

pEC244_2_intI1 <- matrix(unlist(pEC244_2_intI1), 2)


fisher.test(pSF_088_nores_intI1)
fisher.test(pUTI89_intI1)
fisher.test(pU1_F51_B10_intI1)
fisher.test(pACN001_B_intI1)
fisher.test(pBCE049_intI1)
fisher.test(pEC244_2_intI1)

fisher.test(colv_intI1)
fisher.test(plas_counts_intI1)

#### Check association between F plasmids and non-F plasmid carriage ####
table_trait_group<- ColV_vs_AMR %>% filter(grepl("[A-Z]", IncF_RST)) %>% group_by(pUTI89, plas_counts) %>% summarise(counts = n())

pUTI89_plas_counts<- cbind(c(table_trait_group[1,3]
                             ,table_trait_group[2,3]),
                           c(table_trait_group[3,3],
                             table_trait_group[4,3])
)

pUTI89_plas_counts <- matrix(unlist(pUTI89_plas_counts), 2)

table_trait_group<- ColV_vs_AMR %>% filter(grepl("[A-Z]", IncF_RST)) %>% group_by(pSF_088_nores, plas_counts) %>% summarise(counts = n())

pSF_088_nores_plas_counts<- cbind(c(table_trait_group[1,3]
                                    ,table_trait_group[2,3]),
                                  c(table_trait_group[3,3],
                                    table_trait_group[4,3])
)

pSF_088_nores_plas_counts <- matrix(unlist(pSF_088_nores_plas_counts), 2)

table_trait_group<- ColV_vs_AMR %>% filter(grepl("[A-Z]", IncF_RST)) %>% group_by(pACN001_B, plas_counts) %>% summarise(counts = n())

pACN001_B_plas_counts<- cbind(c(table_trait_group[1,3]
                                ,table_trait_group[2,3]),
                              c(table_trait_group[3,3],
                                table_trait_group[4,3])
)

pACN001_B_plas_counts <- matrix(unlist(pACN001_B_plas_counts), 2)

table_trait_group<- ColV_vs_AMR %>% filter(grepl("[A-Z]", IncF_RST)) %>% group_by(pBCE049_1, plas_counts) %>% summarise(counts = n())

pBCE049_plas_counts<- cbind(c(table_trait_group[1,3]
                              ,table_trait_group[2,3]),
                            c(table_trait_group[3,3],
                              table_trait_group[4,3])
)

pBCE049_plas_counts <- matrix(unlist(pBCE049_plas_counts), 2)

table_trait_group<- ColV_vs_AMR %>% filter(grepl("[A-Z]", IncF_RST)) %>% group_by(pU1_F51_B10, plas_counts) %>% summarise(counts = n())

pU1_F51_B10_plas_counts <- cbind(c(table_trait_group[1,3]
                                   ,table_trait_group[2,3]),
                                 c(table_trait_group[3,3],
                                   table_trait_group[4,3])
)

pU1_F51_B10_plas_counts <- matrix(unlist(pU1_F51_B10_plas_counts), 2)

table_trait_group<- ColV_vs_AMR %>% filter(grepl("[A-Z]", IncF_RST)) %>% group_by(ColV, plas_counts) %>% summarise(counts = n())

colv_plas_counts<- cbind(c(table_trait_group[1,3]
                           ,table_trait_group[2,3]),
                         c(table_trait_group[3,3],
                           table_trait_group[4,3]))

colv_plas_counts <- matrix(unlist(colv_plas_counts), 2)


table_trait_group<- ColV_vs_AMR %>% filter(grepl("[A-Z]", IncF_RST)) %>% group_by(pEC244_2, plas_counts) %>% summarise(counts = n())

pEC244_2_plas_counts <- cbind(c(table_trait_group[1,3]
                                ,table_trait_group[2,3]),
                              c(table_trait_group[3,3],
                                table_trait_group[4,3])
)

pEC244_2_plas_counts <- matrix(unlist(pEC244_2_plas_counts), 2)



fisher.test(pSF_088_nores_plas_counts)
fisher.test(pUTI89_plas_counts)
fisher.test(pU1_F51_B10_plas_counts)
fisher.test(pACN001_B_plas_counts)
fisher.test(pBCE049_plas_counts)
fisher.test(pEC244_2_plas_counts)

fisher.test(colv_plas_counts)


#### Check association between IncF_RSTs and fluoroquinolone resistance ####





