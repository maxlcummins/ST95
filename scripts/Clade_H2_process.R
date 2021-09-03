library(ggtree)
library(phytools)
library(tidyr)
library(readr)
library(dplyr)

#Set working directory
working_dir <- "."

#Load paths to files needed for the script
tree_path <- "analysis/output/accessory_overlap_group.tree"
abricate_path <- "analysis/output/genotype.txt"
pointfinder_path <- "analysis/output/pointfinder.txt"
pMLST_data <- "analysis/output/pMLST.txt"
cgMLST_path <- "analysis/data/enterobase/Curated_metadata_all.txt"
phylotype_path <- "analysis/data/enterobase/ST95_phylotypes.txt"
serotype_path <- "analysis/data/enterobase/ST95_serotypes.txt"

#Provide output names
output_name <- "ST95_all"
output_dir <- "processed"

#Set working directory
setwd(working_dir)

#Load abricateR2 script
source("scripts/abricateR2.R")

#Run abricateR
abricateR(
        abricate_in = abricate_path,
        output = output_name,
        output_directory = output_dir,
        writecsv = TRUE,
        pointfinder_data = pointfinder_path,
        pMLST_data = pMLST_data
)

#### Read in and clean tree file ####
tree <-
        read.tree(file = tree_path)

#Trim the names of the assemblies in the tree tip labels
#Note: This deletes everything after a "."
tree$tip.label <- gsub("\\..*", "", tree$tip.label)

#Reassign Reference tip label to the name of the reference, if provided.
#This won't do anything if there isn't a strain called 'Reference' in the tree.
if(exists("refname")){
        tree$tip.label <- gsub("Reference", refname, tree$tip.label)
}

#### Read in and process cgMLST, O/H type and fimH type data ####
#Read in cgMLST data
#Note: This also contains manually curated metadata
Metadata <- read_delim(cgMLST_path,
                       "\t",
                       escape_double = FALSE,
                       trim_ws = TRUE)

##Designate strains as HC50-1106 or Other
# Note: This was a line used to classify as ingroups or outgroups. HC50:1106
#       was originally a group of interest, however this has changed since
# Flag: delete
Metadata$HC50_or_other <-
        gsub("^1106$*", "HC50-1106", Metadata$HC50)
Metadata$HC50_or_other <-
        gsub("^[0-9].*", "Other", Metadata$HC50_or_other)

## Designate strains as being "Other" for infrequent HC groups
# The following lines assign HC groups as  "Other" if they are too 'infrequent'.
# Note: What was defined as "infrequent" was chosen based on how many levels of
#       each factor there were upon manual inspection of groups. To inspect this
#       you can do something like the following, changing HCX to your desired
#       HC group:
#       Metadata %>% 
#               group_by(HCX) %>% 
#               summarise(counts = n()) %>% 
#               arrange(desc(counts))

#HC200
#Identify HC200 groups with < 10 representatives
HC200_others <- Metadata$HC200 %>% table() %>% as.data.frame() %>% dplyr::filter(Freq < 10)
#Duplicate HC200 columnto be our new HC200_other column
Metadata$HC200_Other <- Metadata$HC200
#Replace HC200 groups with < 10 representatives with "Other"
Metadata$HC200_Other[Metadata$HC200_Other %in% HC200_others$.] <- "Other"

#HC100
#Identify HC100 groups with < 3 representatives
HC100_others <- Metadata$HC100 %>% table() %>% as.data.frame() %>% filter(Freq < 10)
#Duplicate HC100 columnto be our new HC100_other column
Metadata$HC100_Other <- Metadata$HC100
#Replace HC100 groups with < 3 representatives with "Other"
Metadata$HC100_Other[Metadata$HC100_Other %in% HC100_others$.] <- "Other"

#HC50
#Identify HC50 groups with < 20 representatives
HC50_others <- Metadata$HC50 %>% table() %>% as.data.frame() %>% filter(Freq <= 20)
#Duplicate HC50 column to be our new HC50_other column
Metadata$HC50_Other <- Metadata$HC50
#Replace HC50 groups with < 3 representatives with "Other"
Metadata$HC50_Other[Metadata$HC50_Other %in% HC50_others$.] <- "Other"


#Read in serotype data
serotype <- read_delim(serotype_path,
                       "\t",
                       escape_double = FALSE,
                       trim_ws = TRUE)

#Remove unwanted columns
serotype <- serotype %>% select(Uberstrain, `O Antigen`, `H Antigen`)

#Read in fimH data
phylotype <- read_delim(phylotype_path,
                        "\t",
                        escape_double = FALSE,
                        trim_ws = TRUE)

#Remove unwanted columns
phylotype <- phylotype %>% select(Uberstrain, `fimH (fimTyper)`)

#Combine Metadata/cgMLST and phylotype data
Metadata <- left_join(Metadata, phylotype)

#Combine Metadata/cgMLST/phylotype data and serotype data
Metadata <- left_join(Metadata, serotype)

##Clean our new column names to get rid of illegal characters and improve
#legibility
colnames(Metadata) <- gsub("fimH \\(fimTyper\\)", "fimH_type", colnames(Metadata))
colnames(Metadata) <- gsub(" Antigen", "_type", colnames(Metadata))

#Combine O50/O2 lipopolysaccharide variants to be O2/O50
Metadata$O_type <- gsub("O50 or O2|O2 or O50","O2/O50", Metadata$O_type)

#Replace uncertain O and H types with 'O?'/'H?'
Metadata$O_type <- gsub("uncertain","O?", Metadata$O_type)
Metadata$H_type <- gsub("uncertain","H?", Metadata$H_type)

#Replace uncertain fimH types with 'fimH?'
Metadata$fimH_type <- gsub("\\*","fimH?", Metadata$fimH_type)

#Replace O and H null strains with a different formatting
Metadata$O_type <- gsub("-","O-", Metadata$O_type)
Metadata$H_type <- gsub("-","H-", Metadata$H_type)


#### General manipulation of Metadata ####

#Add new column 'working_name'
Metadata$working_name <- Metadata$`Assembly barcode`

#Replace the assembly barcode with strain name
#Flag: delete
#Metadata$working_name <-
#        gsub("ESC_SA8243AA_AS", refname, Metadata$working_name)

#Filter metadata table to only contain strains from the tree
Metadata <- Metadata %>% filter(working_name %in% tree$tip.label)

#Change the column order to put working_name first
Metadata <- Metadata %>% select(working_name, everything())

#set rowname to be equal to working name
rownames(Metadata) <- Metadata$working_name

#Remove spaces from column names for metadata - this usually causes issues
colnames(Metadata) <- gsub(" ", "_", colnames(Metadata))


##Fix some strain origins which were incorrectly identified by enterobase
# Note: 
Metadata[Metadata$Source_Details == "Parent; Bone marrow" & Metadata$Revised_Source_Niche == "Canine", "Revised_Source_Niche"] <- "Poultry"

#### Adding in of URLs for icons ####
#Generate a new column called Flag within which we will give URLs based on
#countries of origin, used by ggimage to render flags in our image.
Metadata$Flag <- Metadata$Country
Metadata$Flag <-
        gsub(
                "Australia",
                "https://raw.githubusercontent.com/maxlcummins/AVC171/master/images/Flags/Australia.jpg",
                Metadata$Flag
        )
Metadata$Flag <-
        gsub(
                "Denmark",
                "https://raw.githubusercontent.com/maxlcummins/AVC171/master/images/Flags/Denmark.jpg",
                Metadata$Flag
        )
Metadata$Flag <-
        gsub(
                "United States",
                "https://raw.githubusercontent.com/maxlcummins/AVC171/master/images/Flags/USA.png",
                Metadata$Flag
        )
Metadata$Flag <-
        gsub(
                "^Ireland$",
                "https://raw.githubusercontent.com/maxlcummins/AVC171/master/images/Flags/Ireland.png",
                Metadata$Flag
        )
Metadata$Flag <-
        gsub(
                "Vietnam",
                "https://raw.githubusercontent.com/maxlcummins/AVC171/master/images/Flags/Vietnam.png",
                Metadata$Flag
        )
Metadata$Flag <-
        gsub(
                "United Kingdom",
                "https://raw.githubusercontent.com/maxlcummins/AVC171/master/images/Flags/United_Kindgom.png",
                Metadata$Flag
        )
Metadata$Flag <-
        gsub(
                "Mexico",
                "https://raw.githubusercontent.com/maxlcummins/AVC171/master/images/Flags/Mexico.png",
                Metadata$Flag
        )
Metadata$Flag <-
        gsub(
                "Ghana",
                "https://raw.githubusercontent.com/maxlcummins/AVC171/master/images/Flags/Ghana.png",
                Metadata$Flag
        )
Metadata$Flag <-
        gsub(
                "United Arab Emirates",
                "https://raw.githubusercontent.com/maxlcummins/AVC171/master/images/Flags/UAE.png",
                Metadata$Flag
        )
Metadata$Flag <-
        gsub(
                "Scotland",
                "https://raw.githubusercontent.com/maxlcummins/AVC171/master/images/Flags/Scotland.png",
                Metadata$Flag
        )
Metadata$Flag <-
        gsub(
                "Germany",
                "https://raw.githubusercontent.com/maxlcummins/AVC171/master/images/Flags/Germany.png",
                Metadata$Flag
        )
Metadata$Flag <-
        gsub(
                "Northern Ireland",
                "https://raw.githubusercontent.com/maxlcummins/AVC171/master/images/Flags/Northern_Ireland.jpg",
                Metadata$Flag
        )
Metadata$Flag <-
        gsub(
                "Japan",
                "https://raw.githubusercontent.com/maxlcummins/AVC171/master/images/Flags/japan.png",
                Metadata$Flag
        )
Metadata$Flag <-
        gsub(
                "Sweden",
                "https://raw.githubusercontent.com/maxlcummins/AVC171/master/images/Flags/sweden.png",
                Metadata$Flag
        )
Metadata$Flag <-
        gsub(
                "Norway",
                "https://raw.githubusercontent.com/maxlcummins/AVC171/master/images/Flags/norway.png",
                Metadata$Flag
        )
Metadata$Flag <-
        gsub(
                "China",
                "https://raw.githubusercontent.com/maxlcummins/AVC171/master/images/Flags/china.png",
                Metadata$Flag
        )
Metadata$Flag <-
        gsub(
                "Nepal",
                "https://raw.githubusercontent.com/maxlcummins/AVC171/master/images/Flags/nepal.png",
                Metadata$Flag
        )
Metadata$Flag <-
        gsub(
                "Saudi Arabia",
                "https://raw.githubusercontent.com/maxlcummins/AVC171/master/images/Flags/saudi-arabia.png",
                Metadata$Flag
        )
Metadata$Flag <-
        gsub(
                "Netherlands",
                "https://raw.githubusercontent.com/maxlcummins/AVC171/master/images/Flags/netherlands.png",
                Metadata$Flag
        )
Metadata$Flag <-
        gsub(
                "New Zealand",
                "https://raw.githubusercontent.com/maxlcummins/AVC171/master/images/Flags/new-zealand.png",
                Metadata$Flag
        )
Metadata$Flag <-
        gsub(
                "Hungary",
                "https://raw.githubusercontent.com/maxlcummins/AVC171/master/images/Flags/hungary.png",
                Metadata$Flag
        )
Metadata$Flag <-
        gsub(
                "France",
                "https://raw.githubusercontent.com/maxlcummins/AVC171/master/images/Flags/france.png",
                Metadata$Flag
        )
Metadata$Flag <-
        gsub(
                "Singapore",
                "https://raw.githubusercontent.com/maxlcummins/AVC171/master/images/Flags/singapore.png",
                Metadata$Flag
        )
Metadata$Flag <-
        gsub(
                "Croatia",
                "https://raw.githubusercontent.com/maxlcummins/AVC171/master/images/Flags/croatia.png",
                Metadata$Flag
        )
Metadata$Flag <-
        gsub(
                "Sri Lanka",
                "https://raw.githubusercontent.com/maxlcummins/AVC171/master/images/Flags/sri-lanka.png",
                Metadata$Flag
        )
Metadata$Flag <-
        gsub(
                "Finland",
                "https://raw.githubusercontent.com/maxlcummins/AVC171/master/images/Flags/finland.png",
                Metadata$Flag
        )
Metadata$Flag <-
        gsub(
                "Canada",
                "https://raw.githubusercontent.com/maxlcummins/AVC171/master/images/Flags/canada.png",
                Metadata$Flag
        )
Metadata$Flag <-
        gsub(
                "India",
                "https://raw.githubusercontent.com/maxlcummins/AVC171/master/images/Flags/india.png",
                Metadata$Flag
        )
Metadata$Flag <-
        gsub(
                "Italy",
                "https://raw.githubusercontent.com/maxlcummins/AVC171/master/images/Flags/italy.png",
                Metadata$Flag
        )
#replace countries of origin that are unknown with a URL for a question mark
Metadata$Flag[is.na(Metadata$Flag)] <-
        "https://raw.githubusercontent.com/maxlcummins/AVC171/master/images/Icons/question.png"

#Generate a new column called Classification within which we will give URLs based on
#types of strain source, used by ggimage to render icons in our image.
Metadata$Classification_img <- Metadata$Classification
Metadata$Classification_img <-
        gsub(
                "SEPEC",
                "https://raw.githubusercontent.com/maxlcummins/AVC171/master/images/Icons/blood.png",
                Metadata$Classification_img
        )
Metadata$Classification_img[is.na(Metadata$Classification_img)] <-
        "https://raw.githubusercontent.com/maxlcummins/AVC171/master/images/Icons/question.png"
Metadata$Classification_img <-
        gsub(
                "Faecal",
                "https://raw.githubusercontent.com/maxlcummins/AVC171/master/images/Icons/colon.png",
                Metadata$Classification_img
        )
Metadata$Classification_img <-
        gsub(
                "RMAE",
                "https://raw.githubusercontent.com/maxlcummins/AVC171/master/images/Icons/meat.png",
                Metadata$Classification_img
        )
Metadata$Classification_img <-
        gsub(
                "APEC",
                "https://raw.githubusercontent.com/maxlcummins/AVC171/master/images/Icons/poultry.png",
                Metadata$Classification_img
        )
Metadata$Classification_img <-
        gsub(
                "UPEC",
                "https://raw.githubusercontent.com/maxlcummins/AVC171/master/images/Icons/urine.png",
                Metadata$Classification_img
        )
Metadata$Classification_img <-
        gsub(
                "ExPEC",
                "https://raw.githubusercontent.com/maxlcummins/AVC171/master/images/Icons/ExPEC.png",
                Metadata$Classification_img
        )
Metadata$Classification_img <-
        gsub(
                "Environmental",
                "https://raw.githubusercontent.com/maxlcummins/AVC171/master/images/Icons/earth.png",
                Metadata$Classification_img
        )
Metadata$Classification_img <-
        gsub(
                "Other",
                "https://raw.githubusercontent.com/maxlcummins/AVC171/master/images/Icons/other.png",
                Metadata$Classification_img
        )

#Generate a new column called Pathogen within which we will use to give URLs based on
#pathotype of strains, used by ggimage to render icons in our image.
Metadata$Pathogen <- Metadata$Classification
Metadata$Pathogen <- gsub("ExPEC", "Systemic", Metadata$Pathogen)
Metadata$Pathogen <- gsub("APEC", "Systemic", Metadata$Pathogen)
Metadata$Pathogen <- gsub("UPEC", "Urine", Metadata$Pathogen)
Metadata$Pathogen <- gsub("SEPEC", "Systemic", Metadata$Pathogen)
Metadata$Pathogen <- gsub("RMAE", "Raw Chicken", Metadata$Pathogen)
Metadata$Pathogen <- gsub("Faecal", "Flora", Metadata$Pathogen)

#Here are the URLs for Pathotype images
Metadata$Pathogen_img <- Metadata$Pathogen
Metadata$Pathogen_img <-
        gsub(
                "Flora",
                "https://raw.githubusercontent.com/maxlcummins/AVC171/master/images/Icons/colon.png",
                Metadata$Pathogen_img
        )
Metadata$Pathogen_img <-
        gsub(
                "Raw Chicken",
                "https://raw.githubusercontent.com/maxlcummins/AVC171/master/images/Icons/raw_chicken2.png",
                Metadata$Pathogen_img
        )
Metadata$Pathogen_img <-
        gsub(
                "Systemic",
                "https://raw.githubusercontent.com/maxlcummins/AVC171/master/images/Icons/blood.png",
                Metadata$Pathogen_img
        )
Metadata$Pathogen_img <-
        gsub(
                "Urine",
                "https://raw.githubusercontent.com/maxlcummins/AVC171/master/images/Icons/urine.png",
                Metadata$Pathogen_img
        )
Metadata$Pathogen_img <-
        gsub(
                "Environmental",
                "https://raw.githubusercontent.com/maxlcummins/AVC171/master/images/Icons/earth.png",
                Metadata$Pathogen_img
        )
Metadata$Pathogen_img <-
        gsub(
                "Other",
                "https://raw.githubusercontent.com/maxlcummins/AVC171/master/images/Icons/other.png",
                Metadata$Pathogen_img
        )
Metadata$Pathogen_img[is.na(Metadata$Pathogen_img)] <-
        "https://raw.githubusercontent.com/maxlcummins/AVC171/master/images/Icons/question.png"

#Same as above for source data
Metadata$Revised_Source_Niche_img <- Metadata$Revised_Source_Niche
Metadata$Revised_Source_Niche_img <-
        gsub(
                "Canine",
                "https://raw.githubusercontent.com/maxlcummins/AVC171/master/images/Icons/dog.png",
                Metadata$Revised_Source_Niche_img
        )
Metadata$Revised_Source_Niche_img <-
        gsub(
                "Poultry",
                "https://raw.githubusercontent.com/maxlcummins/AVC171/master/images/Icons/poultry.png",
                Metadata$Revised_Source_Niche_img
        )
Metadata$Revised_Source_Niche_img <-
        gsub(
                "Human",
                "https://raw.githubusercontent.com/maxlcummins/AVC171/master/images/Icons/human.png",
                Metadata$Revised_Source_Niche_img
        )
Metadata$Revised_Source_Niche_img <-
        gsub(
                "Bovine",
                "https://raw.githubusercontent.com/maxlcummins/AVC171/master/images/Icons/cow.png",
                Metadata$Revised_Source_Niche_img
        )
Metadata$Revised_Source_Niche_img <-
        gsub(
                "Environment",
                "https://raw.githubusercontent.com/maxlcummins/AVC171/master/images/Icons/earth.png",
                Metadata$Revised_Source_Niche_img
        )
Metadata$Revised_Source_Niche_img <-
        gsub(
                "Other",
                "https://raw.githubusercontent.com/maxlcummins/AVC171/master/images/Icons/other.png",
                Metadata$Revised_Source_Niche_img
        )
Metadata$Revised_Source_Niche_img[is.na(Metadata$Revised_Source_Niche_img)] <-
        "https://raw.githubusercontent.com/maxlcummins/AVC171/master/images/Icons/question.png"


#### Processing of Genotypic data ####
#Filter the genotypic data to only contain strains in our tree
#Note - non-dynamic variable based on your output names and other variables at
#       the top of this script
df <-
        ST95_all_simple_summary_N90L90 %>% filter(name %in% tree$tip.label)


ST95_all_simple_summary_N90L90 <- ST95_all_simple_summary_N90L90 %>% filter(name %in% tree$tip.label)

#Reorder columns so that name is removed, ColV carriage (Liu) comes first, then
#every other column follows. Also remove pMLST columns
# Note - This may break on inclusion of samples with plasmids in the pMLST db
#       which do not start with "Inc"
# Flag: fix (potential)
df <- df %>% select(ColV, everything(), -name, -starts_with("Inc"))

#Change ColV column to numeric
df$ColV <- as.numeric(df$ColV)

#Replace multiple hits for a given gene with a 1
df[df > 1] <- 1

#create a data frame for generating colsums for given genes (to determine prevelance of a given gene/trait)
colsum <- cbind(colnames(df), colSums(df)) %>% as.data.frame()

#convert our colsums column to be numeric
colsum$V2 <- as.numeric(colsum$V2)

#change all df columns to be numeric
df <- data.frame(lapply(df, as.numeric), stringsAsFactors = FALSE)

#remove columns where colsum is zero - this gets rid of blank colums that may be
#introduced from subsetting by tree tip labels earlier on
df <- df[, colSums(df != 0) > 0]

#Remove unwanted AMR genes from card. 
# I determined which genes to remove by working out which were common to the majority of ST95 isolates:
#df %>% select(starts_with("card_")) %>% 
#               colSums() %>%
#               View()
# You can then determine which genes you consider too common and extract their column names after you filter out genes above the threshold
#df %>% select(starts_with("card_")) %>% 
#        colSums() %>% as.data.frame() %>%
#        rename('count_of_gene' = ".") %>%
#        filter(count_of_gene >= 314) %>%
#        rownames()
#Note - I have decided to hard code the resulting gene names here.
#       If you are replicating this analysis on a different cohort you will need
#       to check if you are happy with these genes being removed.
#       Perhaps some you will want to retain - perhaps there are others you need
#       to add.
df <- df %>% select(
        -starts_with("card_acrB"),
        -starts_with("card_acrD"),
        -starts_with("card_acrE"),
        -starts_with("card_acrF"),
        -starts_with("card_acrS"),
        -starts_with("card_bacA"),
        -starts_with("card_baeR"),
        -starts_with("card_baeS"),
        -starts_with("card_cpxA"),
        -starts_with("card_CRP"),
        -starts_with("card_emrA"),
        -starts_with("card_emrB"),
        -starts_with("card_emrK"),
        -starts_with("card_emrR"),
        -starts_with("card_emrY"),                 
        -starts_with("card_eptA"),
        -starts_with("card_Escherichia_coli_acrA"),
        -starts_with("card_Escherichia_coli_ampC"),
        -starts_with("card_Escherichia_coli_ampH"),
        -starts_with("card_Escherichia_coli_emrE"),
        -starts_with("card_Escherichia_coli_mdfA"),
        -starts_with("card_evgA"),
        -starts_with("card_evgS"),
        -starts_with("card_gadW"),
        -starts_with("card_gadX"),               
        -starts_with("card_H.NS"),
        -starts_with("card_kdpE"),
        -starts_with("card_marA"),
        -starts_with("card_mdtA"),
        -starts_with("card_mdtB"),               
        -starts_with("card_mdtC"),
        -starts_with("card_mdtE"),
        -starts_with("card_mdtF"),
        -starts_with("card_mdtG"),
        -starts_with("card_mdtH"),
        -starts_with("card_mdtM"),
        -starts_with("card_mdtN"),
        -starts_with("card_mdtO"),
        -starts_with("card_mdtP"),
        -starts_with("card_msbA"),
        -starts_with("card_pmrF"),
        -starts_with("card_tolC"),
        -starts_with("card_ugd"),
        -starts_with("card_yojI"),
        #Also got rid of card_SAT.1 and determinant_of_bleomycin_resistance as they encode antiviral or anticancer agents
        #It appears they may also have an antibiotic capacity but they were removed nonetheless
        -starts_with("card_SAT.1"),
        -starts_with("card_determinant_")
)

#Remove unwanted virulence genes - these are in our custom DB or members of operons we dont want to count additively
df <- df %>% select(
        -starts_with("EC_custom_malX"),
        -matches("EC_custom_eit[B-D]"),
        -starts_with("EC_custom_VGI"),
        -starts_with("EC_custom_malX"),
        -starts_with("EC_custom_yeeT"),
        -starts_with("EC_custom_iucD"),
        -starts_with("EC_custom_irp"),
        -starts_with("EC_custom_fim"),
        -starts_with("EC_custom_fyuA")
)

#Remove IS elements
#df <- df %>% select(-starts_with("ISfinder_Feb_2020"))

#Remove unwanted virulence genes - these are in our custom DB or members of operons we dont want to count additively
df <- df %>% select(
        -starts_with("vfdb_chu"),
        -starts_with("vfdb_ent"),
        -starts_with("vfdb_chu"),
        -starts_with("vfdb_fepB|C|D"),
        -starts_with("vfdb_chu"),
        -starts_with("vfdb_fimA|B|C|D|E|F|G|I"),
        -starts_with("vfdb_gsp[D-M]"),
        -starts_with("vfdb_iro[B-D]"),
        -starts_with("vfdb_iuc[ABC]"),
        -starts_with("vfdb_iro[B-D]"),
        -starts_with("vfdb_pap[DEFHIJKX]"),
        -starts_with("vfdb_sfa[ABCDEFGHXY]"),
        -starts_with("vfdb_yag[WXYZ]"),
        -starts_with("vfdb_ybt[EPQSTUX]"),
        -starts_with("vfdb_ykgK")
)

#Remove unwanted virulence genes - these are in our custom DB or members of operons we dont want to count additively
df <- df %>% select(
        -matches("vfdb_chu"),
        -matches("vfdb_ent"),
        -matches("vfdb_chu"),
        -matches("vfdb_fep[BCD]"),
        -matches("vfdb_chu"),
        -matches("vfdb_fim[ABCDEFGI]"),
        -matches("vfdb_gsp[D-M]"),
        -matches("vfdb_iro[B-D]"),
        -matches("vfdb_iuc[ABC]"),
        -matches("vfdb_iro[B-D]"),
        -matches("vfdb_pap[BDEFHIJKX]"),
        -matches("vfdb_sfa[ABCDEFGHXY]"),
        -matches("vfdb_yag[WXYZ]"),
        -matches("vfdb_ybt[EPQSTUX]"),
        -matches("vfdb_ykgK")
)

#Save a copy of this df we use later
db_df <- df

#Change colnames to reflect gene function rather than database source
colnames(df) <- gsub("EC_custom_merA.*", "res_merA", colnames(df))
colnames(df) <- gsub("EC_custom_terA.*", "res_terA", colnames(df))
colnames(df) <- gsub("EC_custom_intI1.*", "res_intI1", colnames(df))
colnames(df) <- gsub("EC_custom_", "vir_", colnames(df))
colnames(df) <- gsub("card_", "res_", colnames(df))
colnames(df) <- gsub("vfdb_", "vir_", colnames(df))
colnames(df) <- gsub("plasmidfinder_", "plas_", colnames(df))
colnames(df) <- gsub("ISfinder_Feb_2020_", "IS_", colnames(df))


#creates a df to be used in generation of new colnames with gene sums
old_new_colnames <- rbind(colnames(df),colnames(df))

#Sum each hit for each gene
genesums <- as.data.frame(colSums(data.matrix(df)))

#rename genesum column to 'sum'
colnames(genesums) <- 'sum'

#paste together the new colnames and assign to our df with old and new names
genesums$rowsumcat <- paste(rownames(genesums), " (", genesums$sum , "/", nrow(df),")", sep="")

#assigns new colnames
colnames(df) <- genesums$rowsumcat

#Separate out our gene hits from different databases to allow us to change their cell contents for dataviz
r <- df %>% select(starts_with("res_"))
p <- df %>% select(starts_with("plas_"))
cv <- df %>% select(starts_with("ColV"))
v <- df %>% select(starts_with("vir_"))
i <- df %>% select(starts_with("IS_"))
po <- df %>% select(starts_with("point"))

#Set plasmid gene hits to be 2, virulence gene hits to be 3, colv genes to be 4
p[p == 1] <- 2
v[v == 1] <- 3
cv[cv == 1] <- 4
i[i == 1] <- 5
po[po == 1] <- 5

#Bind these columns back togeteer
df <- cbind(cv, r, p, v, i, po)

#replace multiple columns for sitA hits to be a column for a single sitA hits
df$vir_sitA <- df %>% select(starts_with("vir_sitA")) %>% rowSums()
df$vir_sitA[df$vir_sitA > 1] <- 3
df <- df %>% select(-starts_with("vir_sitA_"))

#replace multiple columns for incBOKZ hits to be a column for a single incBOKZ hit
df$plas_IncB.O.K.Z <-
        df %>% select(starts_with("plas_IncB")) %>% rowSums()
df$plas_IncB.O.K.Z[df$plas_IncB.O.K.Z > 1] <- 2
df <- df %>% select(-starts_with("plas_IncB.O.K.Z_"))

#Order the genes alphabetically (as they are prefixed by DB name they will be sorted alphabetically by DB name then by gene name)
df <- df[, order(names(df))]

#Trim database name from column names
colnames_save <- gsub("(res|vir|plas)_", "", colnames(df))

#Convert all columns to character type
df <- data.frame(lapply(df, as.character), stringsAsFactors = FALSE)

#Reassign column names
colnames(df) <- colnames_save

#reassign row names.
#Note - This necessitates the row order being the same!
rownames(df) <- ST95_all_simple_summary_N90L90$name

## Combine Data and metadata/cgMLST/phylotype/serotype
# Duplicate our dataframe for the purpose of joining
df_join <- df

# Assign rownames
df_join$working_name <- rownames(df_join)

#Join data and metadata etc
Metadata <- left_join(Metadata, df_join)

#Similar df to Metadata later used, also containst pMLST etc
metadata <- left_join(Metadata, ST95_all_simple_summary_N90L90, by = c("working_name" = "name"))
