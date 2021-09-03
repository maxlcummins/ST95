library(ggtree)
library(phytools)
library(tidyr)
library(readr)
library(dplyr)

working_dir <- "/Users/131785/Dropbox/Doctorate/Manuscripts/AVC171/ST95"

#Set working directory
setwd(working_dir)

#Load paths to files needed for the script
tree_path <- "analysis/output/accessory_ST95_all.tree"
abricate_path <- "analysis/output/genotype.txt"
pointfinder_path <- "analysis/output/pointfinder.txt"
pMLST_data <- "analysis/output/pMLST.txt"
cgMLST_path <- "analysis/data/enterobase/Curated_metadata_all.txt"
phylotype_path <- "analysis/data/enterobase/ST95_phylotypes.txt"
serotype_path <- "analysis/data/enterobase/ST95_serotypes.txt"

#Provide output names
output_name <- "ST95_all"
output_dir <- "processed"

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
Metadata$Pathogen <- gsub("Faecal", "Gastrointestinal", Metadata$Pathogen)

#Here are the URLs for Pathotype images
Metadata$Pathogen_img <- Metadata$Pathogen
Metadata$Pathogen_img <-
        gsub(
                "Gastrointestinal",
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

#Set plasmid gene hits to be 2, virulence gene hits to be 3, colv genes to be 4
p[p == 1] <- 2
v[v == 1] <- 3
cv[cv == 1] <- 4
i[i == 1] <- 5

#Bind these columns back together
df <- cbind(cv, r, p, v, i)

#Keep a copy of this df for future analyses
df_bk <- df
df_bk$working_name <- ST95_all_simple_summary_N90L90$name

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

#### Reading and processing of plasmid mapping data ####

#Read in plasmid map data
source("scripts/plasmid_mapR.R")

#Run plasmid mapR if the pBCE049_1 percent coverage file doesnt exist
if(file.exists("analysis/delims/pBCE049_1_plasmid_coverage_percentage.csv") == FALSE){
        plasmid_mapR(path_to_abricate = "/Users/131785/Dropbox/Doctorate/Manuscripts/AVC171/ST95/analysis/output/abricate/pBCE049_1.txt",
                     plasmid_reference_name = "pBCE049_1",
                     output_directory = "analysis/delims",
                     min_hit_id = 95,
                     min_hit_length = 0.5,
                     writecsv = TRUE,
                     path_to_tree = tree_path,
                     tree_reference_name = "Placeholder")
}

#Read in plasmid percent coveraage dataframe
plas_perc_cov <- read_csv("analysis/delims/pBCE049_1_plasmid_coverage_percentage.csv")

plas_perc_cov <- plas_perc_cov[,3:ncol(plas_perc_cov)]

colnames(plas_perc_cov) <- c("working_name","pBCE049_1")

Metadata <- left_join(Metadata, plas_perc_cov)

#Same again for pUTI89
if(file.exists("analysis/delims/pUTI89_plasmid_coverage_percentage.csv") == FALSE){
        plasmid_mapR(path_to_abricate = "/Users/131785/Dropbox/Doctorate/Manuscripts/AVC171/ST95/analysis/output/abricate/pUTI89.txt",
                     plasmid_reference_name = "pUTI89",
                     output_directory = "analysis/delims",
                     min_hit_id = 95,
                     min_hit_length = 0.5,
                     writecsv = TRUE,
                     path_to_tree = tree_path,
                     tree_reference_name = "Placeholder")
}

#Read in plasmid percent coveraage dataframe
plas_perc_cov <- read_csv("analysis/delims/pUTI89_plasmid_coverage_percentage.csv")

plas_perc_cov <- plas_perc_cov[,3:ncol(plas_perc_cov)]

colnames(plas_perc_cov) <- c("working_name","pUTI89")

Metadata <- left_join(Metadata, plas_perc_cov)

#Same again for pSF_088 (minus the resistance region)
if(file.exists("analysis/delims/pSF_088_nores_plasmid_coverage_percentage.csv") == FALSE){
        plasmid_mapR(path_to_abricate = "/Users/131785/Dropbox/Doctorate/Manuscripts/AVC171/ST95/analysis/output/abricate/pSF_088_nores.txt",
                     plasmid_reference_name = "pSF_088_nores",
                     output_directory = "analysis/delims",
                     min_hit_id = 95,
                     min_hit_length = 0.5,
                     writecsv = TRUE,
                     path_to_tree = tree_path,
                     tree_reference_name = "Placeholder")
}

plas_perc_cov <- read_csv("analysis/delims/pSF_088_nores_plasmid_coverage_percentage.csv")

plas_perc_cov <- plas_perc_cov[,3:ncol(plas_perc_cov)]

colnames(plas_perc_cov) <- c("working_name","pSF_088_nores")

Metadata <- left_join(Metadata, plas_perc_cov)

#Same again for pAPEC_O2_ColV (minus the resistance region)
if(file.exists("analysis/delims/pAPEC_O2_ColV_plasmid_coverage_percentage.csv") == FALSE){
        plasmid_mapR(path_to_abricate = "/Users/131785/Dropbox/Doctorate/Manuscripts/AVC171/ST95/analysis/output/abricate/pAPEC_O2_ColV.txt",
                     plasmid_reference_name = "pAPEC_O2_ColV",
                     output_directory = "analysis/delims",
                     min_hit_id = 95,
                     min_hit_length = 0.5,
                     writecsv = TRUE,
                     path_to_tree = tree_path,
                     tree_reference_name = "Placeholder")
}


plas_perc_cov <- read_csv("analysis/delims/pAPEC_O2_ColV_plasmid_coverage_percentage.csv")

plas_perc_cov <- plas_perc_cov[,3:ncol(plas_perc_cov)]

colnames(plas_perc_cov) <- c("working_name","pAPEC_O2_ColV")

Metadata <- left_join(Metadata, plas_perc_cov)

#Same again for pU1_F51_B10 (minus the resistance region)
if(file.exists("analysis/delims/pU1_F51_B10_plasmid_coverage_percentage.csv") == FALSE){
        plasmid_mapR(path_to_abricate = "/Users/131785/Dropbox/Doctorate/Manuscripts/AVC171/ST95/analysis/output/abricate/pU1_F51_B10.txt",
                     plasmid_reference_name = "pU1_F51_B10",
                     output_directory = "analysis/delims",
                     min_hit_id = 95,
                     min_hit_length = 0.5,
                     writecsv = TRUE,
                     path_to_tree = tree_path,
                     tree_reference_name = "Placeholder")
}



plas_perc_cov <- read_csv("analysis/delims/pU1_F51_B10_plasmid_coverage_percentage.csv")

plas_perc_cov <- plas_perc_cov[,3:ncol(plas_perc_cov)]

colnames(plas_perc_cov) <- c("working_name","pU1_F51_B10")

Metadata <- left_join(Metadata, plas_perc_cov)


#Same again for pAPEC_O1_ColBM (minus the resistance region)
if(file.exists("analysis/delims/pAPEC_O1_ColBM_plasmid_coverage_percentage.csv") == FALSE){
        plasmid_mapR(path_to_abricate = "/Users/131785/Dropbox/Doctorate/Manuscripts/AVC171/ST95/analysis/output/abricate/pAPEC_O1_ColBM.txt",
                     plasmid_reference_name = "pAPEC_O1_ColBM",
                     output_directory = "analysis/delims",
                     min_hit_id = 95,
                     min_hit_length = 0.5,
                     writecsv = TRUE,
                     path_to_tree = tree_path,
                     tree_reference_name = "Placeholder")
}


plas_perc_cov <- read_csv("analysis/delims/pAPEC_O1_ColBM_plasmid_coverage_percentage.csv")

plas_perc_cov <- plas_perc_cov[,3:ncol(plas_perc_cov)]

colnames(plas_perc_cov) <- c("working_name","pAPEC_O1_ColBM")

Metadata <- left_join(Metadata, plas_perc_cov)

#Same again for pAPEC_O78 (minus the resistance region)
if(file.exists("analysis/delims/pAPEC_O78_ColV_plasmid_coverage_percentage.csv") == FALSE){
        plasmid_mapR(path_to_abricate = "/Users/131785/Dropbox/Doctorate/Manuscripts/AVC171/ST95/analysis/output/abricate/pAPEC_O78_ColV.txt",
                     plasmid_reference_name = "pAPEC_O78_ColV",
                     output_directory = "analysis/delims",
                     min_hit_id = 95,
                     min_hit_length = 0.5,
                     writecsv = TRUE,
                     path_to_tree = tree_path,
                     tree_reference_name = "Placeholder")
}


plas_perc_cov <- read_csv("analysis/delims/pAPEC_O78_ColV_plasmid_coverage_percentage.csv")

plas_perc_cov <- plas_perc_cov[,3:ncol(plas_perc_cov)]

colnames(plas_perc_cov) <- c("working_name","pAPEC_O78_ColV")

Metadata <- left_join(Metadata, plas_perc_cov)


#Same again for pEC244_2 (minus the resistance region)
if(file.exists("analysis/delims/pEC244_2_plasmid_coverage_percentage.csv") == FALSE){
        plasmid_mapR(path_to_abricate = "/Users/131785/Dropbox/Doctorate/Manuscripts/AVC171/ST95/analysis/output/abricate/pEC244_2.txt",
                     plasmid_reference_name = "pEC244_2",
                     output_directory = "analysis/delims",
                     min_hit_id = 95,
                     min_hit_length = 0.5,
                     writecsv = TRUE,
                     path_to_tree = tree_path,
                     tree_reference_name = "Placeholder")
}


plas_perc_cov <- read_csv("analysis/delims/pEC244_2_plasmid_coverage_percentage.csv")

plas_perc_cov <- plas_perc_cov[,3:ncol(plas_perc_cov)]

colnames(plas_perc_cov) <- c("working_name","pEC244_2")

Metadata <- left_join(Metadata, plas_perc_cov)

#Same again for pAMSC2 (minus the resistance region)
if(file.exists("analysis/delims/pAMSC2_plasmid_coverage_percentage.csv") == FALSE){
        plasmid_mapR(path_to_abricate = "/Users/131785/Dropbox/Doctorate/Manuscripts/AVC171/ST95/analysis/output/abricate/pAMSC2.txt",
                     plasmid_reference_name = "pAMSC2",
                     output_directory = "analysis/delims",
                     min_hit_id = 95,
                     min_hit_length = 0.5,
                     writecsv = TRUE,
                     path_to_tree = tree_path,
                     tree_reference_name = "Placeholder")
}


plas_perc_cov <- read_csv("analysis/delims/pAMSC2_plasmid_coverage_percentage.csv")

plas_perc_cov <- plas_perc_cov[,3:ncol(plas_perc_cov)]

colnames(plas_perc_cov) <- c("working_name","pAMSC2")

Metadata <- left_join(Metadata, plas_perc_cov)


#Same again for pAMSC2 (minus the resistance region)
if(file.exists("analysis/delims/pACN001_B_plasmid_coverage_percentage.csv") == FALSE){
        plasmid_mapR(path_to_abricate = "/Users/131785/Dropbox/Doctorate/Manuscripts/AVC171/ST95/analysis/output/abricate/pACN001_B.txt",
                     plasmid_reference_name = "pACN001_B",
                     output_directory = "analysis/delims",
                     min_hit_id = 95,
                     min_hit_length = 0.5,
                     writecsv = TRUE,
                     path_to_tree = tree_path,
                     tree_reference_name = "Placeholder")
}


plas_perc_cov <- read_csv("analysis/delims/pACN001_B_plasmid_coverage_percentage.csv")

plas_perc_cov <- plas_perc_cov[,3:ncol(plas_perc_cov)]

colnames(plas_perc_cov) <- c("working_name","pACN001_B")

Metadata <- left_join(Metadata, plas_perc_cov)

#### Designations of strains as APEC or ExPEC ####

# Select columns that are associated with APEC status (Minimal predictors of APEC status doi: 10.1128/JCM.00816-08)
APEC <- df %>% select(`iutA_KU578032.1 (330/668)`, `hylF_CP000836.1 (336/668)`, `iss_type1_pAPEC.O2.ColV (4/668)`, `iroN (388/668)`, `ompT_episomal_HM210637.1 (324/668)`)

#Convert these columns to APEC
APEC <- sapply(APEC, as.numeric) %>% as.data.frame()

#Reassign gene counts greater than one back to one
APEC[APEC > 1] <- 1

#Generate rowsums to develop the number of these minimal predictors per strain
APEC$APEC_val <- rowSums(APEC)

#Reassign rownames
rownames(APEC) <- rownames(df)

#Assign rownames to working_names column
APEC$working_name <- rownames(APEC)

#Generate two columns, strain name and count of APEC genotypic predictors
APEC <- APEC %>% select(working_name, APEC_val)

#Assign APEC status to samples with >= 3 minimal APEC predictors
APEC$APEC_val[APEC$APEC_val < 3] <- 0
APEC$APEC_val[APEC$APEC_val >= 3] <- 1

#Join APEC status back to metadata
Metadata <- left_join(Metadata,APEC)

#As above, except for ExPEC status - Isolation and molecular characterization of
#       nalidixic acid-resistant extraintestinal pathogenic Escherichia coli from
#       retail chicken products

#Isolate ExPEC predictors
ExPEC <- df %>% select(`papC (585/668)`, `papA_CU928161 (504/668)`, `sfaS (101/668)`, `focG (1/668)`, starts_with("afaA"), starts_with("dra"), `kpsMT.II._K2.CP000468.1..APEC.O1.3376659.3378102 (655/668)`, `iutA_KU578032.1 (330/668)`)

#Reassign as type numeric
ExPEC <- sapply(ExPEC, as.numeric) %>% as.data.frame()

#Reassign multiple gene counts to single gene counts
ExPEC[ExPEC > 1] <- 1

#Cluster into groups as per Johnson criteria
ExPEC$V1 <- ExPEC$`papC (585/668)` + ExPEC$`papA_CU928161 (504/668)`
ExPEC$V2 <- ExPEC$`sfaS (101/668)` + ExPEC$`focG (1/668)`
ExPEC$V4 <- ExPEC$`kpsMT.II._K2.CP000468.1..APEC.O1.3376659.3378102 (655/668)`
ExPEC$V5 <- ExPEC$`iutA_KU578032.1 (330/668)`

#Reassign multiple gene counts to single gene counts
ExPEC[ExPEC > 1] <- 1

#Remove columns other than those with grouped genes
ExPEC <- ExPEC %>% select(starts_with("V"))

#Count the number of genes
ExPEC$ExPEC_val <- rowSums(ExPEC)

#Reassign row names
rownames(ExPEC) <- rownames(df)

#Assign rownames to working_names
ExPEC$working_name <- rownames(ExPEC)

#Select only strain names and ExPEC gene group count sum
ExPEC <- ExPEC %>% select(working_name, ExPEC_val)

#Classify strains as ExPEC if it has < 2 gene groups represented, vice-versa
ExPEC$ExPEC_val[ExPEC$ExPEC_val < 2] <- 0
ExPEC$ExPEC_val[ExPEC$ExPEC_val >= 2] <- 1

#Bind metadata and ExPEC status
Metadata <- left_join(Metadata,ExPEC)

#### Counting of AMR and non-IncF repA genes (Inc only) ####

#Create a data frame with only genes from card
#Note: this requires your gene database used for AMR to be card!
amr_df <- db_df %>% select(starts_with("card_"))

#Generate the count of amr genes
#Note: we have already plucked out genes we are uninterested in reporting...
amr_count <- as.data.frame(rowSums(amr_df))

#Reassign AMR_count rownames
amr_count$working_name <- ST95_all_simple_summary_N90L90$name

#Reassign column names
colnames(amr_count) <- c("amr_counts","working_name")

#Bind Metadata and amr_count
Metadata <- left_join(Metadata, amr_count)

## Pull out a list of phenotypic associations with AMR genes
amr_gene_phenotype <- ST95_all.N90.L90.PASS %>% select(GENE, RESISTANCE)
amr_gene_phenotype <- amr_gene_phenotype %>% unique()
amr_gene_phenotype <- amr_gene_phenotype %>% filter(grepl("card", GENE))

amr_gene_phenotype$GENE <- gsub("\\(|\\)", ".", amr_gene_phenotype$GENE)
amr_gene_phenotype$GENE <- gsub("\\'", ".", amr_gene_phenotype$GENE)
amr_gene_phenotype$GENE <- gsub("\\-", ".", amr_gene_phenotype$GENE)

amr_gene_phenotype <- amr_gene_phenotype %>% filter(GENE %in% colnames(db_df))

#Create a column for whether or not a strain is resistant to aminoglycosides
aminoglycosides <- amr_gene_phenotype %>% filter(grepl("aminoglycoside",RESISTANCE)) %>% select(GENE)
aminoglycosides <- aminoglycosides[[1]]
aminoglycosides <- aminoglycosides[aminoglycosides %in% colnames(amr_df)]

aminoglycoside <- amr_df %>% select(aminoglycosides)

aminoglycoside_res <- rowSums(aminoglycoside)

#Create a column for whether or not a strain is resistant to trimethoprim
trimethoprims <- amr_gene_phenotype %>% filter(grepl("diaminopyrimidine",RESISTANCE)) %>% select(GENE)
trimethoprims <- trimethoprims[[1]]
trimethoprims <- trimethoprims[trimethoprims %in% colnames(amr_df)]

trimethoprim <- amr_df %>% select(trimethoprims) %>% select(starts_with("card_dfrA"))

trimethoprim_res <- rowSums(trimethoprim)

#Create a column for whether or not a strain is resistant to fluoroquinolones
fluoroquinolones <- amr_gene_phenotype %>% filter(grepl("fluoroquinolone",RESISTANCE)) %>% select(GENE)
fluoroquinolones <- fluoroquinolones[[1]]
fluoroquinolones <- fluoroquinolones[fluoroquinolones %in% colnames(amr_df)]

fluoroquinolone <- amr_df %>% select(fluoroquinolones) %>% select(starts_with("card_Qnr"))

fluoroquinolone_snps <- ST95_all_simple_summary_N90L90 %>% select(starts_with("pointfinder_gyrA"), starts_with("pointfinder_parC"))

fluoroquinolone <- cbind(fluoroquinolone, fluoroquinolone_snps)

fluoroquinolone_res <- rowSums(fluoroquinolone)

#Create a column for whether or not a strain is resistant to sulfonamides
sulfonamides <- amr_gene_phenotype %>% filter(grepl("sulfonamide",RESISTANCE)) %>% select(GENE)
sulfonamides <- sulfonamides[[1]]
sulfonamides <- sulfonamides[sulfonamides %in% colnames(amr_df)]

sulfonamide <- amr_df %>% select(sulfonamides)

sulfonamide_snps <- ST95_all_simple_summary_N90L90 %>% select(starts_with("pointfinder_folP"))

sulfonamide <- cbind(sulfonamide, sulfonamide_snps)

sulfonamide_res <- rowSums(sulfonamide)

#sul_tri <- cbind(sulfonamide_res, trimethoprim_res)

#sul_tri_res <- rowSums(sulfonamide)


#Create a column for whether or not a strain is resistant to tetracyclines
tetracyclines <- amr_gene_phenotype %>% filter(grepl("tetracycline",RESISTANCE)) %>% select(GENE)
tetracyclines <- tetracyclines[[1]]
tetracyclines <- tetracyclines[tetracyclines %in% colnames(amr_df)]

tetracycline <- amr_df %>% select(tetracyclines) %>% select(starts_with("card_tet"))

tetracycline_res <- rowSums(tetracycline)

#Create a column for whether or not a strain is resistant to polymyxins
polymyxins <- amr_gene_phenotype %>% filter(grepl("peptide",RESISTANCE)) %>% select(GENE)
polymyxins <- polymyxins[[1]]
polymyxins <- polymyxins[polymyxins %in% colnames(amr_df)]

polymyxin <- amr_df %>% select(polymyxins) %>% select(starts_with("card_MCR"))

polymyxin_res <- rowSums(polymyxin)

#Create a column for whether or not a strain is resistant to amphenicol
amphenicols <- amr_gene_phenotype %>% filter(grepl("phenicol",RESISTANCE)) %>% select(GENE)
amphenicols <- amphenicols[[1]]
amphenicols <- amphenicols[amphenicols %in% colnames(amr_df)]

amphenicol <- amr_df %>% select(amphenicols) %>% select(-starts_with("card_tet"), -starts_with("card_msrA"))

amphenicol_res <- rowSums(amphenicol)

#Create a column for whether or not a strain is resistant to macrolides-lincosamides
macrolide_lincosamides <- amr_gene_phenotype %>% filter(grepl("macrolide|lincosamide",RESISTANCE)) %>% select(GENE)
macrolide_lincosamides <- macrolide_lincosamides[[1]]
macrolide_lincosamides <- macrolide_lincosamides[macrolide_lincosamides %in% colnames(amr_df)]

macrolide_lincosamide <- amr_df %>% select(macrolide_lincosamides) %>% select(-starts_with("card_tet"),-starts_with("card_mdtM"), -starts_with("card_floR"), -starts_with("card_qacH"), -starts_with("card_cmlA1"))

macrolide_lincosamide_res <- rowSums(macrolide_lincosamide)

#Create a column for whether or not a strain is resistant to third_gen_cephalosporin
beta_lactams <- amr_gene_phenotype %>% filter(grepl("ceph|penam|penem",RESISTANCE)) %>% select(GENE)
beta_lactams <- beta_lactams[[1]]
beta_lactams <- beta_lactams[beta_lactams %in% colnames(amr_df)]

beta_lactamase <- amr_df %>% select(beta_lactams) %>% select(-starts_with("card_tet"),
                                                             -starts_with("card_cmlA1"),
                                                             -starts_with("card_floR"))

ESBL <- amr_df %>% select(beta_lactams) %>% select(-starts_with("card_tet"),
                                                   -starts_with("card_cmlA1"),
                                                   -starts_with("card_floR"),
                                                   -starts_with("card_TEM"),
                                                   -starts_with("card_SHV"),
                                                   -starts_with("card_OXA"),
                                                   -starts_with("card_DHA"))


beta_lactam_res <- rowSums(beta_lactamase)
ESBL_res <- rowSums(ESBL)

ESBL_res[ESBL_res > 1] <- 1

ESBL_res <- as.data.frame(ESBL_res)

ESBL_res$working_name <- ST95_all_simple_summary_N90L90$name

#Generate a data frame and subsequent count of different antimicrobial classes which a strain is resistant to

class_res <- cbind(fluoroquinolone_res, aminoglycoside_res, sulfonamide_res, tetracycline_res, polymyxin_res, trimethoprim_res, amphenicol_res, beta_lactam_res, macrolide_lincosamide_res)

#Collapse cases where a strain has multiple genes conferring resitance to the
#same antibiotic class
class_res[class_res > 1] <- 1

#Sum the number of antibiotic classes each sample is resistant too
class_res_count <- rowSums(class_res)

#Convert this to a dataframe
class_res <- as.data.frame(class_res)

#Reassign names
class_res$working_name <- ST95_all_simple_summary_N90L90$name

#Create a dataframe with count of antimicrobial AMR class resistances
class_res_count <- as.data.frame(class_res_count)

#Assign names
class_res_count$working_name <- ST95_all_simple_summary_N90L90$name

#Generate a data frame and subsequent count of different antimicrobial classes which a strain is resistant to
ESBL_res_column <- ESBL_res$ESBL_res
CIA_class_res <- cbind(fluoroquinolone_res, macrolide_lincosamide_res, ESBL_res_column)

CIA_class_res <- as.data.frame(CIA_class_res)

#Sum the number of antibiotic classes each sample is resistant too
CIA_class_res$sum <- rowSums(CIA_class_res)

#Collapse cases where a strain has multiple genes conferring resitance to the
#same antibiotic class
CIA_class_res$sum[CIA_class_res$sum > 1] <- 1

#Reassign names
CIA_class_res$working_name <- ST95_all_simple_summary_N90L90$name

CIA_res <- CIA_class_res %>% select(working_name, sum) %>% rename("CIA_resistance" = "sum")

#Rename columns so that functions dont confuse df names and column names
colnames(class_res_count)[1] <- "class_res_counts"
colnames(ESBL_res)[1] <- "ESBL_ress"

#Bind metadata and new AMR phenotypic resistance data
Metadata <- Metadata %>% left_join(class_res_count)
Metadata <- Metadata %>% left_join(class_res)
Metadata <- Metadata %>% left_join(ESBL_res)
Metadata <- Metadata %>% left_join(CIA_res)

# Generates a table comparing CIA resistance by source niche
CIA_vs_source <- Metadata %>% group_by(CIA_resistance, Revised_Source_Niche) %>% summarise(count = n())
CIA_vs_source$total <- (CIA_vs_source$count/668)*100
CIA_vs_source %>% group_by(CIA_resistance) %>% summarise(proportion = sum(total), count = sum(count))

#Generate a new column of genotypic MDR (resistant to >= 3 classes)
Metadata$MDR <- Metadata$class_res_counts
Metadata$MDR[Metadata$MDR < 3] <- 0
Metadata$MDR[Metadata$MDR >= 3] <- 1

Metadata$class_res_counts <- gsub("^","classes_res_0",Metadata$class_res_counts)
Metadata$ESBL_ress <- gsub("1","ESBL_pos",Metadata$ESBL_ress)
Metadata$ESBL_ress <- gsub("0","ESBL_neg",Metadata$ESBL_ress)


## Generate a count of number of non-IncF incompatibility types detected per strain
plas_df <- db_df %>% select(starts_with("plasmidfinder_Inc")) %>% select(-starts_with("plasmidfinder_IncF"))
plas_count <- as.data.frame(rowSums(plas_df))

#Assign names for trackin purposes
plas_count$working_name <- ST95_all_simple_summary_N90L90$name

colnames(plas_count) <- c("plas_counts","working_name")

#Bind tables
Metadata <- left_join(Metadata, plas_count)

#This gsub is performed to allow correct sorting of resistance counts
# (otherwise count_res_2 is considered > count_res_10)
Metadata$amr_counts <- gsub("^","count_res_0",Metadata$amr_counts)
Metadata$amr_counts <- gsub("count_res_010","count_res_10",Metadata$amr_counts)
Metadata$amr_counts <- gsub("count_res_011","count_res_11",Metadata$amr_counts)
Metadata$amr_counts <- gsub("count_res_012","count_res_12",Metadata$amr_counts)
Metadata$amr_counts <- gsub("count_res_013","count_res_13",Metadata$amr_counts)

#Change to char type and add prefix (this helps later with assigning colours to our plots)
Metadata$plas_counts <- gsub("^","count_plas_",Metadata$plas_counts)


#### Generation of small dataframe for Figure 1 ####
#Change to char type and add prefix (this helps later with assigning colours to our plots)
#df$`ColV (328/668)`<- gsub("4","ColV_pos",df$`ColV (328/668)`)
#df$`ColV (328/668)` <- gsub("0","ColV_neg",df$`ColV (328/668)`)

#Create a table using the ColV column of df as a base
plasmid_coverage_table <- as.data.frame(df$`ColV (328/668)`)

#Fix the ugly resultant colnames
colnames(plasmid_coverage_table) <- gsub("^df.*", "ColV", colnames(plasmid_coverage_table))

#Reassign strain names
plasmid_coverage_table$working_name <- rownames(df)

#Join this with our metadata table
plasmid_coverage_table <- left_join(plasmid_coverage_table, Metadata)

##Aprefix (this helps later with assigning colours to our plots)
plasmid_coverage_table$Revised_Source_Niche <- gsub("^", "Source_", plasmid_coverage_table$Revised_Source_Niche)

#Change to df
plasmid_coverage_table <- as.data.frame(plasmid_coverage_table)


#Clean colnames
colnames(plasmid_coverage_table) <- gsub("intI1.*", "intI1", colnames(plasmid_coverage_table))

#Select columns of interest
plasmid_coverage_table <- plasmid_coverage_table %>% select(Revised_Source_Niche, pU1_F51_B10, pBCE049_1, pAPEC_O2_ColV, pSF_088_nores, pUTI89, pACN001_B, pEC244_2, intI1, class_res_counts, ESBL_ress, plas_counts)

#Note - the below script doesn't work for cases where a plasmid has two equally tied best matches!
#Flag: Fix me
#For every row in our table generated above
for(i in 1:nrow(plasmid_coverage_table)){
        #find the plasmid with the highest percent coverage among the plasmids
        #screened (these columns are 2:6) and store this coverage as best_match_val
        best_match_val <- max(plasmid_coverage_table[i,2:8])
        #If a hit for a given plasmid is less than the best_match_val, set it to 0
        #This will leave us with only best matches in the columns 2:6
        if(plasmid_coverage_table[i,2] < best_match_val){
                plasmid_coverage_table[i,2] <- 0
        }
        if(plasmid_coverage_table[i,3] < best_match_val){
                plasmid_coverage_table[i,3] <- 0
        }
        if(plasmid_coverage_table[i,4] < best_match_val){
                plasmid_coverage_table[i,4] <- 0
        }
        if(plasmid_coverage_table[i,5] < best_match_val){
                plasmid_coverage_table[i,5] <- 0
        }
        if(plasmid_coverage_table[i,6] < best_match_val){
                plasmid_coverage_table[i,6] <- 0
        }
        if(plasmid_coverage_table[i,7] < best_match_val){
                plasmid_coverage_table[i,7] <- 0
        }
        if(plasmid_coverage_table[i,8] < best_match_val){
                plasmid_coverage_table[i,8] <- 0
        }
}


df_small2 <- plasmid_coverage_table

#### Binning of percentage hits for plasmids into ranges #### 
plasmid_coverage_table$pBCE049_1 <- paste0("pBCE049_1_",
                                           as.character(
                                                   cut(
                                                           plasmid_coverage_table$pBCE049_1,
                                                           breaks = c(0, 49, 59, 69, 79, 89, 100),
                                                           labels = c("0-49", "50-59", "60-69", "70-79", "80-89", "90-100"),
                                                           include.lowest = TRUE
                                                   )
                                           ))

plasmid_coverage_table$pSF_088_nores <- paste0("pSF_088_nores_",
                                               as.character(
                                                       cut(
                                                               plasmid_coverage_table$pSF_088_nores,
                                                               breaks = c(0, 49, 59, 69, 79, 89, 100),
                                                               labels = c("0-49", "50-59", "60-69", "70-79", "80-89", "90-100"),
                                                               include.lowest = TRUE
                                                       )
                                               ))

plasmid_coverage_table$pAPEC_O2_ColV <- paste0("pAPEC_O2_ColV_",
                                               as.character(
                                                       cut(
                                                               plasmid_coverage_table$pAPEC_O2_ColV,
                                                               breaks = c(0, 49, 59, 69, 79, 89, 100),
                                                               labels = c("0-49", "50-59", "60-69", "70-79", "80-89", "90-100"),
                                                               include.lowest = TRUE
                                                       )
                                               ))

plasmid_coverage_table$pUTI89 <- paste0("pUTI89_",
                                        as.character(
                                                cut(
                                                        plasmid_coverage_table$pUTI89,
                                                        breaks = c(0, 49, 59, 69, 79, 89, 100),
                                                        labels = c("0-49", "50-59", "60-69", "70-79", "80-89", "90-100"),
                                                        include.lowest = TRUE
                                                )
                                        ))

plasmid_coverage_table$pU1_F51_B10 <- paste0("pU1_F51_B10_",
                                             as.character(
                                                     cut(
                                                             plasmid_coverage_table$pU1_F51_B10,
                                                             breaks = c(0, 49, 59, 69, 79, 89, 100),
                                                             labels = c("0-49", "50-59", "60-69", "70-79", "80-89", "90-100"),
                                                             include.lowest = TRUE
                                                     )
                                             ))


plasmid_coverage_table$pACN001_B <- paste0("pACN001_B_",
                                             as.character(
                                                     cut(
                                                             plasmid_coverage_table$pACN001_B,
                                                             breaks = c(0, 49, 59, 69, 79, 89, 100),
                                                             labels = c("0-49", "50-59", "60-69", "70-79", "80-89", "90-100"),
                                                             include.lowest = TRUE
                                                     )
                                             ))


plasmid_coverage_table$pEC244_2 <- paste0("pEC244_2_",
                                             as.character(
                                                     cut(
                                                             plasmid_coverage_table$pEC244_2,
                                                             breaks = c(0, 49, 59, 69, 79, 89, 100),
                                                             labels = c("0-49", "50-59", "60-69", "70-79", "80-89", "90-100"),
                                                             include.lowest = TRUE
                                                     )
                                             ))

rownames(plasmid_coverage_table) <- rownames(df)

#### Determination of reference plasmids as present or absent ####
plasmid_data <- df_small2

plasmid_data$working_name <- rownames(plasmid_coverage_table)

ColV_only <- Metadata %>% select(working_name, `ColV (328/668)`)

colV_comparisons_ <- plasmid_data %>% select(working_name, pSF_088_nores, pAPEC_O2_ColV, pUTI89, pU1_F51_B10, pBCE049_1, pACN001_B , pEC244_2, Revised_Source_Niche)

colV_comparisons <- colV_comparisons_

colV_comparisons$sums<- rowSums(colV_comparisons[2:8])

colV_comparisons_ <- left_join(colV_comparisons_, ColV_only)

colV_comparisons_$pSF_088_nores[colV_comparisons_$pSF_088_nores < 80] <- 0
colV_comparisons_$pSF_088_nores[colV_comparisons_$pSF_088_nores >= 80] <- 1

colV_comparisons_$pAPEC_O2_ColV[colV_comparisons_$pAPEC_O2_ColV < 80] <- 0
colV_comparisons_$pAPEC_O2_ColV[colV_comparisons_$pAPEC_O2_ColV >= 80] <- 1

colV_comparisons_$pUTI89[colV_comparisons_$pUTI89 < 80] <- 0
colV_comparisons_$pUTI89[colV_comparisons_$pUTI89 >= 80] <- 1

colV_comparisons_$pU1_F51_B10[colV_comparisons_$pU1_F51_B10 < 80] <- 0
colV_comparisons_$pU1_F51_B10[colV_comparisons_$pU1_F51_B10 >= 80] <- 1

colV_comparisons_$pBCE049_1[colV_comparisons_$pBCE049_1 < 80] <- 0
colV_comparisons_$pBCE049_1[colV_comparisons_$pBCE049_1 >= 80] <- 1

colV_comparisons_$pACN001_B[colV_comparisons_$pACN001_B < 80] <- 0
colV_comparisons_$pACN001_B[colV_comparisons_$pACN001_B >= 80] <- 1

colV_comparisons_$pEC244_2[colV_comparisons_$pEC244_2 < 80] <- 0
colV_comparisons_$pEC244_2[colV_comparisons_$pEC244_2 >= 80] <- 1

colV_pos <- colV_comparisons_ %>% filter(`ColV (328/668)` == 4) %>% group_by(pSF_088_nores, pAPEC_O2_ColV, pEC244_2, pACN001_B) %>% summarise(counts = n()) 

colv_neg <- colV_comparisons_ %>% filter(`ColV (328/668)` == 0) %>% group_by(pUTI89, pU1_F51_B10, pBCE049_1) %>% summarise(counts = n()) 

all_plas <- colV_comparisons_ %>% group_by(`ColV (328/668)`, pSF_088_nores, pAPEC_O2_ColV, pUTI89, pU1_F51_B10, pBCE049_1) %>% summarise(counts = n()) 

plas_neg <- colV_comparisons %>% filter(sums == 0) %>% select(working_name) %>% as.data.frame()

plas_neg_pMLSTs <- metadata %>% filter(working_name %in% plas_neg$working_name) %>% group_by(IncF_RST) %>% summarise(counts = n())

colV_comparisons_ %>% filter(pAPEC_O2_ColV == 1) %>% group_by(Revised_Source_Niche) %>% summarise(counts = n())

colV_comparisons_ %>% filter(pSF_088_nores == 1) %>% group_by(Revised_Source_Niche) %>% summarise(counts = n())

colV_comparisons_ %>% filter(pU1_F51_B10 == 1) %>% group_by(Revised_Source_Niche) %>% summarise(counts = n())

colV_comparisons_ %>% filter(pUTI89 == 1) %>% group_by(Revised_Source_Niche) %>% summarise(counts = n())

colV_comparisons_ %>% filter(pBCE049_1 == 1) %>% group_by(Revised_Source_Niche) %>% summarise(counts = n())

colV_comparisons_ %>% filter(pEC244_2 == 1) %>% group_by(Revised_Source_Niche) %>% summarise(counts = n())

colV_comparisons_ %>% filter(pACN001_B == 1) %>% group_by(Revised_Source_Niche) %>% summarise(counts = n())

#### Generation of a separate plasmid_coverage for Figure 1 which collapses plasmids into a single column ####
rownames(df_small2) <- rownames(plasmid_coverage_table)

df_small3 <- plasmid_coverage_table %>% select(starts_with("p"), -plas_counts)

df_small3$working_name <- rownames(df_small2)

df_small3 <- df_small3 %>% melt(id.vars = "working_name") %>% filter(!grepl("0-49$", value))
df_small3 <- df_small3 %>% melt(id.vars = "working_name") %>% filter(!grepl("50-59$", value))

df_small3 <- df_small3 %>% filter(!grepl("variable", variable)) %>% select(working_name, value)

#df_small3 <- df_small3 %>% filter(working_name != "ESC_IB7597AA_AS" & value  != "pSF_088_nores_81")

df_small3 <- df_small3 %>% rename(plasmid_map = value)

df_small4 <- plasmid_coverage_table

df_small4$working_name <- rownames(plasmid_coverage_table)

nems <- df_small4$working_name

df_small4 <- left_join(df_small4, df_small3)

extra_data <- metadata %>% select(working_name, Pathogen, IncF_RST, Continent, ColV)

df_small4 <- left_join(df_small4, extra_data, by = "working_name")

df_small5 <- left_join(df_small4, Metadata, by = "working_name")

df_small4 <- df_small4 %>% select(Pathogen, Revised_Source_Niche, plasmid_map, intI1, class_res_counts, ESBL_ress, plas_counts, Continent, ColV, IncF_RST)

#df_small4$Pathogen <- gsub("^","_pathogen_status_", df_small4$Pathogen)

df_small4$Pathogen <- gsub("Flora","Gastrointestinal", df_small4$Pathogen)

rownames(df_small4) <- nems

metadata$ColV_pos_neg<- gsub("4","ColV_pos",metadata$`ColV (328/668)`)
metadata$ColV_pos_neg <- gsub("0","ColV_neg",metadata$ColV_pos_neg)



#### Generation of Scoary trait table ####

scoary_traits <- Metadata %>% select(working_name, Revised_Source_Niche, Pathogen, HC50, HC200)

scoary_traits$Revised_Source_Niche <- replace_na(scoary_traits$Revised_Source_Niche, replace = "Unknown")

scoary_traits$HC50 <- gsub("^","HC50_",scoary_traits$HC50)
scoary_traits$HC200 <- gsub("^","HC200_",scoary_traits$HC200)

scoary_traits$present <- rep(1, nrow(scoary_traits))

scoary_traits_2 <- dcast(scoary_traits, working_name ~ Revised_Source_Niche)

scoary_traits_2 <- dcast(scoary_traits, working_name ~ Pathogen) %>% left_join(scoary_traits_2, by = "working_name")

scoary_traits_2 <- dcast(scoary_traits, working_name ~ HC50) %>% left_join(scoary_traits_2, by = "working_name")

scoary_traits_2 <- dcast(scoary_traits, working_name ~ HC200) %>% left_join(scoary_traits_2, by = "working_name")

scoary_traits_2[is.na(scoary_traits_2)] <- 0

colnames(scoary_traits_2)[1] <- ""

#write.csv(scoary_traits_2, "delims/scoary_traits.csv", row.names = FALSE)

