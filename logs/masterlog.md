

___
# Collection of strains - information and access

Genomic assemblies of ST95 were downloaded from Enterobase on access date was 20-7-20. Some strains were removed due to poor metadata or sequencing statistics (<= 50000 'Low Quality Bases'), leaving a total of 668 ST95 isolates.

Note that the download was split into two chunks due to the job failing due to a connection issue partway through.

Downloading these assemblies was facilitated using a [script](!https://github.com/C-Connor/EnterobaseGenomeAssemblyDownload) modified from that of GitHub user C-Conner. Thanks C-Conner! Note that if you want to run this script you will need your own API key to access Enterobase.
```
cd ~/Data/tools/EnterobaseGenomeAssemblyDownload
nohup python2 EnterobaseGenomeAssemblyDownload.py \
    -d ecoli \
    -l ~/Data/Manuscripts/2020/AVC171/AVC171/delims/metadata_subset.txt \
    -o ~/Data/Manuscripts/2020/AVC171/AVC171/assemblies/ST95_enterobase_subset > \
    ~/Data/Manuscripts/2020/AVC171/AVC171/logs/Enterobase_dl_ST95_20-7-20.out 2> \
    ~/Data/Manuscripts/2020/AVC171/AVC171/logs/Enterobase_dl_ST95_20-7-20.err &


nohup python2 EnterobaseGenomeAssemblyDownload.py \
    -d ecoli -l ~/Data/Manuscripts/2020/AVC171/AVC171/delims/metadata_subset_minus_first_400.txt \
    -o ~/Data/Manuscripts/2020/AVC171/AVC171/assemblies/ST95_enterobase_subset > \
    ~/Data/Manuscripts/2020/AVC171/AVC171/logs/Enterobase_dl_ST95_20-7-20_minus_first_400.out 2> \
    ~/Data/Manuscripts/2020/AVC171/AVC171/logs/Enterobase_dl_ST95_20-7-20_minus_first_400.err &
```
___

---

# Configuration
Here we setup some variables which will be inherited by future steps.

```
# Define the pipelord directory where our pipelines are kept
PIPELORD_DIR=~/Data/pipelord/

# Define the project directory where we will send our outputs
PROJ_DIR=/projects/AusGEM/Users/Max/Manuscripts/ST95

# Make the appropriate subdirectories for our analysis
cd $PROJ_DIR
mkdir data
mkdir data/assemblies
mkdir output
mkdir logs
mkdir snakemake
mkdir scripts
mkdir misc
mkdir delims


#
dt=$(date '+%d/%m/%Y %H:%M:%S')
MASTER_CONF=/projects/AusGEM/Users/Max/Manuscripts/ST95/snakemake/masterconfigs/masterconfig_ST95_all.yaml
#MASTER_CONF=/projects/AusGEM/Users/Max/Manuscripts/ST95/snakemake/masterconfigs/masterconfig_ST95_HC50.yaml
#MASTER_CONF=/projects/AusGEM/Users/Max/Manuscripts/ST95/snakemake/masterconfigs/masterconfig_ST95_HC50_close.yaml

```

You must also modify the master_config.yaml **configuration** file.




___

# Job Running

## Genotypic analysis and plasmid screening with ABRicate


ABRicate version: 2.1.1
Snakefile:

```
# Set the task variable for our output names
TASK=abricate

# Change to the pipelord directory
cd ${PIPELORD_DIR}/${TASK}lord

# Create a directory for our task output log
mkdir ${PROJ_DIR}/snakemake/${TASK}

#Change Snakefile config variable to our master config
perl -p -i -e "s@^configfile.*@configfile: \"${MASTER_CONF}\"@g" Snakefile

# Copy the Snakefile and Environment yaml/s to our project directory
cp Snakefile config/abricate.yaml ${PROJ_DIR}/snakemake/${TASK}

# Create a directory for our task output log
mkdir ${PROJ_DIR}/logs/${TASK}

# Run task
nohup snakemake -j --use-conda  -p > ${PROJ_DIR}/logs/${TASK}/nohup_${TASK}.err 2> ${PROJ_DIR}/logs/${TASK}/nohup_${TASK}.out &

# Save our job ID
PROCESS_ID=$!
echo "$dt" "$PWD" "JOB_ID =" "$PROCESS_ID" >> ${PROJ_DIR}/logs/${TASK}/JOB_IDs
```

___

## Pangenomic analysis

Prokka version: ___

Roary version: ____

Snakefile:
```
# Set the task variable for our output names
TASK=pangenome

# Change to the pipelord directory
cd ${PIPELORD_DIR}/${TASK}lord

# Create a directory for our task output log
mkdir ${PROJ_DIR}/snakemake/${TASK}

#Change Snakefile config variable to our master config
perl -p -i -e "s@^configfile.*@configfile: \"${MASTER_CONF}\"@g" Snakefile

# Copy the Snakefile and Environment yaml/s to our project directory
cp Snakefile config/prokka.yaml config/roary.yaml ${PROJ_DIR}/snakemake/${TASK}

# Create a directory for our task output log
mkdir ${PROJ_DIR}/logs/${TASK}

# Run task
nohup snakemake -j --use-conda  -p > ${PROJ_DIR}/logs/${TASK}/nohup_${TASK}.err 2> ${PROJ_DIR}/logs/${TASK}/nohup_${TASK}.out &

# Save our job ID
PROCESS_ID=$!
echo "$dt" "$PWD" "JOB_ID =" "$PROCESS_ID" >> ${PROJ_DIR}/logs/${TASK}/JOB_IDs
```


## AMR-associated SNP variant analysis

Pointfinder version: ____

Roary version: ____

Snakefile:
```
# Set the task variable for our output names
TASK=pointfinder

# Change to the pipelord directory
cd ${PIPELORD_DIR}/${TASK}lord

# Create a directory for our task output log
mkdir ${PROJ_DIR}/snakemake/${TASK}

#Change Snakefile config variable to our master config
perl -p -i -e "s@^configfile.*@configfile: \"${MASTER_CONF}\"@g" Snakefile

# Copy the Snakefile and Environment yaml/s to our project directory
cp Snakefile config/pointfinder.yaml ${PROJ_DIR}/snakemake/${TASK}

# Create a directory for our task output log
mkdir ${PROJ_DIR}/logs/${TASK}

# Run task
nohup snakemake -j --use-conda  -p > ${PROJ_DIR}/logs/${TASK}/nohup_${TASK}.err 2> ${PROJ_DIR}/logs/${TASK}/nohup_${TASK}.out &

# Save our job ID
PROCESS_ID=$!
echo "$dt" "$PWD" "JOB_ID =" "$PROCESS_ID" >> ${PROJ_DIR}/logs/${TASK}/JOB_IDs
```

## Phylogenetic SNP analysis

Snippy version: ____

Gubbins version: ____

SNP-sites version: ____

SNP-dists version: ____

Snakefile:
```
# Set the task variable for our output names
TASK=snippy

# Change to the pipelord directory
cd ${PIPELORD_DIR}/${TASK}lord

mv /projects/AusGEM/Users/Max/Manuscripts/ST95/data/assemblies/ST95_all/AVC171.fasta /projects/AusGEM/Users/Max/Manuscripts/ST95/data/assemblies/hide

# Create a directory for our task output log
mkdir ${PROJ_DIR}/snakemake/${TASK}

#Change Snakefile config variable to our master config
perl -p -i -e "s@^configfile.*@configfile: \"${MASTER_CONF}\"@g" Snakefile_assemblies

# Copy the Snakefile and Environment yaml/s to our project directory
cp Snakefile config/*.yaml ${PROJ_DIR}/snakemake/${TASK}

# Create a directory for our task output log
mkdir ${PROJ_DIR}/logs/${TASK}

# Run task
nohup snakemake -j -s Snakefile_assemblies --use-conda  -p > ${PROJ_DIR}/logs/${TASK}/nohup_${TASK}.err 2> ${PROJ_DIR}/logs/${TASK}/nohup_${TASK}.out &

# Save our job ID
PROCESS_ID=$!
echo "$dt" "$PWD" "JOB_ID =" "$PROCESS_ID" >> ${PROJ_DIR}/logs/${TASK}/JOB_IDs
```

## BIGSI analysis

Snippy version: ____


Snakefile:
```
# Set the task variable for our output names
TASK=BIGSI

# Change to the pipelord directory
cd ${PIPELORD_DIR}/${TASK}lord

# Create a directory for our task output log
mkdir ${PROJ_DIR}/snakemake/${TASK}

#Change Snakefile config variable to our master config
perl -p -i -e "s@^configfile.*@configfile: \"${MASTER_CONF}\"@g" Snakefile_assemblies

# Copy the Snakefile and Environment yaml/s to our project directory
cp Snakefile config/*.yaml ${PROJ_DIR}/snakemake/${TASK}

# Create a directory for our task output log
mkdir ${PROJ_DIR}/logs/${TASK}

# Run task
nohup snakemake -j --use-conda  -p > ${PROJ_DIR}/logs/${TASK}/nohup_${TASK}.err 2> ${PROJ_DIR}/logs/${TASK}/nohup_${TASK}.out &

# Save our job ID
PROCESS_ID=$!
echo "$dt" "$PWD" "JOB_ID =" "$PROCESS_ID" >> ${PROJ_DIR}/logs/${TASK}/JOB_IDs

```

# Downstream Analysis

## Abricate

Output files must be concatenated from abricate. Reference plasmid databases are separated into their own files for downstream processing with plasmid_mapR.R.

```
PROJ_DIR=/projects/AusGEM/Users/Max/Manuscripts/ST95

cd ${PROJ_DIR}/output/ST95_all/abricate/

cat card/*.tab colV_zoetis/*.tab dfrA5_848/*.tab EC_custom/*.tab ecoh/*.tab ISfinder_Feb_2020/*.tab plasmidfinder/*.tab vfdb/*.tab > genotype.txt

mkdir plasmids

mv pAPEC_O1_ColBM  pAPEC_O2_ColV  pBCE049_1  pCERC4 pSF_088_nores  pU1_F51_B10  pUTI89 plasmids

for f in plasmids/*; do cat ${f}/* > ${f}.txt; done

mv plasmids/*.txt .

```

## Pangenomic tree building

Next we need to build a tree from our core genome

```
source deactivate
source activate /home/malcummi/Data/pipelord/snippylord/.snakemake/conda/88bf0609 # snp_sites from snplord
snp-sites -c core_gene_alignment.aln > snp_sites/core_gene_alignment_snp_sites.aln
source deactivate

source activate /home/malcummi/Data/pipelord/snippylord/.snakemake/conda/4ead7ed0 # snp_dists from snplord
snp-dists -c full_aln/core_gene_alignment.aln > full_aln/core_gene_alignment.csv
snp-dists -c snp_sites/core_gene_alignment_snp_sites.aln > snp_sites/core_gene_alignment_snp_sites.csv
source deactivate



source activate iqtree
cd ../output

#This was run again on the snp_sites
nohup iqtree -s core_gene_alignment_snp_sites.aln -m MFP -bb 1000 -nt AUTO >iqtree_core_genome_aln.out 2>iqtree_core_genome_aln.err &
#JOBID=202926

#This was run again on the full alignment (rather than the snp_sites)
nohup iqtree -s ../Roary.out/core_gene_alignment.aln -m MFP -bb 1000 -nt AUTO >iqtree_core_genome_aln.out 2>iqtree_core_genome_aln.err &
#JOBID=257443
```

# Extra Analysis

## Plasmid database

Next we screened Sistrom et al (2018) plasmid database for pMLST and abricate genes in a hunt for IS, AMR gene and pMLST associations to see if we can find pUTI-like plasmids with AMR genes

### Database setup

```
# Download Sistrom et al (2018) plasmid database
wget https://datadryad.org/stash/downloads/file_stream/147140 -o plasmid_db.fasta

# Split the multifasta into single fastas
# Credit: https://gist.github.com/astatham/621901 - Thanks Astatham!
cat plasmid_db.fasta | awk '{
        if (substr($0, 1, 1)==">") {filename=(substr($0,2) ".fasta")}
        print $0 > filename
}'
```

### Abricate screening

```
# Call wrapper script to allow conda calling
my_python_paths

# Activate snakemake environment
source activate snakemake

# Reset working directory
PROJ_DIR=/projects/AusGEM/Users/Max/Manuscripts/ST95

# Define the pipelord directory where our pipelines are kept
PIPELORD_DIR=~/Data/pipelord/

# Change the config file
MASTER_CONF=/projects/AusGEM/Users/Max/Manuscripts/ST95/snakemake/masterconfigs/masterconfig_plasdb.yaml

# Set the task variable for our output names
TASK=abricate

# Change to the pipelord directory
cd ${PIPELORD_DIR}/${TASK}lord

# Create a directory for our task output log
mkdir ${PROJ_DIR}/snakemake/${TASK}_plasdb

#Change Snakefile config variable to our master config
perl -p -i -e "s@^configfile.*@configfile: \"${MASTER_CONF}\"@g" Snakefile

# Copy the Snakefile and Environment yaml/s to our project directory
cp Snakefile config/abricate.yaml ${PROJ_DIR}/snakemake/${TASK}

# Create a directory for our task output log
mkdir ${PROJ_DIR}/logs/${TASK}

# Run task
nohup snakemake -j --use-conda  -p > ${PROJ_DIR}/logs/${TASK}/nohup_${TASK}_plasdb.err 2> ${PROJ_DIR}/logs/${TASK}/nohup_${TASK}_plasdb.out &

# Save our job ID
PROCESS_ID=$!
echo "$dt" "$PWD" "JOB_ID =" "$PROCESS_ID - plasdb" >> ${PROJ_DIR}/logs/${TASK}/JOB_IDs
```

### Abricate post processing

```
cd ${PROJ_DIR}/output/plasdb/abricate/

cat card/*.tab colV_zoetis/*.tab dfrA5_848/*.tab EC_custom/*.tab ecoh/*.tab ISfinder_Feb_2020/*.tab plasmidfinder/*.tab vfdb/*.tab > plas_genotype.txt

mkdir plasmids

mv pAPEC_O1_ColBM  pAPEC_O2_ColV  pBCE049_1  pCERC4 pSF_088_nores  pU1_F51_B10  pUTI89 plasmids

for f in plasmids/*; do cat ${f}/* > ${f}.txt; done

mv plasmids/*.txt .

```

### pMLST

## pMLST analysis
Note that -j is disabled!

pMLST version:

pMLST yaml: /projects/AusGEM/Users/Max/Manuscripts/ST95/snakemake/pMLST/pMLST.yaml

```
# Call wrapper script to allow conda calling
my_python_paths

# Activate snakemake environment
source activate snakemake

# Change the config file
MASTER_CONF=/projects/AusGEM/Users/Max/Manuscripts/ST95/snakemake/masterconfigs/masterconfig_plasdb.yaml

# Set the task variable for our output names
TASK=pMLST

# Change to the pipelord directory
cd ${PIPELORD_DIR}/${TASK}lord

# Create a directory for our task output log
mkdir ${PROJ_DIR}/snakemake/${TASK}

#Change Snakefile config variable to our master config
perl -p -i -e "s@^configfile.*@configfile: \"${MASTER_CONF}\"@g" Snakefile

# Copy the Snakefile and Environment yaml/s to our project directory
cp Snakefile config/*.yaml ${PROJ_DIR}/snakemake/${TASK}

# Create a directory for our task output log
mkdir ${PROJ_DIR}/logs/${TASK}

# Run Task -j disabled!
#snakemake -np
nohup snakemake --use-conda -j -p > ${PROJ_DIR}/logs/${TASK}/nohup_${TASK}_pMLST_plasdb.err 2> ${PROJ_DIR}/logs/${TASK}/nohup_${TASK}_pMLST_plasdb.out &

# Save our job ID
PROCESS_ID=$!
echo "$dt" "$PWD" "JOB_ID =" "$PROCESS_ID plasdb pMLST" >> ${PROJ_DIR}/logs/${TASK}/JOB_IDs
```
