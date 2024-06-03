## Assembly with spades
```
cd /hpcstor6/scratch01/p/patrick.kearns/C_use_project
$HOME/SPAdes-3.15.5-Linux/bin/spades.py -t 48 --isolate -1 2023-08-15_256samples_130D.bam-1.fq -p 130D -2 2023-08-15_256samples_130D.bam-2.fq 
$HOME/SPAdes-3.15.5-Linux/bin/spades.py -t 48 --isolate -1 2023-08-15_256samples_P0601142931.bam-1.fq -o P0601142931 -2 2023-08-15_256samples_P0601142931.bam-2.fq
$HOME/SPAdes-3.15.5-Linux/bin/spades.py -t 48 --isolate -1 2023-08-15_256samples_RSM3.2.bam-1.fq -2 2023-08-15_256samples_RSM3.2.bam-1.fq -o RSM32
$HOME/SPAdes-3.15.5-Linux/bin/spades.py -t 48 --isolate -1  2023-08-15_256samples_66B.bam-1.fq -s 2023-08-15_256samples_66B.bam-2.fq -o 66B
$HOME/SPAdes-3.15.5-Linux/bin/spades.py -t 48 --isolate -1  2023-08-15_256samples_HP1C.bam-1.fq -2 2023-08-15_256samples_HP1C.bam-2.fq -o HP1C
$HOME/SPAdes-3.15.5-Linux/bin/spades.py -t 48 --isolate -1 2023-08-15_256samples_THA3.2.bam-1.fq -2 2023-08-15_256samples_THA3.2.bam-2.fq -o THA32
``

## Move files around
```
#rename the bins by folder name, assumes your in the folder where all the genomes are located
## Specify the extension 
target_extension=".fasta"

# Function to process files
process_file() {
    local filepath="$1"
    local foldername="$2"
    local filename=$(basename "$filepath")
    
    # Append folder name to the file
    new_filename="${foldername}_${filename}"
    
    # Rename the file
    mv "$filepath" "$(dirname $filepath)/$new_filename"
}

# Main loop to find and process files
find . -type f -name "*$target_extension" | while read filepath; do
    foldername=$(basename "$(dirname "$filepath")")
    process_file "$filepath" "$foldername"
done

#######find and copy all files into a new folder, assumes your are in the folder where all the genomes are located
#tell it the extension for each genome
mkdir aggregated_genomes
target_extension="_scaffolds.fasta"

# Specify the destination folder where the files will be copied
destination_folder="./aggregated_genomes"

# Create the destination folder if it doesn't exist
mkdir -p "$destination_folder"

# Main loop to find and copy files
find . -type f -name "*$target_extension" -exec cp {} "$destination_folder" \;

#remove extra files
cd aggregated_genomes/
rm *final*
 rm K77*
rm  misc_broken_scaffolds.fasta
cd /hpcstor6/scratch01/p/patrick.kearns/C_use_project
```

## CheckM
```
#run checkM
source activate checkm
conda activate checkm
export CHECKM_DATA_PATH=/home/patrick.kearns/checkm
mkdir checkM_results
checkm lineage_wf ./aggregated_genomes ./checkM_results -x fasta -t 48 --pplacer_threads 24 --tab_table 
conda deactivate
```
## Annotate with Prokka
```
#annotate with prokka
conda activate prokka
mkdir annotations
cd aggregated_genomes
for i in *.fasta
do
prokka $i --cpus 24 --outdir $i --force --outdir /hpcstor6/scratch01/p/patrick.kearns/genomes/Assemblies_fastq_checkM/fastq/Chryseo/annotations/$i
done
```
