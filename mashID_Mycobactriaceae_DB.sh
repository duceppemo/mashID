# Mycobactriaceae- 2025-02-20 - 26,101 genomes

# TaxID = 1762

# Base directory
export baseDir="/mnt/50tb/db/mashID/Mycobactriaceae_2025-02-20"
[ -d "$baseDir" ] || mkdir -p "$baseDir"


# Download the table with accessions from
# https://www.ncbi.nlm.nih.gov/datasets/genome/?taxon=1762

# Remove first line (header) and only keep the first colomn (accession numbers)
cat ~/Downloads/ncbi_dataset.tsv | tail -n+2 | cut -f 1 > "${baseDir}"/acc_list.txt

# Download directory for ncbi zip files
[ -d "${baseDir}"/ncbi_dataset ] || mkdir -p "${baseDir}"/ncbi_dataset

conda activate ncbi

# Download all the genomes
function get_genome()
{
    datasets download genome accession \
        "$1" \
        --api-key c26302fcf16ce79bd64d23a39473e3e21e08 \
        --filename  "${baseDir}"/ncbi_dataset/"$1".zip
}

export -f get_genome

cat "${baseDir}"/acc_list.txt \
| xargs -P $(nproc) -I {} bash -c 'get_genome "$@"' _ {}
# | parallel --bar --env baseDir --env get_genome 'get_genome {}'


# Unzip the files
function decompress(){
    sample=$(basename "$1" ".zip")
    unzip "$1" -d "${baseDir}"/ncbi_dataset/"$sample"
}

export -f decompress

find "${baseDir}"/ncbi_dataset -type f -name "*.zip" -print0 \
| xargs -0 -P $(nproc) -I {} bash -c 'decompress "$@"' _ {}
# | parallel --bar --env baseDir --env decompress 'decompress {}'


# Move the fasta files to a single folder
mkdir "${baseDir}"/fasta
find "${baseDir}"/ncbi_dataset -type f -name "*.fna" \
| parallel --bar "mv {} "${baseDir}"/fasta"

# Delete zip files
rm -rf "${baseDir}"/ncbi_dataset


# Rename the fasta files based on species info in the first header
# Bin the fasta by species
python ~/scripts/rename_and_bin_fasta.py \
    -i "${baseDir}"/fasta \
    -o "${baseDir}"/binned_by_species

rm -rf "${baseDir}"/fasta

conda deactivate  # ncbi

# Derep by species
conda activate mashID

# Install derep tool
# cd ~/prog
# git clone https://github.com/rrwick/Assembly-dereplicator


function derep_assebmlies() {
    species=$(basename "$1")
    output_forlder="$2"/"${species}"

    # v0.3.2
    # Derep at 99.9% homology
    python ~/prog/Assembly-dereplicator/dereplicator.py \
        --distance 0.001 \
        --sketch_size 10000 \
        --threads $(nproc) \
        "$1" "$output_forlder"
}

export -f derep_assebmlies

find "${baseDir}"/binned_by_species -mindepth 1 -type d -print0 \
    | xargs -0 -P $(nproc) -I {} bash -c 'derep_assebmlies "$@" "${baseDir}"/derep0.001' _ {}

# Create mashID database
python ~/prog/mashID/make_mashID_db.py \
    -i "${baseDir}"/derep0.001 \
    -o "$baseDir" \
    -p "mycobacteriaceae_2025-02-20" \
    -s 10000 \
    -k 21 \
    -t $(nproc)

