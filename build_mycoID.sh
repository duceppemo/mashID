db_fodler=/media/36tb/db/mashID/mycoID

# Mbovis mashID db update
https://github.com/duceppemo/mashID

# Download RefSeq -> 10719
conda activate ncbi
bash ~/prog/ncbi/refseqDownloader.sh \
    -d bacteria \
    -q "Mycobacterium,Mycolicibacterium,Mycolicibacter,Mycolicibacillus,Mycobacteroides" \
    -u \
    -f fna \
    -o "$db_fodler" \
    -r

conda deactivate

# Remove bad downloads -> 19
rm "${db_fodler}"/*.gz

# Remove Mycobacterium sp. -> 290
find "$db_fodler" -name "*_sp_*" -exec rm {} \;

# Derep
conda activate mashID
python ~/prog/Assembly-Dereplicator/dereplicator.py \
    --distance 0.001 \
    --threads $(nproc) \
    "${db_fodler}"/fna \
    "${db_fodler}"/fna_derep_0.001

conda deactivate


# Add bovis ref genome
conda activate ncbi
bash get_assemblies_v2.sh \
    -q NC_002945.4 \
    -t refseq \
    -o "${db_fodler}"/fna_derep_0.001 \
    -r \
    -u

conda deactivate

conda activate mashID
python ~/prog/mashID/make_mashID_db.py \
    -i "${db_fodler}"/fna_derep_0.001 \
    -o /media/36tb/db/mashID \
    -p mycobacteria_2024-02-23.msh \
    -t $(nproc) \
    -s 10000 \
    -k 21


baseDir=/home/bioinfo/analyses/MBWGS051

python ~/prog/mashID/mashID.py \
    -i "${baseDir}"/fastq \
    -o "${baseDir}"mashID \
    -d /media/36tb/db/mashID/mycobacteria_2024-02-23.msh \
    -p 3 

