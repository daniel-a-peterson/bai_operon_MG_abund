# Step 1 - Map your reads to the bai reference sequences:

Here is a simple bash script to map all the `.fasta` files in your working directory to the provided reference base sequences (`derep_bai_genes_prot.fasta`) using diamond blastx. You do not need to use this script exactly, but you need to create one file for each sample that gives the diamond output from mapping that metagenome sample's nucleotide reads against the reference bai protein sequences.

```
diamond makedb --in derep_bai_genes_prot_renamed.fasta -d bai_PROT

for CUR_FILE in *.fasta;
do
    BASENAME=$(basename "$CUR_FILE")
    SAMPLE_NAME="${BASENAME%.*}"
    diamond blastx -d bai_PROT -q $CUR_FILE -o ${SAMPLE_NAME}_bai_hits.tsv --id 90 --query-cover 50 --threads 8 --max-target-seqs 10 -b 8.0
done
```
 
# Step 2 - Run bai_predictor.py script to estimate overall bai operon abundance and individual bai operon cluster abundance in each sample:

In addition to a diamond mapping result file for each sample generated in step 1, you will need to create a comma-delimited text file with metadata for each sample. This table must include a column named 'filename' that gives the relevant diamond mapping results file for each sample and a column named 'spots' that gives the total number of raw reads in each sample.

Running the script requires an environment with python3, numpy and pandas.

```
python bai_predictor.py --input_metadata YOUR_INPUT_METADATA.csv --input_diamond *_bai_hits.tsv 
```

# Outputs:

BA_transformation_capacity.csv: Proportion of bai operon genes detected overall and for each cluster
BA_transformation_abun.csv: Abundance of bai operon genes overall and for each cluster