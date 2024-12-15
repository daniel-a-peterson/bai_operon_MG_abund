"""bai abundance predictor tool

Predicts bile acid transformation capacity of a stool sample from shotgun metagenomic sequencing reads

Parameters
----------
input_metadata : Path
    Path to the input metadata, in CSV
input_diamond : Path(s)
    Paths to the input diamond mapping results files; can include director with ./*
"""

import math, re, sys
import numpy as np
import pandas as pd
import logging

from argparse import ArgumentParser, RawDescriptionHelpFormatter

__version__ = '0.1'


logger = logging.getLogger('predictor')
logger.setLevel(logging.INFO)
handler = logging.StreamHandler(sys.stdout)
formatter = logging.Formatter(
    "%(asctime)s - %(levelname)s  %(lineno)s - %(message)s"
)
handler.setFormatter(formatter)
logger.addHandler(handler)

# parameters

constant1 = 1e7
constant2 = 25
constant3 = 90

gene_list = ["BaiA", "BaiB", "BaiCD", "BaiE", "BaiF", "BaiG", "BaiH", "BaiI"]

#

def main(*args):

    parser = ArgumentParser(usage = "script",
                            description=__doc__,
                            formatter_class=RawDescriptionHelpFormatter)
    parser.add_argument("--input_metadata", action="store",
                        dest="input_metadata", help="input metadata csv. This table must include a column named 'filename' that gives the relevant diamond mapping results file for each sample and a column named 'spots' that gives the total number of raw reads in each sample", required=True)
    parser.add_argument("--input_diamond", action="store", nargs='+', 
                        dest="input_diamond", help="list of all mapping result files (eg *_mapping_hits.tsv)", required=True)
    parser.add_argument("-verbose", action='store_true')
    options = parser.parse_args(*args)

    logger.info(f"Running bai-estimator on {options.input_metadata} and {options.input_diamond}")

    if options.verbose == True:
        logger.setLevel(logging.DEBUG)
        logger.debug(f"Version is {__version__}")
        
    # 

    read_num_df = pd.read_csv(options.input_metadata)
    
    if "filename" not in read_num_df.columns.values or "spots" not in read_num_df.columns.values:
        parser.print_help()
        sys.exit("\nERROR in metadata file, must have 'filename' and 'spots' columns")
    
    logger.debug(read_num_df.shape)
    
    all_samples_list = read_num_df["filename"].unique().tolist()

    #
    
    best_hits_list = []
    
    file_list = options.input_diamond

    for cur_file in file_list:
        short_filename = re.sub(".*/","",cur_file)
        cur_hits_df = pd.read_csv(
            cur_file,
            sep="\t",
            header=None,
            names=[
                "qseqid",
                "sseqid",
                "pident",
                "length",
                "mismatch",
                "gapopen",
                "qstart",
                "qend",
                "sstart",
                "send",
                "evalue",
                "bitscore",
            ],
        )
        cur_hits_df = cur_hits_df[cur_hits_df["length"] >= constant2].copy()
        best_hits = cur_hits_df.sort_values("bitscore", ascending=False).drop_duplicates(
            "qseqid"
        )
        best_hits["gene"] = best_hits["sseqid"].apply(lambda x: x.split("%")[0])
        best_hits["filename"] = short_filename
        best_hits_list.append(
            best_hits[["filename", "gene", "qseqid", "sseqid", "pident", "length"]]
        )


    all_best_hits_df = pd.concat(best_hits_list, ignore_index=True)
    all_best_hits_df = all_best_hits_df[all_best_hits_df["gene"].isin(gene_list)].copy()
    logger.debug(all_best_hits_df.shape)
    
    # 

    placeholder_df = pd.DataFrame(
        {
            "filename": all_samples_list,
            "gene": "placeholder",
            "length": 0,
            "pident": 0,
            "qseqid": "placeholder",
            "sseqid": "placeholder%placeholder",
        }
    )
    
    # 

    all_best_hits_df = pd.concat([all_best_hits_df,placeholder_df])
    logger.debug(all_best_hits_df.shape)

    all_best_hits_df = all_best_hits_df[
        (all_best_hits_df["pident"] >= constant3) | (all_best_hits_df["gene"] == "placeholder")
    ].copy()
    logger.debug(all_best_hits_df["filename"].nunique())
    logger.debug(all_best_hits_df.shape)

    abun_df = all_best_hits_df.groupby(["filename", "gene"])["qseqid"].count().reset_index()
    abun_df = abun_df.merge(read_num_df[["filename", "spots"]])
    abun_df = abun_df[abun_df["spots"] >= constant1].copy()
    
    abun_df["norm_cov"] = abun_df["qseqid"] / abun_df["spots"].apply(float)
    abun_df["norm_cov_log10"] = abun_df["norm_cov"].apply(math.log10)
    abun_df["reads_log10"] = abun_df["qseqid"].apply(math.log10)
    logger.debug(abun_df.shape)

    abun_pivot_df = abun_df.pivot_table(
        index=["filename"], columns="gene", values="norm_cov_log10",
    ).fillna(-9)
    
    for cur_gene in gene_list:
        if cur_gene not in abun_pivot_df.columns.values:
            abun_pivot_df[cur_gene] = -9
    
    abun_pivot_df["abun_mean"] = abun_pivot_df.drop(columns="placeholder").apply(
        np.mean, axis="columns"
    )
    
    logger.debug(abun_pivot_df.shape)

    abun_total_dict = dict(abun_pivot_df['abun_mean'])
    
    # 
    
    gene_df = all_best_hits_df.groupby(["filename"])["gene"].nunique().reset_index()
    gene_df['total_genes'] = gene_df['gene'] - 1 # account for placeholder gene
    gene_df['total_prop_capacity'] = gene_df['total_genes'] / float(len(gene_list))
    gene_df = gene_df.merge(read_num_df[["filename", "spots"]])
    gene_df = gene_df[gene_df["spots"] >= constant1].copy()
    
    logger.debug(gene_df.shape)
    
    # 
    
    all_best_hits_df["cluster_name"] = all_best_hits_df["sseqid"].apply(
        lambda x: x.split("%")[1]
    )

    logger.debug(all_best_hits_df.shape)


    tax_abun_df = (
        all_best_hits_df.groupby(["filename", "cluster_name", "gene"])["qseqid"]
        .count()
        .reset_index()
    )

    logger.debug(tax_abun_df.shape)

    # 

    placeholder_list = []

    for cur_sample in all_samples_list:
        new_df = pd.DataFrame(
            {
                "filename": cur_sample,
                "cluster_name": "placeholder",
                "gene": gene_list,
                "qseqid": 0,
            }
        )
        placeholder_list.append(new_df)

    placeholder_df = pd.concat(placeholder_list)
    logger.debug(placeholder_df.shape)
    
    # 

    tax_abun_pivot_df = (
        pd.concat([tax_abun_df,pd.concat(placeholder_list)]).pivot_table(index=["filename", "gene"], columns="cluster_name", values="qseqid").fillna(0)
    )
    tax_abun_pivot_df["gene_total"] = tax_abun_pivot_df.apply(sum, axis="columns")
    tax_abun_pivot_df = tax_abun_pivot_df.reset_index()
    logger.debug(tax_abun_pivot_df.shape)

    tax_abun_melted_df = tax_abun_pivot_df.melt(
        id_vars=["filename", "gene", "gene_total"],
        var_name="cluster_name",
        value_name="reads",
    )
    tax_abun_melted_df = tax_abun_melted_df.merge(read_num_df[["filename", "spots"]])
    tax_abun_melted_df["norm_cov"] = tax_abun_melted_df["reads"] / tax_abun_melted_df["spots"].apply(
        float
    )
    tax_abun_melted_df.loc[tax_abun_melted_df["norm_cov"] == 0, "norm_cov"] = 1e-9
    tax_abun_melted_df["norm_cov_log10"] = tax_abun_melted_df["norm_cov"].apply(math.log10)
    logger.debug(tax_abun_melted_df.shape)

    tax_abun_melted_df = tax_abun_melted_df[tax_abun_melted_df["spots"] >= constant1].copy()

    logger.debug(tax_abun_melted_df.shape)

    tax_abun_mean_df = (
        tax_abun_melted_df.groupby(["filename", "cluster_name"])["norm_cov_log10"]
        .mean()
        .reset_index()
    )
    
    tax_abun_pivot2_df = tax_abun_mean_df.pivot_table(
        index=["filename"], columns="cluster_name", values="norm_cov_log10",
    ).fillna(-9).drop(columns=['placeholder']).reset_index()
    logger.debug(tax_abun_pivot2_df.shape)
    
    tax_abun_pivot2_df['total_abun'] = tax_abun_pivot2_df['filename'].map(abun_total_dict)
    
    tax_abun_pivot2_df.to_csv(
        "BA_transformation_abun.csv",
        index=False,
    )
    
    # 

    tax_gene_df = all_best_hits_df.groupby(["filename",'cluster_name'])["gene"].nunique().reset_index()
    tax_gene_df['prop_capacity'] = tax_gene_df['gene'] / float(len(gene_list))
    tax_gene_df = tax_gene_df.merge(read_num_df[["filename", "spots"]])
    tax_gene_df = tax_gene_df[tax_gene_df["spots"] >= constant1].copy()
    logger.debug(tax_gene_df.shape)
    
    tax_gene_pivot_df = tax_gene_df.pivot_table(
        index=["filename"], columns="cluster_name", values="prop_capacity",
    ).fillna(0).reset_index()
    logger.debug(tax_gene_pivot_df.shape)
    
    all_genes_df = gene_df[['filename','total_prop_capacity']].merge(tax_gene_pivot_df, how='outer').fillna(0)
    
    all_genes_df.to_csv(
        "BA_transformation_capacity.csv",
        index=False,
    )


if __name__ == '__main__':
	main(sys.argv[1:])

