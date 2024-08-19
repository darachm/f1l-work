import marimo

__generated_with = "0.8.0"
app = marimo.App(width="medium")


@app.cell
def __():
    import marimo as mo
    return mo,


@app.cell
def __(mo):
    mo.md(
        r"""
        To reproduce this:

        1. setup a virtual environment in a convienent directory
        1. install `requirements.txt` (the packages used, plus `marimo`)
        1. launch out `marimo edit` from the repo's base directory
        1. edit the below variables if you put data and metadata in a different place
        """
    )
    return


@app.cell
def __():
    data_dir = "./data/kinker2020/"
    scratch_dir = "./tmp/"
    return data_dir, scratch_dir


@app.cell
def __(data_dir, mo, scratch_dir):
    mo.md(
        f"""## Read in the data

        Gonna look in directory `{data_dir}` for the data,
        saving scratch work in `{scratch_dir}`.
        """
    )
    return


@app.cell
def __():
    import os
    import numpy as np
    import pandas as pd
    import scipy
    import anndata
    import scanpy as sc
    import seaborn as sns
    import re
    import itertools
    import csv
    return anndata, csv, itertools, np, os, pd, re, sc, scipy, sns


@app.cell
def __(mo):
    mo.md(
        r"""
        First, the metadata file. This appears to have the cell ID, the cell line, the pool, 
        a cancer type, some statistics, and some computed scores. 
        Here's the head of it:
        """
    )
    return


@app.cell
def __(data_dir, pd):
    meta = (
        pd.read_csv(data_dir + "Metadata.txt", sep="\t", low_memory=False)
        .drop([0], axis=0)
        .rename(
            columns={
                "NAME": "CellID",
                "Cell_line": "CellLine",
                "Pool_ID": "Pool",
                "Cancer_type": "Indication",
            }
        )
        .set_index("CellID")
        .convert_dtypes()
    )
    meta
    return meta,


@app.cell
def __(mo):
    mo.md(
        r"""
        Then for reading the counts file, here we have a little function just to
        read the top three lines in as lists, then to read each line as a split up
        list of ints:
        """
    )
    return


@app.cell
def __(pd, re, scipy):
    def proc_a_chunk(counts_chunks, keep_index, verbose=True):
        """
        We take the next chunk, just keep the columns as indexed in keep_index,
        then return the gene IDs and a CSR matrix of it.

        Params:
            counts_chunks: pandas read_csv iterable returned by specifying chunksize
            keep_index: logical list of which columns (cells) to keep
            verbose: do you want it to report progress to STDOUT?
        """
        if verbose:
            print(f"Handing chunk at {counts_chunks._currow}")
        try:
            chunk = counts_chunks.__next__().loc[:, keep_index]
            return (list(chunk.index), scipy.sparse.csr_matrix(chunk.transpose()))
        except StopIteration:
            if verbose:
                print(f"Finished taking {counts_chunks._currow} rows.")
            return (None, None)


    def read_in_umi_counts(filepath, meta, nrows=None, chunksize=100):
        """
        This is for reading that awful format that the TSV is in.
        We expect the first line is the IDs for each cell,
        the second line is each cell's assigned cell line,
        the third is the pool it is in, then every subsequent line
        is a gene ID and then the UMI counts for each cell.

        **shudder**
        """

        with open(filepath, "r") as f:
            # First we read in the first three lines and drop the
            # leading empty space
            cell_id = re.split("\t", f.readline())[1:]
            cell_line = re.split("\t", f.readline())[1:]
            cell_pool = re.split("\t", f.readline())[1:]

            # Then we ask if each cell ID is in the index of the
            # metadata, ie if we should keep it
            keep_index = [x in meta.index for x in cell_id]

            # We use pandas read_csv to open up the remainder of the
            # TSV as a chunksize chunk reader, with first column as
            # index (which is the gene IDs).
            counts_chunks = pd.read_csv(
                f,
                sep="\t",
                header=None,
                index_col=0,
                chunksize=chunksize,
                nrows=nrows,
            )

            # We grab the first chunk
            (return_gene_ids, return_object) = proc_a_chunk(
                counts_chunks, keep_index
            )

            for chunk in counts_chunks:
                (more_return_gene_ids, more_return_object) = proc_a_chunk(
                    counts_chunks, keep_index
                )

                if more_return_gene_ids is not None:
                    return_gene_ids.extend(more_return_gene_ids)
                    return_object = scipy.sparse.hstack(
                        [return_object, more_return_object]
                    )

        return (cell_id, return_gene_ids, return_object)
    return proc_a_chunk, read_in_umi_counts


@app.cell
def __(data_dir, meta, read_in_umi_counts):
    (cell_id, gene_ids, gene_counts) = read_in_umi_counts(
        data_dir + "UMIcount_data.txt", meta=meta, chunksize=1000
    )
    return cell_id, gene_counts, gene_ids


@app.cell
def __(cell_id, gene_counts, gene_ids, mo):
    mo.md(f"""
        ...and report that we:

        Read in {len(cell_id)} cell ids that look like: [ {", ".join(cell_id[0:3])}, ... ]

        Read in gene count matrix as sparse matrix of shape {gene_counts.shape},
        with counts for {len(gene_ids)} different gene IDs.
        """)
    return


@app.cell
def __(mo):
    mo.md(
        r"""
        ## Generate the AnnData object

        First, we take that `gene_counts` dict, we tack on the cell IDs list,
        and we turn that into a pandas DataFrame, set the index, and `reindex`
        to just use the ones where it matches the metadata.
        """
    )
    return


@app.cell
def __(gene_counts):
    gene_counts[0:5, :]
    return


@app.cell
def __(mo):
    mo.md(r"""Then we can convert it to the AnnData type""")
    return


@app.cell
def __(anndata, gene_counts, gene_ids, meta, np):
    ann_datar = anndata.AnnData(
        X=gene_counts,
        obs=meta.astype(
            {
                "CellLine": np.str_,
                "Pool": np.str_,
                "Indication": np.str_,
                "Genes_expressed": np.int32,
                "Discrete_cluster_minpts5_eps1.8": np.str_,
                "Discrete_cluster_minpts5_eps1.5": np.str_,
                "Discrete_cluster_minpts5_eps1.2": np.str_,
                "CNA_subclone": np.str_,
                "SkinPig_score": np.float32,
                "EMTI_score": np.float32,
                "EMTII_score": np.float32,
                "EMTIII_score": np.float32,
                "IFNResp_score": np.float32,
                "p53Sen_score": np.float32,
                "EpiSen_score": np.float32,
                "StressResp_score": np.float32,
                "ProtMatu_score": np.float32,
                "ProtDegra_score": np.float32,
                "G1/S_score": np.float32,
                "G2/M_score": np.float32,
            }
        ),
        var=gene_ids,
    )
    ann_datar.var.drop(columns=[0], inplace=True)
    ann_datar
    return ann_datar,


@app.cell
def __(mo):
    mo.md(r"""Then we apply the same filters that Dean Lee had in his book:""")
    return


@app.cell
def __(ann_datar, sc):
    sc.pp.filter_genes(ann_datar, min_cells=10)
    sc.pp.filter_cells(ann_datar, min_genes=200)
    ann_datar
    return


@app.cell
def __(mo):
    mo.md(r"""...and save it:""")
    return


@app.cell
def __(ann_datar, scratch_dir):
    ann_datar.write(scratch_dir + "kinker_anndata.h5ad")
    return


@app.cell
def __():
    return


if __name__ == "__main__":
    app.run()
