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
        I've adapted the jupyter notebook to a marimo notebook. That 
        format/runtime/environment has several advantages (and limitations!) and 
        I'm using this project as an opportunity to learn about using it.

        ## What else I modified

        Since I'm a bit limited on memory (16GB RAM on my cruncher computer), I tweaked
        the workflow to read in the TSV as chunks, converted each chunk to a sparse
        matrix in `scipy` to save on space, and then serially concatenated each chunk
        to slowly build the full matrix object (as well as the list of gene IDs).

        I've left the filters applied at this stage as the same, as I don't have an
        intuition for where to place those. Perhaps I ought to plot and take a look
        at how these are distributed, but I'll do that later if I do.

        ## How to run this

        1. setup a virtual environment in a conveinent directory
        1. install `requirements.txt` (the packages used, plus `marimo`)
        1. launch out `marimo edit` from the repo's base directory, likely with
            appropriate arguments for port and password
        1. select the notebook to open it, edit the directory variables below,
            and I believe you can run it all with the button in the lower right
            corner
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
    def proc_a_chunk(chunk, keep_index, verbose=True):
        """
        We take the next chunk, just keep the columns as indexed in keep_index,
        then return the gene IDs and a CSR matrix of it.

        Params:
            chunk: chunk from a pandas read_csv read chunkwise
            keep_index: logical list of which columns (cells) to keep
            verbose: do you want it to report progress to STDOUT?
        """
        try:
            chunk = chunk.loc[:, keep_index]
            return (list(chunk.index), scipy.sparse.csr_matrix(chunk.transpose()))
        except StopIteration:
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
            cell_id = re.split("\t", f.readline().rstrip())[1:]
            cell_line = re.split("\t", f.readline().rstrip())[1:]
            cell_pool = re.split("\t", f.readline().rstrip())[1:]

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
            chunk = counts_chunks.__next__()
            print(f"Handing chunk at {counts_chunks._currow}")
            (return_gene_ids, return_object) = proc_a_chunk(chunk, keep_index)

            for chunk in counts_chunks:
                print(f"Handing chunk at {counts_chunks._currow}")
                (more_return_gene_ids, more_return_object) = proc_a_chunk(
                    chunk, keep_index
                )

                if more_return_gene_ids is not None:
                    return_gene_ids.extend(more_return_gene_ids)
                    return_object = scipy.sparse.hstack(
                        [return_object, more_return_object]
                    )

            print(f"Finished taking {counts_chunks._currow} rows.")

        return (cell_id, return_gene_ids, return_object)
    return proc_a_chunk, read_in_umi_counts


@app.cell
def __(data_dir, meta, read_in_umi_counts):
    (cell_id, gene_ids, gene_counts) = read_in_umi_counts(
        data_dir + "UMIcount_data.txt", meta=meta, chunksize=1000
    )
    return cell_id, gene_counts, gene_ids


@app.cell
def __(cell_id, gene_counts, meta, mo):
    mo.md(f"""
        ...and report that we:

        Read in {len(cell_id)} cell ids that look like: [ {", ".join(cell_id[0:3])}, ... ],
        and filtered to keep the {len(meta.index)} cells for which we have the metadata.

        Read in gene count matrix as sparse matrix of shape {gene_counts.shape},
        with {gene_counts.shape[0]} counts for {gene_counts.shape[1]} different gene IDs.
        """)
    return


@app.cell
def __(cell_id, meta, pd):
    cell_id_index = pd.DataFrame({'CellID':[x for x in cell_id if x in meta.index]}).set_index('CellID').index
    meta_re = meta.reindex(index=cell_id_index)
    meta_re
    return cell_id_index, meta_re


@app.cell
def __(cell_id_index, meta_re):
    meta_re.index.equals(cell_id_index)
    # obviously
    return


@app.cell
def __(mo):
    mo.md(r"""Then we can convert it to the AnnData type""")
    return


@app.cell
def __(anndata, gene_counts, gene_ids, meta_re, np, pd):
    ann_datar = anndata.AnnData(
        X=gene_counts,
        obs=meta_re.astype(
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
        var=pd.DataFrame({'gene_names':gene_ids}).set_index('gene_names'),
    )
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


if __name__ == "__main__":
    app.run()
