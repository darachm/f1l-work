import marimo

__generated_with = "0.8.0"
app = marimo.App(width="medium")


@app.cell
def __():
    import marimo as mo
    import os
    import numpy as np
    import pandas as pd
    import scipy
    import anndata
    import scanpy as sc
    import seaborn as sns
    import matplotlib.pyplot as plt
    return anndata, mo, np, os, pd, plt, sc, scipy, sns


@app.cell
def __():
    scratch_dir = "./tmp/"
    return scratch_dir,


@app.cell
def __(sc, scratch_dir):
    adata = sc.read_h5ad(scratch_dir+'kinker_anndata.h5ad')
    return adata,


@app.cell
def __(adata):
    adata
    return


@app.cell
def __(sc):
    def BasicScanpyPreprocessing(adata, n_top_genes=2000, n_neighbors=10, n_pcs=40, random_state=20):
        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)
        sc.pp.highly_variable_genes(adata, n_top_genes=n_top_genes, flavor='seurat')
        #sc.pp.scale(adata, max_value=10)
        sc.tl.pca(adata, svd_solver='arpack')
        sc.pp.neighbors(adata, n_neighbors=n_neighbors, n_pcs=n_pcs, random_state=random_state)
        sc.tl.umap(adata, random_state=random_state)
        sc.tl.leiden(adata, random_state=random_state)
        return adata
    return BasicScanpyPreprocessing,


@app.cell
def __(BasicScanpyPreprocessing, adata):
    adata_pre = BasicScanpyPreprocessing(adata, n_top_genes=2000, n_neighbors=10, n_pcs=40, random_state=20)
    return adata_pre,


@app.cell
def __(adata_pre, sc):
    sc.pl.umap(adata_pre, color=['Indication'])
    return


if __name__ == "__main__":
    app.run()
