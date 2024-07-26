---
author: Darach Miller
date: 2024-07-26
---

Dean Lee has tasked each of us with addressing the Key Scientific Question (KSQ)
of:

    The KSQ: Using available scRNA-seq data from cancer cell lines, how would
    you explore the use of the following FDA-approved antibody therapies in
    additional cancers?
    
        Trastuzumab: Targets HER2 and is used in the treatment of HER2-positive
            breast and gastric cancers.
        Bevacizumab: Targets VEGF and is used for a variety of cancers, 
            including colorectal, lung, glioblastoma, breast, liver, and 
            kidney cancer.
    
    Take the KSQ and break it down into its constituent components. Test your
    own understanding of each component. Find additional resources to bolster 
    your understanding.

[^1]

Below I'm going to break out those components as sections.

# Background

First I'm going to go into the therapies, then into cancer, then models/assays.

[^1]: This seems like the kind of question that a firm would like to know if 
they had rights to using this therapy, to open up new markets for selling a 
de-risked treatment that they already have experience with.
However, this strategy of generalization also has broader implications for
cancer biology. 
A purely empirical science is limited to only the experiments that
have been done, but the appropriate dynamic of theory and empiricism allows
scientists to take steps into the unknown. If we're clever we may know where
our foot will land, but surprises are important to build our models and allow
us to take better steps.


## Therapies listed

### Trastuzumab

The [NCI page for trastuzumab](https://www.cancer.gov/about-cancer/treatment/drugs/trastuzumab)
states that trastuzumab is a monoclonal antibody therapy, where that antibody
binds to HER2.
It further states that this blocks the growth signal, but also "stimulates the
immune system to kill cells that have a lot of HER2".
It is available alone for intravenous administration, or in combinations with 
other drugs or adjuvants.

shoulda read this first:
[The StatPearls for trastuzumab](https://www.ncbi.nlm.nih.gov/books/NBK532246/)
Trastuzumab deruxtecan is an antibody-drug conjugate to gain selective entry 
to cells expressing high HER2.


Trastuzumab is used to treat:

1. HER2 positive breast cancer that is high-risk or hormone receptor negative 
    (grows independent of hormones), in combination with other drugs that 
    inhibit cell growth.
1. Metastatic HER2 positive breast cancer along with another growth-inhibiting
    drug (a taxane), or alone after previous chemotherapy.
1. Adenocarcinoma of the stomach or gastroesophageal junction that has
    metastitized, in combination with other drugs that inhibit cell growth.

From its 
[DailyMed page](https://dailymed.nlm.nih.gov/dailymed/drugInfo.cfm?setid=492dbdb2-077e-4064-bff3-372d6af0a7a2),
we see that trastuzumab is a humanized IgG1 kappa monoclonal antibody that
"selectively binds with high affinity to the extracellular domain of" HER2,
and that this has been shown to mediate antibody-dependent cellular 
cytotoxicity (ADCC) on "HER2 overexpressing cancer cells".
For tested dosing regimes, available concentration in serum at steady-state 
should be [roughly 33-179 ug/mL](https://dailymed.nlm.nih.gov/dailymed/drugInfo.cfm?setid=492dbdb2-077e-4064-bff3-372d6af0a7a2#Table8).

Drugs used or tested along with trastuzumab include:
- paclitaxel
- doxorubicin
- docetaxel
- carboplatin
- cisplatin
- capecitabine

"HER2 positive" appears to be assessed by immunohistochemistry or
gene copy amplification using DNA FISH.

* what are knock-on effects of HER2? what does binding do? does it block 
    something, or stimulate a halt signal? does basal just HER2 alone stimulate
    growth or does it require a signal external of the cell?


### Bevacizumab

[StatsPearls for bevzcizumab](https://www.ncbi.nlm.nih.gov/books/NBK482126/)



#### HER2 and HER2-positive cancers

The HER2 (or c-erbB2) proto-oncogene encodes a transmembrane receptor protein
of 185 kDa, which is structurally related to the epidermal growth factor
receptor. Herceptin has been shown, in both in vitro assays and in animals, to
inhibit the proliferation of human tumor cells that overexpress HER2.

Herceptin is a mediator of antibody-dependent cellular cytotoxicity (ADCC). In
vitro, Herceptin-mediated ADCC has been shown to be preferentially exerted on
HER2 overexpressing cancer cells compared with cancer cells that do not
overexpress HER2.

What is

c-erbB2

What is HER2 postiive

#### Usage

Breast and gastric cancers

### Bevacizumab

#### VEGF

#### Usage

"including colorectal, lung, glioblastoma, breast, liver, and kidney cancer."

## Generalization of drugs

FDA-approval, usage

to different cancers

## Cancer(s)

What it do?

### Cancer cell lines as a model

Dean Lee provides this link to [the Broad's CCLE](https://sites.broadinstitute.org/ccle/#:~:text=Cancer%20cell%20lines%20are%20the,and%20for%20defining%20drug%20efficacy.).


## Available datasets

Dean Lee specifies scRNAseq data from cancer cell lines, so we'll box it there.
But, we should be mindful to use an architecture that could expand to 
scRNAseq data from whole-animal/tissue collections, or other multi-modal 
datasets.

## Possible next steps

# What is to be done

## Process

Dean Lee provides links to 
[an article about Aviv Regev's presentation at AACR 2024](https://www.genengnews.com/topics/cancer/aacr-2024-aviv-regev-shows-how-single-cell-atlases-foster-new-axis-to-genentechs-drug-discovery/)
and
[a recording of a talk by Regev about cell atlases as roadmaps for cancer treatments](https://www.youtube.com/watch?v=Wk5QHySlMXU).
To summarize, Regev (amongst many many others) has been pushing Cell Atlas 
projects for at least a decade. 
Now we can enjoy the _fruits_ of this labor. 




### Methods

Dean Lee linked to 
[10x's website talking about their single-cell platform](https://www.10xgenomics.com/support/single-cell-gene-expression).
10x is the major vendor of scRNAseq library prep.

Issues to pay attention to

- When applied to tissue this technology critically relies on dissociation. 
    There will be physiological effects on the cells that may affect 
    gene expression, interact with lysis protocols, or affect availability 
    of mRNA for profiling. 
    Validation of the dissociation method's effects of gene expression patterns
    may or may not be available, and would be complex and biased to apply.
    Nevertheless, cancer cell lines will hopefully be less effected by these
    issues.
- It appears that most of their tech relies on 3' sequencing
    (3' end seq, so variation
    (what was that paper about stutter starting, umi like?
    (no 5' intronic variants
- It's possible to also use nuclei with ATACseq or the like. Perhaps pre-mRNA?

Dean Lee provides this link to 
[an online book about single-cell best practices](https://www.sc-best-practices.org/preamble.html),
as well as youtube presentations about processing scRNAseq 
[part1](https://www.youtube.com/watch?v=cmOlCTGX4Ik)
and
[part2](https://www.youtube.com/watch?v=FqG_O12oWR4).


### Available datasets

### Analyses

### Validation

## Questions to ask



<style>
body { max-width: 600px; margin: 1em;}
</style>
