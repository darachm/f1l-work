---
author: Darach Miller
date: 2024-07-31
---

Dean Lee has asked each of us to address this Key Scientific Question (KSQ):

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

- What are these treatments? How do they work? When do we think they'll 
    be effective?
- What kind of data can we expect from cancer cell lines, and in what ways might
    analysis of scRNAseq data be able to predict the efficacy of these drugs?
- A proposed analysis strategy to address the question

Some parts are summarized from sources, and some parts are direct quotes.
This document is not an academic contribution, it's a rote summary.

# The two treatments

## Trastuzumab and HER2

Trastuzumab is a monoclonal antibody that binds HER2 (human epidermal growth
factor receptor 2) extracellularly [^statsPearlsTrastuzumab].
This is thought to inhibit HER2 homodimerization, HER2 cytoplasmic
autophorsphorylation, and thus downstream HER2-mediated signals promoting 
tumorigenesis and cell proliferation [^iqbal2014]. 
HER2 also binds in heterodimers with EGFR, HER3, HER4 [^iqbal2014], and
signals downstream most to "particularly the PI3K/Akt" pathway.
"HER2 can also be activated by complexing with other membrane receptors such 
as insulin-like growth factor receptor 1"[^iqbal2014].

Trastuzumab is thought to also stimulate "the immune system to kill cells 
that have a lot of HER2"[^nciTrastuzumab], which is to say it has been shown to 
mediate antibody-dependent cellular cytotoxicity (ADCC) on 
"HER2 overexpressing cancer cells"[^dailyMedTras].
"Proposed mechanisms of trastuzumab actions include (1) inhibition of HER2 
shedding, (2) inhibition of PI3K-AKT pathway, (3) attenuation of cell
signalling, (4) antibody-dependent cellular cytotoxicity, and (5) inhibition of
tumor angiogenesis"[^iqbal2014]. 
It is also used to target antibody-drug conjugates,
such as in the treatment called trastuzumab deruxtecan[^statsPearlsTrastuzumab].

Trastuzumab is often used in combination with chemotherapies that inhibit cell
proliferation or damage DNA. Some use cases:

1. HER2 positive breast cancer that is high-risk or hormone receptor negative 
    (grows independent of hormones)
1. Metastatic HER2 positive breast cancer
1. Adenocarcinoma of the stomach or gastroesophageal junction that has
    metastitized

Regarding its use as a marker for generalization to other cancers,
"clinical trials using HER2 directed therapies in lung and bladder cancers 
reported disappointing clinical benefits."[^iqbal2014]
However, looking at the two cited studies [^langer2004][^hussain2007] show
trends that the study authors take to _suggest_ a relationship.
"HER2 positive" appears to be assessed by immunohistochemistry or
gene copy amplification using DNA FISH[^dailyMedTras][^iqbal2014].
DNA FISH is an indirect method, and the link between gene copy and expression
is not 1:1. For example, "[u]nlike breast cancer, overexpression without gene 
amplification is more common in [urothelial cancer]"[^hussain2007].
In addition,
"ASCO/CAP guidelines state that intratumoral heterogeneity may contribute to 
HER2 testing inaccuracy."[^iqbal2014]
Thus, the efficacy of applying a HER2 -targeting treatment to any cancer will
likely depend on the fidelity of the measurements of the marker being targeted.

Finally,
resistance to trastuzumab can be mediated by a HER2 mutation lacking the 
extracellular domain (p95)[^iqbal2014] and 
increased activity due to point mutations can cause 
constitutive activation[^brandt-rauf1990].

[^statsPearlsTrastuzumab]: https://www.ncbi.nlm.nih.gov/books/NBK532246/
[^nciTrastuzumab]: https://www.cancer.gov/about-cancer/treatment/drugs/trastuzumab
[^iqbal2014]: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4170925/
[^dailyMedTras]: https://dailymed.nlm.nih.gov/dailymed/drugInfo.cfm?setid=492dbdb2-077e-4064-bff3-372d6af0a7a2
[^langer2004]: https://sci-hub.st/10.1200/JCO.2004.04.105
[^hussain2007]: https://sci-hub.st/10.1200/JCO.2006.08.0994
[^brandt-rauf1990]: doi.org/10.1073/pnas.87.21.8660

## Bevacizumab and VEGF-A

Bevacizumab is a VEGF-A-targeting humanized monoclonal antibody and the first 
approved angiogenesis inhibitor[^garcia2020].
VEGF-A is vascular endothelial growth factor A, and bevacizumaub binds to all
known VEGF-A isoforms to block its interaction with the VEGF receptors (VEGFR)
VEGFR-1 (fit-1) and VEGFR-2 (KDRflk-1)[^statsPearlsBevacizumab].
Thus bevzcizumab is a "vascular endothelial growth factor 
inhibitor"[^dailyMedBeva] that accordingly reduces microvascularization, 
but may also "modulate tumor-induced immunosupression"[^garcia2020].

Bevzcizumab is often used in combination with chemotherapies, and is used on a
wide range of cancers:

1. Cervical cancer that has not responded to other treatment, or is spreading
1. Colorectal cancer that has spread
1. Recurrent Glioblastoma
1. Hepatocellular carcinoma that has spread
1. Nonsquamous non-small cell lung cancer that has spread or is recurrent
1. Ovarian epithelial, fallopian tube, or primary peritoneal cancers that have
    recurred, not responded, or failed earlier treatment
1. Renal cell carcinoma that has spread, used in combination with interferon alpha

[^garcia2020]: https://pubmed.ncbi.nlm.nih.gov/32335505/
[^dailyMedBeva]: https://dailymed.nlm.nih.gov/dailymed/drugInfo.cfm?setid=939b5d1f-9fb2-4499-80ef-0607aa6b114e
[^statsPearlsBevacizumab]: https://www.ncbi.nlm.nih.gov/books/NBK482126/

## Generalizing these treatments

Both of these treatments bind the extracellular domain of a transmembrane 
receptor to presumably inhibit its interaction with other receptors.
To identify other cancers where these treatments may exert their mechanistic
effects, we should identify cancers that are likely to express the targeted
marker. 
Additionally, looking for expression of the other receptors that these targets
interact with may also be useful, although other unknown pathways and effects
may be at play.
An additional proposed mechanism is that these antibodies trigger immunological
effects, and so the interaction with the known partners may not be relevant.

Third, these antibodies could be used within an antibody drug conjugate (ADC)
to target chemotherapies, and so simple surface expression of the target is
a likely prerequisite of it. 

# Cancer cell lines as a model

Cancer cell lines are populations of human cells that can be cultured _ex vivo_
in labs (ie "in vitro"), where the cells are derived are thought to share 
some of the cell biology important for studying a particular type of human 
cancer [^mirabelli2019]. Usually they are isolated from a patient's tumor.
Thus, they are a model system that makes a wider range of experiment modalities
and designs than would be possible _in vivo_.

## Assays to profile these models

For example,
Dean Lee provided this link to [the Broad's CCLE](https://sites.broadinstitute.org/ccle/#:~:text=Cancer%20cell%20lines%20are%20the,and%20for%20defining%20drug%20efficacy.).
This initiative characterized a large number (~1000) of human cancer cell lines
to determine genetic variants, measure RNA and protein expression, and profile
sensitivity to pharmacological panels (amongst other assays[^depmap]).
In addition, other resources include Cell Model Passports[^cellModelPassports],
the TCGA [^tcgaPortal][^li2017], and the CCLE's DepMap data repository[^depmap].

Thus, we can use these and other data to identify cancer cell lines that
express the treatment targets or other possible markers of efficacy. 
This analysis can then inform later experiments testing how the treatments
may work in _ex vivo_ cell culture, thus building evidence towards clinical 
application.

[^cellModelPassports]: https://cellmodelpassports.sanger.ac.uk/
[^tcgaPortal]: https://tcpaportal.org/mclp/#/
[^li2017]: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5501076/
[^mirabelli2019]: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6721418/
[^depmap]: https://depmap.org/portal/data_page/?tab=allData

## scRNA-seq as an assay

Single-cell RNAseq is a technology wherein microfluidics and molecular barcoding
are used to multiplex RNAseq library preparation from many individual cells.
Then, these RNAseq libraries can sequenced on the appropriate platform, commonly
Illumina short read. 
The most common vendor/platform for the library preparation is 10x, which
as far as I can tell is based on a poly-dT primed library reverse transcription
step for first strand. 
Second strand synthesis is primed from a template switching oligo strategy.
These libraries can then be subject to sequencing-platform-specific strategies
of fragmentation or whole-template sequencing.
Commonly I see researchers using 3' data, but there is also a 5' option.

Cancer cell lines should be tractable for assays by this method, as they should
be straightforward to dissociate (if need be) into suspension to feed onto
the library preparation chip. It is also possible to use nuclei for DNA genome
assays or to profile pre-mRNA.

The obvious advantage of scRNAseq is that it is single-cell.
Thus this assay may be able to uncover heterogeneity in marker expression
that may complicate the application of targeted treatments like 
trastuzumab and bevacizumab.
Detecting that a subset of cells express these markers may explain why bulk
assays may have not detected such expression, but treatments that deplete these
target-expression subpopulations may shift the population structure to disfavor
expression of the target.
Let's see.

<!--
Dean Lee provides this link to 
[an online book about single-cell best practices](https://www.sc-best-practices.org/preamble.html),
as well as youtube presentations about processing scRNAseq 
[part1](https://www.youtube.com/watch?v=cmOlCTGX4Ik)
and
[part2](https://www.youtube.com/watch?v=FqG_O12oWR4).
-->


# A possible strategy - what is to be done

These therapies target extracellular components of pro-growth signalling 
pathways, either the angiogenic VEGF-A or the epidermal HER2. 

Thus, it may be a good idea to just find cancers where tumor cells are
typically expressing these genes. 

scRNAseq would indicate which cell lines have the potential to express the gene 
from it's mRNA, and so may be a good start.


Dean Lee provides links to 
[an article about Aviv Regev's presentation at AACR 2024](https://www.genengnews.com/topics/cancer/aacr-2024-aviv-regev-shows-how-single-cell-atlases-foster-new-axis-to-genentechs-drug-discovery/)
and
[a recording of a talk by Regev about cell atlases as roadmaps for cancer treatments](https://www.youtube.com/watch?v=Wk5QHySlMXU).
To summarize, Regev (amongst many many others) has been pushing Cell Atlas 
projects for at least a decade. 
Now we can enjoy the _fruits_ of this labor. 




# Unanswered questions for me 

- Is HER2 or VEGF-A expression post-transcriptionally regulated?
- HER2 overexpression is sometimes determined via DNA FISH, so is actually a
    measure of genomic copy. How tightly does that determine gene expression?
    Could silencing mechansisms complicate that linkage?
- What's a good sign for an ADC? Specificity, internalization, ... what else?
- For HER2, what happens when you treat in a non-amplified tumor? What about
    just low expression? 
- PI3K-AKT

- I do not know how accessible various tissues are. 
    For example, I saw it mentioned that these may not cross the blood-brain 
    barrier, so of course that may limit efficacy and wouldn't be apparent in
    cell line studies.

<style>
body { max-width: 600px; margin: 1em;}
</style>
