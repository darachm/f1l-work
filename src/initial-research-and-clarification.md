---
author: Darach Miller
date: 2024-08-02
---

Dean Lee has asked each of us to address this Key Scientific Question (KSQ):

    The KSQ: Using available scRNA-seq data from cancer cell lines, how would
    you explore the use of the following FDA-approved antibody therapies in
    additional cancers?
    
The therapies of interest are trastuzumab and bevacizumab.

# Summary

Trastuzumab and bevacizumab are antibodies that bind HER2 and VEGF-A, 
respectively.
I'm going to look in cancer cell line scRNAseq datasets to see if particular
lines are expressing these, expressing them highly, if associated interacting
partners are expressed, and if downstream signalling pathways are showing
expression that might signify ongoing activation of the pathways.

# Details

Below I'm going to discuss:

- What are these treatments? 
    How do they work? 
    In what situations do we think they'll be effective?
- What kind of data can we expect from cancer cell lines, and in what ways might
    analysis of scRNAseq data be able to predict the efficacy of these drugs?
- A proposed analysis strategy to address the question
- Outstanding questions that may be useful later.

Some parts are summarized from sources, and some parts are direct quotes copied
from the indicated sources.

## The two treatments

### Trastuzumab and HER2

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
such as in the therapy called trastuzumab deruxtecan[^statsPearlsTrastuzumab].

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
likely depend on the fidelity of the measurements of the marker being targeted -
and potentially the heterogeneity of expression.

Finally, to keep in mind the challenge of tumors escaping the therapy,
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

### Bevacizumab and VEGF-A

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

### Generalizing these treatments

Both of these treatments bind the extracellular domain of a transmembrane 
receptor to presumably inhibit its interaction with other receptors.
An alternative mechanism proposed for each is to modulate the immune-system 
by binding the antibody to the target, and trastuzumab has been used as an
ADC.

So, to identify other cancers where these treatments may exert their mechanistic
effects, we should identify cancers that are likely to express the target
protein.
Noting the expression of the other receptors that these targets
interact with may also be useful, but the expression of the 
target (HER2 or VEGF-A) is a good place to start.

## Cancer cell lines as a model

Cancer cell lines are populations of human cells that can be cultured _ex vivo_
in labs (ie "in vitro"), where the cells are derived are thought to share 
some of the cell biology important for studying a particular type of human 
cancer [^mirabelli2019]. Usually they are isolated from a patient's tumor.
Thus, they are a model system that makes a wider range of experiment modalities
and designs than would be possible _in vivo_.

As we'd expect, this means that there are large amounts of data and studies
reported using these _ex vivo_ cell lines.

### Assays to profile these models

For example,
Dean Lee provided this link to [the Broad's CCLE](https://sites.broadinstitute.org/ccle/#:~:text=Cancer%20cell%20lines%20are%20the,and%20for%20defining%20drug%20efficacy.).
This initiative characterized a large number (~1000) of human cancer cell lines
to determine genetic variants, measure RNA and protein expression, and profile
sensitivity to pharmacological panels (amongst other assays[^depmap]).
In addition, other resources include Cell Model Passports[^cellModelPassports],
the TCGA [^tcgaPortal] [^li2017], and the CCLE's DepMap data repository[^depmap].

Thus, we can use these and other data to identify cancer cell lines that
express the treatment targets or other possible markers of efficacy. 
They would also be useful models for testing our predictions empirically,
were someone to be so inclined.

[^cellModelPassports]: https://cellmodelpassports.sanger.ac.uk/
[^tcgaPortal]: https://tcpaportal.org/mclp/#/
[^li2017]: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5501076/
[^mirabelli2019]: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6721418/
[^depmap]: https://depmap.org/portal/data_page/?tab=allData

### What does scRNA-seq data represent?

Single-cell RNAseq is a technology wherein microfluidics and molecular barcoding
are used to multiplex RNAseq library preparation from many individual cells.
Then, these RNAseq libraries can sequenced on the appropriate platform, commonly
Illumina short read. 

The most common vendor/platform for the library preparation is 10x, which
as far as I can tell is based on a poly-dT primed library reverse transcription
step for first strand. 
Second strand synthesis is primed from a template switching oligo.
These libraries can then be subject to sequencing-platform-specific strategies
of fragmentation or whole-template sequencing.
Commonly I see researchers using 3' data, but there is also a 5' option.
So we should be able to get mRNA abundance estimates for these cells, although
I believe it will depend on the specific protocol if we can pull out pre-mRNA
or different isoforms.

Cancer cell lines should be tractable for assays by this method, as they should
be straightforward to dissociate (if need be) into suspension to feed onto
the library preparation chip. It is also possible to use nuclei for DNA genome
assays or to profile pre-mRNA (presumably, using non-polyA selection methods),
so these datasets may also be available.

### What can scRNA-seq data tell us?

RNAseq will tell us what genes are transcribed in the target cell/tissue.
This is a good screening technology to indicate potential gene expression
of the target protein.
I assume there are also going to be transcriptional programs associated with
the activation of the HER2 and/or VEGF-A signalling, and thus if I can find
a clear annotation of this then I can use expression of these to determine if
the pathway appears to be active (and thus capable of being shutoff!).

The obvious distinction of scRNAseq is that it is single-cell.
While this is more noisy than bulk, this means we will be able to look for
heterogeneity of marker expression within populations.
Heterogeneity may complicate the application of targeted treatments like 
trastuzumab and bevacizumab by offering a "bet-hedging" phenomena of lineages 
escaping selection, so this may complicate the therapy.
But, any heterogeneity ay unveil more of the underlying
biology, if the model well-represents the disease it models, and 
could even conceivably also be an opportunity to target the
target-expressing population subset.

## A possible strategy - what is to be done with this analysis

My plan is to look in the proposed cancer cell line scRNAseq data to look
for expression of VEGF-A and/or HER2. 
For signs that they are active, I'll look for if there are downstream 
transcriptional programs that may indicate activation of these receptors.
Later, I will expand this to also look for expression of proteins that are
known to interact with HER2 or VEGF-A.

Another useful analysis would be to see if there's datasets where cell lines 
are treated with trastuzumab or bevacizumab (or similar antibodies), and if
these are indeed affecting any pathway activation or any 
innate-cellular-immunity or apoptotic markers (regarding the 
antibody-dependent cellular cytotoxicity hypothesis). 
Other data modalities could be useful, like proteomics or bulk RNAseq.

# Unanswered questions - Darach's to-do/to-ask list

- Is HER2 or VEGF-A expression post-transcriptionally regulated?
    How do people think about that in using RNAseq for screening, in humans?
- I saw that HER2 signals significantly into the PI3K-AKT pathway. It might be 
    good to collect some KEGG-style pathway annotations, to look for 
    pathway
- For HER2, what happens when you treat in a non-amplified tumor, ie with
    "normal" expression levels? What about just low expression? 
    Has the focus on HER2 overexpressing cancers been a byproduct of simply
    prioritizing investigation of cases where it's highly expressed?
    Could high-expressed tumors also have compensatory mutaitons to reduce
    the protein expression to buffer the expression level lower 
    (not to baseline, but below where the gene copy number would predict)?
- Look up more about how "antibody-dependent cellular cytotoxicity" is 
    triggered, molecular signatures of that, how to assay if that's a thing here
- HER2 overexpression is sometimes determined via DNA FISH, so is actually a
    measure of genomic copy. How tightly does that determine gene expression?
    Could silencing mechansisms complicate that linkage?
- What's a good sign for using one of these as a drug-conjugate? 
    Specificity, internalization, ... what else?
- I do not know how accessible various tissues are. 
    For example, I saw it mentioned that these may not cross the blood-brain 
    barrier, so of course that may limit efficacy and wouldn't be apparent in
    cell line studies.
- I need to watch Regev's talk again!

<style>
body { max-width: 600px; margin: 1em;}
</style>
