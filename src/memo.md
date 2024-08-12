---
author: Darach Miller
date: 2024-08-08
---

# This is "The Memo"

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


## The Paper

Dean Lee offered us the paper Kinker et al 2020 
([doi 10.1038/s41588-020-00726-6](doi.org/10.1038/s41588-020-00726-6))
as a dataset to use. 

The main focus of the work is to evaluate cancer cell lines using scRNAseq
to get a better look at them as models of tumors. A key aspect of this is if
the heterogeneity of a tumor in tissue can be reflected in a theoretically
more homogenous lab culturing environment.

Technical note - they multiplexed 24-27 cell lines (from CCLE) where each pool
had lines with "comparable proliferation rates", then used 10x Chromium for 
scRNAseq on eight pools of these cells. Then they used clustered to deconvolve
it, based on clustering genetic and expression profiles. Both of these were 
based on comparison to bulk RNAseq, and the two approaches "were consistent
for 98% of cells". They excluded inconsistent assigned cells, and lines with
fewer than 50 cells, "low-quality data" cells, and suspected doublets.
They got 53,513 cells from 198 cell lines. Average of 19,264 UMIs and 3,802
genes per cell.

They do point out that they co-cultured for 3 days before collection, so
that could affect expression, but they point to a control of six pure cell
lines as co-culturing having a "modest effect". 
I think they should have just mixed the cells immediately before prepping 
libraries. If it's possible to methanol/ethanol-fix cells right before prep
then I'd prefer if they did that to minimize co-culturing bias! I don't see
why they would have to grow them multiplexed, except to save on 
flasks/media/effort.

They used tSNE then DBSCAN to cluster cells. So I've heard that kind of 
approach has a lot of caveats, and I need to read up on that more later.
They then did "non-negative matrix factorization" to "identify both
continuous and discrete variability". I don't know how that works, but they
"detected 1,445 robust expression programs across all cell lines",
with 4-9 in individual cell lines. They then filter these and say that
there's about 800 programs that are not similar to other programs or technical
cofounders, and only 4.75% correspond to discrete subpopulations 
(DBSCAN clusters I presume).

They define term RHP for recurrent heterogeneous programs, recurrent across 
multiple cell lines. Most prominent two programs were associated with cell 
cycle (by inspection) for G1/S and G2/M phases, but G1/S was more "context 
dependent". They also point towards some signs of ex-vivo culture adaptations
in rapid growth and loss of G1 checkpoint.
Other RHPs are identified, they picked ten there that were detected across
at least eight cell lines from at least four different pools. 
Seven of these ten were "highly similar to the _in vivo_ RHPs" (they reanalyzed
some _in vivo_ datasets), shown by overlap of signature genes, and 
"high correlation of cell scores". 
They look at some of them for function, some are DNA damage, interferon, two
(protein folding/maturation and proteosomal) are not mapping to _in vivo_ 
heterogeneity.

There's similarity in three "distinct" RHPs to the epithelial-mesenchymal
transition (EMT). They also dig into "RHPs 6 and 7".
They annotate RHP 6 as being "a 'classical' p53-dependent senescence program"
and "RHP 7" as some other senescence similar to keratinocytes, lung bronchial
cells, and other epithelial cells - so they name it the "EpiSen" program.
They explore a bit of that by using some cell lines that are high/low in
that program. 

They then looked for factors associated with heterogeneity. 
They did copy-number analysis to identify subclones,
and then found 39% of expression-based clusters were associated with 
heterogeneity. 

They turned to some perturbations, which was nice, so transformed 
TGF-$\beta1$ and TGF-$\beta3$ and saw changes in some of the RHP programs
"EMT-II" and EpiSen. Cool. 

Then they screened drugs for their ability to affect EpiSen-high and EpiSen-low
subpopulations, sorted out from two selected model cell lines. They screened
2,198 compounds for affecting cell viability, then repeated that on 248 hits,
and found 113 with differential killing of the subpopulations of at least one
of the two cell lines.
>40% of hits were shared by both cell lines, and that was 71.4% when 
"considering the targets of the compounds rather than the exact compounds".
They validated dose-response for fourteen compounds.

Relevant to our original question for this project, they found that EpiSen-high
cells were sensitive to the "senolytic" (new word for me, senescence lytic)
ABT-737 (that's the Bcl-2 inhibitor). They also found them sensitive to
"multiple inhibitors of epidermal growth factor receptor (EGFR), AKT, 
phosphatidylinositol-3-OH kinase (PI(3)K), DNA-dependent protein kinase (DNA-PK),
insulin-like growth factor 1 receptor (IGF1R) and Janus kinase
(JAK)". "Accordingly, EpiSen-high cells (which are defined by
low AXL expression) were more sensitive to inhibitors of PI(3)K
and AKT, as well as those of EGFR and IGF1R that signal via the
PI(3)K–AKT axis."
So linking back to [trastuzumab](#trastuzumab-and-her2), that one's supposed to
inhibit the PI3K-AKT pathway (to which HER2 is thought to signal) and HER2
is also thought to be activated by complexing with other receptors including
that IGFR1 [^iqbal2014]. Also EGFR.
Thus, EpiSen-high cells might respond well to trastuzumab because that's another
way of hitting those pathways.

Then they looked at patients, and EpiSen is predictive of survival I believe.

In conclusion,
the found some recurrent heterogenous programs (RHPs).
Genetic subclones identified by inferred copy number changes don't explain
much of those, so non-genetic variability is at play here too 
(ie epigenetics).
Critically, they're looking at these "programs" not as being clean and 
distinct developmental programs but as being partial and limited trends.
They "hypothesize that cancer cells often activate partial or distorted 
programs, possibly not through canonical developmental mechanisms, in a
context-dependent manner", but I would say that it's rather that 
developmental programs are clean, distinct, and canonical by human 
interpretation and it's likely that even development has this quantitative
squishiness that is on display here. It's true here, and probably true in 
"normal" development as well.


Questions I have:
- What is doublet rate?
- What is cell count?
- What do the pure cultured controls look like?
- Does any of these RHPs correlate with proliferation/growth rate?

Technical details from the methods:
scRNAseq libraries were generated using 10x Chromium 3' kit v2
(or v3.1 for the co-culture experiment). 
"Cell barcode filtering, alignment of reads and UMI
counting were performed using Cell Ranger 3.0.1 (10x Genomics)."
They only considered cells with 2-9k genes detected, 
"yielding an average of 6,758 genes analyzed per cell line."
They do some log2 ( 1 + CPM/10 ) for expression, then 
when "analyzing cell lines individually, we only considered genes expressed 
at high or intermediate levels (E i,j > 3.5) in at least 2% of cells".
Huh.
They later say values were centered by cell line for relative expression,
so that assumes no total-count scaling of the transcriptome.
"When analyzing cell lines collectively, we selected the 7,000 most highly 
expressed genes across all cell lines".

### Questions from Dean Lee

    How did the authors handle the potential caveat of co-culturing cell lines 
    before profiling by scRNA-seq? Why do you think that caveat was or was not 
    adequately addressed?  

They made a analysis argument and did a control experiment.
The analysis argument was strangely phrased in the paper,
but they looked at pairwise correlations between NMF programs and used an
ANOVA to see if correlation looks similar within the pool compared to between
pools. This is a good way to detect uncommon and distinct effects, where some
cell line would have a consistent pattern on multiple lines within one pool
that you only see in that pool, but if there's a common mechanisms of crosstalk
then you'd see it in multiple pools so this doesn't rule out the existence of
common effects.
They also break apart the analysis by cancer type modeled by the cell line,
and it looks like there's greater effects for skin cancer, ovarian, types?
Then the ANOVA test, they just put the proportion of variance, but they don't
assess the significance of this. 
I've forgotten any context of how small this proportion of variances listed
here are, so I don't really know how to interpret that, but it looks to me like
some selective presentation.

So then they do a control experiment where they picked six cell lines 
(varied, but no discussion of why these six were picked) and 
did co-culturing or individual culturing.
They don't discuss which are pooled, so presumably they are all in one pool,
which has the flaw that now you can't use the same analysis as before to detect
pool-specific effects, so that's not the right way to detect this effect 
in my opinion.
The different expression is then detected by t-test, and while I'm new to
single-cell I do believe that you ought to be using some sort of 
variance-moderating approach like limma-voom or DESeq2 before you go doing 
these direct comparisons. 
Anyways, they go for it and then just use those differentially expressed genes
to look for hypergenometic test enriched GO terms. 
They find enrichment of processes, and state that it looks like the response
to co-culturing is strain specific.
Thus, this means that the analysis they do before on the actual data shouldn't
really pick up the kinds of effects they see here - the previous analysis would
require seeing an effect shared by members of a pool, not a heterogenous
response that is instead shared between similar lines in different pools.

I don't mean to harsh on the authors, but I just don't buy what they're selling
here.
I don't think it scuttles the paper, but it would have been better to not
pool cell lines in co-culture!
For example, they show how if you sort out the high/low EpiSen program, 
the culture 
re-establishes itself to a similar program. I don't think this culturing 
involves a lot of death and clearance of members of the population, so either 
the EpiSen program is just reverting to a mean level of
activation (like with cell cycle) or there could inter-subpopulation 
communication.
I don't think it was adequately addressed, and I don't think it was technically 
necessary, so I think it's a major shortcoming of the paper.

    The authors identified discrete subpopulations of cells within a subset of 
    individual cell lines (Fig. 2A-B). What might be the reason why some 
    cell lines have these discrete subpopulations while others do not?

Technical reasons such as depth of sampling could confound this.
Cells may be switching states, and this may reflect some developmental
instability - thus perhaps that part of their conclusion may depend on which
cell lines they focused on?
They may reflect some differentiated (probably not terminally!) populations.

    What are Recurrent Heterogeneous Programs (RHPs) and how were they defined?

These are patterns of heterogeneity that are shared across multiple cell lines.
They did their NMF to identify "programs", then excluded "those associated
with technical confounders" and "those with limited similarity to all other
programs". They were looking for programs that seem similar across cell lines.
Then only a few dozen programs remained, and only a few "corresponded to the 
discrete subpopulations described above".

    How do the identified RHPs relate to in vivo programs of heterogeneity in 
    tumors, and what evidence supports this relationship? 

They re-analyze some _in vivo_ data and see a lot of similarity between these
_ex vivo_ cell lines and _in vivo_ tumors.
A lot of the detected "programs" are similar, but there are still differences.

    Where can you download the scRNA-seq data as shown in Figure 1B?

The authors point to two data repositories. They make the most useful
data form available via Broad Institute's Single Cell Portal that requires
a license agreement, but only make the raw sequencing reads available in GEO.
So it depends on how much analysis you want to do, the easiest is the SCP - 
perhaps by design?

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
