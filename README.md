# CXCL-RCC
Analysis of the role of CXCL9/10/11 and CXCR3 in renal cancer 

## Terms of use

The pipeline results will be included in a future publication by Renate Pichler et al (submitted). The analysis concerns solely publicly available data cited below. To reference and use analysis results, please cite our GitHub repository; the data sources listed below and, if available, the publication. In any questions, pelase contact [Renate Pichler](mailto:renate.pichler@i-med.ac.at) or [Piotr Tymoszuk](mailto:piotr.s.tymoszuk@gmail.com).

## Data sources

This is a complete R analysis pipeline of 8 publicly available renal whole-transcriptome datasets:

* __TCGA KIRC project__: clinical information and normalized RNAseq data can be obtained from the [GDC Data Portal of the National Cancer Insitute](https://portal.gdc.cancer.gov/). the data set was first analyzed by Creighton and colleagues [^1]. Data extraction was done with the TCGA-Assembler-2 script available from [@compgenome365/TCGA-Assembler-2](https://github.com/compgenome365/TCGA-Assembler-2/tree/master/TCGA-Assembler)

* __CheckMate 010, CheckMate 025 everolimus, CheckMate 025 nivolumab datasets__: clinical information and normalized RNAseq data available as Supplementary Tables S1 and S4 accompanying the publication of Braun et al. [^2]

* __RECA-EU project__: clinical information and normalized RNAseq data can be downloaded from the [ICGC Data Portal](https://dcc.icgc.org/projects/RECA-EU)

* __GSE73731 and GSE167093 datasets__: normalized microarray expression data and clinical information [^3][^4] can be accessed from [Gene Expression Ominibus](https://www.ncbi.nlm.nih.gov/geo/), e.g. with the _GEOquery_ R package [^5]

* __E-MTAB 1980 dataset__: normalized microarray expression data and clinical information [^6] can be obtained from [ArrayExpressData](https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-MTAB-1980)

## Usage

To make sure to install required development packages prior to runung the pipeline:

```r

devtools::install_github('PiotrTymoszuk/ExDA')
devtools::install_github('PiotrTymoszuk/microViz')
devtools::install_github('PiotrTymoszuk/kmOptimizer')
devtools::install_github('PiotrTymoszuk/coxExtensions')
devtools::install_github('PiotrTymoszuk/trafo')
devtools::install_github('PiotrTymoszuk/gseaTools')
devtools::install_github('PiotrTymoszuk/biggrExtra')
devtools::install_github('PiotrTymoszuk/clustTools')
devtools::install_github('PiotrTymoszuk/soucer')
devtools::install_github('PiotrTymoszuk/somKernels')
devtools::install_github('PiotrTymoszuk/figur')

```
To launch the entire pipeline, source the `exec.R` file:

```r

source('exec.R')

```

[^1]: Creighton CJ, Morgan M, Gunaratne PH, Wheeler DA, Gibbs RA, Robertson G, Chu A, Beroukhim R, Cibulskis K, Signoretti S, et al. Comprehensive molecular characterization of clear cell renal cell carcinoma. Nature 2013 499:7456 (2013) 499:43–49. doi: 10.1038/nature12222

[^2]: Braun DA, Hou Y, Bakouny Z, Ficial M, Sant’ Angelo M, Forman J, Ross-Macdonald P, Berger AC, Jegede OA, Elagina L, et al. Interplay of somatic alterations and immune infiltration modulates response to PD-1 blockade in advanced clear cell renal cell carcinoma. Nature medicine (2020) 26:909. doi: 10.1038/S41591-020-0839-Y

[^3]: Wei X, Choudhury Y, Lim WK, Anema J, Kahnoski RJ, Lane B, Ludlow J, Takahashi M, Kanayama HO, Belldegrun A, et al. Recognizing the Continuous Nature of Expression Heterogeneity and Clinical Outcomes in Clear Cell Renal Cell Carcinoma. Scientific reports (2017) 7: doi: 10.1038/S41598-017-07191-Y

[^4]: Laskar RS, Li P, Ecsedi S, Abedi-Ardekani B, Durand G, Robinot N, Hubert JN, Janout V, Zaridze D, Mukeria A, et al. Sexual dimorphism in cancer: insights from transcriptional signatures in kidney tissue and renal cell carcinoma. Human molecular genetics (2021) 30:343–355. doi: 10.1093/HMG/DDAB031

[^5]: Sean D, Meltzer PS. GEOquery: a bridge between the Gene Expression Omnibus (GEO) and BioConductor. Bioinformatics (Oxford, England) (2007) 23:1846–1847. doi: 10.1093/BIOINFORMATICS/BTM254

[^6]: Sato Y, Yoshizato T, Shiraishi Y, Maekawa S, Okuno Y, Kamura T, Shimamura T, Sato-Otsubo A, Nagae G, Suzuki H, et al. Integrated molecular analysis of clear-cell renal cell carcinoma. Nature Genetics 2013 45:8 (2013) 45:860–867. doi: 10.1038/ng.2699
