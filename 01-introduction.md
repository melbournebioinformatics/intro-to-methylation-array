---
title: "Introduction: loading packages and downloading data"
teaching: 15
exercises: 5
---

:::::::::::::::::::::::::::::::::::::: questions

- How is whole genome methylation interrogated?
- How are methylation values commonly reported?
- What preliminary R packages and data is necessary to begin methylation analysis?

::::::::::::::::::::::::::::::::::::::::::::::::

::::::::::::::::::::::::::::::::::::: objectives

- Setting up the R environment and required  packages
- Create a sample table and define the variable of interest
- Download the six publicly available sample data from GEO
- Import the IDAT files into a `minfi` object

::::::::::::::::::::::::::::::::::::::::::::::::
## Introduction

Methylation in the human genome is known to be associated with development and disease. The Illumina Infinium methylation arrays are by far the most common way to interrogate methylation across the human genome. This workshop provides an easy and reproducible workflow using multiple packages for the analysis of methylation array data. Specifically, we demonstrate the steps involved in a typical differential methylation analysis pipeline including: quality control, filtering, normalization, data exploration and statistical testing for probe-wise differential methylation.

DNA methylation, the addition of a methyl group to a CG dinucleotide of the DNA, is the most extensively studied epigenetic mark due to its role in both development and disease. Although DNA methylation can be measured in several ways, the epigenetics community has enthusiastically embraced the Illumina HumanMethylation450 (450k) array as a cost-effective way to assay methylation across the human genome. Illumina has increased the genomic coverage of the platform with their MethylationEPIC (850k) array, and more recently MethylationEPICv2 (935k) covering to >935,000 sites. As methylation arrays are likely to remain popular for measuring methylation for the foreseeable future, it is necessary to provide robust workflows for methylation array analysis.

For each CpG, there are two measurements: a methylated intensity (denoted by 𝑀) and an unmethylated intensity (denoted by 𝑈). These intensity values can be used to determine the proportion of methylation at each CpG locus. Methylation levels are commonly reported as either beta values (𝛽=𝑀/(𝑀+𝑈)) or M-values (𝑀𝑣𝑎𝑙𝑢𝑒=𝑙𝑜𝑔2(𝑀/𝑈)).

Beta values (proportion of methylation at a particular CpG site) are generally preferable for describing the level of methylation at a locus or for graphical presentation because percentage methylation is easily interpretable. M-values (log2 ratio of methylated to unmethylated signal intensity) are more appropriate for statistical testing.
## Load packages

Complete the [setup instructions](../learners/setup.md) first. Open RStudio and
create a project using **File > New Project > New Directory > New Project**,
or open your existing workshop project. Create an R script, and run the chunks below in order.

Installing a package downloads it to your computer. Loading it with `library()`
makes its functions available in the current session. Load the packages again
whenever you restart R.


``` r
# Data handling
library(dplyr)

# Methylation analysis and annotation
library(minfi)
library(GEOquery)
library(limma)
library(IlluminaHumanMethylationEPICv2anno.20a1.hg38)
library(IlluminaHumanMethylationEPICv2manifest)
library(AnnotationHub)
library(DMRcate)
library(missMethyl)

# Plotting
library(ComplexHeatmap)
```
## Load annotation data

This is information that describes each CpG probe on the Illumina EPICv2 methylation array. It tells you where each probe is located in the genome and which genes or genomic features it is associated with.


``` r
epic <- getAnnotation(IlluminaHumanMethylationEPICv2anno.20a1.hg38)
head(epic)
```

``` output
DataFrame with 6 rows and 42 columns
                        chr       pos      strand            Name    AddressA
                <character> <integer> <character>     <character> <character>
cg25324105_BC11       chr19  37692358           + cg25324105_BC11     1754126
cg25383568_TC11       chr19  38727081           - cg25383568_TC11    79792482
cg25455143_BC11       chr19   1591515           - cg25455143_BC11    80699190
cg25459778_BC11       chr16   1614581           + cg25459778_BC11    60797262
cg25487775_BC11        chr2 161237458           + cg25487775_BC11     5799427
cg25595446_BC11       chr19  49677322           + cg25595446_BC11    65640459
                   AddressB              ProbeSeqA              ProbeSeqB
                <character>            <character>            <character>
cg25324105_BC11    99753217 ATTTATAAACTAATAACCCA.. GTTTATAAACTAATAACCCG..
cg25383568_TC11    69667133 AAACCAAAAAAATAACAAAC.. AAACCGAAAAAATAACAAAC..
cg25455143_BC11     7659147 ATAAAAAAAAATATACAACT.. ATAAAAAAAAATATACGACT..
cg25459778_BC11    65710482 AAAAATTTAAAACAAACAAC.. AAAAATTTAAAACAAACAAC..
cg25487775_BC11    89606481 AAAAACAACCTAAAAAACAA.. AAAAACAACCTAAAAAACAA..
cg25595446_BC11    39619855 AATAAAAATAACAACAACCA.. AATAAAAATAACGACAACCG..
                       Type    NextBase       Color    Probe_rs Probe_maf
                <character> <character> <character> <character> <numeric>
cg25324105_BC11           I           A         Red          NA        NA
cg25383568_TC11           I           C         Grn          NA        NA
cg25455143_BC11           I           T         Red          NA        NA
cg25459778_BC11           I           C         Grn          NA        NA
cg25487775_BC11           I           A         Red          NA        NA
cg25595446_BC11           I           C         Grn          NA        NA
                     CpG_rs   CpG_maf      SBE_rs   SBE_maf
                <character> <numeric> <character> <numeric>
cg25324105_BC11          NA        NA          NA        NA
cg25383568_TC11          NA        NA          NA        NA
cg25455143_BC11          NA        NA          NA        NA
cg25459778_BC11          NA        NA          NA        NA
cg25487775_BC11          NA        NA          NA        NA
cg25595446_BC11          NA        NA          NA        NA
                          Islands_Name Relation_to_Island         col
                           <character>        <character> <character>
cg25324105_BC11 chr19:37691892-37692..             Island           R
cg25383568_TC11 chr19:38726890-38727..             Island           G
cg25455143_BC11 chr19:1591428-159162..             Island           R
cg25459778_BC11  chr16:1610053-1615094             Island           G
cg25487775_BC11 chr2:161237949-16123..              Shore           R
cg25595446_BC11 chr19:49676753-49677..             Island           G
                 Probe_Type   Strand_FR   Strand_TB   Strand_CO Infinium_Design
                <character> <character> <character> <character>     <character>
cg25324105_BC11          cg           F           B           C               1
cg25383568_TC11          cg           R           T           C               1
cg25455143_BC11          cg           R           B           C               1
cg25459778_BC11          cg           F           B           C               1
cg25487775_BC11          cg           F           B           C               1
cg25595446_BC11          cg           F           B           C               1
                    Species     Rep_Num     UCSC_RefGene_Group
                <character> <character>            <character>
cg25324105_BC11       Human           1 TSS200;TSS200;TSS200..
cg25383568_TC11       Human           1        exon_18;exon_18
cg25455143_BC11       Human           1                       
cg25459778_BC11       Human           1                       
cg25487775_BC11       Human           1                       
cg25595446_BC11       Human           1                       
                     UCSC_RefGene_Name UCSC_RefGene_Accession
                           <character>            <character>
cg25324105_BC11 ZNF781;ZNF781;ZNF781.. NR_173332.1;NR_17333..
cg25383568_TC11            ACTN4;ACTN4 NM_001322033.2_2;NM_..
cg25455143_BC11                                              
cg25459778_BC11                                              
cg25487775_BC11                                              
cg25595446_BC11                                              
                      GencodeV41_Group        GencodeV41_Name
                           <character>            <character>
cg25324105_BC11 TSS200;TSS200;TSS200.. ENSG00000120784.17;E..
cg25383568_TC11 exon_18;exon_18;exon.. ACTN4;ACTN4;ACTN4;AC..
cg25455143_BC11                                              
cg25459778_BC11                 TSS200     ENSG00000007545.16
cg25487775_BC11                                              
cg25595446_BC11  exon_1;TSS200;TSS1500 PRMT1;ENSG0000012645..
                  GencodeV41_Accession Phantom5_Enhancers  HMM_Island
                           <character>        <character> <character>
cg25324105_BC11 ENST00000587199.5;EN..                               
cg25383568_TC11 ENST00000252699.7;EN..                               
cg25455143_BC11                                                      
cg25459778_BC11      ENST00000293925.9                               
cg25487775_BC11                                                      
cg25595446_BC11 ENST00000524771.5;EN..                               
                Regulatory_Feature_Group Methyl450_Enhancer         DMR
                             <character>        <character> <character>
cg25324105_BC11   Unclassified_Cell_ty..               TRUE         DMR
cg25383568_TC11   Gene_Associated_Cell..               TRUE            
cg25455143_BC11                                       FALSE            
cg25459778_BC11                                       FALSE            
cg25487775_BC11                                       FALSE            
cg25595446_BC11      Promoter_Associated               TRUE            
                Methyl450_Loci Methyl27_Loci EPICv1_Loci Manifest_probe_match
                   <character>   <character> <character>          <character>
cg25324105_BC11     cg25324105                cg25324105                 TRUE
cg25383568_TC11     cg25383568                cg25383568                 TRUE
cg25455143_BC11     cg25455143                cg25455143                 TRUE
cg25459778_BC11     cg25459778    cg25459778  cg25459778                 TRUE
cg25487775_BC11     cg25487775                cg25487775                 TRUE
cg25595446_BC11     cg25595446    cg25595446  cg25595446                 TRUE
```
## Download data

We use six EPICv2 arrays from [GEO GSE240469](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE240469):
three technical replicates each of LNCaP and PrEC, all at 500 ng DNA input.

- **LNCaP** is a prostate cancer cell line: cells grown in culture that model
  aspects of prostate cancer.
- **PrEC** is a primary culture of normal prostate epithelial cells. It provides
  a non-cancer comparison, but is not matched normal tissue from the same patient.

These data come from Peters *et al.* (2024),
[Characterisation and reproducibility of the HumanMethylationEPIC v2.0 BeadChip for DNA methylation profiling](https://link.springer.com/article/10.1186/s12864-024-10027-5).
The study evaluated the array's performance, including reproducibility and the
impact of DNA input, across cell cultures and tumour-derived samples. Our
workshop uses only the six 500 ng prostate culture arrays. 

Because these are technical replicates of the same culture, we expect them to have similar methylation profiles and cluster together. Whilst methylation studies can be performed on cultured cells, it is more traditionally used for primary human tissue (e.g. tumour tissue). 

### Create the sample table

Each row describes one sample. `Accession` identifies the GEO sample, `Name`
matches its filename label, and `Type` identifies the cell line to compare.


``` r
# Create the sample table for IDAT downloads. Replace this table with your own sample information when using your own data.

targets <- data.frame(
  Accession = c(
    "GSM7698438",
    "GSM7698446",
    "GSM7698462",
    "GSM7698435",
    "GSM7698443",
    "GSM7698459"
  ),
  Name = c(
    "LNCAP_500_1",
    "LNCAP_500_2",
    "LNCAP_500_3",
    "PREC_500_1",
    "PREC_500_2",
    "PREC_500_3"
  ),
  Type = c("LNCaP", "LNCaP", "LNCaP", "PREC", "PREC", "PREC")
)

targets$Name2 <- paste(
  targets$Accession, targets$Name, sep = "_"
)

# Define the variable of interest and convert it to a factor.
variable.of.choice <- as.factor(targets$Type)

targets
```

``` output
   Accession        Name  Type                  Name2
1 GSM7698438 LNCAP_500_1 LNCaP GSM7698438_LNCAP_500_1
2 GSM7698446 LNCAP_500_2 LNCaP GSM7698446_LNCAP_500_2
3 GSM7698462 LNCAP_500_3 LNCaP GSM7698462_LNCAP_500_3
4 GSM7698435  PREC_500_1  PREC  GSM7698435_PREC_500_1
5 GSM7698443  PREC_500_2  PREC  GSM7698443_PREC_500_2
6 GSM7698459  PREC_500_3  PREC  GSM7698459_PREC_500_3
```

### Download the IDAT files

The following chunk downloads the IDAT files for the six accessions in
`targets`. Each sample is saved in its own accession-named folder within
your current project directory.


``` r
# Download only IDAT files for the selected tutorial samples.
gsmlist <- lapply(
  targets$Accession,
  GEOquery::getGEOSuppFiles,
  makeDirectory = TRUE,
  baseDir = getwd(),
  filter_regex = "[.]idat([.]gz)?$"
)
names(gsmlist) <- targets$Accession
```
## Read in .idats

Add the sample names and file paths to `targets`, preserving its sample order.
`Basename` is the path to each file pair.


``` r
# Link the sample information to the downloaded IDAT files.
targets$Sample <- targets$Name
targets$Basename <- file.path(
  getwd(), targets$Accession, targets$Name2
)

# Read in raw data from the red/green IDAT files.
rgSet <- read.metharray.exp(targets = targets)
sampleNames(rgSet) <- targets$Sample
rgSet@annotation <- c(
  array = "IlluminaHumanMethylationEPICv2",
  annotation = "20a1.hg38"
)
rgSet
```

``` output
class: RGChannelSet 
dim: 1105209 6 
metadata(0):
assays(2): Green Red
rownames(1105209): 1600157 1600179 ... 99810982 99810990
rowData names(0):
colnames(6): LNCAP_500_1 LNCAP_500_2 ... PREC_500_2 PREC_500_3
colData names(7): Accession Name ... Basename filenames
Annotation
  array: IlluminaHumanMethylationEPICv2
  annotation: 20a1.hg38
```

`rgSet` contains the raw intensities for the six samples. For your own data,
set `targets$Basename` to your file paths and use the annotation matching your
array type. Quality control, filtering and normalisation follow in the next
episode.
## Save the input objects for reuse

`targets` contains sample information; `rgSet` contains the raw measurements.


``` r
dir.create("data", showWarnings = FALSE)
saveRDS(targets, "data/targets.rds")
saveRDS(rgSet, "data/rgSet.rds")
```

::::::::::::::::::::::::::::::::::::: keypoints

- Raw Illumina methylation array data are stored in `.idat` files, with one red-channel and one green-channel file per sample.
- The `targets` data frame holds the sample metadata, linking each sample's name and cell-line group to its IDAT file paths.
- `read.metharray.exp(targets = targets)` reads the raw intensities into `rgSet`, an `RGChannelSet` object, ready for quality control and preprocessing.

::::::::::::::::::::::::::::::::::::::::::::::::
