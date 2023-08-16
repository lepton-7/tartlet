# Transcriptionally-active Riboswitch Tracer (TaRT)

Finds transcriptional riboswitch activity in RNA-seq data.

## Validation

TaRT is validated on published transcriptomic and metatranscriptomic datasets. Datasets are downloaded with  `validation/rna_seq/download_SRAs.sh` using compilation lists of SRA run IDs.

### RNA-Seq data description

**e_coli_sastry2019_acc.txt**: Run accession list for *E. coli* RNA-seq datasets from [Sastry 2019](https://www.nature.com/articles/s41467-019-13483-w). This list contains original runs from the paper (GSE122211, GSE122295, GSE122296, and GSE122320) as well as datasets from [Guzman2018](https://www.biorxiv.org/content/10.1101/310946v2.full) (GSE114358) and [Anand2019](https://www.nature.com/articles/s41564-018-0340-2) (GSE122779).

**b_sub_guo2023_acc.txt**: Run accession list for *B. subtilis* str. 168 from [Guo 2023](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC9879117/) ([GSE219221](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE219221)).

**b_sub_tilburg2022_acc.txt**: Run accession list for *B. subtilis* str. 168 from [Tilburg 2022](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC9049611/) ([GSE169409](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE169409)). *These reads may not be paired since the fasterq-dump did not produce paired files*.

---
*<u>Turns out HISAT-2 could not process these because they seem to be single-end reads packed as paired-end. Need to be thrown out</u>*

**p_fluor_wiesmann2022_acc.txt**: Run accession list for *P. fluorescens* from [Wiesmann 2022](https://www.nature.com/articles/s41396-022-01343-3) ([GSE190448](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE190448)).

**p_aeru_wiesmann2022_acc.txt**: Run accession list for *P. aeruginosa* from [Wiesmann 2022](https://www.nature.com/articles/s41396-022-01343-3) ([GSE190448](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE190448)).

---
**m_tuber_H37Rv_martini2023_acc.txt**: Run accession list for *M. tuberculosis* str H37Rv from [Martini 2023](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC10178195/) ([GSE218354](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE218354))

**s_elon_lu2023_acc.txt**: Run accession list for *S. elongatus* from [Lu 2023](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC10365998/) ([GSE227397](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE227397))
