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

~~**p_fluor_wiesmann2022_acc.txt**: Run accession list for *P. fluorescens* from [Wiesmann 2022](https://www.nature.com/articles/s41396-022-01343-3) ([GSE190448](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE190448)).~~

~~**p_aeru_wiesmann2022_acc.txt**: Run accession list for *P. aeruginosa* from [Wiesmann 2022](https://www.nature.com/articles/s41396-022-01343-3) ([GSE190448](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE190448)).~~

---
**m_tuber_H37Rv_martini2023_acc.txt**: Run accession list for *M. tuberculosis* str H37Rv from [Martini 2023](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC10178195/) ([GSE218354](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE218354)).

**s_elon_lu2023_acc.txt**: Run accession list for *S. elongatus* from [Lu 2023](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC10365998/) ([GSE227397](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE227397)).

**c_vibrioides_mclaughlin2023_acc**: Run accession list for *Caulobacter vibrioides* from [McLaughlin 2023](https://pubmed.ncbi.nlm.nih.gov/37645952/) ([GSE241057](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE241057)).

~~**a_fischeri_ES114_thompson2016_acc**: Run accession list for *Aliivibrio fischeri ES114* from [Thompson 2016](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5409853/) ([GSE80607](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE80607)).~~ *These ended up being single ended reads.*

**a_fischeri_ES114_griend2023_acc**: Run accession list for *Aliivibrio fischeri ES114* from [Griend 2023](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC10441365/) ([GSE237189](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE237189)).

**b_theta_bolam2019_acc**: Run accession list for *Bacteroides thetaiotaomicron* from [Bolam 2019](https://www.nature.com/articles/s41564-019-0466-x) ([GSE129572](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE129572)).

**b_theta_lewis2023_acc**: Run accession list for *Bacteroides thetaiotaomicron* from [GSE169260](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE169260).

**s_spcc6803_wilde2022_wt_acc**: Run accession list for *Synechocystis sp. PCC 6803* from [GSE132275](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE132275). This list only captures wild-type runs.

**p_fluor_hasnain2023_acc**: Run accession list for *Pseudomonas fluorescens* from [Hasnain 2023](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC10229592/) ([GSE200822](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE200822)).

**s_coelicolor_wezel2023_wt_acc**: Run accesion list for *Streptomyces coelicolor* from [GSE234439](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE234439). This list only captures wild-type runs.

~~**a_rabiei_snrc2015_acc**: Run accession list for the fungus *Ascochyta rabiei* from [PRJNA288273](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA288273); [SRX1080581](https://www.ncbi.nlm.nih.gov/sra/SRX1080581[accn]). Of the 31 runs, only the ones in paired read layout are captured in the list. *Seems like the paired reads are split into individual SRA runs. The second in pair file was renamed according to first in pair's SRA ID*~~. *Did not work; not worth troubleshooting*.

**g_hirsutum_sun2023_acc**: Run accession list for C O T T O N *Gossypium hirsutum* from [GSE182982](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE182982).

**g_hirsutum_mei2022_acc**: Run accession list for C O T T O N *Gossypium hirsutum* from [GSE206663](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE206663). >36Gbyte total; do not use unless the one above is insufficient.

## Miscellaneous

Start an interactive shell on pitzer:

```bash
srun -n 1 --cpus-per-task=48 --account=PDS0325 --time=1:00:00 --pty /bin/bash
```
