# GP2 Local Ancestry

*This repository contains the code and overall pipeline that was designed in order to implement local ancestry inference in GP2's release 10 data and to run a local-ancestry-aware genome-wide association study (plus downstream analysis).*

---
## Reference panel assembly 

---

## QC of non-European populations altogether
We must first know which samples are labeled to belong to the EUR label assigned by genotools.

```
##get the list of EUR samples 
awk 'NR>1 {print $2}' /config/workspace/gp2_tier2_eu_release10/raw_genotypes/EUR/EUR_release10_vwb.psam > fullEUR_POPlist.txt


##make sure the LA ref samples are in the all EUR list 
grep -Ff /config/workspace/ws_files/LA_noAJ/Best1/reference4ProjectedPCA/runProjectedPCA/plotProjectedPCA/parentals/EUR.txt \
fullEUR_POPlist.txt | wc -l 
```
And then get the list of samples from EUR that do not belong to the reference panel. Because the one that belong to the ref panel must go inside the QC with all the other non-EUR samples. Because ultimately they need to go into topmed for phasing. Since the reference panel for local ancestry inference, and the target data must be in the same phase. 

```
##get the list of EUR samples WITHOUT the members of the LA reference panel 
awk 'NR==FNR {remove[$1]; next} !($1 in remove)' \
/config/workspace/ws_files/LA_noAJ/Best1/reference4ProjectedPCA/runProjectedPCA/plotProjectedPCA/parentals/EUR.txt \
fullEUR_POPlist.txt \
> fullEUR_POPlist_noLAref.txt
```
Now, the exclusion of the EUR samples (that are also not part of the reference sample), from the merged dataset
```
plink2   --exclude ultimate_list_snps2remove_nba_qc.txt \
  --make-pgen \
  --memory 50000 \
  --out GP2_r10_final_with_gp2id_ready_genotools_noEUR_exceptLAref \
  --pfile /config/workspace/gp2_makariousmb/R10/GP2_r10_final_with_gp2id_ready_genotools \
  --remove fullEUR_POPlist_noLAref.txt \
  --threads 8
```
The plink2 files named GP2_r10_final_with_gp2id_ready_genotools* are the ones that normally would be QCed with genotools (to partitioning them to the 11 genotools's populations). However, for this particular analysis, the working dataset must have enough representation of the ancestries that are included in the association model. And if we use the 11 populations derived from genotools, we cannot achieve enough representation of other ancestries within each cluster (not even at continental level). 

That is why we must perform the QC with the samples altogether, and thereafter --> imputation with topmed --> local ancestry inference --> GWAS --> post GWAS.

Since we are excluding the EUR samples (as already classified by genotools for GP2 default QC pipeline), the analysis would be, in principle, more computationally affordable. 

For the QC we are going to use a pipeline designed for admixed populations

---
