# Barley Phenology 2023

Brief Introduction


# Blast result collecting and tpm data extraction  
## 1. Blast GEMOMA result with GsRTD protein database  

Protein sequences from GEMOMA result were BLAST against GsRTD protein database with corresponding variety using NCBI blastp

```
# Example Morex
blastp -db Morex.GsRTD.db -query Morex.gene.1.fa -out Morex.gene.1.txt -outfmt 6
```

## 2. Merge blast result

```
merge_blast_result.py
```

Reads the BLAST result from files and write them to a merge file  

## 3. Match TPM by GeneID  

```
match_tpm.py
```

Reads cleaned merged BLAST results from a table and match them with the TPM by geneIDs.  
