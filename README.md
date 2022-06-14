# RNA-seq processing for Chris

## Procedure (on TSCC)
1. Transfer IGM FTP directory to sequencingDataStorage on scratch

2. Run Salmon indexing for the Ensembl hg38 transcriptome

3. Run Salmon quantification for each sample in paired-end read mode

    - Files containing '100_u_IFNa' were replaced with 'induced'
    - Files containing 'Untreated' were replaced with 'uninduced'

4. Transfer the Salmon quantified reads to local

## Procedure (on local)
5. Run DESeq2 tx2gene with the Ensembl hg38 + alt contigs GFF

6. 