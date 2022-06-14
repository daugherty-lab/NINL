# RNA-seq processing for Chris

## Protocol

### RNA-seq
- RNA sequencing libraries were prepared with ??? (likely an Illumina TruSeq Stranded mRNA kit)
- Paired-end RNA-seq was performed on an Illumina ??? (likely Illumina NovaSeq 6000)


### Transcript Quantification and Differential Expression Analysis
- Sequencing reads were quantified with Salmon in a quasi-mapping-based mode to their corresponding reference genomes **- Run Salmon indexing for the Ensembl hg38 transcriptome

3. Run Salmon quantification for each sample in paired-end read mode

    - Files containing '100_u_IFNa' were replaced with 'induced'
    - Files containing 'Untreated' were replaced with 'uninduced'

4. Transfer the Salmon quantified reads to local

## Procedure (on local)
5. Run DESeq2 tx2gene with the Ensembl hg38 + alt contigs GFF

6. 
