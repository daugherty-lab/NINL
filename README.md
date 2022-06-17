# RNA-seq processing for Chris

## Protocol

### RNA-seq
- RNA sequencing libraries were prepared with ??? (likely an Illumina TruSeq Stranded mRNA kit)
- Paired-end RNA-seq was performed on an Illumina ??? (likely Illumina NovaSeq 6000)


### Transcript Quantification and Differential Expression Analysis
- Sequencing reads were quantified with Salmon in a quasi-mapping-based mode to the reference genome
- For Chris:
    - Files containing '100_u_IFNa' were replaced with 'induced'
    - Files containing 'Untreated' were replaced with 'uninduced'
- Read quantifications were imported and differentially expressed genes across experimental conditions were identified using the R package DESeq2

### Confirming NIN Mutation
- FASTQs were converted to SAMs using Bowtie2
- SAMs were imported into Geneious and reads were aligned to a subset of NIN transcripts (ENSG00000100503)
