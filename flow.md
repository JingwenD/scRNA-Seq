# scRNASeq workflow




# 1.Cellranger count

## 1.1 Cellranger count for 5' transcriptome
```{bash}
#!/bin/bash
#
#$ -cwd
#$ -N count_5p
#$ -pe smp 16
#$ -q psoriasis.q
#$ -l h_vmem = 128G


cellranger count --id Y1_5p --fastqs /mnt/lustre/user/jdeng/PBMCsc/5p --sample Y1 --transcriptome /mnt/lustre/user/jdeng/PBMCsc/refdata-gex-GRCh38-2020-A --chemistry fiveprime

```

## 1.2 Cellranger count for 3' transcriptome
```{bash}
#!/bin/bash
#
#$ -cwd
#$ -N count_3p
#$ -pe smp 16
#$ -q psoriasis.q
#$ -l h_vmem = 128G


cellranger count --id Y1_3p --fastqs /mnt/lustre/user/jdeng/lesion/3p --sample Y1 --transcriptome /mnt/lustre/user/jdeng/PBMCsc/refdata-gex-GRCh38-2020-A --chemistry threeprime

```
--fastqs	Path of the FASTQ folder 
e.g. /home/jdoe/runs/HAWT7ADXX/outs/fastq_path
Can take multiple comma-separated paths, which is helpful if the same library was sequenced on multiple flowcells.
Doing this will treat all reads from the library, across flowcells, as one sample

## 1.3 Cellranger count for TCR
```{bash}
#!/bin/bash
#
#$ -cwd
#$ -N count_TCR
#$ -pe smp 8
#$ -q psoriasis.q
#$ -l h_vmem = 64G


cellranger vdj --id Y1_TCR --fastqs /mnt/lustre/user/jdeng/PBMCsc/TCR/Y1 --sample Y1_TCR --reference /mnt/lustre/user/jdeng/PBMCsc/ref/refdata-cellranger-vdj-GRCh38-alts-ensembl-5.0.0
 --localcores 8 --localmem 64
 
cellranger vdj --id Y2_TCR --fastqs /mnt/lustre/user/jdeng/PBMCsc/TCR/Y2 --sample Y2_TCR --reference /mnt/lustre/user/jdeng/PBMCsc/ref/refdata-cellranger-vdj-GRCh38-alts-ensembl-5.0.0
 --localcores 8 --localmem 64
```

## 1.4 Cellranger count for BCR
```{bash}
#!/bin/bash
#
#$ -cwd
#$ -N count_BCR
#$ -pe smp 8
#$ -q psoriasis.q
#$ -l h_vmem = 64G


cellranger vdj --id Y1_BCR --fastqs /mnt/lustre/user/jdeng/PBMCsc/BCR/Y1 --sample Y1_BCR --reference /mnt/lustre/user/jdeng/PBMCsc/ref/refdata-cellranger-vdj-GRCh38-alts-ensembl-5.0.0
 --localcores 8 --localmem 64
 
cellranger vdj --id Y2_BCR --fastqs /mnt/lustre/user/jdeng/PBMCsc/BCR/Y2 --sample Y2_BCR --reference /mnt/lustre/user/jdeng/PBMCsc/ref/refdata-cellranger-vdj-GRCh38-alts-ensembl-5.0.0
 --localcores 8 --localmem 64
```
