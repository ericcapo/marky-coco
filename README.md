# marky-coco
Marky-coco is a ready-to-use pipeline to detect and identify hgcAB genes from raw paired-end fastq files. This pipeline is a collaborative project from researchers of the <a href="https://ercapo.wixsite.com/meta-hg" target="_blank"><b>Meta-Hg working group</b></a> and is paired to the hgcAB gene catalogue <a href="https://smithsonian.figshare.com/articles/dataset/Hg-MATE-Db_v1_01142021/13105370/1?file=26193689" target="_blank"><b>Hg-MATE database</b></a>. Marky-coco works currently only with paired-end fastq files (sample_1.fastq and sample_2.fastq).


## INSTALL

```
git clone https://github.com/ericcapo/marky-coco.git
cd marky-coco
chmod +x workflow/genesearch.sh
conda env create -f environment.yml
```
To install conda, read instructions here https://docs.conda.io/projects/conda/en/latest/user-guide/getting-started.html


## BASIC USAGE
* Copy your files (sample_1.fastq and sample_2.fastq) in the folder marky. Use 'gunzip' if your fastq are in fastq.gz format.
```
cp /remote/folder/sample_1.fastq .
cp /remote/folder/sample_2.fastq .
```
* Activate the conda environment
```
conda activate coco
```
* Run the marky function
```
bash marky.sh sample
```

## SLURM USAGE
* Run the marky_to_slurm file. You can modify the requested amount of time in the file.
```
sbatch marky_to_slurm.sh sample
```

## WORKFLOW
The raw paired-end fastq files are processed with a suite of sofware. 
* fastp # trim and clean the raw reads
* megahit # de-novo assembly of cleaned reads
* bowtie2 # mapp the cleaned reads to the de-novo assembly
* prodigal # predict protein-coding genes from the de-novo assembly
* featureCounts # count reads associated to each gene
* workflow/genesearch.sh. # custom script detecting hgc gene homologs and extract their features
You can modify parameters for each step in the file workflow/Snakefile or workflow/genesearch.sh.

## INPUTS
Standards input files are
* sample_1.fastq & sample_2.fastq

Alternately, you can use intermediate files but you need to copy then in marky folder as following:
* sample_tmp/sample_P1.fastq & sample_tmp/sample_P2.fastq # cleaned fastq files
* sample_tmp/sample_megahit/final.contigs.fa # megahit outputs
* sample_tmp/sample.bam # bowtie2 outputs
* sample_tmp/sample_proteins.faa # prodigal (or prokka) outputs
* sample_tmp/sample_counts.tsv # featureCounts outputs

## OUTPUTS
This software will produce a folder {sample}_outputs including:
* a file hgcA_final.txt with: the gene id, the number of reads, gene length (bp), coverage values (nb of read/bp), taxonomic identification (txid) and the amino acid sequences. See <b>IMPORTANT NOTES</b> for data intepretation.
* a file hgcB_final.txt with: the gene id, the number of reads, gene length (bp), coverage values (nb of read/bp) and the amino acid sequences. 
* a file rpoBb_final.txt (bacterial rpoB genes) and a file rpoBa_final.txt (archaeal rpoB genes) including coverage values (nb of read/bp).
* a file fastp.html and and a file fastp.html that are outputs for the fastp step.
* a file bowtie2.log that is output for the bowtie2 step.

## IMPORTANT NOTES
* <b>True hgcA genes</b> include one of the amino acids motifs: NVWCAAGK, NVWCASGK, NVWCAGGK, NIWCAAGK, NIWCAGGK or NVWCSAGK
* <b>True hgcB genes</b> include one of the amino acids motifs: CMECGA and CIEGCA
* To <b>find hgcB genes side-by-side with hgcA genes</b> in the same contig (so co-located in a microbial genome), look the 3nd number in their gene_id (corresponding to contigs id). Note that the 3rd number of the gene_id corresponds to the number of the genes on contigs. Ex k141_6000_1 and k141_6000_2 would be co-located genes.
* To <b>assign NCBI txid to the corresponding taxonomy</b>, you can use the R script below  or do a manual assignment with db/db_txid_2202220  if you have only few  hgcA gene homologs. If the database do not include the txid you found in your sample, check the identity here https://www.ncbi.nlm.nih.gov/taxonomy/
```
R
> db <- read.table("db/db_txid_220220.txt",h=T)
> a <- read.table("{sample}_outputs/{sample}_hgcA_final.txt", h=T)
> b <- merge(db, a, by="txid", all.x=F, all.y=T)
> write.table(b, file="{sample}_outputs/{sample}_hgcA_final2.txt", sep="\t", row.names=F)
> quit()
```
* To <b>normalize hgc coverage values</b>, sum the coverage values obtained from bacterial and archaeal rpoB genes.

## METHODS
The detection, counting and taxonomic identification of hgcAB genes was done with <a href="https://academic.oup.com/bioinformatics/article/34/17/i884/5093234" target="_blank"><b>marky-coco</b></a>. The metagenomes were trimmed and cleaned using <a href="https://academic.oup.com/bioinformatics/article/34/17/i884/5093234" target="_blank"><b>fastp</b></a> (Chen et al. 2018) with following parameters: -q 30 -l 25--detect_adapter_for_pe --trim_poly_g --trim_poly_x. A de novo single assembly approach was applied using the assembler <a href="https://github.com/voutcn/megahit" target="_blank"><b>megahit</b></a> 1.1.2 (Li et al 2016) with default settings. The annotation of the contigs for prokaryotic protein-coding gene prediction was done with the software <a href="https://github.com/hyattpd/Prodigal" target="_blank"><b>prodigal</b></a> 2.6.3 (Hyatt et al 2010) (Hyatt et al., 2010). The DNA reads were mapped against the contigs with <a href="http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml" target="_blank"><b>bowtie2</b></a> (Langdmead and Salzberg 2012), and the resulting .sam files were converted to .bam files using <a href="http://www.htslib.org/" target="_blank"><b>samtools</b></a> 1.9 (Li et al 2009). The .bam files and the prodigal output .gff file were used to estimate read counts by using <a href="https://rnnh.github.io/bioinfo-notebook/docs/featureCounts.html" target="_blank"><b>featureCounts</b></a>  (Liao et al 2014). In order to detect hgc homologs, HMM profiles derived from the Hg-MATE db v1. were applied to the amino acid FASTA file generated from each assembly with the function hmmsearch from <a href="http://hmmer.org/" target="_blank"><b>hmmer</b></a> 3.2.1 (Finn et al 2011). Once done, the counts (read numbers), length (bp) and amino acid sequences of each hgc homologs was extracted from Prodigal and featureCounts outputs. We define genes as hgcA if their amino acid sequence includes one of the following motifs: NVWCAAGK, NVWCASGK, NVWCAGGK, NIWCAAGK, NIWCAGGK or NVWCSAGK and as hgcB if their amino acid sequence includes one of the following motifs CMECGA and CIEGCA. The reference package ‘hgcA’ from Hg-MATE.db.v1 was used for phylogenetic analysis of the HgcA amino acid sequences. Briefly, amino acid sequences from gene identified as hgcA gene homolog were (i) compiled in a FASTA file, (ii) aligned to Stockholm formatted alignment of hgcA sequences from the reference package with the function hmmalign from hmmer 3.2.2, (iii) placed onto the HgcA reference tree with the function pplacer and (iv) classified using the functions rppr and guppy_classify from the program <a href="https://matsen.fhcrc.org/pplacer/" target="_blank"><b>pplacer</b></a> (Masten et al. 2010). For more details, see in the README.txt of the <a href="https://smithsonian.figshare.com/articles/dataset/Hg-MATE-Db_v1_01142021/13105370/1?file=26193689" target="_blank"><b>Hg-MATE db</b></a>. Counts of hgc genes were calculated, for each gene and each sample, as the number of reads divided by the length of detected genes (read/bp). 

