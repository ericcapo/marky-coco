# marky-coco
<p align="justify">
Marky-coco is a ready-to-use pipeline to detect and identify hgcAB genes from raw paired-end and single end fastq files. This pipeline is a collaborative project from researchers of the <a href="https://ercapo.wixsite.com/meta-hg" target="_blank"><b>Meta-Hg working group</b></a> and is paired to the hgcAB gene catalogue <a href="https://smithsonian.figshare.com/articles/dataset/Hg-MATE-Db_v1_01142021/13105370/1?file=26193689" target="_blank"><b>Hg-MATE database</b></a>. The metagenomes are processed with a suite of sofware: fastp to trim and clean the raw reads, megahit for de-novo assembly, bowtie2 to map the cleaned reads to the de-novo assembly, prodigal to predict protein-coding genes, featureCounts to count the number of reads associated to each gene. Finally workflow/genesearch.sh is a custom bash script allowing to detect hgc gene homologs and extract their features (coverage values, taxonomy, amino acid sequences).  In the current version of marky-coco, the script is also providing outputs with detected merA and merB gene homologs (but with no tips yet for manual inspection and taxonomic identification). A step-by-step <b>tutorial</b> is included in the folder tutorial => download the .html and open it on google chrome or firefox.

![Marky-coco_workflow](https://user-images.githubusercontent.com/10795529/213127826-77844383-3a59-41b3-80f6-b7e3ab6b2ae9.png)
</p>

## INSTALL

```
git clone https://github.com/ericcapo/marky-coco.git
cd marky-coco
conda env create -f environment.yml
```
To install conda, read instructions here https://docs.conda.io/projects/conda/en/latest/user-guide/getting-started.html

<br>
<br>
<br>

## BASIC USAGE WITH PAIRED END METAGENOMES
* Copy your fastq files (sample_1.fastq and sample_2.fastq) in the folder marky-coco.
```
cp /remote/folder/sample_1.fastq .
cp /remote/folder/sample_2.fastq .
```
* Activate the conda environment
```
conda activate coco
```
* Run the marky script
```
bash marky_pe.sh sample
```

## BASIC USAGE WITH SINGLE END METAGENOMES
* Copy your fastq file (sample.fastq) in the folder marky-coco. 
```
cp /remote/folder/sample.fastq .
```
* Activate the conda environment
```
conda activate coco
```
* Run the marky script
```
bash marky_se.sh sample
```

## RUN MARKY-COCO WITH TEST METAGENOMES
Download test files MG01 (paired-end fastq files) and MG02 (single-end fastq file)
```
wget https://figshare.com/ndownloader/articles/19221213/versions/2
unzip 2
```
* Run the marky script for the test paired-end metagenome MG01
```
bash marky_pe.sh MG01
```
* Run the marky script for the test single-end metagenome MG02
```
bash marky_se.sh MG02
```

## SLURM USAGE
<p align="justify">
Run the marky_to_slurm file. See how slurm work on your computer/servor here https://blog.ronin.cloud/slurm-intro/
</p>
* For paired-end metagenomes
```
sbatch marky_pe_to_slurm.sh sample
```
* For single-end metagenomes
```
sbatch marky_se_to_slurm.sh sample
```

## ADVANCED USAGE WITH INTERMEDIATE FILES
Standards input files are fastq files but intermediate files can be used because the pipeline in based on a snakemake structure. For marky-coco to work this way, you would need to put your files in the sample_tmp folder as following:
* sample_tmp/sample_P1.fastq & sample_tmp/sample_P2.fastq # cleaned fastq files
* sample_tmp/sample_megahit/final.contigs.fa # megahit outputs
* sample_tmp/sample.bam # bowtie2 outputs
* sample_tmp/sample_proteins.faa # prodigal outputs
* sample_tmp/sample_counts.tsv # featureCounts outputs

## OUTPUTS FROM METAGENOMES
* hgcA_final.txt with columns: gene id, number of reads, gene length (bp), coverage values (nb of read/bp), taxonomic identification (txid) and amino acid sequences. 
* hgcB_final.txt with columnds : gene id, number of reads, gene length (bp), coverage values (nb of read/bp), and amino acid sequences.  
* merA_final.txt with columnds : gene id, number of reads, gene length (bp), coverage values (nb of read/bp), and amino acid sequences. 
* merB_final.txt with columnds : gene id, number of reads, gene length (bp), coverage values (nb of read/bp), and amino acid sequences. 
* rpoBb_final.txt (bacterial rpoB genes) with coverage values (nb of read/bp).
* rpoBa_final.txt (archaeal rpoB genes) with coverage values (nb of read/bp).
* fastp.html with metrics about the cleaning.
* fastp.json with metrics about the cleaning.
* bowtie2.log with metrics about the read mapping.

## RECOVER HGC FROM GENOMES (ISOLATED, SAGS AND MAGS)
<b>"NEW SINCE 30 JUNE 2023</b>To detect the presence of hgc genes in your genomes, you only need the script detect_hgc_from_fna.sh, the db folder of marky-coco and a folder with all your genomes in fna format.

* Copy your fastq file (sample.fastq) in the folder marky-coco. 
```
cp -r /remote/folder .
```
* Activate the conda environment
```
conda activate coco
```
* Run the script
```
bash detect_hgc_from_fna.sh folder
```

The output file is called "detected_hgc_homologs" provided you the list hgcB homologs found in your genomes of interest. Columns includes the id of the genes, the amino acid sequences, genome_id and gene information (hgcA homolog or hgcB homolog). You still have to conifmr if there are true hgcA and hgcB genes using the criteria described in the section below "IMPORTANT NOTES FOR DATA INTERPRETATION". This script do not provide any information about the abundance/coverage of this gene in the genome, neither a specific taxonomy against Hg-MATE database. This is just the presence of this gene pair in your genome(s).

## IMPORTANT NOTES FOR DATA INTERPRETATION
This software and insights into data interpretation are presented in the paper "A consensus protocol for the recovery of mercury methylation genes from metagenomes" <i>Molecular Ecology Resources</i> <a href="https://doi.org/10.1111/1755-0998.13687" target="_blank"><u>doi: 10.1111/1755-0998.13687</u></a>.  
* <b>True hgcA genes</b> are those with amino acids motifs: NVWCAAGK, NVWCASGK, NVWCAGGK, NIWCAAGK, NIWCAGGK or NVWCSAGK
* <b>True hgcB genes</b> are those with the amino acids motifs: CMECGA and CIECGA + colocated with true hgcA genes.
* To detected <b>co-located hgcB genes</b>, examine gene_ids. Ex k141_6000_1 and k141_6000_2 would be co-located genes.
* <b>merA and merB gene´s</b> homologs are detected in this pipeline using HMM profiles created from the amazing database of Christakis, Boyd and Barkay et al. (2021) <a href="https://doi.org/10.3389/fmicb.2021.682605" target="_blank"><u>doi: 10.3389/fmicb.2021.682605</u></a>. Manual inspection of each homolog is required to identify "true" merA and "true" merB genes but tjhe strategy is not described (yet) here. DO NOT USE this output saying you detect merA and merB genes because it will be clearly wrong.
* To <b>normalize hgc coverage values</b>, sum the coverage values obtained from bacterial and archaeal rpoB genes. Read Capo et al. 2023 MER for other strategies.
* To <b>assign NCBI txid to the corresponding taxonomy</b>, you can use the R script below  or do a manual assignment with db/db_txid_2202220  if you have only few  hgcA gene homologs. If the database do not include the txid you found in your sample, check the identity here https://www.ncbi.nlm.nih.gov/taxonomy/

```
R
> db <- read.table("db/db_txid_220220.txt",h=T)
> a <- read.table("{sample}_outputs/{sample}_hgcA_final.txt", h=T)
> b <- merge(db, a, by="txid", all.x=F, all.y=T)
> write.table(b, file="{sample}_outputs/{sample}_hgcA_final2.txt", sep="\t", row.names=F)
> quit()
```

## METHODS
<p align="justify">
The detection, counting and taxonomic identification of hgcAB genes was done with <a href="https://academic.oup.com/bioinformatics/article/34/17/i884/5093234" target="_blank"><b>marky-coco</b></a>. The metagenomes were trimmed and cleaned using <a href="https://academic.oup.com/bioinformatics/article/34/17/i884/5093234" target="_blank"><b>fastp</b></a> (Chen et al. 2018) with following parameters: -q 30 -l 25 --detect_adapter_for_pe --trim_poly_g --trim_poly_x. A de novo single assembly approach was applied using the assembler <a href="https://github.com/voutcn/megahit" target="_blank"><b>megahit</b></a> 1.1.2 (Li et al 2016) with default settings. The annotation of the contigs for prokaryotic protein-coding gene prediction was done with the software <a href="https://github.com/hyattpd/Prodigal" target="_blank"><b>prodigal</b></a> 2.6.3 (Hyatt et al 2010) (Hyatt et al., 2010). The DNA reads were mapped against the contigs with <a href="http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml" target="_blank"><b>bowtie2</b></a> (Langdmead and Salzberg 2012), and the resulting .sam files were converted to .bam files using <a href="http://www.htslib.org/" target="_blank"><b>samtools</b></a> 1.9 (Li et al 2009). The .bam files and the prodigal output .gff file were used to estimate read counts by using <a href="https://rnnh.github.io/bioinfo-notebook/docs/featureCounts.html" target="_blank"><b>featureCounts</b></a>  (Liao et al 2014). In order to detect hgc homologs, HMM profiles derived from the <a href="https://smithsonian.figshare.com/articles/dataset/Hg-MATE-Db_v1_01142021/13105370/1?file=26193689" target="_blank"><b>Hg-MATE.db.v1</b></a> were applied to the amino acid FASTA file generated from each assembly with the function hmmsearch from <a href="http://hmmer.org/" target="_blank"><b>hmmer</b></a> 3.2.1 (Finn et al 2011). The reference package ‘hgcA’ from Hg-MATE.db.v1 was used for phylogenetic analysis of the HgcA amino acid sequences. Briefly, amino acid sequences from gene identified as hgcA gene homolog were (i) compiled in a FASTA file, (ii) aligned to Stockholm formatted alignment of hgcA sequences from the reference package with the function hmmalign from hmmer 3.2.1 (iii) placed onto the HgcA reference tree with the function pplacer and (iv) classified using the functions rppr and guppy_classify from the program <a href="https://matsen.fhcrc.org/pplacer/" target="_blank"><b>pplacer</b></a> (Masten et al. 2010).</p>
