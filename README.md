# marky-coco
This software is a ready-to-use pipeline to detect and identify hgcAB genes from raw paired-end fastq files. This pipeline is a collaborative project from researchers of the <a href="https://ercapo.wixsite.com/meta-hg" target="_blank"><b>Meta-Hg working group</b></a> and is associated to the hgcAB gene catalogue <a href="https://smithsonian.figshare.com/articles/dataset/Hg-MATE-Db_v1_01142021/13105370/1?file=26193689" target="_blank"><b>Hg-MATE database</b></a>. This software works currently only with paired-end fastq files (sample_1.fastq and sample_2.fastq).


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
* Run the marky_to_slurm function
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


## OUTPUTS
This software will produce a folder {sample}_outputs including:
* a file hgcA_final.txt with: the gene id, the number of reads, gene length (bp), coverage values (nb of read/bp), taxonomic identification (txid) and the amino acid sequences. See <b>IMPORTANT NOTES</b> for data intepretation.
* a file hgcB_final.txt with: the gene id, the number of reads, gene length (bp), coverage values (nb of read/bp) and the amino acid sequences. 
* a file rpoBb_final.txt and a file rpoBa_final.txt with: gene id, the number of reads, gene length (bp), coverage values (nb of read/bp).
* a file fastp.html and and a file fastp.html that are outputs for the fastp step.
* a file bowtie2.log that is output for the bowtie2 step.


## IMPORTANT NOTES
* <b>True hgcA genes</b> include one of the amino acids motifs: NVWCAAGK, NVWCASGK, NVWCAGGK, NIWCAAGK, NIWCAGGK or NVWCSAGK
* <b>True hgcB genes</b> include one of the amino acids motifs: CMECGA and CIEGCA
* To <b>find hgcB genes side-by-side with hgcA genes</b> in the same contig (so co-located in a microbial genome), look the 3nd number in their gene_id (corresponding to contigs id). Note that the 3rd number of the gene_id corresponds to the number of the genes on contigs. Ex k141_6000_1 and k141_6000_2 would be co-located genes.
* To <b>assign NCBI txid to the corresponding taxonomy</b>, you can use the R script below  or use the file db/db_txid_2202220 for manual assignment can be done if you have only few targeted genes.
```
R
db <- read.table("db/db_txid_220220.txt",h=T)
a <- read.table("{sample}_outputs/{sample}_hgcA_final.txt", h=T)
b <- merge(db, a, by="txid", all.x=F, all.y=T)
write.table(b, file="{sample}_outputs/{sample}_hgcA_final2.txt", sep="\t", row.names=F)
```
* To <b>normalize hgc coverage values</b>, sum the coverage values obtained from bacterial and archaeal rpoB genes.
