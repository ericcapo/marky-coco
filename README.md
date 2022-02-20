# marky-coco
This software is a ready-to-use pipeline to detect and identify hgcAB genes from raw paired-end fastq files. This pipeline is a collaborative project from researchers of the <a href="https://ercapo.wixsite.com/meta-hg" target="_blank"><b>Meta-Hg working group</b></a> and is associated to the <a href="https://smithsonian.figshare.com/articles/dataset/Hg-MATE-Db_v1_01142021/13105370/1?file=26193689" target="_blank"><b>Hg-MATE database</b></a>, a hgcAB gene catalogue.


## INSTALL

```
git clone https://github.com/ericcapo/marky-coco.git
cd marky
source conda env create -f environment.yml
```

To install & activate conda, read instructions here https://docs.conda.io/projects/conda/en/latest/user-guide/getting-started.html

## ACTIVATE THE CONDA ENVIRONMENT

```
conda activate environment.yml
```

## USAGE
- Copy your metagenomes (sample_1.fastq and sample_2.fastq) in the folder marky. 
- Run the command below

```
bash marky.sh sample
```

If doesn´t work:
- Only paired-end data are supported right now.
- Input files should look like that sample_1.fastq and sample_2.fastq. 
- Do "chmod +x marky.sh" and re-try. In the following line below,
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
* <b>Not all detected hgcA gene homologs are true hgcA genes</b>. True hgcA genes are only the ones with the following amino acids motifs: NVWCAAGK, NVWCASGK, NVWCAGGK, NIWCAAGK, NIWCAGGK or NVWCSAGK
* <b>Not all detected hgcB gene homologs are true hgB genes</b>. True hgcA genes are only the ones with the following amino acids motifs: CMECGA and CIEGCA
* It is possible to identify hgcB genes found side-by-side with hgcA genes by looking at their gene_id (the second number corresponds to the contigs numbers. If similar between a hgcA and a hgcB gene, that means they are colocated on a microbial genome, the third number of the gene_id corresponds to the number of the genes on the contigs so hgcA and hgcB that would be pairs could be like : k141_6000_1 and k141_6000_2. 
* To <b>assign NCBI txid to the corresponding taxonomy</b>, use the file "db_txid_2202220" located in the subfolder marky_db. A simple R function "merge(db, a, by="txid", all.x=F, all.y=T)" can be used to automatically generated a new column with the corresponding taxonomy. Otherwise, a manual assignment can be done if you have only few targeted genes.
* To normalize your data, rpoB coverage values calculated from both bacteria and archael genomes can be used. To do so, sum the coverage values obtained in the files {sample}_rpoBb_final.txt and {sample}_rpoBa_final.txt and use the obtained value (the sum give you only one value) to normalize your hgcA gene coverage values.
