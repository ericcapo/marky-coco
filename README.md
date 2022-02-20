# marky
Ready-to-use pipeline to detect and identify hgcAB genes


## INSTALL

manually:
```
git clone https://github.com/ericcapo/marky.git
cd marky
Source conda env create -f Coco.environment.yml
```

You need to activate conda on your computer prior to running the "Source" line before. 
Go here for instructions https://docs.conda.io/projects/conda/en/latest/user-guide/getting-started.html

## ACTIVATE THE CONDA ENVIRONMENT

```
Source activate Coco_environment.yml
```

## USAGE

Copy your metagenomes in the folder Marky. Input files should look like that {sample}_1.fastq and {sample}_2.fastq. Only paired-end data are supported right now.

```
bash marky.sh {sample}
```

If it doesn´t work, ensure that the marky.sh can be run. Do
```chmod +x marky.sh```

## OUTPUTS
This software will produce a folder {sample}_outputs that include outputs
- {sample}_hgcA_final.txt includes, for each detected gene, the gene id, the number of reads, gene length (bp), coverage values (nb of read/bp), taxonomic identification (txid) and the amino acid sequences. NOTE THAT NOT ALL HGCA_HOM ARE TRUE HGCA GENES. See the section "IMPORTANT NOTES" BELOW TO DETECT TRUE GENES.
- {sample}_hgcB_final.txt includes, for each detected gene, the number of reads, gene length (bp), coverage values (nb of read/bp) and the amino acid sequences. NOTE THAT NOT ALL HGCB_HOM ARE TRUE HGCA GENES. See the section "IMPORTANT NOTES" BELOW TO DETECT TRUE GENES.
- {sample}_fastp.html and {sample}_fastp.html are outputs files for the fastp step
- {sample}_bowtie2.log provide informations about the number of mapped reads

## IMPORTANT NOTES FOR THE INTERPRETATION OF THE OUTPUTS
* Not all detected hgcA genes are  true hgcA genes. True hgcA genes are only the ones with the following amino acids motifs: NVWCAAGK, NVWCASGK, NVWCAGGK, NIWCAAGK, NIWCAGGK or NVWCSAGK
* Not all detected hgcB genes are true hgB genes. True hgcA genes are only the ones with the following amino acids motifs: CMECGA and CIEGCA
* It is possible to identify hgcB genes found side-by-side with hgcA genes by looking at their gene_id (the second number corresponds to the contigs numbers. If similar between a hgcA and a hgcB gene, that means they are colocated on a microbial genome, the third number of the gene_id corresponds to the number of the genes on the contigs so hgcA and hgcB that would be pairs could be like : k141_6000_1 and k141_6000_2. 
* To <b>assign NCBI txid to the corresponding taxonomy</b>, use the file "db_txid_2202220" located in the subfolder marky_db. A simple R function "merge(db, a, by="txid", all.x=F, all.y=T)" can be used to automatically generated a new column with the corresponding taxonomy. Otherwise, a manual assignment can be done if you have only few targeted genes.
