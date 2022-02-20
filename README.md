# marky
Ready-to-use pipeline to detect and identify hgcAB genes


0. INSTALL

manually:
```
git clone https://github.com/ericcapo/marky.git
cd marky
```

1. ACTIVATE THE CONDA ENVIRONMENT
```
Source activate Coco_environment.yml
```

USAGE

Copy your metagenomes in the folder Marky. Input files should look like that {sample}_1.fastq and {sample}_2.fastq. Only paired-end data are supported right now.

```
bash marky.sh {sample}
```

If it doesn´t work, ensure that the marky.sh can be run. Do
```chmod +x marky.sh```





