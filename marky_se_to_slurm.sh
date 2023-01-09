#!/bin/bash				
#SBATCH -A snic2022-5-51				
#SBATCH -n 12				
#SBATCH --time=10:00:00	

module load conda
source conda_init.sh
conda activate coco

INPUT=$1

mkdir ${INPUT}_outputs

fastp -i ${INPUT}.fastq -q 30 -l 25 -o ${INPUT}_cleaned.fastq -w 6 -h ${INPUT}_outputs/${INPUT}_fastp.html -j ${INPUT}_outputs/${INPUT}_fastp.json --detect_adapter_for_pe --trim_poly_g --trim_poly_x  
megahit -r ${INPUT}_cleaned.fastq -t 6 -o ${INPUT}_tmp
bowtie2-build ${INPUT}_tmp/final.contigs.fa ${INPUT}_tmp/${INPUT}.index
bowtie2 -q ${INPUT}_cleaned.fastq -x ${INPUT}_tmp/${INPUT}.index -p 6 | samtools view -Sb | samtools sort > ${INPUT}_tmp/${INPUT}_automapped_sorted_bam

prodigal -i ${INPUT}_tmp/final.contigs.fa -o ${INPUT}_tmp/${INPUT}_genes.gff -f gff -a ${INPUT}_tmp/${INPUT}_proteins.faa
awk '/^>/ {printf("%s%s\t",(N>0?"\n":""),$0);N++;next;} {printf("%s",$0);} END {printf("\n");}' < ${INPUT}_tmp/${INPUT}_proteins.faa > ${INPUT}_tmp/${INPUT}_proteins2.faa
awk '{print $NF}' ${INPUT}_tmp/${INPUT}_proteins2.faa > ${INPUT}_tmp/${INPUT}_protseq.txt
awk -F '#' '{print $1 }' ${INPUT}_tmp/${INPUT}_proteins2.faa  > ${INPUT}_tmp/${INPUT}_geneid.txt
paste ${INPUT}_tmp/${INPUT}_geneid.txt ${INPUT}_tmp/${INPUT}_protseq.txt | cut -f 1,2 > ${INPUT}_tmp/${INPUT}_prot1.txt
sed 's/>//' ${INPUT}_tmp/${INPUT}_prot1.txt > ${INPUT}_tmp/${INPUT}_prot2.txt
awk '{sub("-", "", $2); print}' < ${INPUT}_tmp/${INPUT}_prot2.txt > ${INPUT}_tmp/${INPUT}_prot2b.txt
awk '{gsub(/*$/,""); print}' ${INPUT}_tmp/${INPUT}_prot2b.txt > ${INPUT}_tmp/${INPUT}_prot2c.txt
echo -e "gene_id prot" | cat - ${INPUT}_tmp/${INPUT}_prot2c.txt  > ${INPUT}_tmp/${INPUT}_prot3.txt

featureCounts -t CDS -o ${INPUT}_tmp/${INPUT}_counts.tsv -g ID -a ${INPUT}_tmp/${INPUT}_genes.gff ${INPUT}_tmp/${INPUT}_automapped_sorted_bam
sed -e '1,2d' ${INPUT}_tmp/${INPUT}_counts.tsv > ${INPUT}_tmp/${INPUT}_counts.txt
cat ${INPUT}_tmp/${INPUT}_counts.txt | awk '{print $2}' > ${INPUT}_tmp/${INPUT}_id1.txt
cat ${INPUT}_tmp/${INPUT}_counts.txt | awk '{print $1}' | cut -d'_' -f2 > ${INPUT}_tmp/${INPUT}_id2.txt
paste -d":" ${INPUT}_tmp/${INPUT}_id1.txt ${INPUT}_tmp/${INPUT}_id2.txt | sed 's/:/_/g' > ${INPUT}_tmp/${INPUT}_id3.txt
cat ${INPUT}_tmp/${INPUT}_counts.txt | awk '{print $7}' > ${INPUT}_tmp/${INPUT}_counts1.txt
cat ${INPUT}_tmp/${INPUT}_counts.txt | awk '{print $6}' > ${INPUT}_tmp/${INPUT}_length.txt
paste ${INPUT}_tmp/${INPUT}_id3.txt ${INPUT}_tmp/${INPUT}_counts1.txt ${INPUT}_tmp/${INPUT}_length.txt > ${INPUT}_tmp/${INPUT}_counts2.txt
echo -e "gene_id seq length" | cat - ${INPUT}_tmp/${INPUT}_counts2.txt  > ${INPUT}_tmp/${INPUT}_counts3.txt

hmmsearch -o ${INPUT}_tmp/${INPUT}_hgcA.txt --tblout ${INPUT}_outputs/${INPUT}_hgcA_hmmer.out db/Hg-MATE-Db.v1/Hg-MATE-Db.v1.01142021_ISOCELMAG_HgcA.hmm ${INPUT}_tmp/${INPUT}_proteins.faa
grep -v '^#' ${INPUT}_outputs/${INPUT}_hgcA_hmmer.out | awk {'{print $1}'} | sort > ${INPUT}_tmp/${INPUT}_hgcA_geneid.txt
seqtk subseq ${INPUT}_tmp/${INPUT}_proteins.faa ${INPUT}_tmp/${INPUT}_hgcA_geneid.txt | awk {'{gsub(/*$/,""); print}'} > ${INPUT}_tmp/${INPUT}_hgcA_proteins.faa
awk '/^>/ {{printf("%s%s\t",(N>0?"\n":""),$0);N++;next;}} {{printf("%s",$0);}} END {{printf("\n");}}' < ${INPUT}_tmp/${INPUT}_hgcA_proteins.faa  |sort > ${INPUT}_tmp/${INPUT}_hgcA_proteins2.faa 
sed 's/>//' ${INPUT}_tmp/${INPUT}_hgcA_proteins2.faa > ${INPUT}_tmp/${INPUT}_hgcA_proteins3.faa
awk {'{print $1, $10}'} ${INPUT}_tmp/${INPUT}_hgcA_proteins3.faa > ${INPUT}_tmp/${INPUT}_hgcA_seq.txt
grep -f ${INPUT}_tmp/${INPUT}_hgcA_geneid.txt ${INPUT}_tmp/${INPUT}_counts3.txt | sort > ${INPUT}_tmp/${INPUT}_hgcA_counts.txt
hmmalign -o ${INPUT}_tmp/${INPUT}_hgcA_proteins.sto --mapali db/Hg-MATE-Db.v1/Hg-MATE-Db.v1.ISOCELMAG_HgcA_full.refpkg/Hg-MATE-Db.v1.ISOCELMAG_HgcA_full.stockholm db/Hg-MATE-Db.v1/Hg-MATE-Db.v1.ISOCELMAG_HgcA_full.refpkg/Hg-MATE-Db.v1.ISOCELMAG_HgcA_full.hmm ${INPUT}_tmp/${INPUT}_hgcA_proteins.faa
pplacer --keep-at-most 1 --max-pend 1 -p -c db/Hg-MATE-Db.v1/Hg-MATE-Db.v1.ISOCELMAG_HgcA_full.refpkg ${INPUT}_tmp/${INPUT}_hgcA_proteins.sto -o ${INPUT}_tmp/${INPUT}_hgcA_proteins.jplace
rppr prep_db --sqlite ${INPUT}_tmp/${INPUT}_hgcA_classify_output -c db/Hg-MATE-Db.v1/Hg-MATE-Db.v1.ISOCELMAG_HgcA_full.refpkg
guppy classify -c db/Hg-MATE-Db.v1/Hg-MATE-Db.v1.ISOCELMAG_HgcA_full.refpkg --pp --sqlite ${INPUT}_tmp/${INPUT}_hgcA_classify_output ${INPUT}_tmp/${INPUT}_hgcA_proteins.jplace
guppy to_csv --point-mass --pp -o ${INPUT}_tmp/${INPUT}_hgcA_classifications.csv ${INPUT}_tmp/${INPUT}_hgcA_proteins.jplace
guppy tog --pp -o ${INPUT}_outputs/${INPUT}_hgcA_tree.nwk ${INPUT}_tmp/${INPUT}_hgcA_proteins.jplace
sed -e '1,1d' ${INPUT}_tmp/${INPUT}_hgcA_classifications.csv > ${INPUT}_tmp/${INPUT}_hgcA_txid.txt
awk -vFPAT='([^,]*)|({"[^"]+"})' '{print $2, $11}' ${INPUT}_tmp/${INPUT}_hgcA_txid.txt  | sort > ${INPUT}_tmp/${INPUT}_hgcA_txid2.txt
paste ${INPUT}_tmp/${INPUT}_hgcA_counts.txt ${INPUT}_tmp/${INPUT}_hgcA_txid2.txt ${INPUT}_tmp/${INPUT}_hgcA_seq.txt | sed 's/$/ hgcA_hom/'  > ${INPUT}_tmp/${INPUT}_hgcA_almost.txt
awk {'{print $1,$2,$3,$5,$7,$8}'} ${INPUT}_tmp/${INPUT}_hgcA_almost.txt > ${INPUT}_tmp/${INPUT}_hgcA_almost2.txt
echo -e 'gene_id length read cov txid sequence gene_type' | cat - ${INPUT}_tmp/${INPUT}_hgcA_almost2.txt > ${INPUT}_outputs/${INPUT}_hgcA_final.txt

hmmsearch -o ${INPUT}_tmp/${INPUT}_hgcB.txt --tblout ${INPUT}_outputs/${INPUT}_hgcB_hmmer.out db/Hg-MATE-Db.v1/Hg-MATE-Db.v1.01142021_ISOCELMAG_HgcB.hmm ${INPUT}_tmp/${INPUT}_proteins.faa
grep -v '^#' ${INPUT}_outputs/${INPUT}_hgcB_hmmer.out | awk '{print $1}' | sort > ${INPUT}_tmp/${INPUT}_hgcB_geneid.txt
seqtk subseq ${INPUT}_tmp/${INPUT}_proteins.faa ${INPUT}_tmp/${INPUT}_hgcB_geneid.txt | awk '{gsub(/*$/,""); print}' > ${INPUT}_tmp/${INPUT}_hgcB_proteins.faa
awk '/^>/ {printf("%s%s\t",(N>0?"\n":""),$0);N++;next;} {printf("%s",$0);} END {printf("\n");}' < ${INPUT}_tmp/${INPUT}_hgcB_proteins.faa  |sort > ${INPUT}_tmp/${INPUT}_hgcB_proteins2.faa 
sed 's/>//' ${INPUT}_tmp/${INPUT}_hgcB_proteins2.faa > ${INPUT}_tmp/${INPUT}_hgcB_proteins3.faa
awk '{print $1, $10}' ${INPUT}_tmp/${INPUT}_hgcB_proteins3.faa > ${INPUT}_tmp/${INPUT}_hgcB_seq.txt
grep -f ${INPUT}_tmp/${INPUT}_hgcB_geneid.txt ${INPUT}_tmp/${INPUT}_counts3.txt | sort > ${INPUT}_tmp/${INPUT}_hgcB_counts.txt
paste ${INPUT}_tmp/${INPUT}_hgcB_counts.txt ${INPUT}_tmp/${INPUT}_hgcB_seq.txt | sed 's/$/ hgcB_hom/'  > ${INPUT}_tmp/${INPUT}_hgcB_almost.txt
awk '{print $1,$2,$3,$5,$6}' ${INPUT}_tmp/${INPUT}_hgcB_almost.txt > ${INPUT}_tmp/${INPUT}_hgcB_almost2.txt
echo -e "gene_id length read cov sequence gene_type" | cat - ${INPUT}_tmp/${INPUT}_hgcB_almost2.txt > ${INPUT}_outputs/${INPUT}_hgcB_final.txt

hmmsearch -o ${INPUT}_tmp/${INPUT}_merA.txt --tblout ${INPUT}_outputs/${INPUT}_merA_hmmer.out db/merAB_Christakis2021/211026_Christakis_reduced_merA_msa.hmm ${INPUT}_tmp/${INPUT}_proteins.faa
grep -v '^#' ${INPUT}_outputs/${INPUT}_merA_hmmer.out | awk {'{print $1}'} | sort > ${INPUT}_tmp/${INPUT}_merA_geneid.txt
seqtk subseq ${INPUT}_tmp/${INPUT}_proteins.faa ${INPUT}_tmp/${INPUT}_merA_geneid.txt | awk {'{gsub(/*$/,""); print}'} > ${INPUT}_tmp/${INPUT}_merA_proteins.faa
awk '/^>/ {{printf("%s%s\t",(N>0?"\n":""),$0);N++;next;}} {{printf("%s",$0);}} END {{printf("\n");}}' < ${INPUT}_tmp/${INPUT}_merA_proteins.faa  |sort > ${INPUT}_tmp/${INPUT}_merA_proteins2.faa 
sed 's/>//' ${INPUT}_tmp/${INPUT}_merA_proteins2.faa > ${INPUT}_tmp/${INPUT}_merA_proteins3.faa
awk {'{print $1, $10}'} ${INPUT}_tmp/${INPUT}_merA_proteins3.faa > ${INPUT}_tmp/${INPUT}_merA_seq.txt
grep -f ${INPUT}_tmp/${INPUT}_merA_geneid.txt ${INPUT}_tmp/${INPUT}_counts3.txt | sort > ${INPUT}_tmp/${INPUT}_merA_counts.txt
paste ${INPUT}_tmp/${INPUT}_merA_counts.txt ${INPUT}_tmp/${INPUT}_merA_seq.txt | sed 's/$/ merA_hom/'  > ${INPUT}_outputs/${INPUT}_merA_homologs.txt

hmmsearch -o ${INPUT}_tmp/${INPUT}_merB.txt --tblout ${INPUT}_outputs/${INPUT}_merB_hmmer.out db/merAB_Christakis2021/211026_Christakis_reduced_merB_msa.hmm ${INPUT}_tmp/${INPUT}_proteins.faa
grep -v '^#' ${INPUT}_outputs/${INPUT}_merB_hmmer.out | awk {'{print $1}'} | sort > ${INPUT}_tmp/${INPUT}_merB_geneid.txt
seqtk subseq ${INPUT}_tmp/${INPUT}_proteins.faa ${INPUT}_tmp/${INPUT}_merB_geneid.txt | awk {'{gsub(/*$/,""); print}'} > ${INPUT}_tmp/${INPUT}_merB_proteins.faa
awk '/^>/ {{printf("%s%s\t",(N>0?"\n":""),$0);N++;next;}} {{printf("%s",$0);}} END {{printf("\n");}}' < ${INPUT}_tmp/${INPUT}_merB_proteins.faa  |sort > ${INPUT}_tmp/${INPUT}_merB_proteins2.faa 
sed 's/>//' ${INPUT}_tmp/${INPUT}_merB_proteins2.faa > ${INPUT}_tmp/${INPUT}_merB_proteins3.faa
awk {'{print $1, $10}'} ${INPUT}_tmp/${INPUT}_merB_proteins3.faa > ${INPUT}_tmp/${INPUT}_merB_seq.txt
grep -f ${INPUT}_tmp/${INPUT}_merB_geneid.txt ${INPUT}_tmp/${INPUT}_counts3.txt | sort > ${INPUT}_tmp/${INPUT}_merB_counts.txt
paste ${INPUT}_tmp/${INPUT}_merB_counts.txt ${INPUT}_tmp/${INPUT}_merB_seq.txt | sed 's/$/ merB_hom/'  > ${INPUT}_outputs/${INPUT}_merB_homologs.txt

hmmsearch --cut_tc db/rpoB/TIGR02013.hmm ${INPUT}_tmp/${INPUT}_proteins.faa > ${INPUT}_tmp/${INPUT}_rpoBb.txt
sed -n -e '18,/Domain/p' ${INPUT}_tmp/${INPUT}_rpoBb.txt | head -n -3  | sed 's/ \+/|/g' | cut -f10 -d"|" > ${INPUT}_tmp/${INPUT}_rpoBb2.txt
grep -f ${INPUT}_tmp/${INPUT}_rpoBb2.txt ${INPUT}_tmp/${INPUT}_counts3.txt | sort > ${INPUT}_tmp/${INPUT}_rpoBb3.txt
echo -e "gene_id length read cov" | cat - ${INPUT}_tmp/${INPUT}_rpoBb3.txt > ${INPUT}_outputs/${INPUT}_rpoBb_final.txt

hmmsearch --cut_tc db/rpoB/TIGR03670.hmm ${INPUT}_tmp/${INPUT}_proteins.faa > ${INPUT}_tmp/${INPUT}_rpoBa.txt
sed -n -e '18,/Domain/p' ${INPUT}_tmp/${INPUT}_rpoBa.txt | head -n -3  | sed 's/ \+/|/g' | cut -f10 -d"|" > ${INPUT}_tmp/${INPUT}_rpoBa2.txt
grep -f ${INPUT}_tmp/${INPUT}_rpoBa2.txt ${INPUT}_tmp/${INPUT}_counts3.txt | sort > ${INPUT}_tmp/${INPUT}_rpoBa3.txt
echo -e "gene_id length read cov" | cat - ${INPUT}_tmp/${INPUT}_rpoBa3.txt > ${INPUT}_outputs/${INPUT}_rpoBa_final.txt

# Cleaning folder from temporary files
#rm -r ${INPUT}_tmp
