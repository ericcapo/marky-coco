#!/usr/bin/env sh
input=$(realpath $1)
output=$(realpath $2)
sample="$3" 

sed -e {'1,2d'} ${input} >  ${sample}_tmp/${sample}_counts.txt
cat ${sample}_tmp/${sample}_counts.txt | awk {'{print $2}'} > ${sample}_tmp/${sample}_counts_sub1.txt
cat ${sample}_tmp/${sample}_counts.txt | awk {'{print $1}'} | cut -d'_' -f2 > ${sample}_tmp/${sample}_counts_sub2.txt
paste -d':' ${sample}_tmp/${sample}_counts_sub1.txt ${sample}_tmp/${sample}_counts_sub2.txt | sed 's/:/_/g' > ${sample}_tmp/${sample}_counts_id.txt
cat ${sample}_tmp/${sample}_counts.txt | awk {'{print $6, $7, $7/$6}'} > ${sample}_tmp/${sample}_counts2.txt
paste ${sample}_tmp/${sample}_counts_id.txt ${sample}_tmp/${sample}_counts2.txt > ${sample}_tmp/${sample}_counts3.txt

hmmsearch -o ${sample}_tmp/${sample}_hgcA.txt --tblout ${sample}_outputs/${sample}_hgcA_hmmer.out db/Hg-MATE-Db.v1/Hg-MATE-Db.v1.01142021_ISOCELMAG_HgcA.hmm ${sample}_tmp/${sample}_proteins.faa
grep -v '^#' ${sample}_outputs/${sample}_hgcA_hmmer.out | awk {'{print $1}'} | sort > ${sample}_tmp/${sample}_hgcA_geneid.txt
seqtk subseq ${sample}_tmp/${sample}_proteins.faa ${sample}_tmp/${sample}_hgcA_geneid.txt | awk {'{gsub(/*$/,""); print}'} > ${sample}_tmp/${sample}_hgcA_proteins.faa
awk '/^>/ {{printf("%s%s\t",(N>0?"\n":""),$0);N++;next;}} {{printf("%s",$0);}} END {{printf("\n");}}' < ${sample}_tmp/${sample}_hgcA_proteins.faa  |sort > ${sample}_tmp/${sample}_hgcA_proteins2.faa 
sed 's/>//' ${sample}_tmp/${sample}_hgcA_proteins2.faa > ${sample}_tmp/${sample}_hgcA_proteins3.faa
awk {'{print $1, $10}'} ${sample}_tmp/${sample}_hgcA_proteins3.faa > ${sample}_tmp/${sample}_hgcA_seq.txt
grep -w -f ${sample}_tmp/${sample}_hgcA_geneid.txt ${sample}_tmp/${sample}_counts3.txt | sort > ${sample}_tmp/${sample}_hgcA_counts.txt
hmmalign -o ${sample}_tmp/${sample}_hgcA_proteins.sto --mapali db/Hg-MATE-Db.v1/Hg-MATE-Db.v1.ISOCELMAG_HgcA_full.refpkg/Hg-MATE-Db.v1.ISOCELMAG_HgcA_full.stockholm db/Hg-MATE-Db.v1/Hg-MATE-Db.v1.ISOCELMAG_HgcA_full.refpkg/Hg-MATE-Db.v1.ISOCELMAG_HgcA_full.hmm ${sample}_tmp/${sample}_hgcA_proteins.faa
pplacer --keep-at-most 1 --max-pend 1 -p -c db/Hg-MATE-Db.v1/Hg-MATE-Db.v1.ISOCELMAG_HgcA_full.refpkg ${sample}_tmp/${sample}_hgcA_proteins.sto -o ${sample}_tmp/${sample}_hgcA_proteins.jplace
rppr prep_db --sqlite ${sample}_tmp/${sample}_hgcA_classify_output -c db/Hg-MATE-Db.v1/Hg-MATE-Db.v1.ISOCELMAG_HgcA_full.refpkg
guppy classify -c db/Hg-MATE-Db.v1/Hg-MATE-Db.v1.ISOCELMAG_HgcA_full.refpkg --pp --sqlite ${sample}_tmp/${sample}_hgcA_classify_output ${sample}_tmp/${sample}_hgcA_proteins.jplace
guppy to_csv --point-mass --pp -o ${sample}_tmp/${sample}_hgcA_classifications.csv ${sample}_tmp/${sample}_hgcA_proteins.jplace
guppy tog --pp -o ${sample}_outputs/${sample}_hgcA_tree.nwk ${sample}_tmp/${sample}_hgcA_proteins.jplace
sed -e '1,1d' ${sample}_tmp/${sample}_hgcA_classifications.csv > ${sample}_tmp/${sample}_hgcA_txid.txt
awk -vFPAT='([^,]*)|({"[^"]+"})' '{print $2, $11}' ${sample}_tmp/${sample}_hgcA_txid.txt  | sort > ${sample}_tmp/${sample}_hgcA_txid2.txt
paste ${sample}_tmp/${sample}_hgcA_counts.txt ${sample}_tmp/${sample}_hgcA_txid2.txt ${sample}_tmp/${sample}_hgcA_seq.txt | sed 's/$/ hgcA_hom/'  > ${sample}_tmp/${sample}_hgcA_almost.txt
awk {'{print $1,$2,$3,$4,$6,$8,$9}'} ${sample}_tmp/${sample}_hgcA_almost.txt > ${sample}_tmp/${sample}_hgcA_almost2.txt
echo 'gene_id length read cov txid sequence gene_type' | cat - ${sample}_tmp/${sample}_hgcA_almost2.txt > ${output}

hmmsearch -o ${sample}_tmp/${sample}_hgcB.txt --tblout ${sample}_outputs/${sample}_hgcB_hmmer.out db/Hg-MATE-Db.v1/Hg-MATE-Db.v1.01142021_ISOCELMAG_HgcB.hmm ${sample}_tmp/${sample}_proteins.faa
grep -v '^#' ${sample}_outputs/${sample}_hgcB_hmmer.out | awk '{print $1}' | sort > ${sample}_tmp/${sample}_hgcB_geneid.txt
seqtk subseq ${sample}_tmp/${sample}_proteins.faa ${sample}_tmp/${sample}_hgcB_geneid.txt | awk '{gsub(/*$/,""); print}' > ${sample}_tmp/${sample}_hgcB_proteins.faa
awk '/^>/ {printf("%s%s\t",(N>0?"\n":""),$0);N++;next;} {printf("%s",$0);} END {printf("\n");}' < ${sample}_tmp/${sample}_hgcB_proteins.faa  |sort > ${sample}_tmp/${sample}_hgcB_proteins2.faa 
sed 's/>//' ${sample}_tmp/${sample}_hgcB_proteins2.faa > ${sample}_tmp/${sample}_hgcB_proteins3.faa
awk '{print $1, $10}' ${sample}_tmp/${sample}_hgcB_proteins3.faa > ${sample}_tmp/${sample}_hgcB_seq.txt
grep -w -f ${sample}_tmp/${sample}_hgcB_geneid.txt ${sample}_tmp/${sample}_counts3.txt | sort > ${sample}_tmp/${sample}_hgcB_counts.txt
paste ${sample}_tmp/${sample}_hgcB_counts.txt ${sample}_tmp/${sample}_hgcB_seq.txt | sed 's/$/ hgcB_hom/'  > ${sample}_tmp/${sample}_hgcB_almost.txt
awk '{print $1,$2,$3,$4,$6,$7}' ${sample}_tmp/${sample}_hgcB_almost.txt > ${sample}_tmp/${sample}_hgcB_almost2.txt
echo "gene_id length read cov sequence gene_type" | cat - ${sample}_tmp/${sample}_hgcB_almost2.txt > ${sample}_outputs/${sample}_hgcB_final.txt

hmmsearch -o ${sample}_tmp/${sample}_merA.txt --tblout ${sample}_outputs/${sample}_merA_hmmer.out db/merAB_Christakis2021/211026_Christakis_reduced_merA_msa.hmm ${sample}_tmp/${sample}_proteins.faa
grep -v '^#' ${sample}_outputs/${sample}_merA_hmmer.out | awk {'{print $1}'} | sort > ${sample}_tmp/${sample}_merA_geneid.txt
seqtk subseq ${sample}_tmp/${sample}_proteins.faa ${sample}_tmp/${sample}_merA_geneid.txt | awk {'{gsub(/*$/,""); print}'} > ${sample}_tmp/${sample}_merA_proteins.faa
awk '/^>/ {{printf("%s%s\t",(N>0?"\n":""),$0);N++;next;}} {{printf("%s",$0);}} END {{printf("\n");}}' < ${sample}_tmp/${sample}_merA_proteins.faa  |sort > ${sample}_tmp/${sample}_merA_proteins2.faa 
sed 's/>//' ${sample}_tmp/${sample}_merA_proteins2.faa > ${sample}_tmp/${sample}_merA_proteins3.faa
awk {'{print $1, $10}'} ${sample}_tmp/${sample}_merA_proteins3.faa > ${sample}_tmp/${sample}_merA_seq.txt
grep -w -f ${sample}_tmp/${sample}_merA_geneid.txt ${sample}_tmp/${sample}_counts3.txt | sort > ${sample}_tmp/${sample}_merA_counts.txt
paste ${sample}_tmp/${sample}_merA_counts.txt ${sample}_tmp/${sample}_merA_seq.txt | sed 's/$/ merA_hom/'  > ${sample}_outputs/${sample}_merA_homologs.txt

hmmsearch -o ${sample}_tmp/${sample}_merB.txt --tblout ${sample}_outputs/${sample}_merB_hmmer.out db/merAB_Christakis2021/211026_Christakis_reduced_merB_msa.hmm ${sample}_tmp/${sample}_proteins.faa
grep -v '^#' ${sample}_outputs/${sample}_merB_hmmer.out | awk {'{print $1}'} | sort > ${sample}_tmp/${sample}_merB_geneid.txt
seqtk subseq ${sample}_tmp/${sample}_proteins.faa ${sample}_tmp/${sample}_merB_geneid.txt | awk {'{gsub(/*$/,""); print}'} > ${sample}_tmp/${sample}_merB_proteins.faa
awk '/^>/ {{printf("%s%s\t",(N>0?"\n":""),$0);N++;next;}} {{printf("%s",$0);}} END {{printf("\n");}}' < ${sample}_tmp/${sample}_merB_proteins.faa  |sort > ${sample}_tmp/${sample}_merB_proteins2.faa 
sed 's/>//' ${sample}_tmp/${sample}_merB_proteins2.faa > ${sample}_tmp/${sample}_merB_proteins3.faa
awk {'{print $1, $10}'} ${sample}_tmp/${sample}_merB_proteins3.faa > ${sample}_tmp/${sample}_merB_seq.txt
grep -w -f ${sample}_tmp/${sample}_merB_geneid.txt ${sample}_tmp/${sample}_counts3.txt | sort > ${sample}_tmp/${sample}_merB_counts.txt
paste ${sample}_tmp/${sample}_merB_counts.txt ${sample}_tmp/${sample}_merB_seq.txt | sed 's/$/ merB_hom/'  > ${sample}_outputs/${sample}_merB_homologs.txt

hmmsearch --cut_tc db/rpoB/TIGR02013.hmm ${sample}_tmp/${sample}_proteins.faa > ${sample}_tmp/${sample}_rpoBb.txt
sed -n -e '18,/Domain/p' ${sample}_tmp/${sample}_rpoBb.txt | head -n -3  | sed 's/ \+/|/g' | cut -f10 -d"|" > ${sample}_tmp/${sample}_rpoBb2.txt
grep -w -f ${sample}_tmp/${sample}_rpoBb2.txt ${sample}_tmp/${sample}_counts3.txt | sort > ${sample}_tmp/${sample}_rpoBb3.txt
echo "gene_id length read cov" | cat - ${sample}_tmp/${sample}_rpoBb3.txt > ${sample}_outputs/${sample}_rpoBb_final.txt

hmmsearch --cut_tc db/rpoB/TIGR03670.hmm ${sample}_tmp/${sample}_proteins.faa > ${sample}_tmp/${sample}_rpoBa.txt
sed -n -e '18,/Domain/p' ${sample}_tmp/${sample}_rpoBa.txt | head -n -3  | sed 's/ \+/|/g' | cut -f10 -d"|" > ${sample}_tmp/${sample}_rpoBa2.txt
grep -w -f ${sample}_tmp/${sample}_rpoBa2.txt ${sample}_tmp/${sample}_counts3.txt | sort > ${sample}_tmp/${sample}_rpoBa3.txt
echo "gene_id length read cov" | cat - ${sample}_tmp/${sample}_rpoBa3.txt > ${sample}_outputs/${sample}_rpoBa_final.txt
