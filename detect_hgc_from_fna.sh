#!/usr/bin/env sh

mkdir tmp tmp/tmp01 tmp/tmp02 tmp/tmp03 tmp/tmp04 tmp/tmp05 tmp/tmp06 tmp/tmp07 tmp/tmp08 tmp/tmp09 tmp/tmp10 tmp/tmp11 tmp/tmp12 tmp/tmp13 tmp/tmp14 tmp/tmp15 tmp/tmp16 outputs

inputs="$1"

FILES=${inputs}/*
for f in $FILES
do
lb1file="tmp/tmp01/${f##*/}.$$"
lb2file="tmp/tmp02/${f##*/}.$$"
lb3file="tmp/tmp03/${f##*/}.$$"
lb4file="tmp/tmp04/${f##*/}.$$"
lb5file="tmp/tmp05${f##*/}.$$"
lb6file="tmp/tmp06${f##*/}.$$"
lb7file="tmp/tmp07/${f##*/}.$$"
lb8file="tmp/tmp08/${f##*/}.$$"
lb9file="tmp/tmp09/${f##*/}.$$"
lb10file="tmp/tmp10/${f##*/}.$$"
lb11file="tmp/tmp11/${f##*/}.$$"
lb12file="tmp/tmp12/${f##*/}.$$"
lb13file="tmp/tmp13/${f##*/}.$$"
lb14file="tmp/tmp14/${f##*/}.$$"
lb15file="tmp/tmp15/${f##*/}.$$"
lb16file="tmp/tmp16/${f##*/}.$$"
lb17file="outputs/${f##*/}.$$"
input_name="${f##*/}"
prodigal -i "$f" -o "${lb1file}" -f gff -a "${lb2file}"
hmmsearch -o "${lb3file}" --tblout "${lb4file}" db/Hg-MATE-Db.v1/Hg-MATE-Db.v1.01142021_ISOCELMAG_HgcA.hmm "${lb2file}"
grep -v '^#' "${lb4file}" | awk {'{print $1}'} | sort > "${lb5file}"
seqtk subseq "${lb2file}" "${lb5file}" | awk {'{gsub(/*$/,""); print}'} > "${lb6file}"
awk '/^>/ {{printf("%s%s\t",(N>0?"\n":""),$0);N++;next;}} {{printf("%s",$0);}} END {{printf("\n");}}' < "${lb6file}"  |sort > "${lb7file}"
sed 's/>//' "${lb7file}" > "${lb8file}"
awk -v input_name="$input_name" '{print $1, $10, "hgcA_homologs", input_name}' "${lb8file}" > "${lb9file}"

hmmsearch -o "${lb10file}" --tblout "${lb11file}" db/Hg-MATE-Db.v1/Hg-MATE-Db.v1.01142021_ISOCELMAG_HgcB.hmm "${lb2file}"
grep -v '^#' "${lb11file}" | awk {'{print $1}'} | sort > "${lb12file}"
seqtk subseq "${lb2file}" "${lb12file}" | awk {'{gsub(/*$/,""); print}'} > "${lb13file}"
awk '/^>/ {{printf("%s%s\t",(N>0?"\n":""),$0);N++;next;}} {{printf("%s",$0);}} END {{printf("\n");}}' < "${lb13file}"  |sort > "${lb14file}"
sed 's/>//' "${lb14file}" > "${lb15file}"
awk -v input_name="$input_name" '{print $1, $10, "hgcB_homologs", input_name}' "${lb15file}" > "${lb16file}"

paste "${lb9file}" "${lb16file}" > "${lb17file}"
done

cat outputs/* > detected_hgc_homologs.txt
rm -r tmp outputs
