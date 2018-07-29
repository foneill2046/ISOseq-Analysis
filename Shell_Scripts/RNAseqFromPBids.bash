#/bin/bash

cat /path/to/file/ISO_PAM/TempPBids_EXACT | xargs -I '{}' gsed -n '/\<{}\>/,/>/{;/\<{}\>/p;/>/d;p;}' /path/to/file/ISO_PAM/annotated.dontScreenMatchAnnot.dontScreenProteinCoding.fa > /path/to/file/ISO_PAM/TempRNA_fromPBid.txt
