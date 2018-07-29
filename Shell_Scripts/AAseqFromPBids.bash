#/bin/bash

cat /path/to/file/ISO_PAM/TempPBids_EXACT | xargs -I '{}' gsed -n '/\<{}\>/,/*/ p' /path/to/file/ISO_PAM/longest_orfs.pep > /path/to/file/ISO_PAM/TempAA_fromPBid.txt
