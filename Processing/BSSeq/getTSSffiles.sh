

grep -f Conv.txt ../TSS_1kb.sh > PromConv_1kb.txt
grep -f Conv.txt ../TSS_5kb.sh > PromConv_5kb.txt


grep -f TSC.txt ../TSS_1kb.sh > TSCConv_1kb.txt
grep -f TSC.txt ../TSS_5kb.sh > TSCConv_5kb.txt

grep -f Niave.txt ../TSS_1kb.sh > NiaveConv_1kb.txt
grep -f Niave.txt ../TSS_5kb.sh > NiaveConv_5kb.txt


grep -f TSS_1kb.sh Conv.txt > PromConv_1kb.txt
Promoter_5kb.bed



../gene_body.bed
awk '{ print $1 }' foo



grep -f Conv.txt ../gene_body_s.bed > ConvGB_1kb.txt
grep -f TSC.txt ../gene_body_s.bed > TSCGB_1kb.txt
grep -f Niave.txt ../gene_body_s.bed > NiaveGB_1kb.txt


#gene_body_s.bed
