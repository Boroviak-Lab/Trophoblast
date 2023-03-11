./cpg_lh /mnt/scratch/gurdon/cap76/Thorsten/Marmoset/all/cjacchus.fa \
 | awk '{$2 = $2 - 1; width = $3 - $2;
   printf("%s\t%d\t%s\t%s %s\t%s\t%s\t%0.0f\t%0.1f\t%s\t%s\n",
    $1, $2, $3, $5, $6, width, $6, width*$7*0.01, 100.0*2*$6/width, $7, $9);}' \
     | sort -k1,1 -k2,2n > CpG.bed


