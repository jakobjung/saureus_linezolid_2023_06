# reference sequences

## sRNA annotation

I blasted the sRNAs from http://srd.genouest.org/ (14-06-2023), and blasted them against my genome: 

```bash
blastn -query sRNA_N315.fasta -subject Staphylococcus_aureus_101588.fasta -max_target_seqs 1 -out blast_srnas.txt -outfmt 6
```

Then I checked for the srna lengths:

```bash
bioawk -c fastx '{ print $name, length($seq) }' < sRNA_N315.fasta > lengths.tsv
```

Then I ran the R-script to merge srna annotations.

```bash
R add_sRNAs_to_gff.R
```



### curate N315 and do proteinortho

```bash
grep "ORFID:" SAN315.gff3 | sed -E "s/(.*\t)CDS(\t.*ORFID:([^;%]+).*)/\\1\\3\\2/" | head -n 100 > SAN315_curated_for_portho.gff

bedtools getfasta -nameOnly -s -fo SAN315_curated_for_portho.fasta -fi genome_N315.fasta -bed SAN315_curated_for_portho.gff
```

run protho:

```bash
proteinortho6.pl --project=portho ./SAN315_curated_for_portho.portho.fasta ./SA101588.ffn --p=blastn 
```



## combine sRNAs:

```
grep -P "\tsRNA" SA101588_prokka.gff  >> sRNAs_combined.gff 
grep -P "\tsRNA" WG_with_sRNAs_Staphylococcus_aureus_101588.gff3  >> sRNAs_combined.gff
```

