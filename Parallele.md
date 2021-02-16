Generating .db and .IR.out.tab in parallel
============
PSI-Sigma is currently a single-thread software. The majority of the computing time was used to create IR.out.tab files (intronic read counts). To generate IR.out.tab files in parallel, here are the instructions:

To generating a universe .db in parallel, you first link all the .SJ.out.tab and .bam files in the `afolder` and run the `PSIsigma-db-v.1.0.pl` script as below. The first parameter of `PSIsigma-db-v.1.0.pl` is the name of Chromosome. `5` is the minimum number of novel junctions to construct an event in .db. The third and fourth parameter of `PSIsigma-db-v.1.0.pl` are `--type` and `--irmode`.
```
mkdir afolder
cd afolder
ln -s bamfolder/*.bam* .
ln -s bamfolder/*.SJ.* .
get ftp://ftp.ensembl.org/pub/release-87/gtf/homo_sapiens//Homo_sapiens.GRCh38.87.gtf.gz
gzip -d Homo_sapiens.GRCh38.87.gtf.gz
(grep "^#" Homo_sapiens.GRCh38.87.gtf; grep -v "^#" Homo_sapiens.GRCh38.87.gtf | sort -k1,1 -k4,4n) > Homo_sapiens.GRCh38.87.sorted.gtf
rm Homo_sapiens.GRCh38.87.gtf

#one file name per line in groupa.txt and groupb.txt
total=$(ls *.SJ.out.tab|wc -l)
counta=$(printf "%.0f" $(($total/2)))
countb=$((total-counta))
ls *.SJ.out.tab|sed 's/.SJ.out.tab//'|head -n $counta > groupa.txt
ls *.SJ.out.tab|sed 's/.SJ.out.tab//'|head -n $countb > groupb.txt

perl ~/tools/PSI-Sigma-1.9j/PSIsigma-db-v.1.0.pl ~/index/Homo_sapiens.GRCh38.100.sorted.gtf 1 5 1 1 &
perl ~/tools/PSI-Sigma-1.9j/PSIsigma-db-v.1.0.pl ~/index/Homo_sapiens.GRCh38.100.sorted.gtf 2 5 1 1 &
perl ~/tools/PSI-Sigma-1.9j/PSIsigma-db-v.1.0.pl ~/index/Homo_sapiens.GRCh38.100.sorted.gtf 3 5 1 1 &
perl ~/tools/PSI-Sigma-1.9j/PSIsigma-db-v.1.0.pl ~/index/Homo_sapiens.GRCh38.100.sorted.gtf 4 5 1 1 &
perl ~/tools/PSI-Sigma-1.9j/PSIsigma-db-v.1.0.pl ~/index/Homo_sapiens.GRCh38.100.sorted.gtf 5 5 1 1 &
perl ~/tools/PSI-Sigma-1.9j/PSIsigma-db-v.1.0.pl ~/index/Homo_sapiens.GRCh38.100.sorted.gtf 6 5 1 1 &
perl ~/tools/PSI-Sigma-1.9j/PSIsigma-db-v.1.0.pl ~/index/Homo_sapiens.GRCh38.100.sorted.gtf 7 5 1 1 &
perl ~/tools/PSI-Sigma-1.9j/PSIsigma-db-v.1.0.pl ~/index/Homo_sapiens.GRCh38.100.sorted.gtf 8 5 1 1 &
perl ~/tools/PSI-Sigma-1.9j/PSIsigma-db-v.1.0.pl ~/index/Homo_sapiens.GRCh38.100.sorted.gtf 9 5 1 1 &
perl ~/tools/PSI-Sigma-1.9j/PSIsigma-db-v.1.0.pl ~/index/Homo_sapiens.GRCh38.100.sorted.gtf 10 5 1 1 &
perl ~/tools/PSI-Sigma-1.9j/PSIsigma-db-v.1.0.pl ~/index/Homo_sapiens.GRCh38.100.sorted.gtf 11 5 1 1 &
perl ~/tools/PSI-Sigma-1.9j/PSIsigma-db-v.1.0.pl ~/index/Homo_sapiens.GRCh38.100.sorted.gtf 12 5 1 1 &
perl ~/tools/PSI-Sigma-1.9j/PSIsigma-db-v.1.0.pl ~/index/Homo_sapiens.GRCh38.100.sorted.gtf 13 5 1 1 &
perl ~/tools/PSI-Sigma-1.9j/PSIsigma-db-v.1.0.pl ~/index/Homo_sapiens.GRCh38.100.sorted.gtf 14 5 1 1 &
perl ~/tools/PSI-Sigma-1.9j/PSIsigma-db-v.1.0.pl ~/index/Homo_sapiens.GRCh38.100.sorted.gtf 15 5 1 1 &
perl ~/tools/PSI-Sigma-1.9j/PSIsigma-db-v.1.0.pl ~/index/Homo_sapiens.GRCh38.100.sorted.gtf 16 5 1 1 &
perl ~/tools/PSI-Sigma-1.9j/PSIsigma-db-v.1.0.pl ~/index/Homo_sapiens.GRCh38.100.sorted.gtf 17 5 1 1 &
perl ~/tools/PSI-Sigma-1.9j/PSIsigma-db-v.1.0.pl ~/index/Homo_sapiens.GRCh38.100.sorted.gtf 18 5 1 1 &
perl ~/tools/PSI-Sigma-1.9j/PSIsigma-db-v.1.0.pl ~/index/Homo_sapiens.GRCh38.100.sorted.gtf 19 5 1 1 &
perl ~/tools/PSI-Sigma-1.9j/PSIsigma-db-v.1.0.pl ~/index/Homo_sapiens.GRCh38.100.sorted.gtf 20 5 1 1 &
perl ~/tools/PSI-Sigma-1.9j/PSIsigma-db-v.1.0.pl ~/index/Homo_sapiens.GRCh38.100.sorted.gtf 21 5 1 1 &
perl ~/tools/PSI-Sigma-1.9j/PSIsigma-db-v.1.0.pl ~/index/Homo_sapiens.GRCh38.100.sorted.gtf 22 5 1 1 &
perl ~/tools/PSI-Sigma-1.9j/PSIsigma-db-v.1.0.pl ~/index/Homo_sapiens.GRCh38.100.sorted.gtf X 5 1 1 &
perl ~/tools/PSI-Sigma-1.9j/PSIsigma-db-v.1.0.pl ~/index/Homo_sapiens.GRCh38.100.sorted.gtf Y 5 1 1 &
perl ~/tools/PSI-Sigma-1.9j/PSIsigma-db-v.1.0.pl ~/index/Homo_sapiens.GRCh38.100.sorted.gtf MT 5 1 1 &
wait
cat *.db.tmp > PSIsigma1d9j.db
cat *.bed.tmp > PSIsigma1d9j.bed
rm *.db.tmp
rm *.bed.tmp
```
Generating IR.out.tab in parallel
============
In the same `afolder`, run `PSIsigma-ir-v.1.2.pl` for each .bam file. The parameter for `PSIsigma-ir-v.1.2.pl` is `--type`.
```
perl ~/tools/PSI-Sigma-1.9j/PSIsigma-ir-v.1.2.pl PSIsigma1d9j.db Sample1.Aligned.sortedByCoord.out.bam 1 &
perl ~/tools/PSI-Sigma-1.9j/PSIsigma-ir-v.1.2.pl PSIsigma1d9j.db Sample2.Aligned.sortedByCoord.out.bam 1 &
perl ~/tools/PSI-Sigma-1.9j/PSIsigma-ir-v.1.2.pl PSIsigma1d9j.db Sample3.Aligned.sortedByCoord.out.bam 1 &
wait
echo "All processes done!"
```


