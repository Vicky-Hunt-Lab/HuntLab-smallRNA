#!/bin/sh
# Copyright 2022 Vicky Hunt Lab Members
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
echo "Starting..."

# echo "Enter Genome file name: "
# read genome
# #"celgenome.fa"

# echo "First sequence file (e.g. 22nt sequences): "
# read ms
# #"Cel_M_col.fa"

# echo "Second sequence file (e.g. 27nt sequences): "
# read ps
# #"Allseq_Cel_male_LF.fa"

genome=$1
ms=$2
ps=$3
longestseq=$4

# Index the genome of interest
bowtie2-build $genome genome

# Map target sequences and passenger sequences to the genome using Bowtie2

bowtie2 -x genome -U $ms -f -a -N 0 --no-1mm-upfront --norc  > Main.sam
bowtie2 -x genome -U $ps -f -a -N 0 --no-1mm-upfront --norc  > Pass.sam

echo "bowtie2 complete"

# revcomp_rna -i $ms -o Rev_main.fa
# revcomp_rna -i $ps -o Rev_pass.fa

# grep -v ">" Rev_main.fa > ms_nohead.txt
# grep -v ">" Rev_pass.fa > ps_nohead.txt

# echo "Reverse complement files created"

# Convert sam output to bam output with samtools

for  f1 in `ls -1 *.sam`
do
  f2="$(basename $f1 .sam)"
  samtools view -S -b $f1 > $f2.bam
done
echo "sam to bam files complete"

# samtools view -h Main.bam | grep -vf ms_nohead.txt | samtools view -bS -o Main_filter.bam -
# samtools view -h Pass.bam | grep -vf ps_nohead.txt | samtools view -bS -o Pass_filter.bam -

cp Main.bam Main_filter.bam
cp Pass.bam Pass_filter.bam

mkdir Originals
mv Main.bam Originals
mv Pass.bam Originals

# Sort the bam output file with samtools
for  f1 in `ls -1 *.bam`
do
  f2="$(basename $f1 .bam)"
  samtools sort $f1 -o $f2.sort.bam
done
echo "bam files sorted"

# Index the sorted bam file using samtools
for  f1 in `ls -1 *.sort.bam`
do
  samtools index $f1
done
echo "bam file index complete"

# Convert indexed bam files to bed files using bedtools Bamtobed
for  f1 in `ls -1 *.sort.bam`
do
  f2="$(basename $f1 .sort.bam)"
  bedtools bamtobed -i $f1 > $f2.bed
done
echo "bed files created"

bedtools intersect -s -a Main_filter.bed -b Pass_filter.bed -wo > intersectOUTSS.txt

echo "same strand intersect file created"

awk '{print $1,$2,$3,$4,$7,$8,$9,$10,$13}' intersectOUTSS.txt > cleanupTEST.txt

echo "File clean up complete"

# Caculate length columns.
# Column 9 will be column 3 minus column 2. Column 10 will be column 7 minus column 6

awk '{$10 = $3-$2; print}' cleanupTEST.txt > length1.txt
awk '{$11 = $7-$6; print}' length1.txt > length2.txt

echo "Lengths calculated"

# Calculate begining and end overhangs.
# Column 11 will be column 6 minus column 2 - Overhang beginning
awk '{$12 = $2-$6; print}' length2.txt > overhang1.txt

# Column 12 will be column 7 minus column 3 - Overhang end
awk '{$13 = $3-$7; print}' overhang1.txt > overhang2.txt

echo "Overhangs calculated"


awk ' $12 < 0 { print $0 ; }' overhang2.txt > beg_less0.txt
awk ' $12 > 0 { print $0 ; }' overhang2.txt > beg_more0.txt

cat beg_less0.txt beg_more0.txt > beg_file.txt

awk ' $13 < 0 { print $0 ; }' overhang2.txt > end_less0.txt
awk ' $13 > 0 { print $0 ; }' overhang2.txt > end_more0.txt

cat end_less0.txt end_more0.txt > end_file.txt

cat beg_file.txt end_file.txt > No_blunts.txt
echo "Blunts removed"

sort -u No_blunts.txt > No_Dups.txt

echo "Duplicates removed"

awk '{print $1,$2,$3,$4}' No_Dups.txt > Cut_no_dup.txt
sort -u Cut_no_dup.txt > No_Dups_sorted.txt
awk '{print $4}' No_Dups_sorted.txt > Read_num.txt
sort -u Read_num.txt > Nodup_read.txt

echo "Read_Number, Occurences" > Read_number_out.txt
for read in `cat Nodup_read.txt`
do
awk -v num=$read '$1 == num { print $0 ; }' Read_num.txt > ReadTemp.txt
count=`wc -l ReadTemp.txt  | cut -d' ' -f1`
echo "$read,$count" >> Read_number_out.txt
rm ReadTemp.txt
done

awk '{print $1,$2,$3,$5,$6,$7,$9}' No_Dups.txt > OverlapOUT.txt

echo "Overlap, Occurrences" > OverlapLengths.txt
for num in $(seq 1 $longestseq);
do
awk -v num=$num '$7 == num { print $0 ; }' OverlapOUT.txt > OverTemp.txt
count=`wc -l OverTemp.txt  | cut -d' ' -f1`
echo "$num,$count" >> OverlapLengths.txt
done
rm OverTemp.txt

echo "Common Overlap Lengths counted"

awk '{print $3}' OverlapOUT.txt > end.txt
sort -u end.txt > unique.txt

echo "Sequence1End, matches" > Distribution.txt
for coord in `cat unique.txt`
do
awk -v num=$coord '$3 == num { print $0 ; }' OverlapOUT.txt > MatchTemp.txt
count=`wc -l MatchTemp.txt  | cut -d' ' -f1`
echo "$coord,$count" >> Distribution.txt
rm MatchTemp.txt
done

tar -cf Ouput_zip.tar unique.txt genome* Main* Pass* beg* end* length* overhang* cleanupTEST.txt No* ps* Rev* ms*
rm  unique.txt genome* Main* Pass* beg* end* length* overhang* cleanup.txt No* ps* Rev* ms*

echo "Distribution of sequences counted"
echo "Complete"