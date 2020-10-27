### THIS IS NOT A SCRIPT TO BE RUN AT ONCE
# script used in Maurin et al., Myrtales phylogenomics, submitted to the American Journal of Botany special issue Angiosperms353 (Nov. 2020)
# sh extension is there to facilitate reading in some text editors
# The commands below allow to:
# 1: calculate some simple stats on sequencing data (Illumina fastq)
# 2: use the results of HybPiper to calculate on and off target coverage and combine with the result of 1



######### 1- MAKE A FILE WITH SUMMARY NUMBERS FOR ALL SAMPLES #########################

### Print species name, total (R1+R2) number of raw reads and number of Gb (can take time to count Gb if many/big files)

for f in *R1.fastq; do (<$f awk -v name=${f/_L001_R1.fastq} '(NR%4)==0 {sum += length($1)} END {reads = NR/2 ;Gb = sum*2/1000000000 ; printf ("%s\t%d\t%.12f\n", name, reads, Gb)}' >> raw_num.txt); done

### Alternative (faster but only valid if always same read length - replace 301 by your own read length):

for f in *R1_001.fastq; do (<$f awk -v name=${f/_L001_R1_001.fastq} '{reads = NR/2 ; Gb = reads*301/1000000000 ; printf ("%s\t%d\t%.12f\n", name, reads, Gb)}' >> raw_num301.txt); done


### Print species name, total (R1+R2) number of trimmed paired reads and of Gb (can take time to count Gb if many/big files)

for f in *R1_Tpaired.fastq; do (<$f awk -v name=${f/_L001_R1_Tpaired.fastq} '(NR%4)==0 {sum += length($1)} END {reads = NR/2 ;Gb = sum*2/1000000000 ; printf ("%s\t%d\t%.12f\n", name, reads, Gb)}' >> Tpaired_num.txt); done

### Alternative only including read number, no Gb counts

for f in *R1_001_Tpaired.fastq; do (<$f awk -v name=${f/_L001_R1_001_Tpaired.fastq} '{reads = NR/2 ; printf ("%s\t%d\n", name, reads)}' >> Tpaired_num_reads.txt); done


### Print species name, total (All) number of trimmed unpaired reads and of Gb (can take time to count Gb if many/big files)
### !! Only working after you concatenated the unpaired reads to generate the *TupairedAll.fastq files!!

for f in *TunpairedAll.fastq; do (<$f awk -v name=${f/_L001_TunpairedAll.fastq} '(NR%4)==0 {sum += length($1)} END {reads = NR/4 ;Gb = sum/1000000000 ; printf ("%s\t%d\t%.12f\n", name, reads, Gb)}' >> Tunpaired_num.txt); done

### Alternative only including read number, no Gb counts

for f in *TunpairedAll.fastq; do (<$f awk -v name=${f/_L001_TunpairedAll.fastq} '{reads = NR/4 ; printf ("%s\t%d\n", name, reads)}' >> Tunpaired_num_reads.txt); done



### put all files where you want, for instance with:
mv trimmed/*num.txt .

### paste all results together by species name:
# (possibly need to check first if same names in both files using diff, if not, add or delete lines appropriately, manually)

# create a header
echo "species\tnumber of raw reads\tnumber of raw Gb\tnumber of trimmed paired reads\tnumber of trimmed paired Gb\tnumber of trimmed unpaired reads\tnumber of trimmed unpaired Gb" > header_tmp.txt
# paste data together
join -t $'\t' raw_num.txt Tpaired_num.txt >> num_summary_1tmp.txt
join -t $'\t' num_summary_1tmp.txt Tunpaired_num.txt >> num_summary_2tmp.txt
# concatenate the header and the data
cat header_tmp.txt num_summary_2tmp.txt >> num_summary.txt
#remove intermediary files
rm *tmp.txt

### The result is num_summary.txt that contains all read and Gb numbers.
# To delete the files corresponding only to raw or trimmed data, one can use:
rm *num.txt




######### 2- CALCULATE COVERAGE FOR ALL SAMPLES #########################

### ONLY WORKS ON SUPERCONTIGS AFTER A HYBPIPER ANALYSIS USING BWA...I think.
### Could use the bam file produced for each species after the bwa mapping but it distinguishes between variants of a same reference
### Instead we remap the reads corresponding to each gene on the final gene sequence obtained for the sample
### We produce one mapping per sample, with all recovered genes in one mapping


# Put all supercontigs in a file for each species separately
# In all following commands XXX is the common point to all gene names (eg. "g" if genes are names g1, g2, g3 etc.)

for f in *L001; do (echo ${f/} >> namelist.txt); done

while read name
do
for f in $name/g*/$name/sequences/intron/*_supercontig.fasta; do (cat $f >> ${name}_Allsupercontigs.fasta); done
done < namelist.txt

# cat all interleaved and unpaired IN THIS ORDER, i.e. PAIRED FIRST

while read name
do
for f in $name/g*/*interleaved.fasta; do (cat $f >> ${name}_AllMappedReads.fasta); done
done < namelist.txt

while read name
do
for f in $name/g*/*unpaired.fasta; do (cat $f >> ${name}_AllMappedReads.fasta); done
done < namelist.txt

# remove duplicated reads (there are some, meaning that some reads mapped 2 different genes): put seq on same line as id, remove duplicated reads, put back seq as normal
for f in *AllMappedReads.fasta; do (<$f awk 1 ORS='______' | sed 's/______>/\n>/g' | awk '!seen[$0]++' | sed 's/______/\n/g' > ${f/.fasta}_u.fasta ); done

# then index the reference files
parallel 'bwa index {}' ::: *_Allsupercontigs.fasta
# and map reads on them with bwa (-p tells that you have paired end, normally bwa will be able to differenciate PE and SE in a same file as long as PE are given first in the interleaved format)
parallel 'bwa mem -a -p {}_Allsupercontigs.fasta {}_AllMappedReads_u.fasta > {}_BWA.sam' :::: namelist.txt

#create a new "parsed" sam without the reads that don't match well using my py script (currently AS> 30; mismatches <3, no clipping - easy to change script to increase or decrease stringency)
parallel 'python /PATH-TO-SCRIPT/parse_sam_bwa.py {}_BWA.sam {}_BWA_parsed.sam' :::: namelist.txt


# get number of reads on target using the parsed sam (i.e. only counting reads that match well):
# combine read name and read seq of the parsed SAM (to avoid counting mates as duplicates), and remove duplicated lines
while read name
do
grep -v '@' ${name}_BWA_parsed.sam | awk 'OFS="______" {print $1, $10}' | awk '!seen[$0]++' >> ${name}_AllMappedReads_parsed_names_u.txt
done < namelist.txt

# count the number of lines and format the counts to be able to join them to the num_summary.txt statistic table
wc -l *AllMappedReads_parsed_names_u.txt | head -n -1 | awk 'OFS="" ; {sub("_L001_AllMappedReads_parsed_names_u.txt", "", $2); print $2, "\t", $1}' >> reads_on_target_tmp.txt


# join this to the data from num_summary.txt (see first section) to make a new num_summary_tmp.txt result file
mv reads_on_target_tmp.txt ../
cd ../
head -n 1 num_summary.txt | awk 'OFS="\t" {print $0, "number of reads on target"}' | cat - < <(tail -n +2 num_summary.txt | join -t $'\t' - reads_on_target_tmp.txt) > num_summary2.txt
mv num_summary2.txt num_summary.txt
# remove tmp files
rm -I *tmp.txt

# use samtools to convert the parsed sam in a sorted indexed bam (speeds up the idxstats step below, and bam necessary for the mpileup)
parallel 'samtools view -b {} -o {.}.bam' ::: *parsed.sam
# sort the bam
parallel 'samtools sort {} -o {.}_sorted.bam' ::: *.bam

# index the sorted bam files
parallel 'samtools index {}' ::: *_sorted.bam

# use samtools to get some basic stats (may not be very interesting but just in case)
parallel 'samtools idxstats {} > {.}_BasicStats.txt' ::: *_sorted.bam


# create a pileup (it will create a file with per base info of how many reads match and how well they match, allowing to calculate a mean read depth for each gene, and even separate by exons and introns)
# use -ff and to filter only the reads following a particular condition (see options)

HERE

parallel 'samtools mpileup -A --ff UNMAP -f {}_Allsupercontigs.fasta {}_BWA_parsed_sorted.bam -o {}_BWA_parsed_sorted_pileup.txt' :::: namelist.txt

# create a file with the intron positions for each gene per species:
while read name
do
for f in $name/g*/$name/intronerate.gff; do (grep '\bintron\b' $f >> ${name}_Allintrons.gff); done
done < namelist.txt
HERE

# run my python script to get read depth stats for each gene, one file per species
# for one sample:
python parse_mpileup_covStat_many_contigs_introns-exons.py Veitchia-metiti-SBL191_S17_L001_BWA_parsed_sorted_pileup.txt Veitchia-metiti-SBL191_S17_L001_Allintrons.gff Veitchia-metiti-SBL191_S17_L001_ReadDepthStats.txt
# on all samples (may take time, use screen, and parallel):
parallel 'python /PATH-TO-SCRIPT/parse_mpileup_covStat_many_contigs_introns-exons.py {}_BWA_parsed_sorted_pileup.txt {}_Allintrons.gff {}_ReadDepthStats.txt' :::: namelist.txt

HERE
# test: check that the exon mean given by awk fit with the exon mean given by python (need to test a contig without introns)
grep 'HEY977' Veitchia-metiti-SBL191_S17_L001_BWA_parsed_sorted_pileup.txt >> test.txt
awk '{sum += $4} END {mean = sum/NR ; print mean}' test.txt


# Combine all Read depth stats in one file, with the species name as first column
HEADER="Species\tContig\tMean exon read depth\tMin exon read depth\tMax exon read depth\tExon read depth std\tMean intron read depth\tMin intron read depth\tMax intron read depth\tIntron read depth std\n"
echo $HEADER > All_ReadDepthStats.txt
for f in *L001_ReadDepthStats.txt; do (awk -v species=${f/_L001_ReadDepthStats.txt} 'NR>2 {RS="\n" ; FS="\t" ; OFS="\t" ; {print species, $0 }}' $f >> All_ReadDepthStats.txt); done

# Calculate mean, min, and max mean exon RD and mean, min, and max mean intron RD for each species and output this in a temporary file Mean_ReadDepths_tmp.txt
for f in *L001_ReadDepthStats.txt; do (awk -v species=${f/_L001_ReadDepthStats.txt} 'NR == 1 {minEx = $2; minInt = $6} NR>2 {sum_ex += $2 ; sum_int += $6 ; if($2 != "na"){RecEx += 1} ; if($2 != "na" && $2 < minEx) {minEx = $2} ; if($2 != "na" && $2 > maxEx) {maxEx = $2} ; if($6 != "na"){RecIn += 1} ; if($6 != "na" && $6 < minInt) {minInt = $6} ; if($6 != "na" && $6 > maxInt) {maxInt = $6} } END {OFS="\t"; if(RecEx != 0) {meanEx = sum_ex/RecEx} else {RecEx = "na" ; meanEx = "na" ; minEx = "na" ; maxEx = "na"}; if(RecIn != 0) {meanInt = sum_int/RecIn} else {RecIn = "na" ; meanInt = "na" ; minInt = "na" ; maxInt = "na"}; print species, RecEx, meanEx, minEx, maxEx, RecIn, meanInt, minInt, maxInt }' $f >> Mean_ReadDepths_tmp.txt); done

# join this to the data from num_summary.txt to make a new num_summary2.txt result file (possibly need to check first if same names in both files using diff, if not, add or delete lines appropriately, manually)
mv Mean_ReadDepths_tmp.txt ../
cd ../
head -n 1 num_summary.txt | awk 'OFS="\t" {print $0, "number of genes with exons", "mean average exon read depth", "minimum average exon read depth", "maximum average exon read depth", "number of genes with introns", "mean average intron read depth", "minimum average intron read depth", "maximum average intron read depth"}' | cat - < <(tail -n +2 num_summary.txt | join -t $'\t' - Mean_ReadDepths_tmp.txt) > num_summary2.txt
mv num_summary2.txt num_summary.txt
# remove tmp files
rm -I *tmp.txt




######### 3- COMBINE THIS WITH OTHER STATS AND MAKE GRAPHS #########################

# You can paste the num_summary.txt table to your own table with library quality, sample type etc., to make a multivariate analysis

# You can get the number of genes retrived at more than xx% of the reference length very easily by running Matt Johnson R script to do the heatmap,
# export the sample_len and percent.len tables created by the script and use conditional sums etc in excel, or do it directly in R.

# Then you can combine this with the other stats.



