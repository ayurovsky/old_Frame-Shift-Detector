# Frame-Shift-Detector

# full pipeline attempt # 1: 
# use reads from 25-32, get biggest drop, get p-values by permutations, distributed, p-value from permutations  
time python get_statistics.py codons.txt reference_files/s_cerevisae_orf_coding_no_Mito_no_Plasmid.fasta reference_files/saccharomyces_cerevisiae_R64-1-1_20110208.fa reference_files/saccharomyces_cerevisiae_R64-1-1_20110208.gtf all_jan2014_sort.sam bla_25_to_32 > blu&
time python distribute_frameshifts.py # calls find_frameshifts 

# get statistics for all_jan2016 - with just the 28-based reads, to get YIL009C-A and YOR239W_28_starts.txt statistics
now, only try 28 length reads, no need to move them out of frame
time python get_statistics.py codons.txt reference_files/s_cerevisae_orf_coding_no_Mito_no_Plasmid.fasta reference_files/saccharomyces_cerevisiae_R64-1-1_20110208.fa reference_files/saccharomyces_cerevisiae_R64-1-1_20110208.gtf all_jan2014_sort.sam bla > blu&
python frameshift_stats_EST3_ABP140.py bla_28 > temp

# these are the original stats before messing with frame-shifted genes
5862 protein coding orfs before any checks
Read in the reference chromosomes
5547 protein coding orfs after checking for introns and other mis-matches
4742 protein coding orfs after intersect removals

Length 21 75059 reads,  1: 0.36 2: 0.41 3: 0.23
Length 22 166755 reads,  1: 0.51 2: 0.16 3: 0.33
Length 23 256091 reads,  1: 0.39 2: 0.16 3: 0.45
Length 24 412034 reads,  1: 0.37 2: 0.39 3: 0.24
Length 25 596749 reads,  1: 0.62 2: 0.16 3: 0.22
Length 26 893131 reads,  1: 0.31 2: 0.15 3: 0.54
Length 27 3165437 reads,  1: 0.26 2: 0.63 3: 0.1
Length 28 19533951 reads,  1: 0.91 2: 0.05 3: 0.04
Length 29 17991094 reads,  1: 0.36 2: 0.04 3: 0.6
Length 30 7704565 reads,  1: 0.37 2: 0.17 3: 0.47
Length 31 1537940 reads,  1: 0.16 2: 0.25 3: 0.59
Length 32 128257 reads,  1: 0.08 2: 0.65 3: 0.27
Length 33 5290 reads,  1: 0.28 2: 0.48 3: 0.24
Length 34 459 reads,  1: 0.46 2: 0.29 3: 0.25
Length 35 84 reads,  1: 0.26 2: 0.21 3: 0.52
Length 36 24 reads,  1: 0.12 2: 0.29 3: 0.58
Length 37 4 reads,  1: 0.25 2: 0.75 3: 0.0


# get statistics and get coverage lines similar to ribosome_footprints_stage1.pl 
time python get_statistics.py codons.txt reference_files/s_cerevisae_orf_coding_no_Mito_no_Plasmid.fasta reference_files/saccharomyces_cerevisiae_R64-1-1_20110208.fa reference_files/saccharomyces_cerevisiae_R64-1-1_20110208.gtf dLys_FT_original.sam bla > blu&

# get statistics for display, similar to ribosome_footprints_stage2.pl
python get_statistics_part2.py dLys_positive_allReadsOverThirty.txt

# eLife article - dLys, dHis, can also get YPD1 and YPD2 from Ying's 2013 paper, get newer data from Ingolia!!!
https://elifesciences.org/articles/03735
# justin's links for dLys, dHis footprint and mrna
https://www.ncbi.nlm.nih.gov/sra/?term=SRP044053


# start with dLys Saccharomyces cerevisia

# get clipped fastq
wget https://sra-download.ncbi.nlm.nih.gov/traces/sra22/SRR/001471/SRR1506633
sratoolkit.2.9.6-1-ubuntu64/bin/fastq-dump SRR1506633&
cat SRR1506633.fastq | fastx_toolkit_0.0.13/bin/fastx_clipper -Q33 -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA -l 20 -n -v -o  SRR_clipped.fastq&
bowtie2-2.3.4-linux-x86_64/bowtie2 --phred33 -p 8 --end-to-end -x reference_files/yeast SRR_clipped.fastq -S dLys_fromSRR.sam&


# get old bowtie2 - no need for that, new one works just dandy
wget https://sourceforge.net/projects/bowtie-bio/files/bowtie2/2.1.0/bowtie2-2.1.0-linux-x86_64.zip/download
# get old reference files from the server, based on yeast genome, NOT! ncbi
bowtie2-2.1.0/bowtie2-build reference_files/saccharomyces_cerevisiae_R64-1-1_20110208.fa s_cerevisiae_bowtie2-2.1.0_index
# build the sam file
#bowtie2-2.3.4-linux-x86_64/bowtie2 --phred33 -p 8 --end-to-end -x bowtie2-indexes/yeast dLys_clipped.fastq -S dLys.sam&
bowtie2-2.1.0/bowtie2 --phred33 -p 8 --end-to-end -x reference_files/s_cerevisiae_bowtie2-2.1.0_index dLys_small.fastq -S dLys_small.sam
bowtie2-2.3.4-linux-x86_64/bowtie2-build s_cerevisiae_R64.fna yeast

TODO HERE : try getting fastqc and getting an over-represented sequence 


# try re-building the sam file with new files
# inside reference files:
../bowtie2-2.3.4-linux-x86_64/bowtie2-build saccharomyces_cerevisiae_R64-1-1_20110208.fa yeast


# OLD NOTES, for reference only
wget ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByExp/sra/SRX/SRX476/SRX476345/SRR1177157/SRR1177157.sra&

/home/alisa/projects/ngs/neiman_ribosome_profiling/sratoolkit.2.3.2-4-ubuntu64/bin/fastq-dump /home/alisa/projects/ribosome_profiling_bias/SRR1177157.sra&

cat SRR1177157.fastq | /home/programs/fastx/fastx_clipper -Q33 -a CTGTAGGCACCATCAATAGATCGGAAGAGCACACGTCTGAA -l 20 -n -v -o /home/alisa/projects/ribosome_profiling_bias/SRR1177157_clipped.fastq&

/home/programs/bowtie2/bowtie2-2.1.0/bowtie2 --phred33 -p 8 --end-to-end -x /home/alisa/projects/ngs/indexes/bowtie2_2.1.0/s_cerevisiae_bowtie2-2.1.0_index SRR1177157_clipped.fastq -S SRR1177157.sam&

perl end_bias_quantification_paper_cerevisiae_v2.pl -n s_cerevisae_orf_coding_no_Mito_no_Plasmid.fasta -c codons.txt -a saccharomyces_cerevisiae_R64-1-1_20110208.gtf -s SRR1177157.sam -r saccharomyces_cerevisiae_R64-1-1_20110208.fa --o5pEnd v2_SRR1177157_5pEnd.txt --oBef5pEnd v2_SRR1177157_before5pEnd.txt --o3pEnd v2_SRR1177157_3pEnd.txt --pAft3pEnd v2_SRR1177157_after3pEnd.txt -min_len=28 -max_len=28 > bla_rib_v2_srr&

perl end_bias_correction_paper_cerevisiae_v2.pl -n s_cerevisae_orf_coding_no_Mito_no_Plasmid.fasta -c codons.txt -a saccharomyces_cerevisiae_R64-1-1_20110208.gtf -s SRR1177157.sam -r saccharomyces_cerevisiae_R64-1-1_20110208.fa -o v2_SRR1177157_self_corrected_stage1.csv --o5pEnd v2_SRR1177157_5pEnd.txt --obefore5pEnd v2_SRR1177157_before5pEnd.txt --o3pEnd v2_SRR1177157_3pEnd.txt --oafter3pEnd v2_SRR1177157_after3pEnd.txt&

perl ribosome_footprints_stage1_using_corrections_paper_cerevisiae_v2.pl -n s_cerevisae_orf_coding_no_Mito_no_Plasmid.fasta -c codons.txt -a saccharomyces_cerevisiae_R64-1-1_20110208.gtf -r saccharomyces_cerevisiae_R64-1-1_20110208.fa -s v2_SRR1177157_self_corrected_stage1.csv -o v2_SRR1177157_self_corrected_stage1_canonical.csv&

perl ribosome_footprints_stage2.pl -e 9 -c codons.txt -s v2_SRR1177157_self_corrected_stage1_canonical.csv -o window_stage2_SRR1177157_self_corrected_v2.csv&
