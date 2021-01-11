# Frame-Shift-Detector
# jupyter nbconvert PeriodicityStats.ipynb --to pdf

# a small change
#############################################
# now, try to combine only the >90% in f1 starts for all datasets - no hanson, no geraschenko, no williams
all_jan2014_bowtie_28_starts
Wu2019_bowtie_unique_28_29_starts
Guydosh2014_bowtie_unique_28_29_starts
Young2015_bowtie_unique_28_starts

# merge the four files
python merge_starts.py


python get_mean_and_std_f1.py output/combined_starts_four_sets_above_90
F1 mean with statistics is 0.8763924094763768
F1 stdev with statistics is 0.22906614440659637
Total f1 is 0.914418446440859
Total f2 is 0.04545015288643027
Total is: 48066725

TODO: which one to use? 0.876 or 0.91444???? probably the 0.91....

python find_frameshifts_approach_5.py output/combined_starts_four_sets_above_90& 

# genes to check YGL105W YGR239C YGL106W YOR123C YNL286W YLR372W YJL057C YGR007W YDL126C YIL161W 

# re-do the simulation on these
# do for shifts of 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 1.0, 2.0, 3.0 
time python simulate_shifts_at_percent.py output/combined_starts_four_sets_above_90 0.05&
python find_frameshifts_approach_5.py output/combined_starts_four_sets_above_90_with_simulated_shifts_at_0.05 &

# First Run ResultsForPaper.ipynb

For the Before gene starts - redo...
time python get_p_value_for_before_gene_starts.py codons.txt reference_files/s_cerevisae_orf_coding_no_Mito_no_Plasmid.fasta reference_files/saccharomyces_cerevisiae_R64-1-1_20110208.fa reference_files/saccharomyces_cerevisiae_R64-1-1_20110208.gtf BeforeGene_FoundMay11.csv BeforeGene_FoundMay11_with_pattern_pvalue.csv BeforeGene_FoundMay11_with_all_pvalues.csv > bla

# for plus one - get p values from stops, then from speeds
time python get_p_value_for_plus_one_gene_starts.py codons.txt reference_files/s_cerevisae_orf_coding_no_Mito_no_Plasmid.fasta reference_files/saccharomyces_cerevisiae_R64-1-1_20110208.fa reference_files/saccharomyces_cerevisiae_R64-1-1_20110208.gtf PlusOneGene_FoundMay15.csv PlusOneGene_FoundMay15_with_pattern_pvalue.csv PlusOneGene_FoundMay15_all_pattern_pvalues.csv > bla

time python get_p_value_for_plus_one_gene_starts_with_codon_speeds.py codons.txt reference_files/s_cerevisae_orf_coding_no_Mito_no_Plasmid.fasta reference_files/saccharomyces_cerevisiae_R64-1-1_20110208.fa reference_files/saccharomyces_cerevisiae_R64-1-1_20110208.gtf PlusOneGene_FoundMay15.csv codon_speed_PlusOneGene_FoundMay15_with_pattern_pvalue.csv codon_speed_PlusOneGene_FoundMay15_all_pattern_pvalues.csv > bla


# then run ResultsForPaper.ipynb again to display tables

# Re-Run with June24 dates to get extra information for PlusOne shifts - may not need it for final tables
time python get_p_value_for_plus_one_gene_starts.py codons.txt reference_files/s_cerevisae_orf_coding_no_Mito_no_Plasmid.fasta reference_files/saccharomyces_cerevisiae_R64-1-1_20110208.fa reference_files/saccharomyces_cerevisiae_R64-1-1_20110208.gtf PlusOneGene_FoundJune24.csv PlusOneGene_FoundJune24_with_pattern_pvalue.csv PlusOneGene_FoundJune24_all_pattern_pvalues.csv > bla

time python get_p_value_for_plus_one_gene_starts_with_codon_speeds.py codons.txt reference_files/s_cerevisae_orf_coding_no_Mito_no_Plasmid.fasta reference_files/saccharomyces_cerevisiae_R64-1-1_20110208.fa reference_files/saccharomyces_cerevisiae_R64-1-1_20110208.gtf PlusOneGene_FoundJune24.csv codon_speed_PlusOneGene_FoundJune24_with_pattern_pvalue.csv codon_speed_PlusOneGene_FoundJune24_all_pattern_speed.csv > bla

# for the Minus direction
time python get_p_value_for_minus_one_gene_starts.py codons.txt reference_files/s_cerevisae_orf_coding_no_Mito_no_Plasmid.fasta reference_files/saccharomyces_cerevisiae_R64-1-1_20110208.fa reference_files/saccharomyces_cerevisiae_R64-1-1_20110208.gtf MinusOneGene_FoundJune30.csv MinusOneGene_FoundJune30_with_pattern_pvalue.csv MinusOneGene_FoundJune30_all_pattern_dists.csv > bla

#############################################
# to combine
Hanson2018_bowtie_unique_25_32_starts
Gerashchenko2014_bowtie_unique_25_32_starts
Guydosh2014_bowtie_unique_25_32_starts
Williams2014_bowtie_unique_25_32_starts
Young2015_bowtie_unique_25_32_starts
Wu2019_bowtie_unique_25_32_starts
all_jan2014_bowtie_25_to_32_starts


# run for two kinds of outputs
python merge_starts.py

# f1 = 0.761
python get_mean_and_std_f1.py output/combined_starts_no_wu_no_jan
python find_frameshifts_approach_5.py output/combined_starts_no_wu_no_jan&

# f1 = 0.759
python get_mean_and_std_f1.py output/combined_starts_all_seven
python find_frameshifts_approach_5.py output/combined_starts_all_seven& # moved output to Approach5_Real_Frameshifts_Info_Jan14.csv

#outFile = open("output/combined_starts_no_wu_no_jan", "w")
outFile = open("output/combined_starts_all_seven", "w")
################################
# Shiber2018 - confusing, maybe do later

################################
# Wu 2019 - excellent dataset, and showing concordance with ours

wget https://sra-download.ncbi.nlm.nih.gov/traces/sra61/SRR/007072/SRR7241903
wget https://sra-download.ncbi.nlm.nih.gov/traces/sra52/SRR/007072/SRR7241904 
sratoolkit.2.9.6-1-ubuntu64/bin/fastq-dump SRR7241903&
sratoolkit.2.9.6-1-ubuntu64/bin/fastq-dump SRR7241904&
cat SRR72419*.fastq > Wu.fastq&

cat Wu.fastq | fastx_toolkit_0.0.13/bin/fastx_clipper -Q33 -l 20 -n -v -o Wu2019.fastq -a CTGTAGGCACCATCAAT&

bowtie2-2.1.0/bowtie2 --phred33 -p 8 --end-to-end -x reference_files/s_cerevisiae_bowtie2-2.1.0_index Wu2019.fastq -S Wu2019_bowtie.sam&

Length 25 359930 reads,  1: 0.82 2: 0.09 3: 0.1
Length 26 508683 reads,  1: 0.19 2: 0.21 3: 0.61
Length 27 2092667 reads,  1: 0.24 2: 0.74 3: 0.03
Length 28 10647372 reads,  1: 0.93 2: 0.02 3: 0.05
Length 29 10397157 reads,  1: 0.08 2: 0.02 3: 0.91
Length 30 3175640 reads,  1: 0.03 2: 0.71 3: 0.26
Length 31 393889 reads,  1: 0.21 2: 0.56 3: 0.23
Length 32 47626 reads,  1: 0.17 2: 0.58 3: 0.25
Length 33 7493 reads,  1: 0.21 2: 0.62 3: 0.17
time python get_statistics_unique.py codons.txt reference_files/s_cerevisae_orf_coding_no_Mito_no_Plasmid.fasta reference_files/saccharomyces_cerevisiae_R64-1-1_20110208.fa reference_files/saccharomyces_cerevisiae_R64-1-1_20110208.gtf Wu2019_bowtie.sam output/Wu2019_bowtie_unique_25_32_starts > bla&

# 0.829 f1
python get_mean_and_std_f1.py output/Wu2019_bowtie_unique_25_32_starts&

# find frameshifts
python find_frameshifts_approach_5.py output/Wu2019_bowtie_unique_25_32_starts&

# now try only very high accuracy reads - 28 and 29
time python get_statistics_unique.py codons.txt reference_files/s_cerevisae_orf_coding_no_Mito_no_Plasmid.fasta reference_files/saccharomyces_cerevisiae_R64-1-1_20110208.fa reference_files/saccharomyces_cerevisiae_R64-1-1_20110208.gtf Wu2019_bowtie.sam output/Wu2019_bowtie_unique_28_29_starts > bla&
#now all reads are 28 length, at 92%
# find frameshifts
python find_frameshifts_approach_5.py output/Wu2019_bowtie_unique_28_29_starts&



################################
# Weinberg2013                  - some problem with this sample, most reads did not align, really nothing to work with 

wget https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos1/sra-pub-run-5/SRR1049521/SRR1049521.2& 
sratoolkit.2.9.6-1-ubuntu64/bin/fastq-dump SRR1049521.2

cat SRR1049521.2.fastq | fastx_toolkit_0.0.13/bin/fastx_clipper -Q33 -l 20 -n -v -o Weinberg2013.fastq -a CTGTAGGCACCATCAAT&
cat Weinberg2013.fastq | fastx_toolkit_0.0.13/bin/fastx_clipper -Q33 -l 20 -n -v -o Weinberg2013_clipped.fastq -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA&

bowtie2-2.1.0/bowtie2 --phred33 -p 8 --end-to-end -x reference_files/s_cerevisiae_bowtie2-2.1.0_index Weinberg2013_clipped.fastq -S Weinberg2013_bowtie.sam&

time python get_statistics_unique.py codons.txt reference_files/s_cerevisae_orf_coding_no_Mito_no_Plasmid.fasta reference_files/saccharomyces_cerevisiae_R64-1-1_20110208.fa reference_files/saccharomyces_cerevisiae_R64-1-1_20110208.gtf Weinberg2013_bowtie.sam output/Weinberg2013_bowtie_unique_25_32_starts > bla&

################################
# Young2015 - very few found frameshifts... 
wget https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos1/sra-pub-run-5/SRR2046309/SRR2046309.1&
wget https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos1/sra-pub-run-5/SRR2046310/SRR2046310.1&
sratoolkit.2.9.6-1-ubuntu64/bin/fastq-dump SRR2046309.1&
sratoolkit.2.9.6-1-ubuntu64/bin/fastq-dump SRR2046310.1&
cat SRR20463*.fastq > Young.fastq&

cat Young.fastq | fastx_toolkit_0.0.13/bin/fastx_clipper -Q33 -l 20 -n -v -o Young2015.fastq -a CTGTAGGCACCATCAAT&
cat Young2015.fastq | fastx_toolkit_0.0.13/bin/fastx_clipper -Q33 -l 20 -n -v -o Young2015_clipped.fastq -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA&

bowtie2-2.1.0/bowtie2 --phred33 -p 8 --end-to-end -x reference_files/s_cerevisiae_bowtie2-2.1.0_index Young2015_clipped.fastq -S Young2015_bowtie.sam&

time python get_statistics_unique.py codons.txt reference_files/s_cerevisae_orf_coding_no_Mito_no_Plasmid.fasta reference_files/saccharomyces_cerevisiae_R64-1-1_20110208.fa reference_files/saccharomyces_cerevisiae_R64-1-1_20110208.gtf Young2015_bowtie.sam output/Young2015_bowtie_unique_25_32_starts > bla&

# 0.788 f1!
python get_mean_and_std_f1.py output/Young2015_bowtie_unique_25_32_starts&


# find frameshifts
python find_frameshifts_approach_5.py output/Young2015_bowtie_unique_25_32_starts&


# do only for very high quality reads - 28 lengths reads only
time python get_statistics_unique.py codons.txt reference_files/s_cerevisae_orf_coding_no_Mito_no_Plasmid.fasta reference_files/saccharomyces_cerevisiae_R64-1-1_20110208.fa reference_files/saccharomyces_cerevisiae_R64-1-1_20110208.gtf Young2015_bowtie.sam output/Young2015_bowtie_unique_28_starts > bla_young&
Length 25 481171 reads,  1: 0.68 2: 0.12 3: 0.2
Length 26 467603 reads,  1: 0.25 2: 0.14 3: 0.61
Length 27 955572 reads,  1: 0.21 2: 0.69 3: 0.09
Length 28 4432507 reads,  1: 0.9 2: 0.04 3: 0.06
Length 29 4785050 reads,  1: 0.13 2: 0.03 3: 0.84
Length 30 2426998 reads,  1: 0.08 2: 0.58 3: 0.33
Length 31 785346 reads,  1: 0.15 2: 0.48 3: 0.37
Length 32 248335 reads,  1: 0.11 2: 0.66 3: 0.23

time python get_statistics_unique.py codons.txt reference_files/s_cerevisae_orf_coding_no_Mito_no_Plasmid.fasta reference_files/saccharomyces_cerevisiae_R64-1-1_20110208.fa reference_files/saccharomyces_cerevisiae_R64-1-1_20110208.gtf Young2015_bowtie.sam output/Young2015_bowtie_unique_28_starts > bla_young&


################################
# get Nedialkova2015 to run - do not include!!! no clear consensus in any of read classes, including 28-len reads;  too few 28-length reads; peak at 31...
wget https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos1/sra-pub-run-5/SRR1944912/SRR1944912.1&
wget https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos1/sra-pub-run-5/SRR1944913/SRR1944913.1&
wget https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos1/sra-pub-run-5/SRR1944914/SRR1944914.1&
wget https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos1/sra-pub-run-5/SRR1944963/SRR1944963.1&
wget https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos1/sra-pub-run-5/SRR1944964/SRR1944964.1&
wget https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos1/sra-pub-run-5/SRR1944965/SRR1944965.1&
sratoolkit.2.9.6-1-ubuntu64/bin/fastq-dump SRR1944912.1&
sratoolkit.2.9.6-1-ubuntu64/bin/fastq-dump SRR1944913.1& 
sratoolkit.2.9.6-1-ubuntu64/bin/fastq-dump SRR1944914.1&
sratoolkit.2.9.6-1-ubuntu64/bin/fastq-dump SRR1944963.1&
sratoolkit.2.9.6-1-ubuntu64/bin/fastq-dump SRR1944964.1&
sratoolkit.2.9.6-1-ubuntu64/bin/fastq-dump SRR1944965.1&
cat SRR19449*.fastq > Nedialkova.fastq&
cat Nedialkova.fastq | fastx_toolkit_0.0.13/bin/fastx_clipper -Q33 -l 20 -n -v -o Nedialkova2015.fastq -a CTGTAGGCACCATCAAT&
cat Nedialkova2015.fastq | fastx_toolkit_0.0.13/bin/fastx_clipper -Q33 -l 20 -n -v -o Nedialkova2015_clipped.fastq -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA& 

bowtie2-2.1.0/bowtie2 --phred33 -p 8 --end-to-end -x reference_files/s_cerevisiae_bowtie2-2.1.0_index Nedialkova2015_clipped.fastq -S Nedialkova2015_bowtie.sam&

time python get_statistics_unique.py codons.txt reference_files/s_cerevisae_orf_coding_no_Mito_no_Plasmid.fasta reference_files/saccharomyces_cerevisiae_R64-1-1_20110208.fa reference_files/saccharomyces_cerevisiae_R64-1-1_20110208.gtf Nedialkova2015_bowtie.sam output/Nedialkova2015_bowtie_unique_25_32_starts > bla&


################################
# get Williams2014 -            maybe include in analysis... 
wget https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos1/sra-pub-run-2/SRR1562873/SRR1562873.1&
wget https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos1/sra-pub-run-1/SRR1562874/SRR1562874.1&
wget https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos1/sra-pub-run-2/SRR1562875/SRR1562875.1&
wget https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos1/sra-pub-run-2/SRR1562876/SRR1562876.1&
wget https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos1/sra-pub-run-1/SRR1562877/SRR1562877.1&
wget https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos1/sra-pub-run-2/SRR1562878/SRR1562878.1& 
sratoolkit.2.9.6-1-ubuntu64/bin/fastq-dump SRR1562873.1&
sratoolkit.2.9.6-1-ubuntu64/bin/fastq-dump SRR1562874.1&
sratoolkit.2.9.6-1-ubuntu64/bin/fastq-dump SRR1562875.1&
sratoolkit.2.9.6-1-ubuntu64/bin/fastq-dump SRR1562876.1&
sratoolkit.2.9.6-1-ubuntu64/bin/fastq-dump SRR1562877.1&
sratoolkit.2.9.6-1-ubuntu64/bin/fastq-dump SRR1562878.1&
cat SRR15628*.fastq > Williams.fastq 
cat Williams.fastq | fastx_toolkit_0.0.13/bin/fastx_clipper -Q33 -l 20 -n -v -o Williams2014.fastq -a CTGTAGGCACCATCAAT&
cat Williams2014.fastq | fastx_toolkit_0.0.13/bin/fastx_clipper -Q33 -l 20 -n -v -o Williams2014_clipped.fastq -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA&

bowtie2-2.1.0/bowtie2 --phred33 -p 8 --end-to-end -x reference_files/s_cerevisiae_bowtie2-2.1.0_index Williams2014_clipped.fastq -S Williams2014_bowtie.sam& 

time python get_statistics_unique.py codons.txt reference_files/s_cerevisae_orf_coding_no_Mito_no_Plasmid.fasta reference_files/saccharomyces_cerevisiae_R64-1-1_20110208.fa reference_files/saccharomyces_cerevisiae_R64-1-1_20110208.gtf Williams2014_bowtie.sam output/Williams2014_bowtie_unique_25_32_starts > bla&

# f1 mean is 0.762
python get_mean_and_std_f1.py output/Williams2014_bowtie_unique_25_32_starts&

# find frameshifts
python find_frameshifts_approach_5.py output/Williams2014_bowtie_unique_25_32_starts&

# do only for high confidence reads - do NOT!!! include williams in high quality reads
time python get_statistics_unique_williams.py codons.txt reference_files/s_cerevisae_orf_coding_no_Mito_no_Plasmid.fasta reference_files/saccharomyces_cerevisiae_R64-1-1_20110208.fa reference_files/saccharomyces_cerevisiae_R64-1-1_20110208.gtf Williams2014_bowtie.sam output/Williams2014_bowtie_unique_28_starts > bla_williams&
Length 25 773515 reads,  1: 0.5 2: 0.14 3: 0.36
Length 26 780876 reads,  1: 0.43 2: 0.14 3: 0.42
Length 27 1433672 reads,  1: 0.44 2: 0.35 3: 0.21
Length 28 7016903 reads,  1: 0.84 2: 0.06 3: 0.09
Length 29 9802900 reads,  1: 0.41 2: 0.04 3: 0.55
Length 30 6152497 reads,  1: 0.41 2: 0.08 3: 0.51
Length 31 2350177 reads,  1: 0.22 2: 0.08 3: 0.7
Length 32 354981 reads,  1: 0.07 2: 0.22 3: 0.72



################################
# get Lareau2014 two cylohexamite replicates (from Hansen supplement) -DO NOT INCLUDE - too few 28-length reads
wget https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos1/sra-pub-run-5/SRR1363415/SRR1363415.1
wget https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos1/sra-pub-run-5/SRR1363416/SRR1363416.1
sratoolkit.2.9.6-1-ubuntu64/bin/fastq-dump SRR1363415.1&
cat SRR1363415.1.fastq | fastx_toolkit_0.0.13/bin/fastx_clipper -Q33 -l 20 -n -v -o Lareau2014_1_a.fastq -a CTGTAGGCACCATCAAT&
cat Lareau2014_1_a.fastq | fastx_toolkit_0.0.13/bin/fastx_clipper -Q33 -l 20 -n -v -o Lareau2014_clipped.fastq -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA&

bowtie2-2.1.0/bowtie2 --phred33 -p 8 --end-to-end -x reference_files/s_cerevisiae_bowtie2-2.1.0_index Lareau2014_clipped.fastq -S Lareau2014_bowtie.sam&

time python get_statistics_unique.py codons.txt reference_files/s_cerevisae_orf_coding_no_Mito_no_Plasmid.fasta reference_files/saccharomyces_cerevisiae_R64-1-1_20110208.fa reference_files/saccharomyces_cerevisiae_R64-1-1_20110208.gtf Lareau2014_bowtie.sam output/Lareau2014_bowtie_unique_25_32_starts > bla&

################################
################################
# get Guydosh2104 sets, all correct ribosome conformation, no strange conditions, cyclohexamite added - other strains and additions looked suspect as far as pausing/speeding up translation, but they have a ton of ribosome profiling datasets (from Hansen supplement)
wget https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos1/sra-pub-run-5/SRR1042853/SRR1042853.1
sratoolkit.2.9.6-1-ubuntu64/bin/fastq-dump SRR1042853.1
cat SRR1042853.1.fastq | fastx_toolkit_0.0.13/bin/fastx_clipper -Q33 -l 20 -n -v -o Guydosh2014_clipped.fastq -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA
cat Guydosh2014_clipped.fastq | fastx_toolkit_0.0.13/bin/fastx_clipper -Q33 -l 20 -n -v -o Guydosh2014_clipped2.fastq -a CTGTAGGCACCATCAAT&

# 11 million unique reads
bowtie2-2.1.0/bowtie2 --phred33 -p 8 --end-to-end -x reference_files/s_cerevisiae_bowtie2-2.1.0_index Guydosh2014_clipped2.fastq -S Guydosh2014_bowtie.sam&

time python get_statistics_unique.py codons.txt reference_files/s_cerevisae_orf_coding_no_Mito_no_Plasmid.fasta reference_files/saccharomyces_cerevisiae_R64-1-1_20110208.fa reference_files/saccharomyces_cerevisiae_R64-1-1_20110208.gtf Guydosh2014_bowtie.sam output/Guydosh2014_bowtie_unique_25_32_starts > bla&

Length 25 186594 reads,  1: 0.72 2: 0.1 3: 0.18
Length 26 359395 reads,  1: 0.09 2: 0.24 3: 0.67
Length 27 1179094 reads,  1: 0.16 2: 0.83 3: 0.02
Length 28 2754334 reads,  1: 0.9 2: 0.02 3: 0.08
Length 29 2075073 reads,  1: 0.03 2: 0.01 3: 0.96
Length 30 248794 reads,  1: 0.02 2: 0.75 3: 0.23
Length 31 10903 reads,  1: 0.39 2: 0.45 3: 0.16
# f1 mean is 0.852!!!! - maybe we should try getting more sets???
python get_mean_and_std_f1.py output/Guydosh2014_bowtie_unique_25_32_starts&

# find frameshifts
python find_frameshifts_approach_5.py output/Guydosh2014_bowtie_unique_25_32_starts&

# redo, take only 28 and 29 footprints
time python get_statistics_unique.py codons.txt reference_files/s_cerevisae_orf_coding_no_Mito_no_Plasmid.fasta reference_files/saccharomyces_cerevisiae_R64-1-1_20110208.fa reference_files/saccharomyces_cerevisiae_R64-1-1_20110208.gtf Guydosh2014_bowtie.sam output/Guydosh2014_bowtie_unique_28_29_starts > bla&
# now f1 mean is 92%
python find_frameshifts_approach_5.py output/Guydosh2014_bowtie_unique_28_29_starts&


################################
get the Gerashchenko 2014 dataset, control riboseq (included in Hansen2018 supplement)
wget https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos1/sra-pub-run-5/SRR1520311/SRR1520311.1& 
sratoolkit.2.9.6-1-ubuntu64/bin/fastq-dump SRR1520311.1&
cat SRR1520311.1.fastq | fastx_toolkit_0.0.13/bin/fastx_clipper -Q33 -l 20 -n -v -o Gerashchenko2014_a.fastq -a CTGTAGGCACCATCAAT
cat Gerashchenko2014_a.fastq | fastx_toolkit_0.0.13/bin/fastx_clipper -Q33 -l 20 -n -v -o Gerashchenko2014_b.fastq -a CAAGCAGAAGACGGCATACGA&
cat Gerashchenko2014_b.fastq | fastx_toolkit_0.0.13/bin/fastx_clipper -Q33 -l 20 -n -v -o Gerashchenko2014_c.fastq -a AATGATACGGCGACCACCGA&
bowtie2-2.1.0/bowtie2 --phred33 -p 8 --end-to-end -x reference_files/s_cerevisiae_bowtie2-2.1.0_index Gerashchenko2014_c.fastq -S Gerashchenko2014_bowtie.sam&

#15 million unique reads
time python get_statistics_unique.py codons.txt reference_files/s_cerevisae_orf_coding_no_Mito_no_Plasmid.fasta reference_files/saccharomyces_cerevisiae_R64-1-1_20110208.fa reference_files/saccharomyces_cerevisiae_R64-1-1_20110208.gtf Gerashchenko2014_bowtie.sam output/Gerashchenko2014_bowtie_unique_25_32_starts > bla&

# get the f1 mean is 0.701
python get_mean_and_std_f1.py output/Gerashchenko2014_bowtie_unique_25_32_starts

# find frameshifts
python find_frameshifts_approach_5.py output/Gerashchenko2014_bowtie_unique_25_32_starts&

# do only for high-quality reads - do NOT!!!! include Gerashenko in high quality set!!!
time python get_statistics_unique_gerashchenko.py codons.txt reference_files/s_cerevisae_orf_coding_no_Mito_no_Plasmid.fasta reference_files/saccharomyces_cerevisiae_R64-1-1_20110208.fa reference_files/saccharomyces_cerevisiae_R64-1-1_20110208.gtf Gerashchenko2014_bowtie.sam output/Gerashchenko2014_bowtie_unique_28_starts > bla_geraschenko&
Length 25 215519 reads,  1: 0.61 2: 0.13 3: 0.26
Length 26 201302 reads,  1: 0.38 2: 0.17 3: 0.45
Length 27 432739 reads,  1: 0.24 2: 0.62 3: 0.14
Length 28 2478978 reads,  1: 0.86 2: 0.08 3: 0.06
Length 29 3178842 reads,  1: 0.25 2: 0.06 3: 0.7
Length 30 1811379 reads,  1: 0.17 2: 0.45 3: 0.38
Length 31 404894 reads,  1: 0.11 2: 0.43 3: 0.45
Length 32 54297 reads,  1: 0.09 2: 0.71 3: 0.2



################################
# get the dataset from Hanson 2018 paper, their own ribosome profiling dataset - conclusion - very few frameshifts found, because low read counts..., but proportion is good - combine??
wget https://sra-download.ncbi.nlm.nih.gov/traces/sra51/SRR/005806/SRR5945809
sratoolkit.2.9.6-1-ubuntu64/bin/fastq-dump SRR5945809&
cat SRR5945809.fastq | fastx_toolkit_0.0.13/bin/fastx_clipper -Q33 -a CTGTAGGCACCATCAAT -l 20 -n -v -o  Hanson2018_clipped.fastq&
bowtie2-2.1.0/bowtie2 --phred33 -p 8 --end-to-end -x reference_files/s_cerevisiae_bowtie2-2.1.0_index Hanson2018_clipped.fastq -S Hanson2018_bowtie.sam& 

# 11 million unique reads, no need to sort sam, will just take unique reads
time python get_statistics_unique.py codons.txt reference_files/s_cerevisae_orf_coding_no_Mito_no_Plasmid.fasta reference_files/saccharomyces_cerevisiae_R64-1-1_20110208.fa reference_files/saccharomyces_cerevisiae_R64-1-1_20110208.gtf Hanson2018_bowtie.sam output/Hanson2018_bowtie_unique_25_32_starts > bla&

# get the f1 mean is 0.702
python get_mean_and_std_f1.py output/Hanson2018_bowtie_unique_25_32_starts

# find frameshifts with above f1 mean
time python find_frameshifts_approach_5.py output/Hanson2018_bowtie_unique_25_32_starts&


# now do only for high-confidence reads  - do NOT!!!! include Hanson in high confidence reads - the proportion is way too low
time python get_statistics_unique_hanson.py codons.txt reference_files/s_cerevisae_orf_coding_no_Mito_no_Plasmid.fasta reference_files/saccharomyces_cerevisiae_R64-1-1_20110208.fa reference_files/saccharomyces_cerevisiae_R64-1-1_20110208.gtf Hanson2018_bowtie.sam output/Hanson2018_bowtie_unique_28_starts > bla_hanson&
Length 25 128129 reads,  1: 0.61 2: 0.23 3: 0.15
Length 26 152566 reads,  1: 0.38 2: 0.19 3: 0.43
Length 27 314624 reads,  1: 0.32 2: 0.53 3: 0.15
Length 28 1575510 reads,  1: 0.82 2: 0.09 3: 0.09
Length 29 1897877 reads,  1: 0.26 2: 0.07 3: 0.67
Length 30 1311167 reads,  1: 0.14 2: 0.49 3: 0.37
Length 31 446566 reads,  1: 0.15 2: 0.5 3: 0.36
Length 32 114788 reads,  1: 0.12 2: 0.62 3: 0.27

################################
# re-do with all_jan with bowtie with only unique reads at 5.5 level!!!!!

# now re-do the mapping with bowtie, to ensure that everything maps the same way again......, and there are no multi-mapping problems
bowtie2-2.1.0/bowtie2 --phred33 -p 8 --end-to-end -x reference_files/s_cerevisiae_bowtie2-2.1.0_index all_jan2014_sorted.fastq -S all_jan2014_bowtie.sam& 
samtools view -S -b all_jan2014_bowtie.sam > all_jan2014_bowtie.bam&
samtools sort -n all_jan2014_bowtie.bam -o all_jan2014_bowtie_sorted.bam&
samtools view -h all_jan2014_bowtie_sorted.bam > all_jan2014_bowtie_sorted.sam&

time python get_statistics_unique.py codons.txt reference_files/s_cerevisae_orf_coding_no_Mito_no_Plasmid.fasta reference_files/saccharomyces_cerevisiae_R64-1-1_20110208.fa reference_files/saccharomyces_cerevisiae_R64-1-1_20110208.gtf all_jan2014_bowtie_sorted.sam output/all_jan2014_bowtie_unique_25_to_32_starts > blu&
time python find_frameshifts_approach_5.py output/all_jan2014_bowtie_unique_25_to_32_starts&

# do for shifts of 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 1.0, 2.0, 3.0 
time python simulate_shifts_at_percent.py output/all_jan2014_bowtie_unique_25_to_32_starts 0.05&
# now run approach 5.7 with the newly-redone bowtie simulated datasets 
time python find_frameshifts_approach_5.py output/all_jan2014_bowtie_unique_25_to_32_starts_with_simulated_shifts_at_0.5&

# old stuff for bowtie re-mapping
time python get_statistics.py codons.txt reference_files/s_cerevisae_orf_coding_no_Mito_no_Plasmid.fasta reference_files/saccharomyces_cerevisiae_R64-1-1_20110208.fa reference_files/saccharomyces_cerevisiae_R64-1-1_20110208.gtf all_jan2014_bowtie_sorted.sam output/all_jan2014_bowtie_25_to_32_starts > blu&
original stats after moving reads:
{28: {1: 30541687, 2: 7888219, 3: 2275609}, 25: {1: 370749, 2: 97483, 3: 128731}, 30: {1: 901298, 2: 248053, 3: 389018}, 33: {1: 82777, 2: 34625, 3: 10896}}
Length 25 596963 reads,  1: 0.62 2: 0.16 3: 0.22
Length 28 40705515 reads,  1: 0.75 2: 0.19 3: 0.06
Length 30 1538369 reads,  1: 0.59 2: 0.16 3: 0.25
Length 33 128298 reads,  1: 0.65 2: 0.27 3: 0.08

# do for shifts of 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 1.0, 2.0, 3.0 
time python simulate_shifts_at_percent.py output/all_jan2014_bowtie_25_to_32_starts 0.05&
# now run approach 5.7 with the newly-redone bowtie simulated datasets 
time python find_frameshifts_approach_5.py output/all_jan2014_bowtie_25_to_32_starts_with_simulated_shifts_at_0.5&

# and real data
time python find_frameshifts_approach_5.py output/all_jan2014_bowtie_25_to_32_starts

# and now re-do for only 28 length reads
time python get_statistics_unique.py codons.txt reference_files/s_cerevisae_orf_coding_no_Mito_no_Plasmid.fasta reference_files/saccharomyces_cerevisiae_R64-1-1_20110208.fa reference_files/saccharomyces_cerevisiae_R64-1-1_20110208.gtf all_jan2014_bowtie_sorted.sam output/all_jan2014_bowtie_25_to_32_starts > blu&
Length 25 534267 reads,  1: 0.63 2: 0.16 3: 0.21
Length 26 801296 reads,  1: 0.32 2: 0.15 3: 0.53
Length 27 2895042 reads,  1: 0.27 2: 0.63 3: 0.1
Length 28 18250353 reads,  1: 0.91 2: 0.05 3: 0.04
Length 29 16545955 reads,  1: 0.37 2: 0.04 3: 0.6
Length 30 6965842 reads,  1: 0.37 2: 0.16 3: 0.47
Length 31 1341356 reads,  1: 0.16 2: 0.24 3: 0.6
Length 32 101425 reads,  1: 0.08 2: 0.64 3: 0.28
# do just 28-length reads
time python get_statistics_unique.py codons.txt reference_files/s_cerevisae_orf_coding_no_Mito_no_Plasmid.fasta reference_files/saccharomyces_cerevisiae_R64-1-1_20110208.fa reference_files/saccharomyces_cerevisiae_R64-1-1_20110208.gtf all_jan2014_bowtie_sorted.sam output/all_jan2014_bowtie_28_starts > blu&
#now all reads are 28 length at 91%
python find_frameshifts_approach_5.py output/all_jan2014_bowtie_28_starts&


################################
# explore the hisat mapping for approach 5.7 - not as good as bowtie!!!!

hisat2 --no-softclip -x reference_files/s_cerevisiae_hisat -U all_jan2014_sorted.fastq -S all_jan2014_hisat.sam&

# run the edited version of the script to get the distribution of reads before shifts - verified no bugs, similar enough to bowtie results
# uncommented the shifting code, and started again

time python get_statistics_hisat_unique.py codons.txt reference_files/s_cerevisae_orf_coding_no_Mito_no_Plasmid.fasta reference_files/saccharomyces_cerevisiae_R64-1-1_20110208.fa reference_files/saccharomyces_cerevisiae_R64-1-1_20110208.gtf all_jan2014_hisat.sam output/all_jan2014_hisat_unique_25_to_32_starts > blue2&

# do for shifts of 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 1.0, 2.0, 3.0 
time python simulate_shifts_at_percent.py output/all_jan2014_hisat_unique_25_to_32_starts 3.0&

# now run approach 5.5 with the hisat simulated datasets 
time python find_frameshifts_approach_5.py output/all_jan2014_hisat_unique_25_to_32_starts_with_simulated_shifts_at_0.5

time python get_statistics_hisat.py codons.txt reference_files/s_cerevisae_orf_coding_no_Mito_no_Plasmid.fasta reference_files/saccharomyces_cerevisiae_R64-1-1_20110208.fa reference_files/saccharomyces_cerevisiae_R64-1-1_20110208.gtf all_jan2014_hisat.sam output/all_jan2014_hisat_25_to_32_starts > blue2&
{28: {1: 30122128, 2: 7908416, 3: 2170787}, 25: {1: 367301, 2: 89545, 3: 118857}, 30: {1: 969105, 2: 266775, 3: 399596}, 33: {1: 84445, 2: 35697, 3: 10314}}
Length 25 575703 reads,  1: 0.64 2: 0.16 3: 0.21
Length 28 40201331 reads,  1: 0.75 2: 0.2 3: 0.05
Length 30 1635476 reads,  1: 0.59 2: 0.16 3: 0.24
Length 33 130456 reads,  1: 0.65 2: 0.27 3: 0.08

# do for shifts of 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 1.0, 2.0, 3.0 
time python simulate_shifts_at_percent.py output/all_jan2014_hisat_25_to_32_starts 3.0

# verify that f1 mean is 0.717 - yes!!!
python get_mean_and_std_f1.py output/all_jan2014_hisat_25_to_32_starts

# now run approach 5.7 with the hisat simulated datasets 
time python find_frameshifts_approach_5.py output/all_jan2014_hisat_25_to_32_starts_with_simulated_shifts_at_0.5

# now try on the real data
time python find_frameshifts_approach_5.py output/all_jan2014_hisat_25_to_32_starts


################################

# now try it on real data
time python find_frameshifts_approach_5.py output/all_jan2014_25_to_32_starts

# now use the approach tried in verify_frameshifts_approach_5.py to implement approach 5 -
# do for shifts of 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 1.0, 2.0, 3.0 to 8.0
time python find_frameshifts_approach_5.py output/all_jan2014_25_to_32_starts_with_simulated_shifts_at_0.5

# now get the probability bound using the binomial formula and 
time python verify_frameshifts_approach_5.py output/all_jan2014_25_to_32_starts

---------------------------------------------------------------------------------------------------------------------------
# now use the approach tried in verify_frameshifts_approach_4.py to implement approach 4 -
# do for shifts of 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 1.0, 2.0, 3.0 to 8.0
time python find_frameshifts_approach_4.py output/all_jan2014_25_to_32_starts output/all_jan2014_25_to_32_starts_with_simulated_shifts_at_0.5

# now get the probability for all the locations were we introduced the frameshift
time python verify_frameshifts_approach_4.py output/all_jan2014_25_to_32_starts& # output in output/all_jan2014_25_to_32_starts_simulated_shifts_approach_4_theoretic_distributions

# now see how we do with a combinatorial model
time python combinatorial_model.py output/all_jan2014_25_to_32_starts

---------------------------------------------------------------------------------------------------------------------------
# now try a different simulation infrastructure do not shift, but add reads at the same distribution as before
# try with approach 3 - get the same lengths window before and after

# get the starts for all genes
time python get_statistics.py codons.txt reference_files/s_cerevisae_orf_coding_no_Mito_no_Plasmid.fasta reference_files/saccharomyces_cerevisiae_R64-1-1_20110208.fa reference_files/saccharomyces_cerevisiae_R64-1-1_20110208.gtf all_jan2014_sort.sam output/all_jan2014_25_to_32_starts > blu&

# do for shifts of 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 1.0, 2.0, 3.0 to 8.0
time python simulate_shifts_at_percent.py output/all_jan2014_25_to_32_starts 3.0

# do for shifts of 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 1.0, 2.0, 3.0 to 8.0
time python find_frameshifts_approach_1.py output/all_jan2014_25_to_32_starts_with_simulated_shifts_at_0.5
time python find_frameshifts_approach_2.py output/all_jan2014_25_to_32_starts_with_simulated_shifts_at_0.5
time python find_frameshifts_approach_3.py output/all_jan2014_25_to_32_starts_with_simulated_shifts_at_0.5


----------------------------------------------------------------------------------------------------------------------------
# full pipeline, approach version 2 

# get the starts for all genes
time python get_statistics.py codons.txt reference_files/s_cerevisae_orf_coding_no_Mito_no_Plasmid.fasta reference_files/saccharomyces_cerevisiae_R64-1-1_20110208.fa reference_files/saccharomyces_cerevisiae_R64-1-1_20110208.gtf all_jan2014_sort.sam output/all_jan2014_25_to_32_starts > blu&

# add the simulated frameshifts
time python simulate_shifts.py output/all_jan2014_25_to_32_starts > junk # output/all_jan2014_25_to_32_starts_with_simulated_shifts and output/all_jan2014_25_to_32_starts_simulated_shifts_info
time python find_frameshifts_approach_2.py output/all_jan2014_25_to_32_starts_with_simulated_shifts

----------------------------------------------------------------------------------------------------------------------------
# full pipeline, approach version 1 of find_frameshifts - get the max drop - % change in frame 1 reads - for a gene, go from there

# first attempt at simulated reads - use real data and introduce shifts at various frequencies
time python simulate_shifts.py bla_25_to_32 > junk # creates bla_25_to_32_with_simulated_shifts and bla_25_to_32_simulated_shifts_info
time python find_frameshifts_approach_1.py bla_25_to_32_with_simulated_shifts
output "*_simulated_shifts_info" and "*_with_simulated_shifts" displayed in SimulatedShiftReport.ipynb


# full pipeline, approach version 1, attempt #2:
# use reads from 25-32, get biggest drop, get p-values by t-test, no permutations 
time python get_statistics.py codons.txt reference_files/s_cerevisae_orf_coding_no_Mito_no_Plasmid.fasta reference_files/saccharomyces_cerevisiae_R64-1-1_20110208.fa reference_files/saccharomyces_cerevisiae_R64-1-1_20110208.gtf all_jan2014_sort.sam bla_25_to_32 > blu&
time python find_frameshifts_approach_1.py bla_25_to_32 
output displayed in FrameShiftingReport.ipynb


# full pipeline, approach version 1, attempt # 1: 
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
