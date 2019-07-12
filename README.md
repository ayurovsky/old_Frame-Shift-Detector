# Frame-Shift-Detector


# eLife article - dLys, dHis, can also get YPD1 and YPD2 from Ying's 2013 paper, get newer data from Ingolia!!!
https://elifesciences.org/articles/03735
# justin's links for dLys, dHis footprint and mrna
https://www.ncbi.nlm.nih.gov/sra/?term=SRP044053


# start with dLys Saccharomyces cerevisiae BY4741

# get clipped fastq
wget https://sra-download.ncbi.nlm.nih.gov/traces/sra22/SRR/001471/SRR1506633
sratoolkit.2.9.6-1-ubuntu64/bin/fastq-dump SRR1506633&
mv SRR1506633.fastq dLys.fastq
cat dLys.fastq | fastx_toolkit_0.0.13/bin/fastx_clipper -Q33 -a CTGTAGGCACCATCAATAGATCGGAAGAGCACACGTCTGAA -l 20 -n -v -o dLys_clipped.fastq&
# get bowtie2 and index the genome #https://www.ncbi.nlm.nih.gov/genome/15?genome_assembly_id=22535
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/146/045/GCF_000146045.2_R64/GCF_000146045.2_R64_genomic.gff.gz
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/146/045/GCF_000146045.2_R64/GCF_000146045.2_R64_genomic.fna.gz
wget https://sourceforge.net/projects/bowtie-bio/files/bowtie2/2.3.4/bowtie2-2.3.4-linux-x86_64.zip/download
bowtie2-2.3.4-linux-x86_64/bowtie2-build s_cerevisiae_R64.fna yeast
# build the sam file
bowtie2-2.3.4-linux-x86_64/bowtie2 --phred33 -p 8 --end-to-end -x bowtie2-indexes/yeast dLys_clipped.fastq -S dLys.sam&


# OLD NOTES, for reference only
wget ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByExp/sra/SRX/SRX476/SRX476345/SRR1177157/SRR1177157.sra&

/home/alisa/projects/ngs/neiman_ribosome_profiling/sratoolkit.2.3.2-4-ubuntu64/bin/fastq-dump /home/alisa/projects/ribosome_profiling_bias/SRR1177157.sra&

cat SRR1177157.fastq | /home/programs/fastx/fastx_clipper -Q33 -a CTGTAGGCACCATCAATAGATCGGAAGAGCACACGTCTGAA -l 20 -n -v -o /home/alisa/projects/ribosome_profiling_bias/SRR1177157_clipped.fastq&

/home/programs/bowtie2/bowtie2-2.1.0/bowtie2 --phred33 -p 8 --end-to-end -x /home/alisa/projects/ngs/indexes/bowtie2_2.1.0/s_cerevisiae_bowtie2-2.1.0_index SRR1177157_clipped.fastq -S SRR1177157.sam&

perl end_bias_quantification_paper_cerevisiae_v2.pl -n s_cerevisae_orf_coding_no_Mito_no_Plasmid.fasta -c codons.txt -a saccharomyces_cerevisiae_R64-1-1_20110208.gtf -s SRR1177157.sam -r saccharomyces_cerevisiae_R64-1-1_20110208.fa --o5pEnd v2_SRR1177157_5pEnd.txt --oBef5pEnd v2_SRR1177157_before5pEnd.txt --o3pEnd v2_SRR1177157_3pEnd.txt --pAft3pEnd v2_SRR1177157_after3pEnd.txt -min_len=28 -max_len=28 > bla_rib_v2_srr&

perl end_bias_correction_paper_cerevisiae_v2.pl -n s_cerevisae_orf_coding_no_Mito_no_Plasmid.fasta -c codons.txt -a saccharomyces_cerevisiae_R64-1-1_20110208.gtf -s SRR1177157.sam -r saccharomyces_cerevisiae_R64-1-1_20110208.fa -o v2_SRR1177157_self_corrected_stage1.csv --o5pEnd v2_SRR1177157_5pEnd.txt --obefore5pEnd v2_SRR1177157_before5pEnd.txt --o3pEnd v2_SRR1177157_3pEnd.txt --oafter3pEnd v2_SRR1177157_after3pEnd.txt&

perl ribosome_footprints_stage1_using_corrections_paper_cerevisiae_v2.pl -n s_cerevisae_orf_coding_no_Mito_no_Plasmid.fasta -c codons.txt -a saccharomyces_cerevisiae_R64-1-1_20110208.gtf -r saccharomyces_cerevisiae_R64-1-1_20110208.fa -s v2_SRR1177157_self_corrected_stage1.csv -o v2_SRR1177157_self_corrected_stage1_canonical.csv&

perl ribosome_footprints_stage2.pl -e 9 -c codons.txt -s v2_SRR1177157_self_corrected_stage1_canonical.csv -o window_stage2_SRR1177157_self_corrected_v2.csv&
