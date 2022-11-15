# PREPROCESSING

# Entering the folder with the raw data
cd /data/sn/HErik/FASTQ_Generation_2022-04-02_19_04_51Z-547863316


# Renaming the files
mv *ds*/*.gz .
rm -rf *ds*


# FASTQC ANALYSIS


# Adding fastqc to the PATH
export PATH=$PATH:/data/tools/FastQC


# Creating the folder for the QC results
mkdir qc01


# Running fastqc
for f in *fastq.gz 
do 
  fastqc -t 38 $f -o qc01
done 



# Running multiqc
multiqc qc01 . 



# BASH


# Downloading the multiqc report
scp -r pmarci@boci.univet.hu:/data/sn/HErik/FASTQ_Generation_2022-04-02_19_04_51Z-547863316/multiqc_report.html /home/marci/Desktop/projects/KOKI/NEW_pipeline/Methodology_article/RESULTS/initial_qc/



# Concatenating the sample files to the folder of the final analysis
zcat 1-* > /data/sn/HErik/KOKI_220402_modszertan_12_minta/smpl01.fastq
zcat 2-* > /data/sn/HErik/KOKI_220402_modszertan_12_minta/smpl02.fastq
zcat 3-* > /data/sn/HErik/KOKI_220402_modszertan_12_minta/smpl03.fastq
zcat 4-* > /data/sn/HErik/KOKI_220402_modszertan_12_minta/smpl04.fastq
zcat 5-* > /data/sn/HErik/KOKI_220402_modszertan_12_minta/smpl05.fastq
zcat 6-* > /data/sn/HErik/KOKI_220402_modszertan_12_minta/smpl06.fastq
zcat 7-* > /data/sn/HErik/KOKI_220402_modszertan_12_minta/smpl07.fastq
zcat 8-* > /data/sn/HErik/KOKI_220402_modszertan_12_minta/smpl08.fastq
zcat 9-* > /data/sn/HErik/KOKI_220402_modszertan_12_minta/smpl09.fastq
zcat 10-* > /data/sn/HErik/KOKI_220402_modszertan_12_minta/smpl10.fastq
zcat 11-* > /data/sn/HErik/KOKI_220402_modszertan_12_minta/smpl11.fastq
zcat 12-* > /data/sn/HErik/KOKI_220402_modszertan_12_minta/smpl12.fastq

# REANALYSIS WITH TRIMMOMATIC SLIDING WINDOW THRESHOLD 15

# Creating a folder for the new results
cd /data/sn/HErik/KOKI_220402_modszertan_12_minta

mkdir new_trimming


# Quality trimming


# Entering the folder of the final analysis
cd /data/sn/HErik/KOKI_220402_modszertan_12_minta



# Creating a folder for the posttrim_qc
mkdir new_trimming/posttrim_qc



# Creating an executable for trimmomatic
TRIM='/data/tools/Trimmomatic-0.38/trimmomatic-0.38.jar'



# Iterating through the fatsq files and doing the trimming
for f in *.fastq
do 
  o='new_trimming/'${f/'.fastq'/'_trimmed.fastq'}
  log='new_trimming/posttrim_qc/'${f/'.fastq'/'_trim.out.log'}
  java -jar $TRIM SE -threads 38 \
    $f $o \
    LEADING:3 \
    TRAILING:3 \
    SLIDINGWINDOW:4:15 \
    MINLEN:36 2> $log
done



# Checking the results
cd new_trimming


ls


# Adapter trimming of the files - Cutadapt


# Entering the folder with the fastq files
cd /data/sn/HErik/KOKI_220402_modszertan_12_minta/new_trimming


# Adapter trimming
for f in *.fastq
do
    cutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA --minimum-length 36 -o ${f/'_trimmed.fastq'/'_trimmed_final.fastq'} $f
done



# Alignment

# Adding star to the PATH
export PATH=/data/tools/STAR-2.7.9a/bin/Linux_x86_64:$PATH
export PATH=/data/tools/subread-2.0.2-Linux-x86_64/bin:$PATH


# PREPARING THE NEWEST REFERENCE GENOME

# Entering the folder
cd /data/pmarci/GENOME_INDEXES


# Creating a folder for the genome
mkdir mouse_GRCm39_107


# Entering the folder
cd mouse_GRCm39_107


# Downloading the genome
wget http://ftp.ensembl.org/pub/release-107/fasta/mus_musculus/dna/Mus_musculus.GRCm39.dna.primary_assembly.fa.gz


gunzip Mus_musculus.GRCm39.dna.primary_assembly.fa.gz



# Next the annotation
wget http://ftp.ensembl.org/pub/release-107/gtf/mus_musculus/Mus_musculus.GRCm39.107.gtf.gz

gunzip Mus_musculus.GRCm39.107.gtf.gz



# Preparing the index
STAR --runThreadN 30 \
    --runMode genomeGenerate \
    --genomeDir /data/pmarci/GENOME_INDEXES/mouse_GRCm39_107 \
    --genomeFastaFiles Mus_musculus.GRCm39.dna.primary_assembly.fa \
    --sjdbGTFfile Mus_musculus.GRCm39.107.gtf


# Creating a variable for the index
idx='/data/pmarci/GENOME_INDEXES/mouse_GRCm39_107'


# Entering the folder with the trimmed fastq files
cd /data/sn/HErik/KOKI_220402_modszertan_12_minta/new_trimming


# FIRST RUN


# Running the tool
for f in *_trimmed_final.fastq
do
  root=${f/'_final.fastq'/''}
  STAR -- genomeDir $idx \
       -- genomeLoad LoadAndKeep \
       -- limitBAMsortRAM 10000000000 \
       -- readFilesIn $f \
       -- outFileNamePrefix 'FIRST_GRCm39_107_'$root \
       -- outReadsUnmapped Fastx \
       -- outSAMtype BAM SortedByCoordinate \
       -- outMultimapperOrder Random \
       -- runThreadN 14 
done


# Unloading shared memory

# Listing all shared memory object
ipcs -m


# Removing segments by ID
ipcrm -m 0
ipcrm -m 32769


# SECOND RUN


# Creating a folder for the final bam files
mkdir final_bams



# Running the tool
for f in *_trimmed_final.fastq
do
  root=${f/'.fastq'/''}
  STAR -- genomeDir $idx \
       -- sjdbFileChrStartEnd *SJ.out.tab \
       -- readFilesIn $f \
       -- outFileNamePrefix 'final_bams/GRCm39_107_'$root \
       -- outReadsUnmapped Fastx \
       -- outSAMtype BAM SortedByCoordinate \
       -- outMultimapperOrder Random \
       -- runThreadN 14 
done

# FEATURECOUNTS


# Entering the folder with the final bams
cd /data/sn/HErik/KOKI_220402_modszertan_12_minta/new_trimming/final_bams


   
# FEATURECOUNTS: without multimappers, with multioverlappers (all counted each times)
# , reads on - strand are considered

# Creating a folder for the results
mkdir feature_counts_2_strand_no_multimappers



# Creating the feature counts object
featureCounts -O -s 2 -T 15 \
   -a /data/pmarci/GENOME_INDEXES/mouse_GRCm39_107/Mus_musculus.GRCm39.107.gtf \
   -o feature_counts_2_strand_no_multimappers/featureCounts_GRCm39.107_s2_no_multi.txt \
   GRCm39_107*bam 


# Downloading the feature_counts
scp -r pmarci@boci.univet.hu:/data/sn/HErik/KOKI_220402_modszertan_12_minta/new_trimming/final_bams/feature_counts_2_strand_no_multimappers/featureCounts_GRCm39.107_s2_no_multi.txt /home/marci/Desktop/projects/KOKI/NEW_pipeline/Methodology_article/RESULTS/new_pipeline_test/feature_counts_runs/s2_multioverlap_no_multimapper/

  
# ALIGNMENT QC


# Entering the folder where the alignment is
cd /data/sn/HErik/KOKI_220402_modszertan_12_minta/new_trimming/final_bams


# Running mutliqc
multiqc . -n alignment_report --ignore feature_counts_1_strand_no_multimappers \
    --ignore feature_counts_2_strand_fraction_multimappers \
    --ignore feature_counts_2_strand_no_multimappers \
    --ignore feature_counts_no_strand_or_multimappers


# Downloading the report
scp -r pmarci@boci.univet.hu:/data/sn/HErik/KOKI_220402_modszertan_12_minta/new_trimming/final_bams/alignment_report.html /home/marci/Desktop/projects/KOKI/NEW_pipeline/Methodology_article/RESULTS/new_pipeline_test/alignment_qc/


# FEATURECOUNTS_QCs

# Entering the folder
cd feature_counts_2_strand_no_multimappers


# Running mutliqc
multiqc . -n feature_counts_2_strand_no_multimappers


# Downloading the report
scp -r pmarci@boci.univet.hu:/data/sn/HErik/KOKI_220402_modszertan_12_minta/new_trimming/final_bams/feature_counts_2_strand_no_multimappers/feature_counts_2_strand_no_multimappers.html /home/marci/Desktop/projects/KOKI/NEW_pipeline/Methodology_article/RESULTS/new_pipeline_test/feature_counts_runs/s2_multioverlap_no_multimapper/