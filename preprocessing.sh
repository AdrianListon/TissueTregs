#!/usr/bin/env bash


###configuration
#setting directories for data and software
wrkspc="/home/tareens/Workspace"

#locaton of the various software
loc_fastqc=$wrkspc"/Software/FastQC"
loc_trimGalore=$wrkspc"/Software/TrimGalore"
loc_STAR=$wrkspc"/Software/STAR"
loc_fCounts=$wrkspc"/Software/subread/bin"

#location of data
loc_rawReads=$wrkspc"/Tregs/01_raw_reads"
dataAnnotation=$wrkspc"/Tregs/01_raw_reads/dataAnnotation.txt"

#location of outputs
loc_fastqc_output_init=$wrkspc"/Tregs/02_qc_output_init"
loc_trimGalore_output=$wrkspc"/Tregs/03_processed_reads"
loc_fastqc_output_final=$wrkspc"/Tregs/03_qc_output_final"
loc_aligned_output=$wrkspc"/Tregs/04_aligned_reads"
loc_counts=$wrkspc"/Tregs/05_counts"

#location of genome info and annotations
loc_genomeIndex=$wrkspc"/Tregs/04_genome_Musmusculus/Index"
loc_genomeGTF=$wrkspc"/Tregs/04_genome_Musmusculus/GTF"




###Initial QC on raw data
#perform initial QC
printf "\n###\nRunning loopFastQC.sh on ${loc_rawReads} directory\n###\n\n"
timestamp=$(date +"%Y-%m-%d_%T")
printf "Started: "${timestamp}"\n\n"

#perform multi-thread FastQC analysis
${loc_fastqc}/fastqc -t 10 ${loc_rawReads}/*.fastq.gz \
-o ${loc_fastqc_output_init}

timestamp=$(date +"%Y-%m-%d_%T")
printf "\nFinished: "${timestamp}"\n\n"



###MultiQC report generation
#using MultiQC to read all the FastQC results - input/output directory and
# output filename
printf "\n###\nRunning multiQC.sh on ${loc_input} directory\n###\n\n"
timestamp=$(date +"%Y-%m-%d_%T")
printf "Started: "${timestamp}"\n\n"

multiqc -f -p -n "QC_initial_results.html" -o $loc_fastqc_output_init \
--no-megaqc-upload $loc_fastqc_output_init

timestamp=$(date +"%Y-%m-%d_%T")
printf "\nFinished: "${timestamp}"\n\n"



###Seqeuence trimming
#perform trimming via trimGalore
printf "\n###\nRunning trimGalore.sh on ${loc_rawReads} directory\n###\n\n"
timestamp=$(date +"%Y-%m-%d_%T")
printf "Started: "${timestamp}"\n\n"

#perform multi-thread TrimGalore analysis
${loc_trimGalore}/trim_galore -j 3 \
-o ${loc_trimGalore_output} ${loc_rawReads}/*.fastq.gz

##core usage from documentation:
## 1. Actual core usage: It should be mentioned that the actual number of 
#cores used is a little convoluted. Assuming that Python 3 is used and pigz is 
#installed, `--cores 2` would use 2 cores to read the input (probably not at a 
#high usage though), 2 cores to write to the output (at moderately high usage), 
#and 2 cores for Cutadapt itself + 2 additional cores for Cutadapt (not sure 
#what they are used for) + 1 core for Trim Galore itself. So this can be up to 
#9 cores, even though most of them won't be used at 100% for most of the time. 
#Paired-end processing uses twice as many cores for the validation (= writing 
#out) step. `--cores 4` would then be: 4 (read) + 4 (write) + 4 (Cutadapt) 
#+ 2 (extra Cutadapt) + 1 (Trim Galore) = 15, and so forth.

## 2. It seems that `--cores 4` could be a sweet spot, anything above has 
#diminishing returns.

timestamp=$(date +"%Y-%m-%d_%T")
printf "\nFinished: "${timestamp}"\n\n"



###Post trimming QC
#perform second QC
printf "\n###\nRunning loopFastQC.sh on ${loc_trimGalore_output} directory\n\
###\n\n"
timestamp=$(date +"%Y-%m-%d_%T")
printf "Started: "${timestamp}"\n\n"

#perform multi-thread FastQC analysis
${loc_fastqc}/fastqc \
-t 10 ${loc_trimGalore_output}/*.fq.gz \
-o ${loc_fastqc_output_final}

timestamp=$(date +"%Y-%m-%d_%T")
printf "\nFinished: "${timestamp}"\n\n"



###Sequence alignment
#using STAR to align reads - grouped by lane
printf "\n###\nRunning alignSTAR.sh via ${loc_genomeIndex} on \
${loc_trimGalore_output} directory\n###\n"
timestamp=$(date +"%Y-%m-%d_%T")
printf "Started: "${timestamp}"\n\n"

#using a loop to run STAR on files in loc_trimGalore_output via the 
#dataAnnotation file containing list of raw lane files for each sample
while IFS= read -r line
do
	words=($line)

	#reformatting file extensions to pick up trimmed files
	#file extensions differ between raw files (.fastq.gz) and trimmed files 
	#(_trimmed.fq.gz)
	for i in $(seq 1 $((${#words[@]}-1)))
	do
		words[$i]=${words[$i]%.*}
		words[$i]=${words[$i]%.*}_trimmed.fq.gz

	done

	#running STAR aligner
	${loc_STAR}/STAR --runThreadN 8 --runMode alignReads \
	--genomeDir ${loc_genomeIndex} \
	--readFilesPrefix ${loc_trimGalore_output}/ \
	--readFilesIn ${words[1]},${words[2]},${words[3]},${words[4]} \
	--readFilesCommand gunzip -c \
	--outFileNamePrefix ${loc_aligned_output}/${words[0]}. \
	--outSAMtype BAM SortedByCoordinate

#loading data annotation file at the end of the loop
done < ${dataAnnotation}

timestamp=$(date +"%Y-%m-%d_%T")
printf "\nFinished: "${timestamp}"\n\n"



###Creating counts matrix
#using featureCounts to create counts matrix
printf "\n###\nRunning fCounts.sh on ${loc_aligned_output} directory\n###\n\n"
timestamp=$(date +"%Y-%m-%d_%T")
printf "Started: "${timestamp}"\n\n"

#perform multi-thread featureCounts analysis
${loc_fCounts}/featureCounts -T 10 -a $loc_genomeGTF/*.gtf \
-o $loc_counts/counts.txt $(find $loc_aligned_output/*.bam ! \
	-name '*Undetermined*')

#post process counts_matrix to remove descriptor line 
sed -i '1d' $loc_counts/counts.txt

#post process counts_matrix to trim column headers using an R script
#R script is in the same location as this BASH script
Rscript trimColHeads.R $loc_counts/counts.txt

timestamp=$(date +"%Y-%m-%d_%T")
printf "\nFinished: "${timestamp}"\n\n"



###MultiQC report generation
#using MultiQC to read the all post filtering results - input directories (last 
#input directory is output directory) and output filename
printf "\n###\nRunning multiQC.sh on ${loc_input} directory\n###\n\n"
timestamp=$(date +"%Y-%m-%d_%T")
printf "Started: "${timestamp}"\n\n"

#running multiqc
multiqc -f -p -n "MultiQc_combined_report.html" -o $loc_counts \
--no-megaqc-upload $loc_fastqc_output_final $loc_trimGalore_output \
$loc_aligned_output $loc_counts

timestamp=$(date +"%Y-%m-%d_%T")
printf "\nFinished: "${timestamp}"\n\n"



###finishing message
printf "\nPreprocessing finished\n\n"


