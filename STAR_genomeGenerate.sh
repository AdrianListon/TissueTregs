
## commands to get the mus musculus GRCm38 assembly GTF and Fasta and use STAR to create the genome index
## run in the directory you want the genome work to be done (it was .../04_musmusculus_genome/... for me)

#getting the files from gencode
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M23/gencode.vM23.primary_assembly.annotation.gtf.gz
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M23/GRCm38.primary_assembly.genome.fa.gz

mkdir Fasta
mkdir GTF
mkdir Index

gunzip gencode.vM23.primary_assembly.annotation.gtf.gz
mv gencode.vM23.primary_assembly.annotation.gtf.gz ./GTF/

gunzip GRCm38.primary_assembly.genome.fa.gz
mv GRCm38.primary_assembly.genome.fa.gz ./Fasta

#loading the job submission and STAR modules
module load ssub
module load STAR

#submit the job
ssub -o output_STAR_makegenome.log -c16 --mem=128GB --email \
STAR --runThreadN 16 \
--runMode genomeGenerate \
--genomeDir Index \
--genomeFastaFiles /bi/group/liston/samar/Tregs/04_genome_Musmusculus/GRCm38.primary_assembly.genome.fa \
--sjdbGTFfile /bi/group/liston/samar/Tregs/04_genome_Musmusculus/gencode.vM23.primary_assembly.annotation.gtf \
--sjdbOverhang 100

#print job queue
printf "Submitted STAR genomeGenerate"
squeue

