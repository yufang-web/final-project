module load trimmomatic/0.39
module load fastqc/0.11.7
for file in $(ls *.gz)
do
gzip -d $file
done
for file in $(ls *.fastq)
do
fastqc $file
done
for file in $(ls *_1.fastq)
do
base=$(basename $file "_1.fastq")
trimmomatic PE -phred33   ${base}_1.fastq  ${base}_2.fastq    ${base}_1P.fastq ${base}_1U.fastq   ${base}_2P.fastq  ${base}_2U.fastq ILLUMINACLIP:${HPC_TRIMMOMATIC_ADAPTER}/TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:36
done
fastqc *P.fastq
module load hisat2/2.2.1
for sample in `ls *_1P.fastq`
#Execute the task
do
base=$(basename $sample "_1P.fastq")
hisat2 -p 8 --rg-id=${base} --rg PL:ILLUMINA \
-x /blue/gms6231/share/Final_Project/TEAM_B/Ref_Dir/chr12_index --dta --rna-strandness RF \
-1 ${base}_1P.fastq \
-2 ${base}_2P.fastq \
-S ${base}.sam
done

module load samtools/1.12
module load htseq/0.11.2
for file in  $(ls *.sam)
do
samtools sort -@ 8 -o  ${file%%.*}.bam  $file
done


for file in $(ls *.bam)
do
htseq-count --format bam --order pos --mode intersection-strict --stranded reverse --minaqual 1 --type exon --idattr gene_id $file /blue/gms6231/share/Final_Project/TEAM_B/Ref_dir/chr12_Homo_sapiens.GRCh38.95.gtf > ${file%%.*}.tsv
done

##union
for file in $(ls *.bam)
do
htseq-count --format bam --order pos --mode union --nonunique all --stranded reverse --minaqual 1 --type exon --idattr gene_id $file /blue/gms6231/share/Final_Project/TEAM_B/Ref_dir/chr12_Homo_sapiens.GRCh38.95.gtf > ${file%%.*}_all.tsv
done

for file in $(ls *.bam)
do
htseq-count --format bam --order pos --mode union --nonunique none --stranded reverse --minaqual 1 --type exon --idattr gene_id $file /blue/gms6231/share/Final_Project/TEAM_B/Ref_dir/chr12_Homo_sapiens.GRCh38.95.gtf > ${file%%.*}_none.tsv
done


join KO1_all.tsv KO2_all.tsv | join - KO3_all.tsv | join - rescue1_all.tsv | join - rescue2_all.tsv | join - rescue3_all.tsv > gene_read_counts_table_all.tsv
echo "GeneID KO1_all KO2_all KO3_all rescue1_all rescue2_all rescue3_all" > header.txt
cat header.txt gene_read_counts_table_all.tsv | grep -v "__" | awk -v OFS="\t" '$1=$1' > gene_read_counts_table_all_final.tsv
rm -f gene_read_counts_table_all.tsv header.txt
head gene_read_counts_table_all_final.tsv | column -t 



join KO1_none.tsv KO2_none.tsv | join - KO3_none.tsv | join - rescue1_none.tsv | join - rescue2_none.tsv | join - rescue3_none.tsv > gene_read_counts_table_none.tsv
echo "GeneID KO1_none KO2_none KO3_none rescue1_none rescue2_none rescue3_none" > header.txt
cat header.txt gene_read_counts_table_all.tsv | grep -v "__" | awk -v OFS="\t" '$1=$1' > gene_read_counts_table_all_final.tsv
rm -f gene_read_counts_table_all.tsv header.txt
head gene_read_counts_table_all_final.tsv | column -t 

perl -ne 'if ($_ =~ /gene_id\s\"(ENSG\S+)\"\;/) { $id = $1; $name = undef; if ($_ =~ /gene_name\s\"(\S+)"\;/) { $name = $1; }; }; if ($id && $name) {print "$id\t$name\n";} if ($_=~/gene_id\s\"(ERCC\S+)\"/){print "$1\t$1\n";}' /blue/gms6231/share/Final_Project/TEAM_B/Ref_dir/chr12_Homo_sapiens.GRCh38.95.gtf | sort | uniq > ENSG_ID2Name.txt






