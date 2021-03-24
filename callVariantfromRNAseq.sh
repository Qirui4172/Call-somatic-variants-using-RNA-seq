#!/usr/bin/bash


#===============================================================================================
genomeDir="/home/qirui/ref/mm10"
refGenome="/home/qirui/ref/mm10/mm10.fa"
gtfDir="/home/qirui/ref/mm10"
knownIndels="/home/qirui/snpIndel/mm10-onlyIndels.vcf"
knownSnpIndels="/home/qirui/snpIndel/mm10-snpIndels.vcf"
workDir="/home/qirui/somaticVariants_mouse"
cd ${workDir}

#===============================================================================================
Samples=("sample1_tumor" "sample2_tumor" "sample3_tumor" "sample1_normal" "sample2_normal" "sample3_normal")

#===============================================================================================
# Step 1: 2-pass mapping to genome

module load GCC/5.4.0-2.26  OpenMPI/1.10.3
module load STAR/2.6.0c
module load SAMtools/1.4

for sample in ${Samples[@]}
do
	time=`date "+%Y-%m-%d %H:%M:%S"`
	echo -e "\t$time\t$sample: Mapping ..."

	read1=${sample}".fp"
	read2=${sample}".rp"
	mkdir -p mapped/${sample} mapped/${sample}/pass1 mapped/${sample}/pass2 mapped/${sample}/index

	# 1-pass mapping
	STAR --genomeDir $genomeDir --readFilesIn trimmed/${read1} trimmed/${read2} --outFileNamePrefix mapped/${sample}/pass1/ --runThreadN 10 --readFilesCommand zcat

	# index
	STAR --runMode genomeGenerate --genomeDir mapped/${sample}/index --genomeFastaFiles $refGenome --sjdbFileChrStartEnd mapped/${sample}/pass1/SJ.out.tab --sjdbOverhang 74 --outFileNamePrefix mapped/${sample}/index/ --runThreadN 10 --readFilesCommand zcat

	# 2-pass mapping
	STAR --genomeDir mapped/${sample}/index --readFilesIn trimmed/${read1} trimmed/${read2} --outFileNamePrefix mapped/${sample}/pass2/ --outSAMtype BAM Unsorted --outSAMattrRGline ID:$sample SM:$sample PL:ILLUMINA --runThreadN 10 --readFilesCommand zcat

	# stastic summary
	samtools sort -@ 10 mapped/${sample}/pass2/Aligned.out.bam > mapped/${sample}/pass2/${sample}.bam
	samtools index -@ 10 mapped/${sample}/pass2/${sample}.bam mapped/${sample}/pass2/${sample}.bai
	samtools flagstat -@ 10 mapped/${sample}/pass2/${sample}.bam > mapped/${sample}/pass2/${sample}.flagstat

done
module purge

#===============================================================================================
# Step 2: Mark duplicate

module load picard/2.6.0-Java-1.8.0_131
module load GCC/5.4.0-2.26  OpenMPI/1.10.3 SAMtools/1.4

for sample in ${Samples[@]}
do
	time=`date "+%Y-%m-%d %H:%M:%S"`
	echo -e "\t$time\t$sample: Marking duplicates ..."

	java -jar $EBROOTPICARD/picard.jar MarkDuplicates I=mapped/${sample}/pass2/${sample}.bam O=markDuplicate/${sample}.bam M=markDuplicate/${sample}.mx
	samtools index -@ 10 markDuplicate/${sample}.bam markDuplicate/${sample}.bai

done
module purge

#===============================================================================================# Split N cigar
# Step 3: Split N Cigar
module load GATK/3.6-Java-1.8.0_92
module load GCC/5.4.0-2.26  OpenMPI/1.10.3 SAMtools/1.4

for sample in ${Samples[@]}
do
	time=`date "+%Y-%m-%d %H:%M:%S"`
	echo -e "\t$time\t$sample: Spliting N Cigar ..."

	java -jar $EBROOTGATK/GenomeAnalysisTK.jar -T SplitNCigarReads -R $refGenome -I markDuplicate/${sample}.bam -o splitCigar/${sample}.bam -rf ReassignOneMappingQuality -RMQF 255 -RMQT 60 -U ALLOW_N_CIGAR_READS
	samtools index -@ 10 splitCigar/${sample}.bam splitCigar/${sample}.bai

done
module purge

#===============================================================================================
# Step 4: Realign to indel regions

module load GATK/3.6-Java-1.8.0_92
module load GCC/5.4.0-2.26 OpenMPI/1.10.3 SAMtools/1.4

for sample in ${Samples[@]}
do
	time=`date "+%Y-%m-%d %H:%M:%S"`
	echo -e "\t$time\t$sample: Realigning to indel regions ..."

	# create intervals
	java -jar $EBROOTGATK/GenomeAnalysisTK.jar -T RealignerTargetCreator -R $refGenome -known $knownIndels -I splitCigar/${sample}.bam --filter_reads_with_N_cigar -o realign/${sample}.intervals

	# realign to intervals
	java -jar $EBROOTGATK/GenomeAnalysisTK.jar -T IndelRealigner -R $refGenome -known $knownIndels -I splitCigar/${sample}.bam -targetIntervals realign/${sample}.intervals -o realign/${sample}.bam

done
module purge

#===============================================================================================
# Step 5: Recalibrate

module load GATK/3.6-Java-1.8.0_92
module load GCC/7.3.0-2.30 OpenMPI/3.1.1 R/3.5.1

for sample in ${Samples[@]}
do
	time=`date "+%Y-%m-%d %H:%M:%S"`
	echo -e "\t$time\t$sample: Recalibrating ..."

	java -jar $EBROOTGATK/GenomeAnalysisTK.jar -T BaseRecalibrator -R $refGenome -knownSites $knownSnpIndels -I realign/${sample}.bam -o recalibrate/${sample}.table

	java -jar $EBROOTGATK/GenomeAnalysisTK.jar -T BaseRecalibrator -R $refGenome -knownSites $knownSnpIndels -I realign/${sample}.bam -BQSR recalibrate/${sample}.table -o recalibrate/${sample}.post.table

	java -jar $EBROOTGATK/GenomeAnalysisTK.jar -T PrintReads -R $refGenome -I realign/${sample}.bam -BQSR recalibrate/${sample}.table -o recalibrate/${sample}.bam

done
module purge

#===============================================================================================
# Step 6: Pile up

module load GCC/5.4.0-2.26 OpenMPI/1.10.3 SAMtools/1.4

for sample in ${Samples[@]}
do
	time=`date "+%Y-%m-%d %H:%M:%S"`
	echo -e "\t$time\t$sample: Piling up ..."

	samtools mpileup -B -d 100000 -f $refGenome -o piledup/${sample}.mpup recalibrate/${sample}.bam

done
module purge

#===============================================================================================
# Step 7: Call variant

module load VarScan/2.4.1-Java-1.8.0_92

for sample in "sample1" "sample2" "sample3"
do
	time=`date "+%Y-%m-%d %H:%M:%S"`
	echo -e "\t$time\t$sample: Calling variants ..."

	java -jar $EBROOTVARSCAN/VarScan.v2.4.1.jar somatic piledup/${sample}_normal.mpup piledup/${sample}_tumor.mpup callVariants/${sample} --min-var-freq 0.001 --strand-filter 0 --output-vcf 1

done
module purge

#===============================================================================================
# Step 8: Annotate variant

module load GCC/6.4.0-2.28  OpenMPI/2.1.1
module load annovar/2017Jul16-Perl-5.26.0

for sample in ${Samples[@]}
do
	time=`date "+%Y-%m-%d %H:%M:%S"`
	echo -e "\t$time\t$sample: Annotating variants ..."

	# extract somatic mutations, merge snv and indel into one file
	grep -E '^#|SS=2' callVariants/${sample}.snp.vcf > callVariants/${sample}.vcf
	rownum=$(grep 'SS=2' callVariants/${sample}.indel.vcf|wc -l)
	if [[ $rownum != 0 ]];then
		grep 'SS=2' callVariants/${sample}.indel.vcf >> callVariants/${sample}.vcf
	fi

	rownum=$(grep 'SS=2' callVariants/${sample}.vcf|wc -l)
	if [[ $rownum == 0 ]];then
		echo "Sample ${sample} has no row of SS=2"
		continue
	fi

	# annotate
	convert2annovar.pl -format vcf4old callVariants/${sample}.vcf > callVariants/${sample}.avi
	table_annovar.pl --thread 10 callVariants/${sample}.avi $gtfDir -buildver mm10 -remove -protocol refGene,snp138 -operation g,f -nastring . -polish -out annotated/${sample}

done
module purge

#------------------------------------------------------------------------
time=`date "+%Y-%m-%d %H:%M:%S"`
echo -e "\t$time\t$sample: Done with all samples!"

