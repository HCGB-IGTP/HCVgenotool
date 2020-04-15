#!/bin/bash

################################################################################
# HCV genotyping pipeline from Illumina short amplicon sequencing.
################################################################################
# Objective: detect and quantify each HCV subtype present in a sample.
#
# Usage: ./genotyping_pipeline.sh \
#        <path/to/FASTQDIR> \
# 		 <path/to/ADAPTERS_FILE> \
#		 <path/to/REFERENCE> \
#        <path/to/WORKINGDIR>
# 	FASTQDIR: Folder where all the fastq files are contained. Only Illumina
# 		pair-end sequencing from PCR amplicons is accepted. One sample pair per
# 		file.
# 		IMPORTANT: *.fastq or *.fastq.gz are accepted. Read 1 must be specified
#		as *R1* and read 2 as *R2* in the filename.	Sample name will be taken
# 		from the string before the first "_" character that should be present in
# 		the filename.
# 		Examples of valid fastq filenames:
# 			2_S1_L001_R1_001.fastq.gz -> sample name: "2", read 1.
# 			sample1_R2.fastq          -> sample name: "sample1", read 2.
# 	ADAPTERS_FILE: path to a fasta file containing adapter sequences used to
# 		create amplicons for sequencing. Only two adapters allowed, named as
# 		">Forward" and ">Reverse" (">F" and ">R" or ">f" and "r" is sufficient).
# 		Adapter sequence has to be in only one line.
# 	REFERENCE: path to a fasta file containing HCV reference sequences.
#		Sequence names should follow the convention TS_XXXXX, where T is HCV
# 		type, S is HCV subtype and XXXXX is the sequence identifier.
# 	WORKINGDIR: path to a folder where output results. If WORKINGDIR does not
# 		exist it will be created. If WORKINGDIR exists, program exits with
#		ERROR message.
#
# HCVgenotool
# Copyright (C) 2018  David Piñeyro
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.
#
# Author: David Piñeyro <igtp.genbioinfo@gmail.com>
# License: GPL-3.0-or-later
# Version: 0.0.1
# Date: February 21st, 2018
################################################################################
# FUNCTIONS
function runtime_print {
	# This function prints a formatted string for a time diff of $2 - $1.
	# Arguments:
	# 	$1: start time from $(date +%s.%N).
	# 	$2: end time from $(date +%s.%N).
	dt=$(echo "$2 - $1" | bc)
	dd=$(echo "$dt/86400" | bc)
	dt2=$(echo "$dt-86400*$dd" | bc)
	dh=$(echo "$dt2/3600" | bc)
	dt3=$(echo "$dt2-3600*$dh" | bc)
	dm=$(echo "$dt3/60" | bc)
	ds=$(echo "$dt3-60*$dm" | bc)
	if [[ $dd -gt 0 ]] ; then
	printf "%d Day/s, %d Hour/s, %d Minute/s, %.4f Seconds\n" $dd $dh $dm $ds
	elif [[ $dh -gt 0 ]] ; then
	printf "%d Hour/s, %d Minute/s, %.4f Seconds\n" $dh $dm $ds
	elif [[ $dm -gt 0 ]] ; then
	printf "%d Minute/s, %.4f Seconds\n" $dm $ds
	else
	printf "%.4f Seconds\n" $ds
	fi
}
# Starting message
TIMESTART=$(date +%s.%N)
echo ""
echo "######################################################################"
echo "#                    HCV Genotyping pipeline                         #"
echo "######################################################################"
echo ""
##########################################################
# Step0: Generating global variables.
##########################################################
# Check for the correct number of command line arguments (4).
if [[ $# -ne 4 ]] ; then
	printf "[ERROR] Incorrect number of command line arguments. Remember, use:\n
	$./genotyping_pipeline.sh \ \n
	  <path/to/FASTQDIR> \ \n
	  <path/to/ADAPTERS_FILE> \ \n
	  <path/to/REFERENCE> \ \n
	  <path/to/WORKINGDIR>\n"
	exit
fi
# Command line arguments.
FASTQDIR=$1
ADAPTERS_FILE=$2
REFERENCE=$3
WORKINGDIR=$4
# Printing info of the arguments.
echo "COMMAND LINE ARGUMENTS."
echo "FASTQDIR: $FASTQDIR"
echo "ADAPTERS_FILE: $ADAPTERS_FILE"
echo "REFERENCE: $REFERENCE"
echo "Results will be written to: $WORKINGDIR"
echo ""
# Creating WORKINGDIR.
if [ ! -d "$WORKINGDIR" ] ; then
	mkdir $WORKINGDIR
else
	echo "[ERROR] Incorrect results directory. Results directory '$WORKINGDIR' already exists."
	exit
fi
# Generating SAMPLES array with samples names.
DATAFOLDER="$WORKINGDIR/data"
mkdir $DATAFOLDER
SOURCEDIR=`pwd`  # Storing previous dir.
cd $FASTQDIR  # Caution, we are now in FASTQDIR.
ls *.fastq* | cut -f1 -d "_" | sort | uniq > "$DATAFOLDER/samples.txt"
readarray -t SAMPLES < "$DATAFOLDER/samples.txt"
# Collecting primers from ADAPTERS_FILE.
write_fwd=false
write_rev=false
while IFS='' read -r line || [[ -n "$line" ]] ; do
	l=$line
if [ $write_fwd = true ]; then
    adapter_F=$line
	write_fwd=false
fi
if [ $write_rev = true ]; then
	adapter_R=$line
	write_rev=false
fi
if [ ${l:0:2} = ">F" ] || [ ${l:0:2} = ">f" ] ; then
	write_fwd=true
fi
if [ ${l:0:2} = ">R" ] || [ ${l:0:2} = ">r" ] ; then
	write_rev=true
fi
done < "$ADAPTERS_FILE"
##########################################################
# Step1: adapter and quality trimming with Cutadapt.
##########################################################
echo "Starting adapter and quality trimming using Cutadapt..."
# Starting time.
TIME_0=$(date +%s.%N)
# Creating dir.
TRIMDIR="$WORKINGDIR/Trimm"
mkdir $TRIMDIR
# Loop over SAMPLES to trimm adapters
for sample in ${SAMPLES[@]} ; do
	echo "	Processing sample $sample"
	# Creating OUT dir
	OUT="$TRIMDIR/$sample"
	mkdir $OUT
	# Assign fastq sample names.
	FASTQ_1=`ls $sample"_"*R1*.fastq*`
	FASTQ_2=`ls $sample"_"*R2*.fastq*`
	# RUN Cutadapt
	cutadapt -g $adapter_F -G $adapter_R \
	-m 35 \
	-q 10,10 \
	-o $OUT/$sample"_trimmed_R1.fastq" \
	-p $OUT/$sample"_trimmed_R2.fastq" \
	$FASTQ_1 \
	$FASTQ_2 &> "$OUT/run.log"
done
cd $SOURCEDIR  # Comming back to original directory.
TIME_1=$(date +%s.%N)
dif_t=$( runtime_print $TIME_0 $TIME_1 )
echo "Adapter and quality trimming completed in $dif_t"
echo ""
##########################################################
# Step2: join read pairs with ea-utils/fastq-join.
##########################################################
echo "Joining read pairs using fastq-join..."
# Starting time.
TIME_0=$(date +%s.%N)
# Creating dir.
JOINDIR="$WORKINGDIR/Joined"
mkdir $JOINDIR
for sample in ${SAMPLES[@]} ; do
	echo "	Processing sample $sample"
	# Creating OUT dir
	OUT="$JOINDIR/$sample"
	mkdir $OUT
	# Assign fastq sample names.
	FASTQ_1=$TRIMDIR/$sample/$sample"_trimmed_R1.fastq"
	FASTQ_2=$TRIMDIR/$sample/$sample"_trimmed_R2.fastq"
	# RUN fastq-join
	../bin/ea-utils/clipper/./fastq-join \
	-v ' ' \
	$FASTQ_1 \
	$FASTQ_2 \
	-o "$OUT/$sample" &> "$OUT/run.log"
done
TIME_1=$(date +%s.%N)
dif_t=$( runtime_print $TIME_0 $TIME_1 )
echo "Joining completed in $dif_t"
echo ""
##########################################################
# Step3: Quality filtering with FASTX Toolkit.
##########################################################
echo "Quality filtering using FASTX..."
# Starting time.
TIME_0=$(date +%s.%N)
# Only reads with >= 95% of their bases over Q30 are kept.
# Quality filtered files in the same join/$sample dirs.
for sample in ${SAMPLES[@]} ; do
	echo "	Processing sample $sample"
	# OUTDIR
	OUT="$JOINDIR/$sample"
	# Joined fastq sample name.
	J_FASTQ=$OUT/$sample"join"
	# RUN fastq-join
	../bin/fastx-toolkit/./fastq_quality_filter \
    -Q 33 \
	-q 30 \
	-p 95 \
	-i "$J_FASTQ" \
	-o $OUT/$sample"join_qcfiltered.fastq"
done
TIME_1=$(date +%s.%N)
dif_t=$( runtime_print $TIME_0 $TIME_1 )
echo "Quality filtering completed in $dif_t"
echo ""
##########################################################
# Step4: calculate reference sequences identities.
##########################################################
echo "Calculating HCV reference identities using optimal global alignment..."
# Starting time.
TIME_0=$(date +%s.%N)
python HCV_amplicon_identity.py \
-p $ADAPTERS_FILE \
-r $REFERENCE \
-o "$DATAFOLDER/HCV"
TIME_1=$(date +%s.%N)
dif_t=$( runtime_print $TIME_0 $TIME_1 )
echo "HCV reference identities calculated in $dif_t"
echo ""
##########################################################
# Step5: BWA-mem alignment against HCV amplicons generated in Step4.
##########################################################
echo "BWA-mem alignment against HCV genotypes..."
echo ""
# Starting time.
TIME_0=$(date +%s.%N)
# Creating dir.
ALIGNDIR="$WORKINGDIR/Align"
mkdir $ALIGNDIR
# Amplicons as reference.
AMPLICONS="$DATAFOLDER/HCV_amplicons.fasta"
# Generate bwa-index files.
echo "Generating bwa index files..."
../bin/bwa/./bwa index $AMPLICONS
echo "bwa index files generated."
echo ""
# Alignment of each sample. AMPLICONS as references.
for sample in ${SAMPLES[@]} ; do
	echo "BWA alignment of sample $sample"
	# Creating OUT dir
	OUT="$ALIGNDIR/$sample"
	mkdir $OUT
	# Assign fastq sample names.
	FASTQjoined=$JOINDIR/$sample/$sample"join_qcfiltered.fastq"
	# RUN BWA mem
	../bin/bwa/./bwa mem \
	-M \
	-t 4 \
	$AMPLICONS \
	$FASTQjoined > "$OUT/$sample.sam"
	echo ""
done
TIME_1=$(date +%s.%N)
dif_t=$( runtime_print $TIME_0 $TIME_1 )
echo "BWA alignment completed in $dif_t"
echo ""
##########################################################
# Step6: calculate samples subtypes.
##########################################################
echo "Genotyping samples..."
# Starting time.
TIME_0=$(date +%s.%N)
# Creating dir.
# Generate a list of sam files.
sams=""
for sample in ${SAMPLES[@]} ; do
	sams="$sams,$ALIGNDIR/$sample/$sample.sam"
done
sams=${sams:1}  # Deleting the initial ","
# Run HCV_genotyping.py tool.
python HCV_genotyping.py \
-i "$DATAFOLDER/HCV_identity_df_per_group.csv" \
-s "$sams" \
-o "$WORKINGDIR/result.csv"
TIME_1=$(date +%s.%N)
dif_t=$( runtime_print $TIME_0 $TIME_1 )
echo "Genotyping completed in $dif_t"
echo ""

# Runtime calculation and exit message
dif_t=$( runtime_print $TIMESTART $TIME_1 )
echo "HCV genotyping pipeline completed successfully."
echo "Total runtime: $dif_t"
