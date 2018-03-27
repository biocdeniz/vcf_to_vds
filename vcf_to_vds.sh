#!/bin/sh


# This script will download genomics files from dbGaP, and process them to create .vds files

# Set some environment variables
export SRA_TOOLKIT=/opt/sratoolkit.2.9.0-centos_linux64
export BCF_TOOLS=/opt/bcftools-1.7

# set the working directory
if [ -z "$WD" ]; then
    echo "Working directory is not set. Using default `/var/tmp`"
    export WD=/var/tmp
fi
cd $WD

# install the sratoolkit

# install aspera

# install bcftools

# install parallel

# get the key
KEY=*.ngc

# Get the project number
PROJECT=`echo $KEY | cut -d "." -f 1 | cut -d "_" -f 2`

# Import the key
$SRA_TOOLKIT/bin/vdb-config --import $KEY

# change the root directory
$SRA_TOOLKIT/bin/vdb-config --set /repository/user/protected/dbGaP-${PROJECT}/root=${WD}

# get the convert_hail file
HAIL_SCRIPT=greg_write_vcf.py

# download the tarred files from dbGaP. There is 1 .krt file (~2Mb) per patient group (between 1 and 4 patient group per study)
# it will download 1 encrypted tar.gz file (~30Gb) per patient group
ls *.krt | /usr/local/bin/parallel $SRA_TOOLKIT/bin/prefetch -X 200G

# decrypt the tar.gz files (still ~30Gb each)
$SRA_TOOLKIT/bin/vdb-decrypt ${WD}/files

# untar them
# each .tar.gz files contains 23 .vcf.gz files (~1.5Gb each)
function untar_files() {
  tar -xzvf $1
	rm -f $1
}
export -f untar_files

find -name *.tar* -type f | /usr/local/bin/parallel untar_file

# this function will unzip the .vcf.gz (~1.5Gb each) to a .vcf file (~100Gb each!!!!!!!!), fix the mistakes inside, and the re-compress them into a .vcf.bgz file (~1.5Gb each)
function gz_to_bgz () {
	VCF=`echo $1 | sed 's/.gz//g'`
  	gunzip -c "$1" | sed --file=VCFConversion.txt > $VCF
	rm -f $1
	$BCF_TOOLS/bin/bcftools view $VCF -Oz -o ${VCF}.bgz
	rm -f $VCF
	$BCF_TOOLS/bin/bcftools index -f ${VCF}.bgz
	
}
export -f gz_to_bgz
	
find -name *.vcf.gz -type f | /usr/local/bin/parallel gz_to_bgz

# this function will merge the vcf from the same chromosome across each consent group, and create 23 merged_VCF files
function merge () {
	VCF_LIST=`find -name cg*chr${1}.*vcf.bgz -type f`
  	$BCF_TOOLS/bin/bcftools merge -m none -O z -o merged_chr${1}.vcf.bgz ${VCF_LIST}
	rm $VCF_LIST
}
export -f merge

/usr/local/bin/parallel merge ::: {1..22}
merge X

# ships everything to s3
aws s3 cp . s3://avl-hail-dev/greg_stuff/ --recursive --exclude "*" --include "merged*"

# cleans up the mess
find -name *.vcf.bgz* -type f -delete

# create the 1 vds file from the 23 merged_files using hail
/opt/hail/HailProxy/src/scripts/spark_submit.sh $HAIL_SCRIPT s3://avl-hail-dev/greg_stuff/merged*.vcf.bgz s3://avl-hail-dev/greg_stuff/${PROJECT}.vds
aws s3 rm s3://avl-hail-dev/greg_stuff/ --recursive --exclude "*" --include "merged*.vcf.bgz"
aws s3 rm s3://avl-hail-dev/greg_stuff/ --recursive --exclude "*" --include "*_*folder*"

