#!/bin/sh

## pre-installed software needed:
# sratoolkit
# bcftools
# parallel

## files needded
# project key
# header.txt
# cart files

function logger() {
	MSG=$1
	TIMESTAMP=`TZ=EST date`
	echo "${TIMESTAMP} ${MSG}"
}
export -f logger

# Set some environment variables
export SRA_TOOLKIT=/opt/sratoolkit.2.9.0-centos_linux64
export WD=`pwd`

# get the key
KEY=`ls *.ngc`

# Get the project number
PROJECT=`echo $KEY | cut -d "." -f 1 | cut -d "_" -f 2`

# Import the key
logger "Import key"
vdb-config --import $KEY

# change the root directory
logger "Running 'vdb-config' to set ROOT directory to ${WD}"
vdb-config --set /repository/user/protected/dbGaP-${PROJECT}/root=${WD}
rm -rf ncbi

# download the tarred files from dbGaP. There is 1 .krt file (~2Mb) per patient group (between 1 and 4 patient group per study)
# it will download 1 encrypted tar.gz file (~30Gb) per patient group
KRT=`ls *.krt`
parallel prefetch -X 200G ::: $KRT

# decrypt the tar.gz files (still ~30Gb each)
parallel --link vdb-decrypt ::: `find ./files -type f` ::: `ls files | sed 's/.ncbi_enc//g'`
rm -rf files


# This function takes a vcf.gz inside a tar.gz file, and create a .vcf.bgz and ship it to s3
function tar_to_bgz {
	logger "  Conversion process for $2"
	
	VCF=`echo "$2" | cut -d "/" -f 2 | sed 's/.vcf.gz//g'`
	CONTIG=`tar -xOzf "$1" "$2" | gunzip -c | sed -n -e 50p | awk '{print $1}'`
	SAMPLES=`tar -xOzf "$1" "$2" | gunzip -c | grep -m 1 "#CHROM"`
	
	logger "    create header file for ${VCF}"
	cat header.txt | sed -e "3i ##contig=<ID=${CONTIG}>" -e "$ a ${SAMPLES}" > ${VCF}_header.txt
	
	logger "    extract $2 from $1, unzip, modify header, bgzip and ship to s3"
  tar -xOzf "$1" "$2" | gunzip -c | bcftools reheader -h ${VCF}_header.txt | bcftools view -Oz | aws s3 cp - s3://versmee-etl/greg_test/${VCF}.vcf.bgz
	
	logger "    remove ${VCF} header file"
  rm ${VCF}_header.txt
}
export -f tar_to_bgz

VCFGZ=`for TARGZ in *.tar.gz;do tar -tf ${TARGZ};done`
TAR=`for TARGZ in *.tar.gz;do for VCF in $(tar -tf ${TARGZ}); do echo $TARGZ ;done; done`

logger "Launch the whole conversion process"
parallel --ungroup --link tar_to_bgz ::: $TAR ::: $VCFGZ
logger "Full conversion process finished"

# mv the tar.gz files to s3
aws s3 mv . s3://versmee-etl/greg-test --recursive --exclude "*" --include "*.tar.gz"

# this function will merge the vcf from the same chromosome across each consent group, and create 23 merged_VCF files
function merge () {
	logger "  Start merging process for chromosome $1"
	
	logger "    download vcf.bgz files from s3 for chromosome $1"
	aws s3 cp s3://versmee-etl/greg_test/ . --recursive --exclude "*" --include "*.chr${1}.*.vcf.bgz"
	
	logger "    create index files from vcf.bgz files for chromosome $1"
	for BGZ in *.chr${1}.*.vcf.bgz
	do
		logger "      create index files for $BGZ"
		bcftools index -f $BGZ
	done
	
	VCFBGZ=`*.chr${1}.*.vcf.bgz`
	logger "    merging vcf.bgz files for chromosome $1 into merged_chr${1}.vcf.bgz and ship it to s3"
  bcftools merge -m none -O z $VCFBGZ | aws s3 cp - s3://versmee-etl/greg_test/merged_chr${1}.vcf.bgz
	
	logger "    remove vcf.bgz files for chromosome $1"
	rm $VCFBGZ
	
	logger "    remove index files for chromosome $1"
	find -name "*.chr${1}.*.vcf.bgz.csi" -delete
}
export -f merge

logger "Launch the whole merging process"
parallel --ungroup merge ::: {1..22} X
logger "Full merging process finished"



