#!/bin/sh

## pre-installed software needed:
# bcftools
# parallel

## files needded
# project key
# cart files

# update yum
sudo yum -y update
sudo yum -y install gcc

# install sratoolkit
wget https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/current/sratoolkit.current-centos_linux64.tar.gz
tar -xf sratoolkit.current-centos_linux64.tar.gz
rm sratoolkit.current-centos_linux64.tar.gz
export PATH=~/sratoolkit.2.9.0-centos_linux64/bin:$PATH

# install bcftools
brew install bcftools

# install parallel
brew install parallel

# get the header file
wget https://raw.githubusercontent.com/gversmee/vcf_to_vds/master/header.txt

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
prefetch -X 200G $KRT

# decrypt the tar.gz.ncbi_enc files (still ~30Gb each)
logger "Decrypting files"
ENC=`find ./files -type f`
parallel --link vdb-decrypt ::: $ENC ::: `echo $ENC | xargs -n1 basename | sed 's/\.ncbi_enc//g'`

# get folder name
export TNAME=`ls *.tar.gz | head -1`
export FNAME=`echo $TNAME| cut -d "." -f 1-5`
export STNAME=`tar -tf ${TNAME} | head -2 | tail -1 | cut -d "/" -f 2 | cut -d "_" -f 3`
export NAME="${STNAME}.${FNAME}"
logger "Will place the files in the folder ${NAME} on s3://versmee-etl"

# mv the cart files and the crypted files to s3
logger "Moving the kart files to s3"
aws s3 mv . s3://versmee-etl/${NAME}/krt --recursive --exclude "*" --include "*.krt"
aws s3 mv files s3://versmee-etl/${NAME}/tar_gz_encrypted --recursive --exclude "*" --include "*"
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
  tar -xOzf "$1" "$2" | gunzip -c | bcftools reheader -h ${VCF}_header.txt | bcftools view -Oz | aws s3 cp - s3://versmee-etl/${NAME}/vcf_bgz/${VCF}.vcf.bgz
	
	logger "    remove ${VCF} header file"
  rm ${VCF}_header.txt
}
export -f tar_to_bgz

VCFGZ=`for TARGZ in *.tar.gz;do tar -tzf ${TARGZ} | grep -v '/$';done`
TAR=`for TARGZ in *.tar.gz;do for VCF in $(tar -tzf ${TARGZ} | grep -v '/$'); do echo $TARGZ ;done; done`

logger "Launch the whole conversion process"
parallel --ungroup --link tar_to_bgz ::: $TAR ::: $VCFGZ
logger "Full conversion process finished"

# mv the tar.gz files to s3
########### define name before
logger "Moving the tar.gz files to s3"
aws s3 mv . s3://versmee-etl/${NAME}/tar_gz --recursive --exclude "*" --include "*.tar.gz"

# this function will merge the vcf from the same chromosome across each consent group, and create 23 merged_VCF files
function merge () {
	logger "  Start merging process for chromosome $1"
	
	logger "    download vcf.bgz files from s3 for chromosome $1"
	aws s3 cp s3://versmee-etl/${NAME}/vcf_bgz/ . --recursive --exclude "*" --include "*.chr${1}.*.vcf.bgz"
	
	logger "    create index files from vcf.bgz files for chromosome $1"
	for BGZ in *.chr${1}.*.vcf.bgz
	do
		logger "      create index files for $BGZ"
		bcftools index -f $BGZ
	done
	
	VCFBGZ=`find -name "*.chr${1}.*.vcf.bgz"`
	logger "    merging vcf.bgz files for chromosome $1 into merged_chr${1}.vcf.bgz and ship it to s3"
  bcftools merge -m none -O z $VCFBGZ | aws s3 cp - s3://versmee-etl/${NAME}/merged/merged_chr${1}.${NAME}.vcf.bgz
	
	logger "    remove vcf.bgz files for chromosome $1"
	rm $VCFBGZ
	
	logger "    remove index files for chromosome $1"
	find -name "*.chr${1}.*.vcf.bgz.csi" -delete
}
export -f merge

logger "Launch the whole merging process"
parallel --ungroup merge ::: 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X
logger "Full merging process finished"

rm dbgap_to_hail.sh header.txt

aws ec2 stop-instances --instance-ids i-0e66363ec98851c6e