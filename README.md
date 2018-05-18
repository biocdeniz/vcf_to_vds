# vcf_to_vds

#need to import the project key and the cart files

wget https://raw.githubusercontent.com/gversmee/vcf_to_vds/master/dbgap_to_hail.sh
wget https://raw.githubusercontent.com/gversmee/vcf_to_vds/master/header.txt
nohup ./dbgap_to_hail.sh > output.txt 2> err.txt &
