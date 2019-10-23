#!/bin/sh

# This script downloads plasmid release from NCBI Refseq
# https://ftp.ncbi.nlm.nih.gov/refseq/release/plasmid/

echo
echo "Downloading current RefSeq plasmid release from:"
echo "https://ftp.ncbi.nlm.nih.gov/refseq/release/plasmid/"

# ASPERA
ascp=/path/to/ascp
ascp_id_dsa=/path/to/asperaweb_id_dsa.openssh
download_directory=/path/to/download/directory/plasmidDB

# FTP download
# if ASCP is not provided the script will attempt to perform FTP download
if [ -z "$ascp" ]; then
	echo " --> Using wget"
	wget -q --recursive -e robots=off --reject "index.html" -nH --cut-dirs=2 --no-parent https://ftp.ncbi.nlm.nih.gov/refseq/release/plasmid/ -P $download_directory
# ASPERA download
else
	echo " --> Using Aspera"
	$ascp -i $ascp_id_dsa -q -k1 -Tr anonftp@ftp.ncbi.nlm.nih.gov:refseq/release/plasmid/ $download_directory
fi

# Checkpoint
if [ -z "$(ls -A $download_directory/plasmid)" ]; then
   echo " --> Download failed!"
   kill $$
else
   echo " --> Download coplete!"
fi

# Unziping the files
echo " --> Gunzipping files"
cd $download_directory/plasmid
gunzip *.gz

# Concatinating files
echo " --> Concatinating files"
cat *.genomic.fna > $download_directory/plasmidDB.fasta
cat *.genomic.gbff > $download_directory/plasmidDB.gbff
cat plasmid.[0-9].protein.faa > $download_directory/pladmid.protein.faa
cat plasmid.[0-9].protein.gpff > $download_directory/pladmid.protein.gpff
cat plasmid.nonredundant_protein.[0-9].protein.faa > $download_directory/plasmid.nonredundant_protein.faa
cat plasmid.nonredundant_protein.[0-9].protein.gpff > $download_directory/plasmid.nonredundant_protein.gpff
cd $download_directory
rm -r $download_directory/plasmid

echo " --> DONE!"
echo
