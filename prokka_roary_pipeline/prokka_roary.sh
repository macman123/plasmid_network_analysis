#!/bin/bash

# Prokka and Roary:
roary=/path/to/roary
prokka=/path/to/prokka
prokka_run=/path/to/prokka_run.sh

# Variables:
# base directory should just contain fasta files to analyze
user=username
baseDirectory=/base/directory/
h_rt=96:00:00
tmem=12G
h_vmem=12G
free_scratch=4000 # in MB
max_iter=10 # max number of re-submissions of failed jobs
refresh_rate=10s #refresh rate for submission status


# RUNNING SCRIPT:
# split fasta file:
# awk '/^>/ {OUT=substr($1,2) ".fasta"}; OUT {print > OUT}' most_abundant_incs.fasta
cd $baseDirectory
#for file in $baseDirectory/*; do
#  echo ${file##*/} >> fasta_files.txt
#done
fasta_files=fasta_files.txt
nr_files=$(wc -l < $fasta_files) 
echo 
echo ">>>> Prokka-Roary PIPLINE <<<<"
date
echo
echo "Author: Mislav Acman"
echo "Base Directory: " $baseDirectory
echo "Fasta files: " $nr_files
echo "Max iterations: " $max_iter

# Modifying prokka_run for submission
echo
echo "    --> Modifying prokka_run.sh for submission"
cp $prokka_run temp_prokka_run.sh
# Line 3:
sed -e "3s@.*@#$ -l h_rt=$h_rt@" temp_prokka_run.sh > temp.sh
mv temp.sh temp_prokka_run.sh
# Line 5:
sed -e "5s@.*@#$ -l tmem=$tmem@" temp_prokka_run.sh > temp.sh
mv temp.sh temp_prokka_run.sh
# Line 6:
sed -e "6s@.*@#$ -l h_vmem=$h_vmem@" temp_prokka_run.sh > temp.sh
mv temp.sh temp_prokka_run.sh
# Line 9:
sed -e "9s@.*@#$ -t 1-$nr_files@" temp_prokka_run.sh > temp.sh
mv temp.sh temp_prokka_run.sh
# Line 12:
sed -e "12s@.*@prokka=$prokka@" temp_prokka_run.sh > temp.sh
mv temp.sh temp_prokka_run.sh
# Line 13:
sed -e "13s@.*@baseDirectory=$baseDirectory@" temp_prokka_run.sh > temp.sh
mv temp.sh temp_prokka_run.sh
# Line 14:
sed -e "14s@.*@fasta_files=$fasta_files@" temp_prokka_run.sh > temp.sh
mv temp.sh temp_prokka_run.sh
# Line 15:
sed -e "15s@.*@free_scratch=$free_scratch@" temp_prokka_run.sh > temp.sh
mv temp.sh temp_prokka_run.sh


# Running prokka
chmod 755 temp_prokka_run.sh
calc(){ awk 'BEGIN { printf "%.2f", '$*' }'; }
iteration=1
while [ "$iteration" -le "$max_iter" -o "$nr_files" -ne 0 ]; do
	echo
	qsub temp_prokka_run.sh
	stat=`qstat | grep $user | grep "prokka_run"`
	jobID=`echo $stat | cut -d ' ' -f 1`
	jobs_complete=1
	perc=0.00
	echo
	while [ -n "$stat" ]; do
                tput cuu1
                tput el
                echo -e "    --> qsub prokka_run.sh\titeration: $iteration\tjobs: $jobs_complete/$nr_files\t$perc%"
		stat=`qstat | grep $user | grep "prokka_run"`
		qw=`qstat | grep $user | grep "prokka_run" | grep "qw"`
		if [ -n "$qw" ]; then
			jobs_complete=`qstat -u $user | grep "prokka_run" | grep "qw" | awk '{ print $9 }' | cut -d'-' -f1 | cut -d' ' -f1`
			perc=`calc $jobs_complete/$nr_files*100`
		fi
		sleep $refresh_rate
	done
	tput cuu1
        tput el
        echo -e "    --> qsub prokka_run.sh\titeration: $iteration\tjobs: $nr_files/$nr_files\t100.00%"

	
	# Checking errors
        awk -F"\t" '{print $1}' errors.txt > temp_fasta_files.txt
        mv errors.txt errors_iter$iteration.txt
	
	# Checking failed prokka runs
	counter=1
	while read p; do
		folder=`echo $p | cut -d'.' -f1`
		if grep -q "Error:" prokka_run.o$jobID.$counter || [ ! -f $folder/*.gff ]; then
			echo $p >> temp_fasta_files.txt
			rm -rf $folder
		fi
		counter=$(($counter+1))
	done <$fasta_files	

	rm $fasta_files
	mv temp_fasta_files.txt $fasta_files
	echo "    --> Failed runs: `wc -l < $fasta_files`"

	# concatinating log file
	cat prokka_run.o$jobID.* > prokka_run_iter$iteration.log
	rm prokka_run.o$jobID.*

	#Updating prokka_run
	iteration=$(($iteration + 1))
	nr_files=`wc -l < $fasta_files`
	sed -e "9s@.*@#$ -t 1-$nr_files@" temp_prokka_run.sh > temp.sh
	mv temp.sh temp_prokka_run.sh
done

echo
echo "    --> Iterations complete!"
# Running roary
echo "    --> Coppying .gff files"
mkdir gff_files
find ./ -name '*.gff' -exec cp -t gff_files {} +
echo 
echo "#$ -S /bin/bash" > roary_run.sh
echo "#$ -N roary_run" >> roary_run.sh
echo "#$ -l h_rt=$h_rt" >> roary_run.sh
echo "#$ -cwd" >> roary_run.sh
echo "#$ -l tmem=16G" >> roary_run.sh
echo "#$ -l h_vmem=16G" >> roary_run.sh
echo "#$ -j y" >> roary_run.sh
echo "#$ -t 1-1" >> roary_run.sh
echo "" >> roary_run.sh
echo "source /share/apps/genomics/Roary-master/roary_env.sh" >> roary_run.sh  # this line is specific for the cluster I was using
echo "cd $baseDirectory" >> roary_run.sh
echo "$roary -e -mafft -p 8 $baseDirectory/gff_files/*.gff -f roary_output -v" >> roary_run.sh
chmod 755 roary_run.sh
qsub roary_run.sh
echo "    --> Running roary"
stat=`qstat | grep $user | grep "roary_run"`
jobID=`echo $stat | cut -d ' ' -f 1`
while [ -n "$stat" ]; do
        sleep $refresh_rate
	stat=`qstat | grep $user | grep "roary_run"`
done

mv roary_run.o$jobID.1 roary_run.log

echo
date
echo ">>>> Prokka-Roary DONE! <<<<"

# Removing files
rm $fasta_files
rm temp_prokka_run.sh
rm roary_run.sh
