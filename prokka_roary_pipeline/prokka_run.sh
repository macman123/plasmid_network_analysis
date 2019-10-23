#$ -S /bin/bash
#$ -N prokka_run
#$ -l h_rt=24:00:00
#$ -cwd
#$ -l tmem=12G
#$ -l h_vmem=12G
#$ -l tscratch=4G
#$ -j y # put output and error into single file
#$ -t 1-1 #number of jobs

# Variables:
prokka=/path/to/prokka
baseDirectory=/base/directory/
fasta_files=fasta_files.txt
free_scratch=4000
fasta_file=$( sed "${SGE_TASK_ID}q;d" $fasta_files )
prefix=$( echo $fasta_file | cut -d'.' -f 1 )
CPUs=$( grep -c ^processor /proc/cpuinfo )
prokka_cmnd="$prokka $fasta_file --prefix $prefix --cpus $CPUs"


echo ">>>> Running Prokka <<<<"
echo
echo date
echo "Node:" `hostname`
echo "Job ID: $JOB_ID"
echo "Process ID: $$"
echo "Base Directory: " $baseDirectory
echo "File: " $fasta_file
echo "Prefix: " $prefix 
echo "CPUs: " $CPUs " - using all"
echo "Prokka command: " $prokka_cmnd
echo

### Checking scratch space
free_space=`df -m /scratch0/ | tail -1 | awk '{print $4}'`
echo -n > errors.txt
if [ "$free_space" -lt "$free_scratch" ]; then
	echo "There is " $free_space "MB scratch space"
	echo "Error: not enough scratch space!"
	cd $baseDirectory
	echo -e "$prefix\t`hostname`\t`date`\tNo scratch space" >> errors.txt
        kill -9 $$
else
	echo "There is " $free_space "MB scratch space"
	echo "Checkpoint passed!"
	echo
fi


### Execute prokka
mkdir /scratch0/prokka_run/
cd $baseDirectory
cp $fasta_file /scratch0/prokka_run/
cd /scratch0/prokka_run/
eval $prokka_cmnd
cp -r $prefix $baseDirectory
rm -rf /scratch0/prokka_run/

