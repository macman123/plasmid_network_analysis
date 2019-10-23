
@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
@ ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ @
@ ~~ Performing clique detection in a dataset of complete bacterial plasmids ~~ @
@ ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ @
@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

This readme file and acompanying scripts are a guide to construct and analyse bacterial plasmid 
networks with OSLOM community detection algorithm. It is directly related to work presented in:

	Acman, M., van Dorp, L., Santini, J. M., & Balloux, F. (2019). Large-scale network analysis
    captures biological features of bacterial plasmids. bioRxiv, 785212.
    doi: https://doi.org/10.1101/785212

Please cite the paper if you are using any of the methods or scripts found in this repository.
For any questions relating the analysis or the results feel free to contact the corresponding 
authors.


###########################
## Required Dependencies ##
###########################

Python3 with the following packages/modules: 
	SeqIO, NCBITaxa, os, ssl, re, math, sys, numpy

BLAST tools (2.6.0 or above) (ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/)
	makeblastdb and blastn

MOBtyping tool (https://github.com/AlexOrlek/MOBtyping)

Prokka (1.12 or above) (https://github.com/tseemann/prokka)

Roary (3.12.0 or above) (https://github.com/sanger-pathogens/Roary)

BinDash (0.2.1 or above) (https://github.com/zhaoxiaofei/bindash)

OSLOM2 (2.5 or above) (http://www.oslom.org)


##############################################
## A Dataset of complete bacterial plasmids ##
##############################################

A dataset of bacterial plasmids was downloaded from NCBI's public RefSeq repository.
(https://ftp.ncbi.nlm.nih.gov/refseq/release/plasmid/ - accessed on 26th September 2018)

The accession numbers of plasmid sequences used in the project are available in 
supplementary_table_1.tsv. The table also contains all the accompanyng metadata associated
with the plasmid sequences.

The custom made script for assembling the dataset is available in building_plasmidDB/ directory.
The directory contains the following scripts:
   --> download_plasmidDB.sh    -  Bash script used to download the dataset
   --> extractmeta_plasmidDB.py -  Python script which extracts the metadata from plasmid GenBank
                                   files. It creates a table with the following information for 
                                   every plasmid sequence: sequence description, given plasmid name,
                                   sequence length, plasmid host taxonomy, collection date and 
                                   location, isolation source, notes
   --> curate_plasmidDB.py      -  Python script used to curate the dataset. The script scans the 
                                   sequence description for the following regular expression: 
                                                   'plasmid.*complete sequence'
                                   Non-matching sequences are discarded. It also requires that the 
                                   sequence host is from Bacteria domain. In the end, the script
                                   plots a plasmid length distribution before and after filtering.
   
   After running the scripts above download the PlasmidFinder database:
         https://bitbucket.org/genomicepidemiology/plasmidfinder_db/src/master/
   This database contains conserved replicon regions of plasmid sequences used in replicon typing
   Blast (blastn) the plasmid sequences (in plasmidDB.fasta) against the PlasmidFinder database.
     NOTE: Use output format 6 with blastn and be careful not to limit the number of outputs 
           (by default is limited to 500). Do not put any other restrictions on blastn.       
   --> append_type_to_meta.py   -  Python script used to analyse PlasmidFinder blast results.
                                   The script is filtering replicon overlaps and calling replicon
                                   types which have 95% coverage and 95% similarity. The results are
                                   added to the metadata table. 

MOB typing was performed using MOBtyping tool (https://github.com/AlexOrlek/MOBtyping) with default
parameters.


###########################
## PROKKA-ROARY PIPELINE ##
###########################

The genetic content of each plasmid sequence was obtained using prokaryotic genome annotation
pipeline (i.e. prokka-roary pipeline). The pipeline uses the following dependencies:
   -->  https://github.com/sanger-pathogens/Roary
   -->  https://github.com/tseemann/prokka

The pipeline is implemented within prokka_roary.sh and prokka_run.sh scripts. As plasmid sequences
are numerous, these scripts were made to be used on a computer cluster as an array job submission.
If you wish to implement this pipeline:
   --> coppy individual plasmid fasta files into an empty directory on your computer cluster
   --> adjust the variables in prokka_roary.sh script pointing to the particular directory
   --> run prokka_roary.sh from the login node
         - the script will adjust the prokka_run.sh and submit it as an array job
         - the script contains several checkpoints to ensure everything is running smoothly
         - the output .gff files are organized in a directory where the plasmid sequences are and
           the script procedes with the roary implementation

Essentially, the scripts are doing the following:
   --> For each plasmid sequence the following prokka command is run:
           prokka plasmid_sequence.fasta --prefix prefix --cpus CPU_number
   --> After all .gff files have been collected and placed in a single directory the following
       roary command is run:
           roary -p 20 *.gff -f roary_output -v
   --> Roary's main purpose is the pan genome analysis with aim of extracting information about core
       and accessory genes within bacterial species/genera. As plasmids are extremely diverse with 
       very different gene content, performing complete roary analysis is likely to be unsuccessful.
       However, partial run is enough to execute Roary's blastp analysis and can thus give us:
          - clustered_proteins        - lists of clustered prokka CDSs
          - gene_presence_absence.csv - table listing association of genes and samples
       which are the only files important for further analysis.

Lastly, the CDSs identified using the pipeline above were assocaiated with Gene Ontology (GO) terms.
To achieve this, the following online resources were used:
       https://www.uniprot.org/uploadlists/
       http://current.geneontology.org/ontology/external2go/hamap2go


########################
## SCORING SIMILARITY ##
########################

The scripts for scoring similarity between pairs od plasmid sequences can be found in, evidently,
scoring_similarity/ directory. The exact Jaccard index (JI) is used to score similarity between 
plasmids. BinDash (https://github.com/zhaoxiaofei/bindash) software was used to compute the JI. 
A requirement to implement the scripts is a text file with paths to individual plasmid sequences 
(see example plasmid_paths.txt).

   --> split_in_batches.sh   - sorts plasmid sequences into batches for computational efficiency and 
                               makes batches_comb.txt file containing all possible combinations of 
                               batch files
   --> pairwise_dist.sh      - takes batches_comb.txt and computes exact pairwise JI
                               The script runs the following BinDash commands:
                  SKETCHING batches:
                  bindash sketch --listfname="$batch" --outfname=exactJ/batch_${counter}_full_sketch 
                  --nthreads=20 --minhashtype=-1
                  COMPUTING JI between two pairs of batches:
                  $bindash dist  exactJ/${batch_1}_full_sketch exactJ/${batch_2}_full_sketch 
                  --nthreads=20 --mthres=1e9 >> exactJ_dist.tsv 

   --> convert_to_DM.py      - converts BinDash output into a distance matrix.
                               The resulting distance matrix is essentially the plasmid network.
                               NOTE: in our research paper, plasmids with less than 100 k-mers in
                               common were considered to have JI=0. This is not implemented in this
                               script.

####################
## OSLOM ANALYSIS ##
####################

The main pupose of using OSLOM (http://www.oslom.org) is to identify significant cliques of plasmids.
In order to perform OSLOM analysis, the JI distance matrix was converted to pajek network format.
A standard OSLOM run would look like this:
      oslom_undir -r 250 -hr 0 -t 0.05 -singlet -seed 42 -cp 0 -f plasmid_net.net -w

However, running OSLOM on a full plasmid network yields inconsistent results with few cliques.
Thus, a 0.3 JI threshold was implemented to thin out the plasmid network. 0.3 JI threshold was chosen
after running several consecutive OSLOM runs over a range of thresholds. A script for iteratively 
runnning OSLOM analysis can be found in: oslom_analysis/OSLOM_iter.sh

OSLOM output contains the "tp1" file which contains information about assigned plasmid 
communities/cliques. This file was used in subsequent analyses performed in R (3.6.1).




