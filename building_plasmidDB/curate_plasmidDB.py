# Dependencies
from Bio import SeqIO
import os
import re
import math

# Variables
download_directory = "/path/to/download/directory/plasmidDB"
gbff = "plasmidDB.gbff"
fasta = "plasmidDB.fasta"
meta = "plasmidDB_meta.tsv"


# ASCII histogram
def ascii_hist(names, values, max_value=None):
    if len(names) != len(values):
        print("Error: Names and Values don't match!")
        return -1

    if max_value is None:
        max_value = max(values)

    # Normalizing values
    values_norm = [(value / max_value)*40 for value in values]
    scale=[str(int(x*max_value/5)) for x in [0,1,2,3,4,5]]

    print("\tCount:")
    print(" ",scale[0],scale[1],scale[2],scale[3],scale[4],scale[5],sep="\t")
    print("\t|-------|-------|-------|-------|-------|")
    print("\t|")
    for i in range(len(names)):
        print(names[i], "\t|", "]" * round(values_norm[i]), "  ", values[i],sep="")
        print("\t|")


## MAIN CODE:
print()
print('>>> Curating plasmid dataset <<<')
print("RefSeq Plasmid release tends to have impurities such as:")
print("  --> gene sequences")
print("  --> partial plasmid sequences")
print("  --> whole genomes")
print("  --> poorly annotated sequences")
print("  ...etc")
print()
print("After filtering there is usually reduction in short and extremely long sequences")
print("The script will remove all entries which:")
print("  --> do not have proper title annotation")
print("  --> are not from Bacterial domain")

os.chdir(download_directory)
title_regex = 'plasmid.*complete sequence'

with open(meta, "r") as f:
    meta_table = f.readlines()

# Locations
table_header = meta_table[0].split(sep="\t")
descript_loc = table_header.index("Description")
domain_loc = table_header.index("Domain")
len_loc = table_header.index("Length")
id_loc = table_header.index("Accession")

# Other variables
lengths_before = []
lengths_after = []
keep = []
meta_table_filtered = open("plasmidDB_meta_filtered.tsv", 'w')
meta_table_filtered.write("\t".join(table_header))

# Filtering
for line in meta_table[1:]:
    entry = line.split(sep="\t")
    lengths_before.append(int(entry[len_loc]))
    if re.search(title_regex, entry[descript_loc]):
        if entry[domain_loc] == "Bacteria":
            meta_table_filtered.write(line)
            lengths_after.append(int(entry[len_loc]))
            keep.append(entry[id_loc])

meta_table_filtered.close()


# Length distributions
powers_max = int(math.log10(max(lengths_before)))
powers_min = int(math.log10(min(lengths_before)))
powers = range(powers_min+1, powers_max+1)

print()
print()
print("Plasmid length distribution BEFORE filtering:")
names = []
values = []
for power in powers:
    names.append("10^"+str(power)+" bp")
    count = sum(( i >= 5*(10**(power-1)) ) & ( i < 0.5*(10**(power+1)) ) for i in lengths_before)
    values.append(count)
ascii_hist(names, values, max_value=math.ceil(max(values)/1000)*1000)

print()
print("Plasmid length distribution AFTER filtering:")
names = []
values = []
for power in powers:
    names.append("10^"+str(power)+" bp")
    count = sum(( i >= 5*(10**(power-1)) ) & ( i < 0.5*(10**(power+1)) ) for i in lengths_after)
    values.append(count)
ascii_hist(names, values, max_value=math.ceil(max(values)/1000)*1000)


# Cleaning up fasta file and GenBank file
print()
print()
print("Cleaning up plasmidDB.fasta")
with open("plasmidDB_filtered.fasta", "w") as f:
    for seq in SeqIO.parse(open(fasta,"r"), "fasta"):
        if seq.id in keep:
            SeqIO.write([seq], f, "fasta")


print("Cleaning up plasmidDB.gbff")
with open("plasmidDB_filtered.gbff", "w") as f:
    for gb_record in SeqIO.parse(open(gbff,"r"), "genbank"):
        if gb_record.id in keep:
            SeqIO.write([gb_record], f, "genbank")

print("Removing unfiltered plasmidDB files")
os.system("mv plasmidDB_filtered.fasta " + fasta)
os.system("mv plasmidDB_filtered.gbff " + gbff)
os.system("mv plasmidDB_meta_filtered.tsv " + meta)
print("DONE!")
