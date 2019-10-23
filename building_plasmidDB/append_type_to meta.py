import re
from Bio import SeqIO


# Overlap function
def getOverlap(a, b):
    return min(a[1], b[1]) - max(a[0], b[0])

# Variables
PlasmidFinderdb = "/path/to/PlasmidFinderdb/PlasmidFinderdb.fasta"
meta = "plasmidDB_meta.tsv"
PlasmidFinder_hits = "/path/to/PlasmidFinder_hits"

print()
print('>>> Appending Plasmid Finder Results to meta <<<')

# Reading meta
with open(meta, "r") as f:
    meta_table = f.readlines()
table_header = meta_table[0].split(sep="\t")
id_loc = table_header.index("Accession")

# Reading PlasmidFinder_hits
with open(PlasmidFinder_hits, "r") as f:
    PlasmidFinder_hits = f.readlines()

# Making dictionary of PlasmidFinder gene lengths
PlasmidFinder_lens = dict()
for gene in SeqIO.parse(open(PlasmidFinderdb,"r"), "fasta"):
    PlasmidFinder_lens[gene.id] = len(gene.seq)

# Appending to meta
new_meta = open("plasmidDB_meta_new.tsv", 'w')
table_header[len(table_header) - 1] = table_header[len(table_header) - 1].rstrip()
table_header.append("Plasmid_type")
table_header.append("Plasmid_type_concise\n")
new_meta.write("\t".join(table_header))

for line in meta_table[1:]:
    entry = line.split(sep="\t")

    plasmid_type = list()
    plasmid_type_concise = list()
    best_hits = list()

    # Getting significant hits
    for hit in PlasmidFinder_hits:
        if re.match(entry[id_loc], hit):
            hit_split = hit.split(sep="\t")
            gene = hit_split[1].split(sep="_")[0]
            gene_cov = abs(int(hit_split[5]) - int(hit_split[6]) + 1) / PlasmidFinder_lens[hit_split[1]]
            gene_loc = [int(hit_split[3]), int(hit_split[4])]
            gene_pident = float(hit_split[2])
            if gene_cov < 0.95:
                continue
            elif gene_pident < 0.95:
                continue
            else:
                best_hits.append((gene,gene_cov,gene_loc,gene_pident))

    # Sorting significant hits
    best_hits = sorted(best_hits, reverse=True, key=lambda best_hits: best_hits[3])
    while len(best_hits) > 0:
        hit = best_hits[0]
        plasmid_type.append(hit[0])
        to_remove = list()
        for hit2 in best_hits:
            if getOverlap(hit[2], hit2[2]) >= -100:
                to_remove.append(best_hits.index(hit2))
        for i in sorted(to_remove, reverse=True):
            del best_hits[i]

    plasmid_type = list(set(plasmid_type))
    for gene in plasmid_type:
        if gene.split(sep="(")[0] not in plasmid_type_concise:
            plasmid_type_concise.append(gene.split(sep="(")[0])

    if len(plasmid_type) == 0:
        plasmid_type = 'NA'
        plasmid_type_concise = 'NA'
    else:
        plasmid_type = "|".join(plasmid_type)
        plasmid_type_concise = "|".join(plasmid_type_concise)

    # Appending to meta table
    entry[len(entry) - 1] = entry[len(entry) - 1].replace("\n", '')
    entry.append(plasmid_type)
    entry.append(plasmid_type_concise + "\n")
    new_meta.write("\t".join(entry))

new_meta.close()
os.system("mv plasmidDB_meta_new.tsv " + meta)

print("Appending compleated!")
