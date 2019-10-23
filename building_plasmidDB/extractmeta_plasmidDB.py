# Dependencies
from Bio import SeqIO
from ete3 import NCBITaxa
import os, ssl
# Bypassing errors with SSL on Mac (for Linux remove this):
if (not os.environ.get('PYTHONHTTPSVERIFY', '') and
    getattr(ssl, '_create_unverified_context', None)):
    ssl._create_default_https_context = ssl._create_unverified_context

# Variables
download_directory = "/path/to/download/directory/plasmidDB"
gbff = "plasmidDB.gbff"
fasta = "plasmidDB.fasta"


# MAIN CODE:
print()
print("Extracting metadata for plasmidDB")
print("Reading GenBank file and storing valuable info in tabulated format")
os.chdir(download_directory)

meta = open('plasmidDB_meta.tsv','w')
meta.write("Accession\tDescription\tPlasmid\tLength\tOrganism\tSpecies\tGenus\tFamily\tOrder\tClass\t"
            "Phylum\tDomain\tCollection_date\tGeo_location\tIsolation_source\tNote\n")

print(' --> Setting up a taxonomy database')
ncbi = NCBITaxa()
#ncbi.update_taxonomy_database()

print(' --> Creating metadata table')
for gb_record in SeqIO.parse(open(gbff,"r"), "genbank"):
    record_meta = ["NA"]*16
    record_meta[0] = gb_record.id                                                   # Accession
    record_meta[1] = gb_record.description                                          # Description
    if "plasmid" in gb_record.features[0].qualifiers.keys():
        record_meta[2] = gb_record.features[0].qualifiers["plasmid"][0]             # Plasmid
    record_meta[3] = str(len(gb_record.seq))                                        # Length

    organism = gb_record.annotations["organism"]                                    # Taxonomy
    record_meta[4] = organism
    taxid = ncbi.get_name_translator([organism])
    if organism in taxid.keys():
        lineage = ncbi.get_lineage(taxid[organism][0])
        ranks = ncbi.get_rank(lineage)
        for rank in ranks:
            if ranks[rank] == 'species':
                record_meta[5] = ncbi.get_taxid_translator([rank])[rank]
            if ranks[rank] == 'genus':
                record_meta[6] = ncbi.get_taxid_translator([rank])[rank]
            if ranks[rank] == 'family':
                 record_meta[7] = ncbi.get_taxid_translator([rank])[rank]
            if ranks[rank] == 'order':
                record_meta[8] = ncbi.get_taxid_translator([rank])[rank]
            if ranks[rank] == 'class':
                record_meta[9] = ncbi.get_taxid_translator([rank])[rank]
            if ranks[rank] == 'phylum':
                record_meta[10] = ncbi.get_taxid_translator([rank])[rank]
            if ranks[rank] == 'superkingdom':
                record_meta[11] = ncbi.get_taxid_translator([rank])[rank]

    if "collection_date" in gb_record.features[0].qualifiers.keys():                    # Collection date
        record_meta[12] = gb_record.features[0].qualifiers["collection_date"][0]
    if "lat_lon" in gb_record.features[0].qualifiers.keys():                            # Geo location
        record_meta[13] = gb_record.features[0].qualifiers["lat_lon"][0]
    elif "country" in gb_record.features[0].qualifiers.keys():
        record_meta[13] = gb_record.features[0].qualifiers["country"][0]
    if "isolation_source" in gb_record.features[0].qualifiers.keys():                   # Isolation source
        record_meta[14] = gb_record.features[0].qualifiers["isolation_source"][0]
    if "note" in gb_record.features[0].qualifiers.keys():                               # Note
        record_meta[15] = gb_record.features[0].qualifiers["note"][0]

    meta.write("\t".join(record_meta) + "\n")                                           # Appending to the table

meta.close()
