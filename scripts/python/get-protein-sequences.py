import pandas as pd
# import sys, time
from Bio import Entrez
Entrez.email = 'james.santangelo37@gmail.com'


def get_protein_sequences(accession_numbers, species):
    """Retrieve protein sequences from NCBI

    Based on filtered balstx hits, retrieves protein sequences of to subject sequences from NCBI and writes sequences to fasta file.

    Args:
      accession_numbers (:obj:'list' of 'str'): List of subject sequence accession number that can be queried against the NCBI Entrez database.

    Returns:
      None: Writes fasta file to disk.
    """

    # Open fasta file.
    protein_sequences = open("../../data-raw/011_retrieve-proteins/{0}/M-truncatula_ProteinSeqs.fasta".format(species), 'w')

    # Retrieve amino acid sequence for each subject sequence in list and write to fasta
    for ID in accession_numbers:
        handle = Entrez.efetch(db="protein", id=ID, rettype="fasta", retmode="text")
        time.sleep(3)
        record = handle.read()
        protein_sequences.write(record)
    protein_sequences.close()


# Which species is being processed? Entered as command line argument.
species = sys.argv[1]

# Load subject ID's (i.e. accession numbers) from dataframe with filtered blastx results
ProteinBlast_filtered = pd.read_csv("../../data-clean/010_blastx/{0}/Nagy_transcriptome_ProteinBlast_Filtered.csv".format(species), usecols = ["sseqid"])

# Extract list of subject accession numbers to query against NCBI Entrez
accession_numbers = list(ProteinBlast_filtered)

# Retrieve sequences from NCBI and write to fasta
get_protein_sequences(accession_numbers, species)
