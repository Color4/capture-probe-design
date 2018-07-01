import pandas as pd
import sys, time
from Bio import Entrez, SeqIO
Entrez.email = 'james.santangelo37@gmail.com'


def get_protein_sequences(qseq_sseq_IDs, species):
    """Retrieve protein sequences from NCBI

    Based on filtered balstx hits, retrieves protein sequences of to subject sequences from NCBI and writes sequences to fasta file.

    Args:
      accession_numbers (:obj:'list' of 'str'): List of subject sequence accession number that can be queried against the NCBI Entrez database.

    Returns:
      None: Writes fasta file to disk.
    """

    # Open fasta file.
    protein_sequences = open("../../data-raw/021_retrieve-proteins/{0}/NagyTranscriptHits_ProteinSeqs.fasta".format(species), 'w')

    # Retrieve amino acid sequence for each subject sequence in list and write to fasta
    for pair in qseq_sseq_IDs:
        qseqid = pair[0]
        sseqid = pair[1]
        handle = Entrez.efetch(db="protein", id=sseqid, rettype="fasta", retmode="text")
        record = SeqIO.read(handle, "fasta")
        protein_sequences.write(">{0};{1}\n{2}\n".format(qseqid,record.description,record.seq))
    protein_sequences.close()


# Which species is being processed? Entered as command line argument.
species = sys.argv[1]

# Load subject ID's (i.e. accession numbers) from dataframe with filtered blastx results
ProteinBlast_filtered = pd.read_csv("../../data-clean/020_blastx/{0}/Nagy_transcriptome_ProteinBlast_Filtered.csv".format(species), usecols = ["qseqid", "sseqid"])

# Extract list of subject accession numbers to query against NCBI Entrez
qseq_sseq_IDs = ProteinBlast_filtered.values.tolist()

# Retrieve sequences from NCBI and write to fasta
get_protein_sequences(qseq_sseq_IDs, species)
