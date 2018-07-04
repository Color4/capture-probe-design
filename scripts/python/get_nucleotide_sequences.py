import pandas as pd
import sys, time
from Bio import Entrez, SeqIO
from tqdm import tqdm
from urllib2 import HTTPError, URLError

Entrez.email = 'james.santangelo37@gmail.com'


def get_nucleotide_sequences(Accession_numbers, fout):
    """Download nucleotide sequences from GenBank

    Args:
        Accession_numbers (:obj:'list' of 'str'): Accession numbers for sequences to download

    Returns:
        None: Writes sequences to disk in FASTA format
    """

    # Open fasta file to which nucleotide sequences should be written
    nucleotide_sequences = open(fout, 'w')

    # Retrieve amino acid sequence for each subject sequence in list and write to fasta
    for num in tqdm(Accession_numbers):
        print(num)  # Print sequence being downloaded to command line

        # Exception handling. Unlimited retries if denied access
        while True:
            try:
                handle = Entrez.efetch(db="nucleotide", id=num, rettype="fasta", retmode="text")
            except HTTPError as err:
                print("Received error from server %s" % err)  # Which error was returned? (e.g. 404, 500).
                time.sleep(3)  # Wait before retrying
                continue  # Retry same Accession number
            except URLError as err:
                print("Received error from server %s" % err)
                time.sleep(3)
                continue
            break  # Break out of loop if sequence successfully retrieved
        record = SeqIO.read(handle, "fasta")
        nucleotide_sequences.write(">{0}\n{1}\n".format(record.description,record.seq))
    nucleotide_sequences.close()


if __name__ == "__main__":
    colnames = ["Accession_number", "Description"]
    HCN_genes = pd.read_table("../../data-raw/T-repens_HCN_AccessionNumbers.txt", sep = " ", quotechar='"', names = colnames)
    Accession_numbers = list(HCN_genes["Accession_number"])
    fout = "../../data-clean/T-repens_HCN_genes.fasta"
    get_nucleotide_sequences(Accession_numbers, fout)
