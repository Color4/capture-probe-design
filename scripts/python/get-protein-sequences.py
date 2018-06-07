import pandas as pd
import sys
from Bio import Entrez
Entrez.email = 'james.santangelo37@gmail.com'

def filter_proteinBlast(blast_results):
    filtered = blast_results[(blast_results['evalue'] <= 1.000000e-10)]
    filtered = filtered.sort_values('length', ascending=False).drop_duplicates(['qseqid'])
    return(filtered)


def get_protein_sequences(ID_list):
    protein_sequences = open('protein-sequences/{0}_ProteinSeqs.fasta'.format(sys.argv[1]), 'w')
    for ID in ID_list:
        handle = Entrez.efetch(db="protein", id=ID, rettype="fasta", retmode="text")
        record = handle.read()
        protein_sequences.write(record)
    protein_sequences.close()


colnames = ["qseqid", "sseqid", "pident",
            "length", "mismatch", "gapopen",
            "qstart", "qend", "qlen", "sstart",
            "send", "slen", "evalue", "bitscore",
            "qcovs", "qcovhsp"]

ProteinBlast = pd.read_table("protein-blastx/Nagy_transcriptome_ProteinBlast_{0}.csv".format(sys.argv[1]), names=colnames)

ProteinBlast_filtered = filter_proteinBlast(ProteinBlast)

Protein_FilteredSseqids = list(ProteinBlast_filtered["sseqid"])

get_protein_sequences(Protein_FilteredSseqids)
