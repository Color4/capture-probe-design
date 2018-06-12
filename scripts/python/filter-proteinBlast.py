import pandas as pd
import sys


def filter_proteinBlast(blast_results):
    """Filters a .csv file of blastx results outputted using -outfmt 6

    Args
        blast_results (:obj:'DataFrame'): pandas DataFrame genererated by loading a .csv file output from blastx using the -outfmt 6 format.

    Returns:
        filtered (:obj:'DataFrame'): pandas DataFrame of filtered blastx results that can be written to csv.
    """

    # Keep only strong hits
    filtered = blast_results[(blast_results['evalue'] <= 1.000000e-10)]

    # If query sequence has multiple hits, keep only the longest alignment
    filtered = filtered.sort_values('length', ascending=False).drop_duplicates(['qseqid'])

    return(filtered)


# Names of columns output from blastx -outfmt 6. Used for adding header to imported csv.
colnames = ["qseqid", "sseqid", "pident",
            "length", "mismatch", "gapopen",
            "qstart", "qend", "qlen", "sstart",
            "send", "slen", "evalue", "bitscore",
            "qcovs", "qcovhsp"]

# Which species is being processed? Entered as command line argument.
species = sys.argv[1]

# Load in blast results as pandas DataFrame.
proteinBlast_results = pd.read_table("../../data-raw/020_blastx/{0}/Nagy_transcriptome_ProteinBlast.csv".format(species), names=colnames)

# Filter blastx results
ProteinBlast_filtered = filter_proteinBlast(proteinBlast_results)

# Write filtered blastx results to clean data folder.
ProteinBlast_filtered.to_csv("../../data-clean/020_blastx/{0}/Nagy_transcriptome_ProteinBlast_Filtered.csv".format(species), header=colnames, index = False)
