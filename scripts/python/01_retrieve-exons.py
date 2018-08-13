import os
import argparse
import sys
import time

# Add paths to gffutils (only for running on HPCNODE)
sys.path.append("/home/santang3/Python-modules/lib/python2.7/site-packages/")
sys.path.append("/home/santang3/.local/lib/python2.7/site-packages/")

import gffutils
from tqdm import tqdm
from Bio import SeqIO


def write_exons(clean_db_in, exon_gff3_out, exon_fasta_out, reference_genome):
    """From a gff3 database generated using gffutils, extract only the coding sequences and writes them to a new gff3 file. Also writes coding sequences
    to fasta file

    Args:
        clean_db_in ('str'): Path to GFF3 database generated using the gffutils module. See 'create_gff3_database.py'
        cds_gff3_out ('str'): Path to GFF3 file to be written with cds annotations.
        cds_fasta_out ('str'): Path to fasta file to be written with cds sequences
        reference_genome ('str'): Path to reference genome from which coding sequences should be extracted.

    Returns:
        None: Writes gff3 and fasta to disk.
    """

    # Open gff3 database
    gff3_database = gffutils.FeatureDB(clean_db_in)

    # Iterate over exons in gff3 annotation file and write to new file.
    with open(exon_gff3_out, "w") as f:
        for exon in gff3_database.features_of_type("exon"):
            f.write(str(exon) + '\n')

    # Length of generator object iterating over cds features in database.
    # Used to print progress bar to shell.
    gen_length = len([exon for exon in gff3_database.features_of_type("exon")])

    # Iterate over exons features in database.
    # Query reference for its sequence and write sequence to fasta.

    with open(exon_fasta_out, "w") as f:
        for exon in tqdm(gff3_database.features_of_type("exon"),
                        total=gen_length):
            # print(exon)
            for parent_gene in gff3_database.parents(id=exon, featuretype="gene"):
                parent_gene_id = parent_gene.id.split("_")[-1]
            exon_id = exon.attributes['exon_id'][0]
            ID = exon_id + "." + parent_gene_id
            sequence = str(exon.sequence(reference_genome, use_strand=True))
            f.write(">{0}\n{1}\n".format(ID, sequence))


def filter_exons(path_to_exon_fasta, fasta_out):
    """Removes coding sequences with "N" and those less than
    200 bp in length

    Args:
        path_to_cds_fasta ('str'): Path to fasta file containing all coding sequences
        fasta_out ('str'): Path to which filtered fasta file should be written.

    Returns:
        None: Writes fasta file to disk.
    """

    # Load in sequence recrods as list
    sequences = list(SeqIO.parse(path_to_exon_fasta, "fasta"))

    seen = []

    # Iterate over sequence records, filter and write to disk.
    with open(fasta_out, "w")as f:
        for record in sequences:
            if record.id not in seen:
                if "N" not in record.seq and "n" not in record.seq:  # Remove sequences with "N"
                    if len(record.seq) >= 200:  # Keep only sequences > 200 bp.
                        ID = record.id
                        seq = record.seq
                        seen.append(ID)
                        f.write(">{0}\n{1}\n".format(ID, seq))


def create_directory(directory):
    """Creates directory if it does not exist

    Args:
        directory ('str'): Directory to be created

    Returns:
        None
    """
    if not os.path.exists(directory):
        os.makedirs(directory)


if __name__ == "__main__":

    # Create directories
    print("Creating output directories")
    raw_data_dir = "../../data-raw/01_reference-exons/"
    clean_data_dir = "../../data-clean/01_reference-exons/"
    create_directory(raw_data_dir)
    create_directory(clean_data_dir)
    time.sleep(1)

    # Define command line arguments
    parser = argparse.ArgumentParser(prog="From a reference geneome and gff3 annotation file, generates a fasta file with all exons greater than 200 bp in length and containing no ambiguous nucleotides.")
    parser.add_argument('genome_dir', help="Path to directory containing reference genome and gff3 annotation file", type=str)
    args = vars(parser.parse_args())

    # Retrieve arguments passed from command line
    genome_dir = args['genome_dir']

    # Retrieve gff3
    gff3_file = [f for f in os.listdir(genome_dir) if f.endswith('.gff3')]
    if len(gff3_file) == 1:
        gff3_file = genome_dir + gff3_file[0]
    else:
        sys.exit("There is more than one gff3 annotation file in the directory. Exiting!")

    # Define filename for database to be written do disk
    gff3_db = genome_dir + gff3_file.split('/')[-1].replace("gff3", "db")
    # print(gff3_file, gff3_db)

    # Define what attributes field in gff3 file should be used as ID in database
    # id_spec = {"gene": "gene_id", "mRNA": "transcript_id", "CDS": "ID", "exon": "Name", "three_prime_UTR": "Parent", "five_prime_utr": "Parent", "contig": "ID", "chromosome": "ID"}
    id_spec = "ID"

    print "Creating SQL database from reference genome's gff3 file now!"
    time.sleep(1)

    # Create database
    gffutils.create_db(gff3_file, dbfn=gff3_db,
                       force=True, keep_order=False,
                       id_spec = id_spec, verbose = True,
                       merge_strategy='merge', sort_attribute_values=False)

    print "Done creating gff3 database!"
    time.sleep(1)

    # Name raw files to be written to disk
    raw_exon_gff3 = raw_data_dir + "referenceGenome_ExonsOnly.gff3"
    raw_exon_fasta = raw_data_dir + "referenceGenome_ExonsOnly.fasta"

    # Retrieve reference genome
    reference_genome = [f for f in os.listdir(genome_dir) if f.endswith('.fa') or f.endswith('.fasta')]
    if len(reference_genome) == 1:
        reference_genome = genome_dir + reference_genome[0]
    else:
        sys.exit("There is more than one genome file in the directory. Exiting!")

    # print(raw_exon_gff3, raw_exon_fasta, reference_genome)

    print "Writing gff3 and fasta files with exons from reference genome!"
    time.sleep(1)

    write_exons(gff3_db, raw_exon_gff3, raw_exon_fasta, reference_genome)

    print "Done writing exons to raw gff3 and fasta file. Both files are located in {0}".format(raw_data_dir)
    time.sleep(1)

    print "Filtering exons now!"
    time.sleep(1)

    # Name final exons file to be written to disk
    filtered_exons_fasta = clean_data_dir + "referenceGenome_ExonsOnly_filtered.fasta"

    filter_exons(raw_exon_fasta, filtered_exons_fasta)

    print "Done filtering exons. The final filtered exons are located in {0}".format(clean_data_dir)

