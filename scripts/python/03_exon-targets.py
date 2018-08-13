import argparse
import re
import sys
import time
import os
import pandas as pd
from Bio import SeqIO
create_directory = __import__("01_retrieve-exons").create_directory


def extract_exons(recipHits_df, transcriptome):
    """Extracts exon sequences from transcripts based on Reciprocal best hits BLAST alignment

    Args:
        recipHits_df (:obj:'Pandas dataframe'): Pandas dataframe with final, filtered reciprocal best hits between transcripts and reference exons.
        transcriptome (:obj:'FASTA'): Fasta file with transcript sequences.

    Retrurns:
        None: Writes fasta to disk
    """

    sequence_counter = 0

    with open(ExonOut_raw, "w") as f:
        for index, row in recipHits_df.iterrows():

            pratense_exon = row['qseqid']
            transcript = row['sseqid']
            start = row['sstart']
            end = row['send']

            Seq_record = transcriptome[transcript]

            # Sequences on forward (+) strand
            if start < end:
                strand = "+"
                sequence = str(Seq_record.seq[start - 1: end])
                length = len(sequence)
                description = Seq_record.description
                pattern = re.compile('\A[a-z]+-[0-9]+\s(.+)\sLength')
                match = pattern.match(description)
                description = match.group(1)
                f.write(">Transcript_id:{0};Strand:{1};PratenseExon_id:{2};Length:{3};Transcript_annot:{4}\n{5}\n".format(transcript, strand, pratense_exon, length, description, sequence))
                sequence_counter += 1

            # Sequences on reverse (-) strand. Reverse complement required.
            elif end < start:
                strand = "-"
                sequence = Seq_record.seq[end - 1: start]
                sequence = str(sequence.reverse_complement())
                length = len(sequence)
                description = Seq_record.description
                pattern = re.compile('\A[a-z]+-[0-9]+\s(.+)\sLength')
                match = pattern.match(description)
                description = match.group(1)
                f.write(">Transcript_id:{0};Strand:{1};PratenseExon_id:{2};Length:{3};Transcript_annot:{4}\n{5}\n".format(transcript, strand, pratense_exon, length, description, sequence))
                sequence_counter += 1

    print "Wrote a total of {0} sequences to file.".format(sequence_counter)
    if sequence_counter != len(recipHits_df):
        print "**WARNING** The number of sequences written to disk does not match the number of rows in the reciprocal hits dataframe. Some sequences were not written."
    time.sleep(2)


def filter_GC_content(exon_fasta, filtered_exons_out):
    """Filters out exons with GC content <30% and > 70%

    Args:
        exon_fasta (:obj:'FASTA'): Fasta file containing the unfiltered exons retrieved from the transcriptome

    Returns:
        None: Writes fasta to disk.
    """

    # Load exons
    exon_targets = list(SeqIO.parse(exon_fasta, "fasta"))

    count_good = 0
    total_length_good = 0

    with open(filtered_exons_out, "w") as f:
        for record in exon_targets:
            sequence = record.seq
            GC_count = sequence.count("G") + sequence.count("C")
            length = len(sequence)
            GC_percent = (float(GC_count) / length) * 100

            if GC_percent >= 30 and GC_percent <= 70:
                f.write(">{0}\n{1}\n".format(record.description, sequence))
                count_good += 1
                total_length_good += length
            else:
                pass

    print "A total of {0} exons were removed from the fasta file due to low or high GC content. There are thus {1} exons remaining, covering a total of {2} Mbp of sequence.".format(len(exon_targets) - count_good, count_good, total_length_good / 10**6.0)
    time.sleep(3)


if __name__ == "__main__":

        # Define command line arguments
    parser = argparse.ArgumentParser(prog="Based on reciprocal best hit blast, retrieve exon sequences from transcriptome using start and end positions from BLAST alignment.")
    parser.add_argument('recipHits_dir', help="Path to directory containing filtered reciprocal hits text file", type=str)
    parser.add_argument('trans_dir', help="Path to directory containing transcriptome fasta file", type=str)
    args = vars(parser.parse_args())

    # Retrieve command line arguments
    recipHits_dir = args['recipHits_dir']
    trans_dir = args['trans_dir']

    # Create output directories
    print("Creating output directories")
    raw_data_dir = "../../data-raw/03_exon-targets/"
    clean_data_dir = "../../data-clean/03_exon-targets/"
    create_directory(raw_data_dir)
    create_directory(clean_data_dir
    time.sleep(1)

    # Path to raw exon sequence fasta file
    ExonOut_raw = raw_data_dir + "exon_targets_raw.fasta"

    # Load text file with reciprocal hits
    print "Loading reciprocal hits dataframe"
    time.sleep(1)
    # colnames = ["qseqid", "sseqid", "pident", "length",
    #             "mismatch", "gapopen", "qstart", "qend",
    #             "qlen", "sstart", "send", "slen", "evalue",
    #             "bitscore", "qcovs", "qcovhsp", "GeneID"]
    recipHits_df = [f for f in os.listdir(recipHits_dir) if f.endswith('.txt')]
    if len(recipHits_df) == 1:
        recipHits_df = recipHits_dir + recipHits_df[0]
    else:
        sys.exit("There is more than one text file in the directory. Exiting!")
    recipHits_df = pd.read_table(recipHits_df,
                                 sep="\t")
    # print recipHits_df
    print "Done loading reciprocal hits dataframe"
    time.sleep(1)

    print "Loading transcriptome into Python dictionary"
    time.sleep(1)
    transcriptome = [f for f in os.listdir(trans_dir) if f.endswith('.fasta')]
    if len(transcriptome) == 1:
        transcriptome = trans_dir + transcriptome[0]
    else:
        sys.exit("There is more than one text file in the directory. Exiting!")
    transcriptome = SeqIO.to_dict(SeqIO.parse(transcriptome, "fasta"))
    print "Done loading transcriptome"
    # print(transcriptome.keys())
    time.sleep(1)

    print "Extracting exons now!"
    time.sleep(1)

    sequence_counter = 0
    extract_exons(recipHits_df, transcriptome)

    print "Filtering exons with low or high GC content."
    filtered_exons_out = clean_data_dir + "exon_targets_filtered.fasta"

    filter_GC_content(ExonOut_raw, filtered_exons_out)
