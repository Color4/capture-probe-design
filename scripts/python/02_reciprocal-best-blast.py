import argparse
import subprocess
import sys
import os
import time
import pandas as pd
create_directory = __import__("01_retrieve-exons").create_directory


def exit_status(exit_code):
    """Prints success/fail of subprocess based on exit code

    Args:
        exit_code: ('int'): Exit code of subprocess

    Returns:
        None: Print to stdout.
    """
    if exit_code == 0:
        print "\nBlast database successfully created/queried \n"
    else:
        print "\nThere was an error with the database creation/querying.\n Double check the command and try again."


def makebBlastDatabase(directory):
    """Creates BLAST database from fasta file.

    Args:
        directory ('str'): Path to directory to fasta file

    Returns"
        SeqsFasta ('str'): Full path to fasta file from which database was created.
        Blastdb ('str'): Full path to BLAST database
    """

    # Get name of fasta with sequences
    SeqsFasta = [file for file in os.listdir(directory) if file.endswith(".fasta")]

    if len(SeqsFasta) == 1:
        SeqsFasta = directory + SeqsFasta[0]
    else:
        sys.exit("There is more than one fasta file in the directory. Exiting!")

    # Name of database
    Blast_db = directory + SeqsFasta.split('/')[-1].replace(".fasta", "_db")

    print(SeqsFasta, Blast_db)

    # Command to make BLAST database
    MakeBlastDb_Command = "makeblastdb -in " + SeqsFasta + " -out " + Blast_db + " -parse_seqids" + " -dbtype nucl"

    exit_code = subprocess.call(MakeBlastDb_Command, shell=True)
    exit_status(exit_code)

    return(SeqsFasta, Blast_db)


def reciprocal_blast(transSeqs, exonSeq, outfile_TransEx, outfile_ExTrans):
    """Performs reciprocal best BLAST

    Args:
        transSeqs (:obj:'Tuple' of :obj:'str'): Tuple containing path to transcriptome fasta file [0] and BLAST database [1].
        exonSeq (:obj:'Tuple' of :obj:'str'): Tuple containing path to filtered exon fasta file [0] and BLAST database [1].
        outfile_TransEx ('str'): Full path to text file containing transcriptome to exon BLAST results.
        outfile_ExTrans ('str'): Full path to text file containing exon to transcriptome BLAST results.

    Returns:
        None: Writes files to disk.
    """

    # Nested list to be expanded for input into BLAST command.
    # First list is input for BLAST of transcriptome to exons
    # Second list in input for BLAST of exons to transcriptome
    Blast_inputs = [[transSeqs[0], exonSeqs[1], outfile_TransEx],
                    [exonSeqs[0], transSeqs[1], outfile_ExTrans]]

    # print(Blast_inputs)

    # Run BLAST command twice: Once for Transcriptome to cds and once
    # for cds to transcriptome. Results written to disk.
    for i in range(2):
        Blast_Command = "blastn -query {0} -db {1} -out {2} -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend qlen sstart send slen evalue bitscore qcovs qcovhsp' -num_threads 2 -evalue 1e-10 -max_target_seqs 1 -max_hsps 1".format(*Blast_inputs[i])
        print(Blast_Command)
        exit_code = subprocess.call(Blast_Command, shell=True)
        exit_status(exit_code)

    return(outfile_TransEx, outfile_ExTrans)

def retrieve_reciprocal_hits(TransEx_df, ExTrans_df):
    """Extract reciprocal hits from reciprocal BLAST results

    Args:
        TransEx_df (:obj:'Pandas datafram'): Pandas dataframe with transcriptome to exons BLAST results
        ExTrans_df (:obj:'Panda datafram'): Pandas dataframe with exon to transcriptome BLAST results

    Returns:
        df_reciprocal_hits (:obj:'Pandas dataframe'): Pandas dataframe with reciprocal hits as rows. Reciprocal hits extacted from exon to transcriptome BLAST results.
    """
    IDs_Trans2Ex = {}  # Qseqids from BLAST of transcriptome to exons
    IDs_Ex2Trans = {}  # Qseqids from BLAST of exons to transcriptomes
    reciprocal_hits = {}  # Store reciprocal hits

    print "Extracting transcriptome to exon BLAST hits"
    time.sleep(1)

    for index, row in TransEx_df.iterrows():
        qseqid_Trans2Ex = row['qseqid']
        sseqid_Trans2Ex = row['sseqid']
        IDs_Trans2Ex[qseqid_Trans2Ex] = sseqid_Trans2Ex

    print "Extracting exon to transcriptome BLAST hits"
    time.sleep(1)

    for index, row in ExTrans_df.iterrows():
        qseqid_Ex2Trans = row['qseqid']
        sseqid_Ex2Trans = row['sseqid']
        IDs_Ex2Trans[qseqid_Ex2Trans] = sseqid_Ex2Trans

    num_IDs_Trans2Ex = len(IDs_Trans2Ex)
    num_IDs_Ex2Trans = len(IDs_Ex2Trans)

    print "There are {0} hits of transcript sequences to reference exons".format(num_IDs_Trans2Ex)
    time.sleep(2)
    print "There are {0} hits of reference exons to transcripts sequences".format(num_IDs_Ex2Trans)
    time.sleep(2)

    print "Retrieving reciprocal hits"
    time.sleep(1)

    # Create sets from BLAST qseqids and sseqids
    # Second set is reversed to facilitate checking intersection
    # (i.e. reciprocal hits) between two sets of BLAST hits.
    Set_IDs_Ex2Trans = set(IDs_Ex2Trans.items())
    Set_IDs_Trans2Ex = [(t[1], t[0]) for t in set(IDs_Trans2Ex.items())]

    for (key, value) in Set_IDs_Ex2Trans.intersection(Set_IDs_Trans2Ex):
        # print '{0}: {1} is present in both IDs_Trans2Ex and IDs_Ex2Trans'.format(key, value)
        reciprocal_hits[key] = value

    num_RecipHits = len(reciprocal_hits)
    print "There are {0} Reciprocal hits".format(num_RecipHits)
    time.sleep(2)

    print "Creating dataframe with reciprocal hits"
    time.sleep(1)

    df_reciprocal_hits = ExTrans_df[ExTrans_df['qseqid'].isin(reciprocal_hits.keys())]

    return(df_reciprocal_hits)


def filter_recipHits(df_reciprocal_hits, qcovs, pident):
    """Filter recirpocal hits dataframe

    Ensures that each gene is represented only by a single exon, that the alignment length is greater than 200 bp, that the percent identity is greater than --pident (80% default) and that the reference exon is covered at least --qcovs by the target transcript (50% default).

    Args:
        df_reciprocal_hits (:obj:'Pandas dataframe'): Pandas dataframe with reciprocal hits as rows.

    Returns:
        df_reciprocal_hits_filtered (:obj:'Pandas dataframe'): Pandas dataframe with filtered reciprocal hits as rows.
    """

    unfiltered_df_length = len(df_reciprocal_hits)

    print "Removing hits where alignment is < 200 bp"
    df_recipHits_filtered = df_reciprocal_hits[df_reciprocal_hits['length'] >= 200]
    time.sleep(1)
    print "Done removing short alignments. Removed {0} alignments.".format(unfiltered_df_length - len(df_recipHits_filtered))
    time.sleep(2)

    print "Identifying duplicate gene IDs"
    time.sleep(1)

    Gene_IDs = [GeneID.split(".")[-1] for GeneID in df_recipHits_filtered['qseqid'].tolist()]
    # print(Gene_IDs)
    # print(len(Gene_IDs))

    unique = []
    duplicates = []

    for ID in Gene_IDs:
        if ID not in unique:
            unique.append(ID)
        else:
            duplicates.append(ID)

    num_unique = len(unique)
    print "There are {0} unique gene IDs in the reciprocal hits database.".format(num_unique)
    time.sleep(2)

    df_recipHits_filtered['GeneID'] = df_recipHits_filtered['qseqid'].str.split(".").str[-1]

    frames = []

    print "Removing rows with duplicate Gene IDs"
    time.sleep(1)

    # Keep only a single alignment for each unique gene ID. First, keep the one with the highest Bitscore. If Bitscrores are equal, keep the one with the highest percent identity.
    for ID in unique:
        df = df_recipHits_filtered[df_recipHits_filtered['GeneID'] == ID]
        if len(df) > 1:
            BitScore = max(df['bitscore'])
            df = df[df['bitscore'] == BitScore]
            if len(df) > 1:
                pident = max(df['pident'])
                df = df[df['pident'] == pident]
                frames.append(df)
            else:
                frames.append(df)
        else:
            frames.append(df)

    num_unique_dfs = len(frames)

    if num_unique_dfs != num_unique:
        print "**WARNING** The number of unique dataframes to concatenate does not equal the number of unique Gene IDs. Something went wrong."

    print "A single alignment for each gene has been kept. Now concatenating dataframes."
    df_recipHits_filtered = pd.concat(frames)
    time.sleep(1)
    print "Done concatenating dataframes."
    time.sleep(1)

    filtered_df_length = len(df_recipHits_filtered)

    print "Removing rows that do not meet qcovs and pident threshholds."
    df_recipHits_filtered_qcovPident = df_recipHits_filtered[(df_recipHits_filtered['pident'] >= pident) & (df_recipHits_filtered['qcovs'] >= qcovs)]
    time.sleep(1)

    total_length = sum([row['length'] for index, row in df_recipHits_filtered_qcovPident.iterrows()])

    print "Done removing rows that do not meet qcovs and pident threshholds. {0} alignments removed. There are {1} unique alignments remaining, covering approximately {2} Mbp in length".format(filtered_df_length - len(df_recipHits_filtered_qcovPident), len(df_recipHits_filtered_qcovPident), (total_length / 10**6.0))
    time.sleep(3)

    return(df_recipHits_filtered_qcovPident)

if __name__ == "__main__":

    # Define command line arguments
    parser = argparse.ArgumentParser(prog="Perform reciprocal best hit blast between transcriptome of focal species and filtered exons from closely related reference species.")
    parser.add_argument('trans_dir', help="Path to directory containing transcriptome fasta file", type=str)
    parser.add_argument('exon_dir', help="Path to directory containing filtered exon fasta file", type=str)
    parser.add_argument('-qcovs', nargs="?", default=50, type=int, help="Minimum coverage of transcript alignment to reference exon to be kept in final, filtered reciprocal hits dataframe. Default 50.")
    parser.add_argument('-pident', nargs="?", default=80, type=int, help="Minimum percent identity between transcript and reference exon alignment to be kept in final, filtered reciprocal hits dataframe. Default 80.")
    args = vars(parser.parse_args())

    # Create output directories
    print("Creating output directories")
    raw_data_dir = "../../data-raw/02_reciprocal-best-blast/"
    clean_data_dir = "../../data-clean/02_reciprocal-best-blast/"
    create_directory(raw_data_dir)
    create_directory(clean_data_dir)
    time.sleep(1)

    # Retrieve arguments passed from command line
    trans_dir = args['trans_dir']
    exon_dir = args['exon_dir']
    qcovs = args['qcovs']
    pident = args['pident']

    # Make BLAST databases for transcriptome and exons
    print "Creating transcriptome BLAST database!"
    time.sleep(1)
    transSeqs = makebBlastDatabase(trans_dir)
    print "Done creating transcriptome BLAST database"
    time.sleep(1)

    print "Creating exons BLAST database!"
    time.sleep(1)
    exonSeqs = makebBlastDatabase(exon_dir)
    print "Done creating exons BLAST database"
    time.sleep(1)

    # Name reciprocal BLAST output files
    outfile_TransEx = raw_data_dir + "transcriptome_to_exons.txt"
    outfile_ExTrans = raw_data_dir + "exons_to_transcriptome.txt"

    # Perform reciprocal best hit blast
    print "Performing reciprocal best BLAST now!"
    time.sleep(1)
    reciprocal_blast(transSeqs, exonSeqs, outfile_TransEx, outfile_ExTrans)
    print "Done with reciprocal BLAST. Results are in {0}".format(raw_data_dir)

    # # Load reciprocal blast results into Pandas dataframes
    print "Loading recirpocal BLAST results intro dataframes"
    time.sleep(1)

    colnames = ["qseqid", "sseqid", "pident", "length",
                "mismatch", "gapopen", "qstart", "qend",
                "qlen", "sstart", "send", "slen", "evalue",
                "bitscore", "qcovs", "qcovhsp"]

    transcriptome_to_exons = pd.read_table(outfile_TransEx,
                                           sep="\t",
                                           names=colnames)

    exons_to_transcriptome = pd.read_table(outfile_ExTrans,
                                           sep="\t",
                                           names=colnames)

    print "Done loading recirpocal BLAST results intro dataframes"
    time.sleep(1)

    # Retrieve dictionary with reciprocal hits.
    reciprocal_hits = retrieve_reciprocal_hits(transcriptome_to_exons, exons_to_transcriptome)

    recipHitsDF_out = raw_data_dir + "reciprocalHits_All.txt"

    print "Writing reciprocal hits datafram to disk. Dataframe is located in {0}".format(recipHitsDF_out)
    reciprocal_hits.to_csv(recipHitsDF_out, sep="\t", index = False, header=colnames)
    time.sleep(1)

    recipHits_filtered = filter_recipHits(reciprocal_hits, qcovs, pident)

    recipHitsFilteredDF_out = clean_data_dir + "reciprocalHits_Filtered.txt"
    print "Writing filtered reciprocal hits dataframe to disk. The dataframe is located in {0}".format(recipHitsFilteredDF_out)

    recipHits_filtered.to_csv(recipHitsFilteredDF_out, sep="\t", index = False, header=colnames + ["GeneID"])
