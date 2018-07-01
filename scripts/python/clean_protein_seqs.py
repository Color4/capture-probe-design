import re
import sys
from Bio import SeqIO


def protein_fasta_replace_J(multi_fasta, species):
    """Replace 'J' in protein sequences with 'I'

    Args:
        multi_fasta (text file): Multi fasta file with protein sequences to modify

    Returns:
        None: Write modified fasta file to disk.
    """

    # Name of new fasta file
    new_file = "../../data-clean/021_retrieve-proteins/{0}/NagyTranscriptHits_ProteinSeqs.fasta".format(species)

    with open(new_file, "w") as f:
        for record in SeqIO.parse(multi_fasta, "fasta"):
            # Write modified record to file only if 'J' is in the sequence.
            if "J" in record.seq:
                sequence = str(record.seq)

                # Get indices of 'J' in sequence
                occurences = [m.start() for m in re.finditer("J", sequence)]

                # Convert sequence to list so index positions from above can be modified.
                recordSeq_list = list(sequence)

                # Change all instances of 'J' to 'I'
                for occurence in occurences:
                    recordSeq_list[occurence] = "I"

                # Convert list back to string and write to file.
                sequence = ''.join(recordSeq_list)
                f.write(">" + str(record.description) + "\n")
                f.write(str(sequence) + "\n")
            else:
                f.write(">" + str(record.description) + "\n")
                f.write(str(record.seq) + "\n")


# Which species' protein sequences are being modified? Input at command line
species = sys.argv[1]

# Patch to raw sequence file
multi_fasta = "../../data-raw/021_retrieve-proteins/{0}/NagyTranscriptHits_ProteinSeqs.fasta".format(species)

protein_fasta_replace_J(multi_fasta, species)
