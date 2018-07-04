import sys

sys.path.append("/home/santang3/Python-modules/lib/python2.7/site-packages/")
sys.path.append("/home/santang3/.local/lib/python2.7/site-packages/")

import gffutils
import glob
import pybedtools


def fields2name(f):
    """replace GFF featuretype field with the attributes field

    Args:
        f (:obj: 'BedTools'): BedTools object with fields representing columns from a gff3 file

    Returns:
        f (:obj: 'BedTools'): Same BedTools object with the featuretype field replaced with the attributes field. Thus, the attributes field will be written as a header when writing out sequences rather than the featuretype field.
    """
    f[2] = f[-1]
    return(f)


def write_cds_to_gff3(species):
    """From a gff3 database generated using gffutils, extract only the coding sequences and writes them to a new gff3 file.

    Args:
        species ('str'): Species for which to extract coding sequences from the gff3 file.

    Returns:
        None: Writes gff3 to disk.
    """

    # Open gff3 database
    gff3_database = gffutils.FeatureDB("../../data-clean/030_exonerate_p2g/{0}/Nagy_transcriptome_p2g.db".format(species))

    # Open multi fasta file to write coding sequence gff3 features
    coding_sequences = "../../data-clean/031_coding-sequences/{0}/NagyTranscriptHits_CodingSequences.gff3".format(species)

    # Iterate over cds in gff3 annotation file and write to new file. Add attributes
    with open(coding_sequences, "w") as f:
        for cds in gff3_database.features_of_type("cds"):
            cds.attributes['length'] = str(len(cds))
            cds.attributes['cds_id'] = cds.id
            cds.attributes['gene_id'] = cds.attributes["Parent"][0] + '_' + cds.id
            f.write(str(cds) + '\n')

def write_cds_to_fasta(species):
    """From a gff3 file with only coding sequences, extract the coding sequences from the species' reference genome and writes seuqences to a FASTA file. Fasta header is attribute column from gff3 file

    Args:
        species ('str'): Species for which to extract coding sequences from the gff3 file.

    Returns:
        None: Writes FASTA to disk.
    """

    # Open gff3 annotation file with coding sequences
    coding_sequences = "../../data-clean/031_coding-sequences/{0}/NagyTranscriptHits_CodingSequences.gff3".format(species)

    # Specify path to FASTA file for reference genome
    genome_fasta = glob.glob("../../../data-reference/reference-genomes/{0}/*.fasta".format(species))

    # Path to file to which sequences should be written
    fout = "../../data-clean/031_coding-sequences/{0}/NagyTranscriptHits_CodingSequences.fasta".format(species)

    # Create BedTools object from gff3 file, specify which field to write as header and extract sequences from reference. Strandedness (i.e + or -) is maintained.
    pybedtools.BedTool(coding_sequences).each(fields2name).sequence(fi=genome_fasta, s = True, name = True, fo = fout)


if __name__ == "__main__":
    species = sys.argv[1]
    write_cds_to_gff3(species)
    write_cds_to_fasta(species)
