import sys
import os

# Add paths to gffutils (only for running on HPCNODE)
sys.path.append("/home/santang3/Python-modules/lib/python2.7/site-packages/")
sys.path.append("/home/santang3/.local/lib/python2.7/site-packages/")

import gffutils


def create_gff3_database(path_to_gff3, data_quality):
    """Create gff3 database from gff3 annotation file

    Args:
        path_to_gff3 ('str'): Path where gff3 file is located. Will also write database to same path.
        data_quality ('str'): One of "raw" or "clean", specifying the state of the gff3 file.

    Returns:
        None: Writes gff3 database to disk.
    """

    # Change path to loction of gff3 file
    os.chdir(path_to_gff3)

    # Specify filenames
    gff3_file = "Nagy_transcriptome_p2g_{0}.gff3".format(data_quality)
    database_file = "Nagy_transcriptome_p2g_{0}.db".format(data_quality)

    # Create database in same directory (default)
    gffutils.create_db(gff3_file, dbfn=database_file,
                       force=True, keep_order=True,
                       merge_strategy='merge', sort_attribute_values=True)


if __name__ == "__main__":
    path_to_gff3 = sys.argv[1]
    data_quality = sys.argv[2]
    create_gff3_database(path_to_gff3, data_quality)
