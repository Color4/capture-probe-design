import sys

sys.path.append("/home/santang3/Python-modules/lib/python2.7/site-packages/")
sys.path.append("/home/santang3/.local/lib/python2.7/site-packages/")

import gffutils


def remove_gff3_duplicates(species):
    """Remove genes and children from gff3

    Removes genes and children if the same protein mapped to multiple
    regions of the genome. Keeps the first one to map since the scores
    do not differ.

    Args:
        species ('str'): Species for which to extract features from the gff3 file.

    Returns:
        None: Writes gff3 to disk.
    """

    # Open raw gff3 database
    gff3_database = gffutils.FeatureDB("../../data-raw/030_exonerate_p2g/{0}/Nagy_transcriptome_p2g.db".format(species))

    # Open gff3 file to which unique feature will be written
    clean_gff3 = "../../data-clean/030_exonerate_p2g/{0}/Nagy_transcriptome_p2g.gff3".format(species)

    # Create list with duplicates
    duplicates = [feature.id for feature in gff3_database.features_of_type("gene") if len(feature.id.split(".")) != 2]

    # Iterate over cds in gff3 annotation file and write to new file. Add attributes
    with open(clean_gff3, "w") as f:
        for feature in gff3_database.all_features(order_by = "start"):
            if feature.featuretype == "gene":
                if feature.id not in duplicates:
                    f.write(str(feature) + '\n')
                    for child in gff3_database.children(id = feature.id, order_by = "start"):
                        f.write(str(child) + '\n')

def create_gff3_database(species, data_state):
    """Create gff3 database from gff3 annotation file

    Args:
        species ('str'): Species for which to extract features from the gff3 file.
        data_state ('str'): One of "raw" or "clean", specifying the state of the gff3 file.

    Returns:
        None: Writes gff3 database to disk.
    """

    # Specify paths
    gff3_file = "../../data-{0}/030_exonerate_p2g/{1}/Nagy_transcriptome_p2g.gff3".format(data_state, species)
    database_file = "../../data-{0}/030_exonerate_p2g/{1}/Nagy_transcriptome_p2g.db".format(data_state, species)

    # Create database
    gffutils.create_db(gff3_file, dbfn=database_file, force=True, keep_order=True, merge_strategy='merge', sort_attribute_values=True)


if __name__ == "__main__":
    species = sys.argv[1]
    print("Creating gff3 database from raw annotation file")
    create_gff3_database(species, "raw")  # Create database from raw gff3

    print("Removing duplicates from raw gff3 annotation file")
    remove_gff3_duplicates(species)  # Remove duplicates from raw database and gff3 file

    print("Creating gff3 database from clean annotation file")
    create_gff3_database(species, "clean")  # Create database from clean fatabase
