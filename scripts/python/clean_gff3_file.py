import sys

sys.path.append("/home/santang3/Python-modules/lib/python2.7/site-packages/")
sys.path.append("/home/santang3/.local/lib/python2.7/site-packages/")

import gffutils


def clean_gff3_file(path_to_raw_gff3, path_to_clean_gff3):
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
    gff3_database = gffutils.FeatureDB(path_to_raw_gff3)

    # Open gff3 file to which unique feature will be written
    clean_gff3 = path_to_clean_gff3

    # Create list with duplicates
    duplicates = [feature.id for feature in gff3_database.features_of_type("gene") if len(feature.id.split(".")) != 2]

    # Iterate over cds in gff3 annotation file and write to new file. Add attributes
    with open(clean_gff3, "w") as f:
        for feature in gff3_database.all_features(order_by = "start"):
            if feature.featuretype == "gene":
                if feature.id not in duplicates:
                    f.write(str(feature) + '\n')
                    i = 1
                    for child in gff3_database.children(id = feature.id,
                                                        order_by = "start", featuretype = "cds"):
                        child.id = "cds" + "_" + str(i)
                        f.write(str(child) + '\n')
                        i += 1
                    i = 1
                    for child in gff3_database.children(id = feature.id,
                                                        order_by = "start", featuretype = "intron"):
                        child.id = "intron" + "_" + str(i)
                        f.write(str(child) + '\n')
                        i += 1
                    i = 1
                    for child in gff3_database.children(id = feature.id,
                                                        order_by = "start", featuretype = "splice_donor"):
                        child.id = "splice_donor" + "_" + str(i)
                        f.write(str(child) + '\n')
                        i += 1
                    i = 1
                    for child in gff3_database.children(id = feature.id,
                                                        order_by = "start", featuretype = "splice_acceptor"):
                        child.id = "splice_acceptor" + "_" + str(i)
                        f.write(str(child) + '\n')
                        i += 1


if __name__ == "__main__":
    path_to_raw_gff3 = sys.argv[1]
    path_to_clean_gff3 = sys.argv[2]
    clean_gff3_file(path_to_raw_gff3, path_to_clean_gff3)

