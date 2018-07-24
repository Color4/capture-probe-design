import sys

sys.path.append("/home/santang3/Python-modules/lib/python2.7/site-packages/")
sys.path.append("/home/santang3/.local/lib/python2.7/site-packages/")

import gffutils


def modify_attributes(feature, feature_type, index):

    Nested_ID = feature_type + "_" + str(index)
    Parent = feature.attributes['Parent'][0]
    ID = Parent + "-" + Nested_ID
    feature.attributes.clear()
    feature.attributes['Nested_ID'] = Nested_ID
    feature.attributes['Parent'] = Parent
    feature.attributes['ID'] = ID

    return(feature)


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
        f.write("##gff-version 3\n")
        for feature in gff3_database.features_of_type(featuretype = "gene", order_by = "start"):
             if feature.id not in duplicates:
                # print(feature.id, feature.attributes['ID'])
                f.write(str(feature) + '\n')
                i = 1
                for child in gff3_database.children(id = feature.id,
                                                    order_by = "start", featuretype = "cds"):
                    child = modify_attributes(child, "cds", i)
                    f.write(str(child) + '\n')
                    i += 1
                i = 1
                for child in gff3_database.children(id = feature.id,
                                                    order_by = "start", featuretype = "intron"):
                    child = modify_attributes(child, "intron", i)
                    f.write(str(child) + '\n')
                    i += 1
                i = 1
                for child in gff3_database.children(id = feature.id,
                                                    order_by = "start", featuretype = "splice_donor"):
                    child = modify_attributes(child, "splice_donor", i)
                    f.write(str(child) + '\n')
                    i += 1
                i = 1
                for child in gff3_database.children(id = feature.id,
                                                    order_by = "start", featuretype = "splice_acceptor"):
                    child = modify_attributes(child, "splice_acceptor", i)
                    f.write(str(child) + '\n')
                    i += 1


if __name__ == "__main__":
    path_to_raw_gff3 = sys.argv[1]
    path_to_clean_gff3 = sys.argv[2]
    clean_gff3_file(path_to_raw_gff3, path_to_clean_gff3)

