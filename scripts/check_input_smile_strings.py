import os
import sys
from bvbrc_docking.utils import validate_smiles

def check_smile_strings(staging_dir, smile_lines):
    inputs = []
    failed_smile_strings = []
    
    for item_clean in smile_lines:
        # item = item.strip("]").strip("[").stri('"')
        item = item_clean.split()
        # flexibiliy with and without names
        if len(item) == 3:
            ident = item[0]
            name = item[1]
            smiles_str =item[2]
        elif len(item) == 2:
            ident = item[0]
            smiles_str =item[1]
        else:
            # if something else comes in just try taking the first and last items
            ident = item[0]
            smiles_str =item[-1]
        # Now validate the strings
        if validate_smiles(smiles_str) == True:
            inputs.append(item_clean)
        else:
            if validate_smiles(ident) == True:
                inputs.append(item_clean)
            else:
                failed_smile_strings.append(item_clean)
    # Write valid inputs to file
    valid_smile_strs = os.path.join(staging_dir,"ligands.smi")
    with open(valid_smile_strs, "w") as file1:
        for item in inputs:
            file1.write(item)
            #file1.write(f"{item}\n")

    # Write failed_smile_strings to file
    failed_smile_strs_txt = os.path.join(staging_dir,"invalid_smile_strings.txt")
    with open(failed_smile_strs_txt, "w") as file2:
        for item in failed_smile_strings:
            file2.write(item)
    return

def open_smile_file(input_smile_file):
    with open(input_smile_file, 'r') as file:
        lines=file.readlines()
    return lines

def main(argv):
    staging_dir = argv[1]
    # input_smile_data = argv[2]
    input_smile_file = argv[2]
    smile_lines = open_smile_file(input_smile_file)
    check_smile_strings(staging_dir, smile_lines)
    print('nicole this is updated 8 6 2024 3')
if __name__ == "__main__":
    main(sys.argv)