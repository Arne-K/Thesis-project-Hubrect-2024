import json
import os
from collections import Counter
import itertools as it
import csv
import argparse
import sys

# Path to the file
Dictionary_FilePath = 'BcCorrect_Dictionary.json'
Barcode_FilePath = "Barcodes_CelSeq1.tsv"

parser = argparse.ArgumentParser(description='This script returns a dictiorary of the barcodes, each having a list of unambiguous permutations within the defined Hamming distance + the barcode itself.')
parser.add_argument('--Barcode', help='the barcode to be corrected')
parser.add_argument('--maxHamming', type=int, help='define a maximum Hamming distance for generating the correction permutations')
args = parser.parse_args()

##################################
# =====  DEFINE FUNCTIONS  ===== #
##################################

def find_compatible_barcodes(barcode, HDmax = 0):
    #"""Given a barcode sequence and a maximum Hammin distance, it returns a list of compatible barcode sequences"""
    nt = ['N'] if HDmax == 0 else ['N','C','T','G','A']
    HDmax = 1 if HDmax == 0 else HDmax

    compatible_barcodes = set([barcode])
    for hd in range(1, HDmax+1): 
        comb = [''.join(l) for l in it.product(nt, repeat = hd)]
        for c in comb:
            for p in it.permutations(range(len(barcode)), hd):
                s0 = barcode
                for x, l in zip(p, c):
                    s0 = s0[:x] + l + s0[x+1:]
                compatible_barcodes.add(s0)
    return list(compatible_barcodes)

def Create_BlackListFiltered_Dictionary(Barcode_File, Defined_HDmax):
    Barcodes = [] 
    with open(Barcode_File, newline='', encoding='utf-8') as tsvfile:
        reader = csv.DictReader(tsvfile, delimiter='\t')  
        for row in reader:
            Barcodes.append(row['Barcode']) # get the barcode column by name
        print(f"Obtained {len(Barcodes)} barcodes from {Barcode_File}", file=sys.stderr)
    
    # Generate the dictionary object of possible permutations within defined Hamming distance for each barcode
    All_Bc_Permutations={}
    for Bc in Barcodes:
        All_Bc_Permutations[Bc] = find_compatible_barcodes(Bc, HDmax = Defined_HDmax) # returns a dictionary where the list of generated permutations (value) is assigned to the corresponding barcode (key)

    # Obtaining permutations that occur more than once across all the barcodes
    Permutations_Collapsed = list(it.chain.from_iterable(All_Bc_Permutations.values())) # collapse the permutation sets into a single list
    Permutation_Counts = Counter(Permutations_Collapsed) #  returns a dictionary-like object where the permutations are the keys, and values are the counts they occur in Permutations_Collapsed list
    Permutation_Counts2 = {permutation: count for permutation, count in Permutation_Counts.items() if count >= 2} # subset key-value pairs where permutation count is at least 2

    # Obtain the blacklist-filtered barcodes
    Permut_Blacklist = list(Permutation_Counts2.keys()) # Define barcode permutations that are not unambiguous
    Bc_Permutations_Filtered = {bc: [p for p in permutations if p not in Permut_Blacklist] for bc, permutations in All_Bc_Permutations.items()} # returns a dictionary where barcodes are the keys and unambiguous permutations are values
    return(Bc_Permutations_Filtered)


###################################
# =====  EXECUTIVE SECTION  ===== #
###################################

# Check if the barcode dictionary with the blacklist-corrected permutations is already created. If not, then create it and save it (serialized with json)
if os.path.exists(Dictionary_FilePath):
    pass
else:
    Dictionary_BlackListFiltered=Create_BlackListFiltered_Dictionary(Barcode_FilePath, Defined_HDmax=args.maxHamming)
    with open('BcCorrect_Dictionary.json', 'w') as f:
        json.dump(Dictionary_BlackListFiltered, f)
    print(f"Unambiguous_BarcodePermutations.py: Created BcCorrect_Dictionary.json", file=sys.stderr)

# Read in the dictionary
with open('BcCorrect_Dictionary.json', 'r') as f:
    Dictionary_BlackListFiltered = json.load(f)
print(f"Unambiguous_BarcodePermutations.py: Opened BcCorrect_Dictionary.json ...", file=sys.stderr)

# Takes a barcode as an input and prints out the unambiguous permutations one by one
# This output is captured in the Barcode_Correction.sh that works with it as an awk array
Unambig_PermN = len(Dictionary_BlackListFiltered[args.Barcode])
print(f"Unambiguous_BarcodePermutations.py: Found {Unambig_PermN} unambiguous permutations for {args.Barcode}", file=sys.stderr)
for permutation in Dictionary_BlackListFiltered[args.Barcode]:
    print(permutation)
