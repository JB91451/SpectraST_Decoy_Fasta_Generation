import argparse
from fuzzywuzzy import fuzz
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import os

# parse arguments from the command line
parser = argparse.ArgumentParser(description="generate a fasta file for correct annotation of decoy hits from spectraST's decoy-generation log file")
parser.add_argument('-f', '--fasta', help = "Target-only fasta file that can be used to map the original peptides", required=True)
parser.add_argument('-l', '--decoy_log', help = "Log file of the spectraST decoy generation step. For each peptide in the spectral library a line similar to 'DECOY: Shuffle PEPTIDE to TIPEPDELK . (2 AAs added randomly.)' is required", required=True)
parser.add_argument('-s', '--separation_string', default = "KRKR" , help = "string to separate newly generated peptides (default = KRKR)")
parser.add_argument('-d', '--decoy_string', default = "Reverse_", help = "decoy string during the original database search, entries in the provided fasta file starting with this string are filtered prior mapping library peptides (default = 'Reverse')")
parser.add_argument('-o', '--out_base', help = "Output file names. Files with the extension '.fasta' and '.tsv' will be generated. Caution: Will be overwritten, if existing!", required=True)
args = parser.parse_args()




# try to open spectraST log file and import all lines starting with "DECOY:" to peptide_list
try:
    log_file = open(args.decoy_log, "r")
except FileNotFoundError:
    print("Error: Log file not found, please provide a valid log file")

peptide_list = list()
for line in log_file:
    if line.startswith("DECOY: "):
        peptide_list.append(line)

    # lists to hold original and modified peptide sequences
original_peptide_list = list()
new_peptide_list = list()
    # loop through peptide_list and obtain whitespace separated words from it
for entry in peptide_list:
    tmp_list = entry.split()

    for word in tmp_list:
        # use fuzzy string comparison to allow some deviation of the original pattern
            # if no 'to' is present (with at least 80 percent similarity) in a certain line, the return would be 'unknown' and ultimately an error
        found_to = "unknown"
        current_ratio = 80
        new_ratio = fuzz.ratio(word, "to")
        if new_ratio > current_ratio:
            current_ratio = new_ratio
            found_to = word

        # locate the position of found_to and identify peptides before and after
        if found_to != "unknown":
            to_position = tmp_list.index(found_to)
            if to_position > 0 & to_position < len(tmp_list):
                original_peptide = tmp_list[to_position - 1]
                new_peptide = tmp_list[to_position + 1]

    # append the peptide pair to the corresponding list after each loop through the current tmp_list
    original_peptide_list.append(original_peptide)
    new_peptide_list.append(new_peptide)


# try to open the fasta file
try:
    fastaFile = SeqIO.parse(args.fasta, "fasta")
except FileNotFoundError:
    print("Error: Fasta file not found, please provide a valid database file")

    # lists to hold identifier and protein sequence of each entry
proteinID = list()
proteinSeq = list()
proteinDescription = list()
    # import fasta and remove entries that contain the Reverse_/Decoy_ Tag provided with -d/--decoy_string on the command line
for record in fastaFile:
    tmpID = record.id
   # if tmpID.startswith(args.decoy_string):
    if args.decoy_string in tmpID:
        print("removed", tmpID, "prior mapping.")
    else:
        tmpSeq = record.seq
        tmpDescr = record.description
        proteinID.append(tmpID)
        proteinSeq.append(tmpSeq)
        proteinDescription.append(tmpDescr)

    # generate a reversed protein sequence list for usage in output fasta (concatenating the fw-sequences with the decoy tags would cause interference during protein mapping)
reverse_proteinSeq = list()
for one_proteinSeq in proteinSeq:
    reverse_proteinSeq.append(one_proteinSeq[::-1])


# get indices of protein sequences that contain each of the original_peptides (multiple matches are possible albeit this should not occur in unique-only libraries)
original_peptide_index_list = list()
for original_peptide in original_peptide_list:
    index = [i for i, s in enumerate(proteinSeq) if original_peptide in s]
    original_peptide_index_list.append(index)

separator = args. separation_string
for new_peptide, index in zip(new_peptide_list, original_peptide_index_list):
    for one_index in index:
        tmpSeq = reverse_proteinSeq[one_index] + separator + new_peptide
        reverse_proteinSeq[one_index] = tmpSeq



# open a new fasta file to write results
Fasta_out = os.path.join(args.out_base + ".fasta")
new_fasta = open(Fasta_out, "w")
records = []
for sequence, identifier, description in zip(reverse_proteinSeq, proteinID, proteinDescription):

    records.append(SeqRecord(seq = sequence, id = identifier, description=description, name=""))

SeqIO.write(records, new_fasta, "fasta")
new_fasta.close()

# open a new tsv file to write overview file
Tsv_out = os.path.join(args.out_base + ".tsv")
Tsv_file = open(Tsv_out, "w")
Tsv_header = str("Protein ID \t Original Peptide Sequence \t Decoy Peptide Sequence \n")
Tsv_file.write(Tsv_header)
for one_originalPep, one_decoyPep, one_index in zip(original_peptide_list, new_peptide_list, original_peptide_index_list):
    for unlist_one_index in one_index:
        Tsv_result = str(proteinID[unlist_one_index]) + "\t" + str(one_originalPep) + "\t" + str(one_decoyPep) +"\n"
        Tsv_file.write(Tsv_result)

Tsv_file.close()

print()
print("DECOY peptides written to", Fasta_out)

