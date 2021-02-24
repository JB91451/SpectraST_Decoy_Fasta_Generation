# SpectraST_Decoy_Fasta_Generation

A small python script that generates a fasta file with correctly annotated decoys from spectraST's decoy-generation log file.

## Installation
- a working python installation (~= 3.5) is required
- ensure that the packages fuzzywuzzy and biopython are available on you computer
- download the file SpectraST_DECOY_Fasta.py

## Application
- open a comand window and navigate to the folder where SpectraST_DECOY_Fasta.py is stored
- type <br>
```SpectraST_DECOY_Fasta.py -f path\to\target_only_fasta -l path\to\decoy_generation_log_file\Decoy_generation_log.txt -o \path\to\outfile\out_basename``` <br>
This comand will load the log file of the decoy generation and create a fasta file with target and decoy entries for all proteins present in the spectral library.

## Options
- -f = Target-only fasta file that can be used to map the original peptides, required

- -l = Log file of the spectraST decoy generation step. For each peptide in the spectral library a line similar to <br>
      'DECOY: Shuffle PEPTIDE to TIPEPDELK . (2 AAs added randomly.)' <br>
has to be present, required

- -s = String to separate newly generated peptides in the DECOY proteins, default = KRKR

- -d = Decoy string during the original database search, entries in the provided fasta file starting with this string are filtered prior mapping library peptides, default = 'Reverse'

- -o = Output file names. Files with the extension '.fasta' and '.tsv' will be generated. *Caution:* Will be overwritten, if existing!, required


