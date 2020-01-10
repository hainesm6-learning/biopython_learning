from Bio import SeqIO
import filecmp

EXAMPLE_SEQRECORD = SeqIO.read("U49845.1.gb", "genbank")

def genbank_details():
    genbank_att = [att for att in dir(EXAMPLE_SEQRECORD) if not att.startswith(
        "__") and not att.startswith("_") and not callable(getattr(EXAMPLE_SEQRECORD, att))]
    for att in genbank_att:
        print(f"{att} of EXAMPLE_SEQRECORD is/are {getattr(EXAMPLE_SEQRECORD, att)}")
        print("\n")

def genbank_round_trip():
    SeqIO.write(EXAMPLE_SEQRECORD, "example_seqrecord.gb", "genbank")
    print(f"It is {filecmp.cmp('U49845.1.gb', 'example_seqrecord.gb')} that round trips of genbank files are possible using Biopython")
    with open("U49845.1.gb", "r") as file1:
        with open("EXAMPLE_SEQRECORD.gb", "r") as file2:
            file_diffs = set(file1).symmetric_difference(file2)
    with open('output_file.txt', 'w') as file_out:
        for line in file_diffs:
            file_out.write(line)

# Although the above suggests round trips of genbank files are not possible using Biopython,
# "output_file.txt" suggests differences are only observed in formatting (e.g. spaces) and 
# in the number of characters per line in Seq objects. The latter is outlined in the Biopython
# tutorial.

# The output of the above also illustrates what attributes are required to describe genbank
# files in biopython so they may accurately be written.

def seqrecord_fasta_format():
    print(EXAMPLE_SEQRECORD.description)
    print(EXAMPLE_SEQRECORD.format("fasta")[:200])

def seqindex_play():
    index = SeqIO.index("U49845.1.gb", "genbank")
    print(type(index))