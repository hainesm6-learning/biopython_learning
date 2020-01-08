from Bio import SeqIO
import filecmp

example_seqrecord = SeqIO.read("U49845.1.gb", "genbank")
genbank_att = [att for att in dir(example_seqrecord) if not att.startswith(
    "__") and not att.startswith("_") and not callable(getattr(example_seqrecord, att))]
for att in genbank_att:
    print(f"{att} of example_seqrecord is/are {getattr(example_seqrecord, att)}")
    print("\n")
SeqIO.write(example_seqrecord, "example_seqrecord.gb", "genbank")
print(f"It is {filecmp.cmp('U49845.1.gb', 'example_seqrecord.gb')} that round trips of genbank files are possible using Biopython")
with open("U49845.1.gb", "r") as file1:
    with open("example_seqrecord.gb", "r") as file2:
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

