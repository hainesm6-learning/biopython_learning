from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from Bio.Seq import MutableSeq

# Seq doesn't inherit from String
print(Seq.__bases__)

# Example of f string formatting using Seq
my_seq = Seq("AGTACACTGGT", IUPAC.unambiguous_dna)
print(f"The value of my_seq is {my_seq}.")

# Can Seq objects be overwritten by reverse complements or are they immutable?
my_seq = my_seq.reverse_complement()
print(f"The value of my_seq is {my_seq}.")

# What methods are available from MutableSeq
print(f"MutableSeq attributes and methods are {dir(MutableSeq)}.")
