from Bio.Seq import Seq
from Bio.Alphabet import IUPAC

# Seq doesn't inherit from String
print(Seq.__bases__)

# Example of f string formatting using Seq
my_seq = Seq("AGTACACTGGT", IUPAC.unambiguous_dna)
print(f"The current sequences is {my_seq}.")