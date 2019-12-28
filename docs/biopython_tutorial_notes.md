# Description
These notes were taken by [hainesm6](https:\\github.com\hainesm6) and are based on the material provided by the [Biopython Tutorial and Cookbook](http://biopython.org/DIST/docs/tutorial/Tutorial.html#htoc15).

# Sections to read
- [ ] Chapter 3:
    - [ ] 3.1 - 3.7
    - [ ] 3.11 - 3.14
- [ ] Chapter 4
- [ ] Chapter 5
- [ ] Chapter 21
- [ ] Chapter 22:
    - [ ] 22.1

# Chapter 1 Introduction
## 1.3 installation
    pip install biopython

# Chapter 3 Sequence Objects
- Sequence objects differ from strings given they have different methods.
- Sequence objects have an alphabet attribute describing what the sequence string characters mean.
## 3.1 Sequences and Alphabets
- Biopython available alphabets are defined in the Bio.Alphabet module.
- IUPAC alphabets are commonly used.
- While a default alphabet argument is provided, alphabets are typically specificied as in the following example:
    ```python
    >>> from Bio.Seq import Seq
    >>> from Bio.Alphabet import IUPAC
    >>> my_seq = Seq("AGTACACTGGT", IUPAC.unambiguous_dna)
    >>> my_seq
    Seq('AGTACACTGGT', IUPACUnambiguousDNA())
    >>> my_seq.alphabet
    IUPACUnambiguousDNA()

## 3.2 Sequences act like strings
- Many string methods are available to Seq even though Seq does not inherit from String e.g. len(), enumerate(), count()
- Unless stated, Seq objects are immutable.

## 3.3 Slicing a sequence
- Slicing is consistent with other python objects e.g. lists.

## 3.4 Turning Seq objects into strings
- The str() constructor can be used to convert a Seq object into a String object:
    ```python
    >>> str(my_seq)
    'GATCGATGGGCCTATATAGGATCGAAAATCGC'
- Feasible to format a string representation of a Seq object using a "%s" placeholder:
    ```python
    >>> fasta_format_string = ">Name\n%s\n" % my_seq
    >>> print(fasta_format_string)
    >Name
    GATCGATGGGCCTATATAGGATCGAAAATCGC
    <BLANKLINE>

## 3.5 Concatenating or adding sequences



