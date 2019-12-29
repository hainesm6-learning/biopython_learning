# Description
These notes were taken by [hainesm6](https:\\github.com\hainesm6) and are based on the material provided by the [Biopython Tutorial and Cookbook](http://biopython.org/DIST/docs/tutorial/Tutorial.html#htoc15).

# Sections to read
- [x] Chapter 3:
    - [x] 3.1 - 3.7
    - [x] 3.11 - 3.14
- [ ] Chapter 4
- [ ] Chapter 5
- [ ] Chapter 21
- [ ] Chapter 22:
    - [ ] 22.1

# Chapter 1 Introduction
## 1.3 installation
    pip install biopython

# Chapter 3 Sequence Objects
- Seq objects have different methods compared to strings.
- Sequence objects have an alphabet attribute describing what the sequence string characters mean.
## 3.1 Sequences and Alphabets
- Biopython available alphabets are defined in the **Bio.Alphabet** module.
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
- Many string methods are available to Seq even though Seq does not inherit from String ([biopython_play.py](biopython_learning/scripts/biopython_play.py)) e.g. **len()**, **enumerate()**, **count()**.
- Unless converted to [MutableSeq](##3.12-MutableSeq-objects), Seq objects are immutable.

## 3.3 Slicing a sequence
- Slicing is consistent with other python objects e.g. lists.

## 3.4 Turning Seq objects into strings
- The **str()** constructor can be used to convert a Seq object into a String object:
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
- It is also possible to format sequences using f strings [biopython_play.py](biopython_learning/scripts/biopython_play.py)

## 3.5 Concatenating or adding sequences
- Seq objects written in the same alphabet can be added together using the **+** operator. 

## 3.6 Changing case
- **upper()** and **lower()** methods change the case of sequences.
- Note lower case sequences are not valid IUPAC:
    ```python
    >>> from Bio.Seq import Seq
    >>> from Bio.Alphabet import IUPAC
    >>> dna_seq = Seq("ACGT", IUPAC.unambiguous_dna)
    >>> dna_seq
    Seq('ACGT', IUPACUnambiguousDNA())
    >>> dna_seq.lower()
    Seq('acgt', DNAAlphabet())

## 3.7 Nucleotide sequences and (reverse) complements
- **complement()** and **reverse_complement()** methods calculate the complement and reverse complement of nucleotide sequences, respectively.
- Seq objects are immutable and as such the result of these methods is not applied by default.

## 3.11 Comparing Seq objects
- Biopython will compare Seq objects based on characters. As such, nucleotide and protein sequences can be evaluated as equivalent. However, a warning is provided.

## 3.12 MutableSeq objects
- MutableSeq objects are created from Seq objects using the **tomutable()** method or by invoking the **MutableSeq()** constructor:
    ```python
    >>> from Bio.Seq import MutableSeq
    >>> from Bio.Alphabet import IUPAC
    >>> mutable_seq = MutableSeq("GCCATTGTAATGGGCCGCTGAAAGGGTGCCCGA", IUPAC.unambiguous_dna)
- In addition to Seq methods, **insert()**, **append()** and indexing are useful operations which can be performed on MutableSeq objects.
- Given like lists, MutableSeq objects are not hashable, they cannot be used as dictionary keys.
- The **toseq()** method converts MutableSeq objects to Seq objects.

## 3.13 UnknownSeq objects
- The UnknownSeq class enables unknown sequences of arbitary length to be defined while minimising memory consumption e.g.
    ```python
    >>> from Bio.Seq import UnknownSeq
    >>> from Bio.Alphabet import IUPAC
    >>> unk_dna = UnknownSeq(20, alphabet=IUPAC.ambiguous_dna)
    >>> unk_dna
    UnknownSeq(20, alphabet=IUPACAmbiguousDNA(), character='N')
    >>> print(unk_dna)
    NNNNNNNNNNNNNNNNNNNN
- Unknown sequences are present in GenBank and EMBL files where only continuous overlapping fragments (contigs) may be given.
## 3.14 Working with strings directly
pass

# Chapter 4 Sequence annotation objects
- [SeqRecord](http://biopython.org/DIST/docs/api/Bio.SeqRecord.SeqRecord-class.html) and [SeqFeature](http://biopython.org/DIST/docs/api/Bio.SeqFeature.SeqFeature-class.html) documentation is available online.
## 4.1 The SeqRecord object







