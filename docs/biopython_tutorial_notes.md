# Biopython tutorial notes

- [Biopython tutorial notes](#biopython-tutorial-notes)
  - [Description](#description)
  - [Sections to read](#sections-to-read)
  - [Chapter 1 Introduction](#chapter-1-introduction)
    - [1.3 installation](#13-installation)
  - [Chapter 3 Sequence Objects](#chapter-3-sequence-objects)
    - [3.1 Sequences and Alphabets](#31-sequences-and-alphabets)
    - [3.2 Sequences act like strings](#32-sequences-act-like-strings)
    - [3.3 Slicing a sequence](#33-slicing-a-sequence)
    - [3.4 Turning Seq objects into strings](#34-turning-seq-objects-into-strings)
    - [3.5 Concatenating or adding sequences](#35-concatenating-or-adding-sequences)
    - [3.6 Changing case](#36-changing-case)
    - [3.7 Nucleotide sequences and (reverse) complements](#37-nucleotide-sequences-and-reverse-complements)
    - [3.11 Comparing Seq objects](#311-comparing-seq-objects)
    - [3.12 MutableSeq objects](#312-mutableseq-objects)
    - [3.13 UnknownSeq objects](#313-unknownseq-objects)
    - [3.14 Working with strings directly](#314-working-with-strings-directly)
  - [Chapter 4 Sequence annotation objects](#chapter-4-sequence-annotation-objects)
    - [4.1 The SeqRecord object](#41-the-seqrecord-object)
    - [4.2 Creating a SeqRecord](#42-creating-a-seqrecord)
      - [4.2.1 SeqRecord objects from scratch](#421-seqrecord-objects-from-scratch)
      - [4.2.2 SeqRecord objects from FASTA files](#422-seqrecord-objects-from-fasta-files)
      - [4.2.3 SeqRecord objects from GenBank files](#423-seqrecord-objects-from-genbank-files)
    - [4.3 Feature, location and position objects](#43-feature-location-and-position-objects)
      - [4.3.1 SeqFeature objects](#431-seqfeature-objects)

## Description

These notes were taken by [hainesm6](https:\\github.com\hainesm6) and are based on the material provided by the [Biopython Tutorial and Cookbook](http://biopython.org/DIST/docs/tutorial/Tutorial.html#htoc15).

## Sections to read

- [x] Chapter 3:
  - [x] 3.1 - 3.7
  - [x] 3.11 - 3.14
- [ ] Chapter 4
- [ ] Chapter 5

## Chapter 1 Introduction

### 1.3 installation

```python
pip install biopython
```

## Chapter 3 Sequence Objects

- Seq objects have different methods compared to strings.
- Sequence objects have an alphabet attribute describing what the sequence string characters mean.
  
### 3.1 Sequences and Alphabets

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
```

### 3.2 Sequences act like strings

- Many string methods are available to Seq even though Seq does not inherit from String ([seq_play.py](/scripts/seq_play.py) e.g. **len()**, **enumerate()**, **count()**.
- Unless converted to [MutableSeq](##3.12-MutableSeq-objects), Seq objects are immutable.

### 3.3 Slicing a sequence

- Slicing is consistent with other python objects e.g. lists.

### 3.4 Turning Seq objects into strings

- The **str()** constructor can be used to convert a Seq object into a String object:

```python
>>> str(my_seq)
'GATCGATGGGCCTATATAGGATCGAAAATCGC'
```

- Feasible to format a string representation of a Seq object using a "%s" placeholder:

```python
>>> fasta_format_string = ">Name\n%s\n" % my_seq
>>> print(fasta_format_string)
>Name
GATCGATGGGCCTATATAGGATCGAAAATCGC
<BLANKLINE>
```

- It is also possible to format sequences using f strings [seq_play.py](/scripts/seq_play.py)

### 3.5 Concatenating or adding sequences

- Seq objects written in the same alphabet can be added together using the **+** operator.

### 3.6 Changing case

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
```

### 3.7 Nucleotide sequences and (reverse) complements

- **complement()** and **reverse_complement()** methods calculate the complement and reverse complement of nucleotide sequences, respectively.
- Seq objects are immutable and as such the result of these methods is not applied by default.

### 3.11 Comparing Seq objects

- Biopython will compare Seq objects based on characters. As such, nucleotide and protein sequences can be evaluated as equivalent. However, a warning is provided.

### 3.12 MutableSeq objects

- MutableSeq objects are created from Seq objects using the **tomutable()** method or by invoking the **MutableSeq()** constructor:

```python
>>> from Bio.Seq import MutableSeq
>>> from Bio.Alphabet import IUPAC
>>> mutable_seq = MutableSeq("GCCATTGTAATGGGCCGCTGAAAGGGTGCCCGA", IUPAC.unambiguous_dna)
```

- In addition to Seq methods, **insert()**, **append()** and indexing are useful operations which can be performed on MutableSeq objects.
- Like lists, MutableSeq objects are not hashable, they cannot be used as dictionary keys.
- The **toseq()** method converts MutableSeq objects to Seq objects.

### 3.13 UnknownSeq objects

- The UnknownSeq class enables unknown sequences of arbitary length to be defined while minimising memory consumption e.g.

```python
>>> from Bio.Seq import UnknownSeq
>>> from Bio.Alphabet import IUPAC
>>> unk_dna = UnknownSeq(20, alphabet=IUPAC.ambiguous_dna)
>>> unk_dna
UnknownSeq(20, alphabet=IUPACAmbiguousDNA(), character='N')
>>> print(unk_dna)
NNNNNNNNNNNNNNNNNNNN
```

- Unknown sequences are present in GenBank and EMBL files where only continuous overlapping fragments (contigs) may be given.

### 3.14 Working with strings directly

pass

## Chapter 4 Sequence annotation objects

- [SeqRecord](http://biopython.org/DIST/docs/api/Bio.SeqRecord.SeqRecord-class.html) and [SeqFeature](http://biopython.org/DIST/docs/api/Bio.SeqFeature.SeqFeature-class.html) documentation is available online.

### 4.1 The SeqRecord object

- The **SeqRecord** class is defined in the **Bio.SeqRecord** module.
- The attributes of a SeqRecord instance are sufficient to completely describe genbank entries (Refer to [seqrecord_play.py](/scripts/seqrecord_play.py))

### 4.2 Creating a SeqRecord

#### 4.2.1 SeqRecord objects from scratch

- When creating a SeqRecord object from scratch, a Seq object is required and an **id** attribute is important if the SeqRecord object will be written to file.
- The following generates a SeqRecord object and assigns an id:

```python
>>> from Bio.Seq import Seq
>>> from Bio.SeqRecord import SeqRecord
>>> simple_seq = Seq("GATC")
>>> simple_seq_r = SeqRecord(simple_seq, id="AC12345")
```

- Several additional attributes are available and can be used to completely describe GenBank Records. Refer to [Sample_GenBank_Record](https://www.ncbi.nlm.nih.gov/Sitemap/samplerecord.html) and the output of [seqrecord_play.py](/scripts/seqrecord_play.py).

#### 4.2.2 SeqRecord objects from FASTA files

- SeqRecord objects are generated from FASTA files using the **SeqIO** module and associated methods [Chapter 5](#Chapter-5-Sequence-Input/Output).

```python
>>> from Bio import SeqIO
>>> record = SeqIO.read("NC_005816.fna", "fasta")
>>> record
SeqRecord(seq=Seq('TGTAACGAACGGTGCAATAGTGATCCACACCCAACGCCTGAAATCAGATCCAGG...CTG',
SingleLetterAlphabet()), id='gi|45478711|ref|NC_005816.1|', name='gi|45478711|ref|NC_005816.1|',
description='gi|45478711|ref|NC_005816.1| Yersinia pestis biovar Microtus ... sequence',
dbxrefs=[])
```

- A default alphabet has been used to for the Seq object. If possible, a more specific alphabet should be utilised (refer to SeqIO documentation).
- The first word of the FASTA record’s title line (after removing the greater than symbol) is used for both the id and name attributes.
- The whole title line is used for the record description.

#### 4.2.3 SeqRecord objects from GenBank files

- SeqRecord objects are generated from GenBank files using a similar syntax to that used above for FASTA files.
- An example is given in the [seqrecord_play.py](/scripts/seqrecord_play.py) script. Of note:
  - A more specific alphabet is given.
  - The **name** is derived from the LOCUS, while the **id** includes the version suffix.
  - Most of the information gets stored as annotations in the **example_seqrecord.annotations** attribute as key, value pairs.

### 4.3 Feature, location and position objects

#### 4.3.1 SeqFeature objects

- **SeqFeature** objects describe [Features](https://www.ncbi.nlm.nih.gov/Sitemap/samplerecord.html#FeaturesA) in SeqRecord objects.
- Three important attributes belong to the SeqFeature class:
  - **SeqFeature.type** which is equivalent to the key of Features in GenBank files e.g. CDS.
  - **SeqFeature.location** which is usually a **SeqFeature.FeatureLocation** object with the following arguments:
  
  ```python
  FeatureLocation(start, end, strand=None, ref=None, ref_db=None)
  ```
  