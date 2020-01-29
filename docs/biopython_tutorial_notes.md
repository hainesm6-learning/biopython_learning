# Biopython tutorial notes

- [Biopython tutorial notes](#biopython-tutorial-notes)
  - [Description](#description)
  - [Additional useful sources](#additional-useful-sources)
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
      - [4.3.2 Positions and locations](#432-positions-and-locations)
      - [4.3.3 Sequence described by a feature or location](#433-sequence-described-by-a-feature-or-location)
    - [4.4 Comparison](#44-comparison)
    - [4.5 References](#45-references)
    - [4.6 The format method](#46-the-format-method)
    - [4.7 Slicing a SeqRecord](#47-slicing-a-seqrecord)
    - [4.8 Adding SeqRecord objects](#48-adding-seqrecord-objects)
    - [4.9 Reverse-complementing SeqRecord objects](#49-reverse-complementing-seqrecord-objects)
  - [Chapter 5 Sequence Input/Output](#chapter-5-sequence-inputoutput)
    - [5.1 Parsing or Reading Sequences](#51-parsing-or-reading-sequences)
    - [5.2 Parsing sequences from compressed files](#52-parsing-sequences-from-compressed-files)
    - [5.3 Parsing Sequences from the net](#53-parsing-sequences-from-the-net)
    - [5.4 Sequence files as Dictionaries](#54-sequence-files-as-dictionaries)
    - [5.5 Writing Sequence Files](#55-writing-sequence-files)
      - [5.5.1 round trips](#551-round-trips)
      - [5.5.2 Converting between sequence file formats](#552-converting-between-sequence-file-formats)
      - [5.5.3 converting a file of sequences to their reverse complements](#553-converting-a-file-of-sequences-to-their-reverse-complements)
    - [.6  Low level FASTA and FASTQ parsers](#6-low-level-fasta-and-fastq-parsers)

## Description

These notes were taken by [hainesm6](https:\\github.com\hainesm6) and are based on the material provided by the [Biopython Tutorial and Cookbook](http://biopython.org/DIST/docs/tutorial/Tutorial.html#htoc15).

## Additional useful sources

- [Sample GenBank Record](https://www.ncbi.nlm.nih.gov/Sitemap/samplerecord.html#CDSB)
- [GenBank Feature Table Definition](http://www.insdc.org/files/feature_table.html#3.3)

## Sections to read

- [x] Chapter 3:
  - [x] 3.1 - 3.7
  - [x] 3.11 - 3.14
- [x] Chapter 4
- [x] Chapter 5

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

- Many string methods are available to Seq even though Seq does not inherit from String ([seq_play.py](/scripts/seq_play.py) e.g. **len()**, **enumerate()**, **count()**).
- Unless converted to [MutableSeq](##3.12-MutableSeq-objects), Seq objects are immutable.

### 3.3 Slicing a sequence

- Slicing is consistent with other python objects e.g. strings and lists.

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

- It is also possible to format sequences using **f strings** and the **format()** method [seq_play.py](/scripts/seq_play.py)

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
  - **SeqFeature.type** *(type = str)* which is equivalent to a [Feature Key](http://www.insdc.org/files/feature_table.html#1) in a GenBank file e.g. CDS.
  - **SeqFeature.location** which is often another object as described in the subsequent section.
  - **SeqFeature.qualifiers** formated as an OrderedDict with values given in lists (below). Note all string values must be given in *list* objects. This is because multiple keys of the same value are not permitted in Python *dicts* and yet multiple qualifier values may be necessary ([feature table refer to 3.3.2](http://www.insdc.org/files/feature_table.html#3.3))
  
  ```python
  qualifiers:
    Key: chromosome, Value: ['IX']
    Key: db_xref, Value: ['taxon:4932']
    Key: mol_type, Value: ['genomic DNA']
    Key: organism, Value: ['Saccharomyces cerevisiae']
  ```

#### 4.3.2 Positions and locations

- Most commonly a SeqFeature.location attribute will contain a **FeatureLocation** instance unless discontinuous locations are required. For further information on the FeatureLocation class the following code should be executed:

```python
from Bio.SeqFeature import FeatureLocation
help(FeatureLocation)
```

- **CompoundLocation** objects hold multiple FeatureLocation objects, enabling for instance eukaryotic gene locations to be described:

```python
 |      >>> from Bio.SeqFeature import FeatureLocation, CompoundLocation
 |      >>> f1 = FeatureLocation(10, 40, strand=+1)
 |      >>> f2 = FeatureLocation(50, 59, strand=+1)
 |      >>> f = CompoundLocation([f1, f2])
 |      >>> len(f) == len(f1) + len(f2) == 39 == len(list(f))
 |      True
 |      >>> print(f.operator)
 |      join
 |      >>> 5 in f
 |      False
 |      >>> 15 in f
 |      True
 |      >>> f.strand
 |      1
 ```

- Passing in *int* values into **start** and **end** FeatureLocation arguments acutally generates **ExactPosition** objects.
- To allow more abstract defitions of **start** and **end** FeatureLocation arguments, four classes are available:
  - **AfterPosition**
  - **BeforePosition**
  - **BetweenPosition**
  - **UnknownPosition**
- As an example:

```python
>>> from Bio import SeqFeature
>>> start_pos = SeqFeature.AfterPosition(5)
>>> end_pos = SeqFeature.BetweenPosition(9, left=8, right=9)
>>> my_location = SeqFeature.FeatureLocation(start_pos, end_pos)
>>> my_location.start
AfterPosition(5)
>>> print(my_location.start)
>5
>>> my_location.end
BetweenPosition(9, left=8, right=9)
>>> print(my_location.end)
(8^9)
```

- It is possible to search for which SeqFeatures occur at a given position e.g.

```python
>>> from Bio import SeqIO
>>> my_snp = 4350
>>> record = SeqIO.read("NC_005816.gb", "genbank")
>>> for feature in record.features:
...     if my_snp in feature:
...         print("%s %s" % (feature.type, feature.qualifiers.get("db_xref")))
...
source ['taxon:229193']
gene ['GeneID:2767712']
CDS ['GI:45478716', 'GeneID:2767712']
```

#### 4.3.3 Sequence described by a feature or location

- To extract the sequence from a SeqFeature object, the **extract()** method can be applied:

> extract(self, parent_sequence)
    Extract the feature's sequence from supplied parent sequence.
    The parent_sequence can be a Seq like object or a string, and will
    generally return an object of the same type. The exception to this is
    a MutableSeq as the parent sequence will return a Seq object.

### 4.4 Comparison

- SeqRecord instances should not be compared directly as is the case for all Python objects. Instead compare attributes of this class between instances.

### 4.5 References

- The Reference class contains several attributes to describe references associated with SeqRecord and SeqFeature objects:

```terminal
class Reference(builtins.object)
 |  Represent a Generic Reference object.
 |
 |  Attributes:
 |   - location - A list of Location objects specifying regions of
 |     the sequence that the references correspond to. If no locations are
 |     specified, the entire sequence is assumed.
 |   - authors - A big old string, or a list split by author, of authors
 |     for the reference.
 |   - title - The title of the reference.
 |   - journal - Journal the reference was published in.
 |   - medline_id - A medline reference for the article.
 |   - pubmed_id - A pubmed reference for the article.
 |   - comment - A place to stick any comments about the reference.
 ```

- The example_seqrecord variable from [seqrecord_play.py](/scripts/seqrecord_play.py) has several Reference objects stored in a list, accessed via the example_seqrecord.annotations['references'] key.
- An example from this *dict* value is:

```pythonprint
location: [0:5028]
authors: Roemer,T., Madden,K., Chang,J. and Snyder,M.
title: Selection of axial growth sites in yeast requires Axl2p, a novel plasma membrane glycoprotein
journal: Genes Dev. 10 (7), 777-793 (1996)
medline id:
pubmed id: 8846915
comment:
```

### 4.6 The format method

- The **format()** method of the SeqRecord class requires an "Bio.SeqIO output format" and returns the SeqRecord formatted in this manner.
- This method appears useful for displaying SeqRecords in the terminal.

### 4.7 Slicing a SeqRecord

- An incredibly useful feature of Biopython is the ability to slice SeqRecord objects returning a new object:

```python
>>> from Bio import SeqIO
>>> record = SeqIO.read("NC_005816.gb", "genbank")
>>> record
SeqRecord(seq=Seq('TGTAACGAACGGTGCAATAGTGATCCACACCCAACGCCTGAAATCAGATCCAGG...CTG',
IUPACAmbiguousDNA()), id='NC_005816.1', name='NC_005816',
description='Yersinia pestis biovar Microtus str. 91001 plasmid pPCP1, complete sequence',
dbxrefs=['Project:58037'])

>>> sub_record = record[4300:4800]
>>> sub_record
SeqRecord(seq=Seq('ATAAATAGATTATTCCAAATAATTTATTTATGTAAGAACAGGATGGGAGGGGGA...TTA',
IUPACAmbiguousDNA()), id='NC_005816.1', name='NC_005816',
description='Yersinia pestis biovar Microtus str. 91001 plasmid pPCP1, complete sequence.',
dbxrefs=[])
```

- Following slicing **SeqFeature.location** attributes are updated to reflect their new position.
- However, in the above example **sub_record.dbxrefs** has been deleted along with any **sub_record.annotations** key, value pairs. This is to reflect the fact the new sequence is no longer equivalent to the previous.
- To complicate things further, the id, name and description have been preserved in sub_record.
- It is advisable to update general attributes of the new SeqRecord following slicing.

### 4.8 Adding SeqRecord objects

- SeqRecord objects can be added to generate new SeqRecords using **+** operator. Note, similar to slicing SeqRecord objects, certain attributes are lost.

### 4.9 Reverse-complementing SeqRecord objects

- The SeqRecord **reverse_complement()** method calculates the reverse complement of a SeqRecord preserving SeqFeatures and modifying them accordingly. Again annotations and dxrefs are dropped.

## Chapter 5 Sequence Input/Output

- Additional help on the SeqIO module can be accessed as below or through the [wiki](http://biopython.org/wiki/SeqIO) etc.

```python
>>> from Bio import SeqIO
>>> help(SeqIO)
...
```

### 5.1 Parsing or Reading Sequences

- Majority of time will use the **SeqIO.parse()** function:

```pythonhelp
Help on function parse in module Bio.SeqIO:

parse(handle, format, alphabet=None)
    Turn a sequence file into an iterator returning SeqRecords.

    Arguments:
     - handle   - handle to the file, or the filename as a string
       (note older versions of Biopython only took a handle).
     - format   - lower case string describing the file format.
     - alphabet - optional Alphabet object, useful when the sequence type
       cannot be automatically inferred from the file itself
       (e.g. format="fasta" or "tab")
```

- To avoid yielding an iterable, the **SeqIO.read()** function will yield a single SeqRecord object when only one entry is present in **handle**.
- Depending on the number and size of the SeqRecord objects associated with handle, it can be desirable to load all SeqRecord objects into a list object simultaneously e.g.

```python
from Bio import SeqIO

records = list(SeqIO.parse("ls_orchid.gbk", "genbank"))

print("Found %i records" % len(records))
```

```pythonprint
Found 94 records
```

- SeqRecord objects are manipulated as described in [Chapter 4](#chapter-4-sequence-annotation-objects).

### 5.2 Parsing sequences from compressed files

- In addition to strings corresponding with file paths, SeqIO.parse() accepts file objects for **handle** values e.g.

```python
>>> from Bio import SeqIO
>>> with open("ls_orchid.gbk") as handle:
...     print(sum(len(r) for r in SeqIO.parse(handle, "gb")))
67518
```

- This facilitates the decompressing of compressed files using other 3rd party libraries:

```python
>>> import gzip
>>> from Bio import SeqIO
>>> with gzip.open("ls_orchid.gbk.gz", "rt") as handle:
...     print(sum(len(r) for r in SeqIO.parse(handle, "gb")))
...
67518
```

### 5.3 Parsing Sequences from the net

- Refer to Chapters 9 & 10 in the [biopython tutorial](http://biopython.org/DIST/docs/tutorial/Tutorial.html#htoc111).

### 5.4 Sequence files as Dictionaries

- Biopython provides 3 functions to convert iterables of SeqRecord objects to dictionary or dictionary-like objects:
  - **Bio.SeqIO.to_dict()**
  - **Bio.SeqIO.index()** 
  - **Bio.SeqIO.index_db()**
- An example of SeqIO.to_dict() is:

```python
>>> from Bio import SeqIO
>>> orchid_dict = SeqIO.to_dict(SeqIO.parse("ls_orchid.gbk", "genbank"))
```

- This is equivalent to:

```python
>>> from Bio import SeqIO
>>> orchid_dict = {}
>>> for rec in SeqIO.parse("ls_orchid.gbk", "genbank"):
>>>     if not orchid_dict.get(rec.id):
>>>         orchid_dict[rec.id] = rec
>>>     else:
>>>         raise...
```

- The if/else clause ensures SeqIO.to_dict() fails when multiple SeqRecord objects with the same id are found within the iterable.
- To configure dict keys, an optional 2nd argument is provided that accepts *function* objects:

```pythonhelp
Help on function to_dict in module Bio.SeqIO:

to_dict(sequences, key_function=None)
    Turn a sequence iterator or list into a dictionary.

    Arguments:
     - sequences  - An iterator that returns SeqRecord objects,
       or simply a list of SeqRecord objects.
     - key_function - Optional callback function which when given a
       SeqRecord should return a unique key for the dictionary.
```

- The following example illustrates this by assigning [SEGUIDs](https://onlinelibrary.wiley.com/doi/abs/10.1002/pmic.200600032) for each SeqRecord in the file handle:

```python
>>> from Bio import SeqIO
>>> from Bio.SeqUtils.CheckSum import seguid
>>> seguid_dict = SeqIO.to_dict(SeqIO.parse("ls_orchid.gbk", "genbank"),
...                             lambda rec : seguid(rec.seq))
>>> record = seguid_dict["MN/s0q9zDoCVEEc+k/IFwCNF2pY"]
>>> print(record.id)
Z78532.1
>>> record.description
'C.californicum 5.8S rRNA gene and ITS1 and ITS2 DNA'
```

> A key limitation of **SeqIO.to_dict()** is that all attributes etc. of SeqRecords are held in RAM. This can slow calculations for files containing large numbers of entries.

- An alternative **SeqIO.index()** which returns a dictionay-like object  where not all SeqRecords are in memory. Rather they are called into memory on demand as required. I suspect that for computers with SSD **SeqIO.index()** will be close to **SeqIO.to_dict()** in terms of speed.
- A limitation of SeqIO.index() is that handles to files cannot be used, only filenames (below). However, bgz compressed files can be used (see [indexing compressed files 5.4.4](http://biopython.org/DIST/docs/tutorial/Tutorial.html#htoc62))

```pythonhelp
Help on function index in module Bio.SeqIO:

index(filename, format, alphabet=None, key_function=None)
    Indexes a sequence file and returns a dictionary like object.

    Arguments:
     - filename - string giving name of file to be indexed...
```

- *The dictionary-like object returned by the SeqIO.index() function is technically an instance of the _IndexedSeqFileDict class. This class has a get_raw() method which can return an array of bytes. However, it doesn't seem clear what encoding is used by this function.*
- If you are working with an incredibly large files, **SeqIO.index_db** will store record information on disk in an SQLite3 database.

### 5.5 Writing Sequence Files

- Use the **SeqIO.write()** function:

```pythonhelp
Help on function write in module Bio.SeqIO:

write(sequences, handle, format)
    Write complete set of sequences to a file.

    Arguments:
     - sequences - A list (or iterator) of SeqRecord objects, or (if using
       Biopython 1.54 or later) a single SeqRecord.
     - handle    - File handle object to write to, or filename as string
       (note older versions of Biopython only took a handle).
     - format    - lower case string describing the file format to write.

    Note if providing a file handle, your code should close the handle
    after calling this function (to ensure the data gets flushed to disk).

    Returns the number of records written (as an integer).
```

#### 5.5.1 round trips

- Where a written file is the same as the unaltered parsed version.
- This may not happen, refer to [seqrecord_play.py](/scripts/seqrecord_play.py).

#### 5.5.2 Converting between sequence file formats

- Rather than 1st parsing and then writing, biopython provides the **SeqIO.convert()** function:

```pythonhelp
Help on function convert in module Bio.SeqIO:

convert(in_file, in_format, out_file, out_format, alphabet=None)
    Convert between two sequence file formats, return number of records.

    Arguments:
     - in_file - an input handle or filename
     - in_format - input file format, lower case string
     - out_file - an output handle or filename
     - out_format - output file format, lower case string
     - alphabet - optional alphabet to assume
```

#### 5.5.3 converting a file of sequences to their reverse complements

- List and generator comprehensions are powerful tools for batch converting SeqRecords. For instance:

```python
>>> from Bio import SeqIO
>>> records = (rec.reverse_complement(id="rc_"+rec.id, description = "reverse complement") \
...            for rec in SeqIO.parse("ls_orchid.fasta", "fasta") if len(rec)<700)
>>> SeqIO.write(records, "rev_comp.fasta", "fasta")
18
```

- Note in the above example the SeqIO.write() function prints the number of entries when writing from generators.
- A further impressive entry is given in section 20.1.3 where SeqRecords are translated and written to a seperate file.

### .6  Low level FASTA and FASTQ parsers

- Alternatives to SeqIO.parse() for FASTA FASTQ files are **SimpleFastaParser** or **FastqGeneralIterator** e.g.

```python
>>> from Bio.SeqIO.FastaIO import SimpleFastaParser
```

```pythonhelp
Help on function SimpleFastaParser in module Bio.SeqIO.FastaIO:

SimpleFastaParser(handle)
    Iterate over Fasta records as string tuples.

    For each record a tuple of two strings is returned, the FASTA title
    line (without the leading '>' character), and the sequence (with any
    whitespace removed). The title line is not divided up into an
    identifier (the first word) and comment or description.
```

-Refer to Chapter 20 in the online tutorial for further information.
