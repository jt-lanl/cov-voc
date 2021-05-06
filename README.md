# BACKGROUND

This repository contains several python routines for analyzing and
tracking variants of concern in the SARS-CoV-2 virus, based on the
sequence alignments that are built at the [LANL
cov.lanl.gov](https://cov.lanl.gov) website, which in turn are based
on the sequences provided by GISAID (<https://www.gisaid.org/>).

The emphasis is on the spike protein, though some of the routines can be applied
as well to other proteins in the SARS-CoV-2 genome.

For most of these routines, the input sequence file (usually, but not necessarily, fasta) file whose
first sequence is taken as the *reference* sequence, typically `NC_045512_spike_surface_glycoprotein` which is the ancestral Wuhan strain.

# MAIN PROGRAMS

## SHIVER

SHIVER (SARS CoV-2 Historically Identified Variants in Epitope Regions) provides a list of variants
calculated to cover as much of the variation as possible, while focusing on just the NTD and RBD
SHIVER identifies sets of variant forms of the SARS CoV-2 virus with a
focus on just the NTD and RBD neutralizing antibody epitope regions of
the Spike protein, chosen to maximize coverage globally and/or on
separate continents[*], depending on which of several strategies is
employed.

## XSPIKE

XSPIKE (eXplore the SPIKE protein) does does three main tasks:

* identifies the sites that have the highest entropy

* computes pairwise correlations among those sites, and produces "circle plots" and "heat maps"

* creates a digest of variants and a continent-wise count for each of those variants

## EMBERS

EMBERS (not an acronym...yet!) creates colorful stacked barplots that show how variant counts vary over time.
These are like the "blue wave" plots in our original D614G paper, but with many more variants and many more colors.

The variants used in `embers` are defined in a "color mutation" file, two typical lines of which look like:

    DarkGreen [N501Y,D614G]! N501Y
    Orange [H69-,V70-,Y144-,N501Y,A570D,D614G,P681H,T716I,S982A,D1118H] B.1.1.7 UK VOC

Each line contains a `ColorName`, a `[MutationString]`, possibly a `!` symbol, and then a `VariantName`.
In the second line above, the "UK VOC" is ignored.

# SECONDARY PROGRAMS

## FIXFASTA

`fixfasta` is a preprocessor routine that takes a fasta file (doesn't
have to be fasta, can be tbl or fasta.gz or several other formats) and
applies various filters to clean up the file.  An important one is to
identify all the columns for which the reference sequence exhibits a
dash (`-`) and to strip these columns from all the sequences.  The
reason for this is so that the `n`'th character in each sequence
corresponds to site number `n`.

## MUT2FASTA

`mut2fasta` takes one or more mutant strings (either from a file or from the
command line), each of which looks something like "`[A222V,A262S,S494P,D614G]`",
and produces a fasta file, each sequence of which corresponds to a mutant
specified by the string, relative to the reference sequence, which is the
first sequence in the specified reference fasta file.

## COMMONTYPES

given an input fasta file, the first of which is a reference sequence,
find the most common sequences, and express them in terms of mutation strings

## USITES

computes the highest-entropy sites for the global sequence set and for each of the
continents's sequence sets (which are subsets of the global set), and then returns
the union of those lists.  you specify the number of sites you want in the final union
list.  this list of sites is useful as input to `xspike` if you want different runs
over different geographical regions to all be using a common set of sites.

___

# SOME USEFUL LIBRARIES

### readseq/sequtil 

for reading fasta sequence files 

* also tbl, and several other formats, are supported;
the routine `readseq.read_seqfile(...)` will determine the format of the file based on the file extension.

* Note that gzip'd files also work so `-i sequencefile.fasta.gz` on the command line will also be automatically understood as a compressed fasta-formatted file, and it will be read without explicitly decompressing the file.  

### mutants/spikevariants 

manipulates single-site and multiple-site mutations and parses mutation strings of the form
`[S13I,W152C,L452R,D614G]`

### covid

contains a lot of hard-coded covid-specific routines and constants (eg, definition of where RBD and NTD regions are)

### colornames
simple module for translating common color names into hash-character names

### intlist
generic utilities for dealing with lists of integers; for instance, can create headers such as the following,
which indicate that the sites being displayed are 19, 20, 142, ..., 501.  If it's not obvious, you read *down* each
column to get the number of the site.
    
      1111111222222344444445
    124455555444556613577890
    902423678234362779278441
    TTGYWMEFRLALDSAVKNLSTESN
___

### The 'help' command-line option

Most of the main python scripts take command line arguments, and
those that do take '-h' to provide a usage message; if nothing else,
it will at least list the names of the various options.  For example,
if you type

    python shiver.py -h

then you'll get something like the following:
    
    usage: shiver.py [-h] [--input INPUT] [--filterbyname FILTERBYNAME]
                     [--xfilterbyname XFILTERBYNAME] [--keepdashcols]
                     [--dates DATES DATES] [--daysago DAYSAGO] [--output OUTPUT]
                     [-n N] [--strategy STRATEGY] [--region REGION] [--verbose]
    
    SHIVER: SARS CoV-2 Historically Identified Variants in Epitope Regions
    
    optional arguments:
      -h, --help            show this help message and exit
      --input INPUT, -i INPUT
                            input fasta file with aligned sequences (first is
                            master)
      --filterbyname FILTERBYNAME, -f FILTERBYNAME
                            Only use sequences whose name matches this pattern
      --xfilterbyname XFILTERBYNAME, -x XFILTERBYNAME
                            Do not use sequences whose name matches this pattern
      --keepdashcols        Do not strip columns with dash in ref sequence
      --dates DATES DATES, -d DATES DATES
                            range of dates (two dates, yyyy-mm-dd format)
      --daysago DAYSAGO     Consider date range from DAYSAGO days ago until today
      --output OUTPUT, -o OUTPUT
                            output fasta file with variant sequences
      -n N                  number of components in cocktail
      --strategy STRATEGY, -s STRATEGY
                            how to pick variants: (T)aketurns, (M)ostimproved,
                            (G)lobalonly
      --region REGION       region of spike sequence over which patterns are
                            defined
      --verbose, -v         verbose
