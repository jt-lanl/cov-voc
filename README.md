# BACKGROUND

This repository contains several python routines for analyzing and
tracking variants of concern in the SARS-CoV-2 virus, based on the
sequence alignments that are built at the [LANL
cov.lanl.gov](https://cov.lanl.gov) website, which in turn are based
on the sequences provided by GISAID (<https://www.gisaid.org/>).

The emphasis is on the spike protein, though some of the routines can be applied
as well to other proteins in the SARS-CoV-2 genome.

For most of these routines, the input is a sequence file (usually, but not necessarily, in fasta format) whose
first sequence is taken as the *reference* sequence, typically `NC_045512_spike_surface_glycoprotein` which is the ancestral Wuhan strain.

# MAIN ANALYSIS PROGRAMS

## SHIVER

SHIVER (SARS CoV-2 Historically Identified Variants in Epitope Regions)
identifies sets of variant forms of the SARS CoV-2 virus with a
focus on just the NTD and RBD neutralizing antibody epitope regions of
the Spike protein, chosen to maximize coverage globally and/or on
separate continents, depending on which of several strategies is
employed.

## XSPIKE

XSPIKE (eXplore the SPIKE protein) does three main tasks:

* identifies the sites that have the highest entropy

* computes pairwise correlations among those sites, and produces "circle plots" and "heat maps"

* creates a digest of variants and a continent-wise count for each of those variants

# SECONDARY PROGRAMS

## EMBERS

EMBERS is a display tool that creates colorful stacked barplots showing how variant counts change over time.
These are like the "blue wave" plots in our original D614G paper, but with many more variants and many more colors.
The routine can also create simple line plots on a logarithmic axis, which is useful for tracking variants that
are much rarer than the dominant variants (but that may be increasing in time).

The variants used in `embers` are defined in a "color mutation" file, described below

## FIXFASTA

`fixfasta` is a preprocessor routine that takes an input sequence file
(usually fasta, but can be tbl or fasta.gz or several other formats) and
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

## VFASTA

for viewing fasta files; you can just look at subsets of sequences and/or subsets of sites

## COMMONTYPES

given an input fasta file, the first of which is a reference sequence,
find the most common sequences, and express them in terms of mutation strings

## USITES

computes the highest-entropy sites for the global sequence set and for each of the
continents's sequence sets (which are subsets of the global set), and then returns
the union of those lists.  you specify the number of sites you want in the final union
list.  this list of sites is useful as input to `xspike` if you want different runs
over different geographical regions to all be using a common set of sites.

## COUNTVARIANTS

counts how many variants of a variant appear in an input set of sequences

## SHIVER-BARPLOT

takes output of `shiver` and produces a bar plot showing coverages for different continents
and different number of components. 

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
simple module for translating common color names into hash-hexcodes

### intlist
generic utilities for dealing with lists of integers; for instance, can create headers such as the following,
which indicate that the sites being displayed are 19, 20, 142, ..., 501.  If it's not obvious, you read *down* each
column to get the number of the site.
    
      1111111222222344444445
    124455555444556613577890
    902423678234362779278441
    TTGYWMEFRLALDSAVKNLSTESN
___

# COLOR MUTATION TABLE

Variants are defined in a color mutation table (eg, `color-mutation-table.txt` file); the
table is a list of variants, with one variant per line.

    Format for this file:
      Comments begin with '#' and are ignored
      Each line has: <color> <mutation> <name> <etc>
      <color> is common color name, with no spaces (eg, DarkGreen not Dark Green), limited to X11_color_names
      <mutation> is of form [<ssm>,<ssm>,...<ssm>], possibly followed by an '!'
          <ssm> is a single site mutation of form <rchar><site><mchar>, where
                <rchar> = character in reference sequence
                <site> = integer site number
                <mchar> = character in mutated sequence ('.' indicates any, '!' indicates any except <rchar>)
          '!' indicates an "exact" match; that means that except for the sites indicated by the mutation, the
              mutant sequence must agree with the reference sequence at every "relevant" site, where a "relevant"
              site is among the union of all the sites in all the mutation strings
      <name> is the name associated with the mutation (eg, Pango lineage)
      <etc>  is any further text; it is treated as a comment and is ignored

For example, two typical lines look like:

    Orange  [H69-,V70-,Y144-,N501Y,A570D,D614G,P681H,T716I,S982A,D1118H] B.1.1.7 UK VOC
    Fuchsia [D614G,Q677!]!                                               Near-Furin
    
Each line contains a `ColorName`, a `[MutationString]`, possibly a `!`
symbol, and then a `VariantName`.  These components are separated by
arbitrary whitespace.  In the first line above, the "UK
VOC" is ignored.  In the second line, Q677! indicates that the
mutation can have any character except Q at site 677.  The second '!'
(after the ']') indicates that only mutations at sites 614 and 677 are
permitted; a sequence that disagrees with the reference sequence at
any other site will not be consistent with this pattern.

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


# COPYRIGHT (for SHIVER and XSPIKE)

(c) 2021. Triad National Security, LLC. All rights reserved.

This program was produced under U.S. Government contract
89233218CNA000001 for Los Alamos National Laboratory (LANL), which is
operated by Triad National Security, LLC for the U.S.  Department of
Energy/National Nuclear Security Administration. All rights in the
program are reserved by Triad National Security, LLC, and the
U.S. Department of Energy/National Nuclear Security
Administration. The Government is granted for itself and others acting
on its behalf a nonexclusive, paid-up, irrevocable worldwide license
in this material to reproduce, prepare derivative works, distribute
copies to the public, perform publicly and display publicly, and to
permit others to do so.

# LICENSE (for SHIVER and XSPIKE)

This program is open source under the BSD-3 License.

Redistribution and use in source and binary forms, with or without modification, are permitted
provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright
   notice, this list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright
   notice, this list of conditions and the following disclaimer in the
   documentation and/or other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its
   contributors may be used to endorse or promote products derived from
   this software without specific prior written permission.


THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
"AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

# LANL C-number

C21022 SHIVER was approved for Open Source Assertion on 03/15/2021

C21029 xSpike was approved for Open Source Assertion on 06/15/2021

