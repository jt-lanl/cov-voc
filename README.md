# BACKGROUND

This repository contains several python routines for analyzing and
tracking variants of concern in the SARS-CoV-2 virus, based on the
sequence alignments that are built at the [LANL
cov.lanl.gov](https://cov.lanl.gov) website, which in turn are based
on the sequences provided by GISAID (<https://www.gisaid.org/>).

The emphasis is on the spike protein, though some of the routines can
be applied as well to other proteins in the SARS-CoV-2 genome.

For most of these routines, the input is an aligned amino-acid
sequence file (usually, but not necessarily, in fasta format) whose
first sequence is taken as the *reference* sequence, typically
`NC_045512_spike_surface_glycoprotein` which is the ancestral Wuhan
strain.


# MAIN ANALYSIS PROGRAMS

## SHIVER

SHIVER (SARS CoV-2 Historically Identified Variants in Epitope Regions)
identifies sets of variant forms of the SARS CoV-2 virus with a
focus on just the NTD and RBD neutralizing antibody epitope regions of
the Spike protein, chosen to maximize coverage globally and/or on
separate continents, depending on which of several strategies is
employed.

The standalone `shiver-barplot` routine takes output of `shiver` and
produces a bar plot showing coverages for different continents and
different number of components.


## XSPIKE

XSPIKE (eXplore the SPIKE protein) does three main tasks:

* identifies the sites that have the highest entropy

* computes pairwise correlations among those sites, and produces "circle plots" and "heat maps"

* creates a digest of variants and a continent-wise count for each of those variants

# SECONDARY PROGRAMS

A small ecosystem of python scripts has been developed to support and
display the analysis done in the main routines `shiver` and `xspike`.

Note that I often refer to sequence files as "fasta files" since that
is the format I most often use. But all of these routines read
sequence files using the `readseq.py` module, and this permits a
number of file types (including fasta, mase, tbl, and raw sequences,
and gzip'd versions of these as well); the file type is inferred from
the extension of the file name.

## EMBERS, SPARKS, and CINDERS

`embers` is a display tool that creates colorful stacked barplots
showing how variant counts change over time.  These are like the "blue
wave" plots in our original D614G paper, but with many more variants
and many more colors.  The routine can also create simple line plots
on a logarithmic axis, which is useful for tracking variants that are
much rarer than the dominant variants (but that may be increasing in
time).

The variants used in `embers` are defined in a "color mutation" file,
described below

`sparks` was formerly called `embers_bynames` -- it makes the same
kind of stacked barplots (or lineplots) that embers makes, but based
on the pango lineage names (see `cov-lineages.org`) that are encoded
in the sequence names.  In fact, it works with both sequence files
(.fasta, .tbl, etc) and with simple name files (.nm indicates a list
of sequence names without the sequences.

Because the pango lineage names are somewhat Byzantine in their
structure, `sparks` reads a "lineage table" which has colors, names,
and regular expressions describing the range of pango names that
correspond to the given name.

`cinders` is an experimental code that uses both pango names and
sequences; currently it only considers how single site mutations (eg,
A222V) affect the various pango lineage evolution over time.

## PANGOCOMMONFORMS

The `pangocommonforms` code produces a report that identifies the
common sequence patterns associated with each pango lineage (see
`cov-lineages.org`). In addition to the most common forms, a consensus
form is also identified.  The consensus is defined site-by-site, and
constructs a sequences by taking the most common amino-acid at each
site. Usually the consensus is also the most common form, but that is
not always the case, and in some cases, the consensus itself does not
even appear in the database.

## FIXFASTA

`fixfasta` is a preprocessor routine that takes an input sequence file
(usually fasta, but can be tbl or fasta.gz or several other formats) and
applies various filters to clean up the file.

One "fix" is to identify all the columns for which the reference
sequence exhibits a dash (`-`) and to strip these columns from all the
sequences.  The reason for this is so that the `n`'th character in
each sequence corresponds to site number `n`. [Note that the current
generation of analysis tools in this package is able to handle
insertions, and does not assume that sequence string index matches
site number.]

`fixfasta` also has options for codon-aligning DNA sequences and
for translating DNA sequences in to amino-acid sequences.

## MATCHFASTA

`matchfasta` identifies sequences from a fasta file that matches a
given pattern (sequence pattern and/or geographical region and/or
pango lineage). With no pattern specified, it provides a conveinent
way to view fasta files (eg, subsets of sequences and/or subsets of
sites).

## MUT2FASTA

`mut2fasta` takes one or more mutant strings (either from a file or from the
command line), each of which looks something like "`[A222V,A262S,S494P,D614G]`",
and produces a fasta file, each sequence of which corresponds to a mutant
specified by the string, relative to the reference sequence, which is the
first sequence in the specified reference fasta file.  Note that there is no
checking that such sequences exist in nature; for that, you'll want to use `matchfasta`.

## FIXALIGN, TWEAKALIGN, and COALIGN

`fixalign` takes an alignment of either DNA or amino acid sequences
and identifies inconsistencies among the sequences.  Two subsequences
are inconsistent if their dash-removed strings are identical, but
their dash-included strings differ.  For amino-acid sequences it also
invokes some heuristics to improve manifestly poor sub-alignments.
The `fixalign` code and provide a list of tweaks that can by used (by
`tweakalign` to actually fix the sequences; or it can just go ahead
and apply those tweaks directly.

*Note that `fixalign` does not make alignments from unaligned
sequences, and poorly aligned sequences (especially if the poor
alignments are consistent) will not be improved.*

`tweakalign` reads alignment tweaks either from the command line or
from a file and applies them to the input sequences.  You can use the
tweaks suggested by fixalign, but you can also edit those tweaks to
taste and/or add some of your own.  A single tweak is a pair of mutation
strings; for example the pair

    [E156G,F157-,R158-] [E156-,F157-,R158G]

changes the string 'G--' at sites 156-158 into '--G'.  A tweak should not
affect the raw sequence, only how it is aligned with respect to the reference
sequence.

`coalign` takes a set of amino-acid sequenes that are assumed to be
well-aligned, and a set of codon-aligned DNA sequences from which the
amino-acid sequences were derived, and it then re-aligns the DNA
sequences to be consistent with their newly-aligned amino-acid
sequence counterpargs.

## APOBEC: APOCOUNT, APOPLOT

The routines `apocount` and `apoplot` analyze and plot the occurrence
of mutational patterns in DNA sequences that are associated with the
apobec enzyme. A surfeit of apobec-style mutations is indicative of
apobec presence. The characteristic apobec mutation is G->A in a
context where 'G' is part of a string '...*G*A...' or '...*G*G...'

For DNA sequences, the reverse-compliment is also indicative; thus:
C->T in a context where 'C' is part of a string '...T*C*...' or
'...C*C*...'.  (Note, these are the "loose" definitions of apobec
context, a tighter definition is also implemented as an option.)

## COMMONTYPES

`commontypes` takes an input fasta file, the first of which is a
reference sequence, and finds the most common sequences and then
express them in terms of mutation strings.  One can also list the
geographical regions where those strains are most commonly seen, and
one can request a few ISL numbers for each strain, so you can find
examples in the GISAID database.  By specifying a mutation pattern,
you can restrict consideration to a mutation pattern, and thereby
obtain counts for "variants of variants"

(Note the routine `countvariants` is now considered obsolete; use
`commontypes` instead.)

## USITES

`usites` is a fairly specialized routine.  It first computes the
highest-entropy sites for the global sequence set and for each of the
continents's sequence sets (which are subsets of the global set), and
then returns the *union* of those lists.  You specify the number of
sites you want in the final union list.  This list of sites is useful
as input to `xspike` if you want different runs over different
geographical regions to all be using a common set of sites.

## MKTABLE, MUTISL

These are specialized routines that are used in the creation of the
tables that appear in the LANL corner of the GISAID website.

## BIGHAMMING

Find sequences that are far (in hamming distance) from a reference sequence

# SOME USEFUL LIBRARIES

### sequtil

for reading (also writing, and performing basic manipulations to)
fasta sequence files

* also tbl, mase, and seq (raw); the routine
  `sequtil.read_seqfile(...)` will determine the format of the file
  based on the file extension.

* Note that gzip'd files also work so `-i sequencefile.fasta.gz` on
  the command line will also be automatically understood as a
  compressed fasta-formatted file, and it will be read without
  explicitly decompressing the file.  (Note that `*.xz` files
  are similarly supported.)

* Output files are also written according to their name, with output
  to `*.gz` and `*.xz` files also implemented.

* A recent addition is the 'pkl' and 'ipkl' filetypes -- this is a
  python pickle (and incremental pickle) file; the 'pkl' is a direct
  serialization of the `SequenceSample` array as one object; the
  'ipkl' serializes each sample separately. These are *not* archival
  formats, but can be used to speed up the file-reading time for the
  routines in this package.

* Another recent addition is the 'nm' filetype -- this is a list of the
  sequence names only, no actual sequences.  Since sequences names
  have a lot of meta-information, you can often do analysis using only
  the names (eg, see the `sparks` routine above).

* The sequences are read into a list of `SequenceSample` data types
  (containing a name and a sequence), but if you want a literal list
  be sure to use a command like `seqs=list(seqs)` because by default
  the reading and filtering routines return python *iterators*, which
  [depending on the usage scenario] can be much more memory efficient,
  but only allow a single pass through the data.

### mutants/spikevariants

These library routines are for managing and manipulating single-site and
multiple-site mutations, and for parsing mutation strings of the form
`[S13I,W152C,L452R,D614G]`.

### covid

The `covid` library contains covid-specific
routines and constants (eg, definition of where RBD and NTD regions
are).  It also contains a lot of miscellaneous routines that happen to
be used by various tools in this package

### colornames
is a simple module for translating common color names
into hash-hexcodes.

### intlist
The `intlist` library provides generic utilities for dealing with
lists of integers; for instance, it can create headers such as the
following, which indicate that the sites being displayed are 19, 20,
142, ..., 501.  If it's not obvious, you read *down* each column to
get the number of the site. (In the example below 'W' is at site 152.)

      1111111222222344444445
    124455555444556613577890
    902423678234362779278441
    TTGYWMEFRLALDSAVKNLSTESN

### verbose
encapsulates vprint functions that write messages based on user-set
level of verbosity; typical use:

       import verbose as v
       ...
       v.verbosity(1)
       ...
       v.vprint("read",n_sequences,"sequences from file")

as an alternative to:

       verbosity_level=1
       ...
       if verbosity_level > 0:
           print("read",n_sequences,"sequences from file",file=sys.stderr)
       ...
___

# COLOR MUTATION TABLE

Variants are defined in a color mutation table (eg,
`color-mutation-table.txt` file); the table is a list of variants,
with one variant per line.

    Format for this file:
      Comments begin with '#' and are ignored
      Each line has: <color> <mutation> <name> <etc>
      <color> is common color name, with no spaces (eg, DarkGreen not Dark Green), limited to X11_color_names
      <mutation> is of form [<ssm>,<ssm>,...<ssm>], possibly followed by an '!'
          <ssm> is a single site mutation of form <rchar><site><mchar>, where
                <rchar> = character in reference sequence
		          ('+' indicates an insertion after the site)
                <site> = integer site number
                <mchar> = character in mutated sequence
		          ('.' indicates any, '*' indicates any except <rchar>,
			   '_' indicates <rchar>, 'x' indicates blank [only allowed if <rchar> is '+'])
          '!' indicates an "exact" match; that means that except for the sites indicated by the mutation, the
              mutant sequence must agree with the reference sequence at every "relevant" site, where a "relevant"
              site is among the union of all the sites in all the mutation strings
      <name> is the name associated with the mutation (eg, Pango lineage)
      <etc>  is any further text; it is treated as a comment and is ignored

For example, two typical lines look like:

    Orange  [H69-,V70-,Y144-,N501Y,A570D,D614G,P681H,T716I,S982A,D1118H] B.1.1.7 UK VOC
    Fuchsia [D614G,Q677*]!                                               Near-Furin
    Green   [+214x]                                                      No-Insertion

Each line contains a `ColorName`, a `[MutationString]`, possibly a `!`
symbol, and then a `VariantName`.  These components are separated by
arbitrary whitespace.  In the first line above, the "UK
VOC" is ignored.  In the second line, `Q677*` indicates that the
mutation can have any character except Q at site 677.  The '!'
(after the ']') indicates that only mutations at sites 614 and 677 are
permitted; a sequence that disagrees with the reference sequence at
any other site will not be consistent with this pattern.
In the Green line, any sequence will match as long as there is no
insertion after site 214.

# LINEAGE TABLE

Typical lines in the lineage table looks like:

    Orange          Alpha   (B\.1\.1\.7)|(Q\..*)
    BlueViolet      Delta   (B\.1\.617\.2)|(AY\..*)

with three columns corresponding to color, name, and a regexp that matches the
various pango names associated with the name.  For instance, B.1.1.7 or Q.anything
will correspond to the Alpha variant, and will be plotted in orange.
___

### The 'help' command-line option

Most of the main python scripts take command line arguments, and
those that do take '-h' to provide a usage message; if nothing else,
it will at least list the names of the various options.  For example,
if you type

    python shiver.py -h

then you'll get something like the following:

    usage: shiver.py [-h] [--input INPUT]
                     [--filterbyname FILTERBYNAME [FILTERBYNAME ...]]
                     [--xfilterbyname XFILTERBYNAME [XFILTERBYNAME ...]]
                     [--keepdashcols] [--keeplastchar] [--dates DATES DATES]
                     [--days DAYS] [--fixsiteseventy] [--keepx] [--output OUTPUT]
                     [-n N] [--strategy STRATEGY] [--region REGION]
                     [--colormut COLORMUT] [--verbose]

    SHIVER: SARS CoV-2 Historically Identified Variants in Epitope Regions

    optional arguments:
      -h, --help            show this help message and exit
      --input INPUT, -i INPUT
                            input fasta file with aligned sequences (first is
                            master)
      --filterbyname FILTERBYNAME [FILTERBYNAME ...], -f FILTERBYNAME [FILTERBYNAME ...]
                            Only use sequences whose name matches this pattern
      --xfilterbyname XFILTERBYNAME [XFILTERBYNAME ...], -x XFILTERBYNAME [XFILTERBYNAME ...]
                            Do not use sequences whose name matches this pattern
      --keepdashcols        Do not strip columns with dash in ref sequence
      --keeplastchar        Do not strip final stop codon from end of sequences
      --dates DATES DATES, -d DATES DATES
                            range of dates (two dates, yyyy-mm-dd format)
      --days DAYS           Consider date range of DAYS days ending on the last
                            sampled date
      --fixsiteseventy      Sites 68-70 should be I--, not -I- or --I
      --keepx               Keep sequences that include bad characters, denoted X
      --output OUTPUT, -o OUTPUT
                            output fasta file with variant sequences
      -n N                  number of components in cocktail
      --strategy STRATEGY, -s STRATEGY
                            how to pick variants: (T)aketurns, (M)ostimproved,
                            (G)lobalonly
      --region REGION       region of spike sequence over which patterns are
                            defined
      --colormut COLORMUT   name of color mutation file
                            (mutation_string,lineage_name) are 2nd,3rd columns
      --verbose, -v         verbose


# REFERENCE

B. Korber, W. M. Fischer, S. Gnanakaran, H. Yoon, J. Theiler, W. Abfalterer, N. Hengartner, E. E. Giorgi, T. Bhattacharya, B. Foley, K. M. Hastie, M. D. Parker, D. G. Partridge, C. M. Evans, T. M. Freeman, T. I. de Silva, A. Angyal, R. L. Brown, L. Carrilero, L. R. Green, D. C. Groves, K. J. Johnson, A. J. Keeley, B. B. Lindsey, P. J. Parsons, M. Raza, S. Rowland-Jones, N. Smith, R. M. Tucker, D. Wang, M. D. Wyles, C. McDanal, L. G. Perez, H. Tang, A. Moon-Walker, S. P. Whelan, C. C. LaBranche, E. O. Saphire, and D. C. Montefiori. "Tracking changes in SARS-CoV-2 Spike: evidence that D614G increases infectivity of the COVID-19 virus." Cell 182 (2020) 812-827. 

# COPYRIGHT

(c) 2021-2022. Triad National Security, LLC. All rights reserved.

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

# LICENSE

This program is open source under the BSD-3 License.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are
met:

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
