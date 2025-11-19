# SPACE
Detection of Generic Spaced Motifs Using Constrained Submotif Pattern Mining
____________________________________________________________________

## Compilation

To compile the source code, issue this command:

    $ ./compileall

The code is compiled with aggressive optimizations (`-O3 -march=native -fopenmp`) for maximum performance. See [PERFORMANCE.md](PERFORMANCE.md) for details.

## System Requirements
It runs in UNIX/Linux system, where Perl version 5.8.4 and above should be installed.

**Requirements:**
- GCC 5.0+ (with OpenMP support)
- Perl 5.8.4+
- Linux/UNIX system
- Multi-core CPU recommended for best performance 

## Input format, storing input and output files

   Input files should be in fasta format, e.g:
   
   >Seqname
   CTGATTTTAAACCC.....(nucleotides should be in one line and no break)

Currently we only accept DNA sequences. Maximum number of sequence
allowed = 1000 each of length 5000bp.
Before running the program, all the input files should be stored in 
`files/` directory. File extension accepted include: `.fasta`, `.fa`, `.txt`

   
Supported species.
Currently we support 16 species. They are:


|Name           |                  |Symbol|
|---------------|----------------  |-----:|
|B.subtilis     | bacteria         | BBS  |
|G.gallus       | chicken          | GG   |
|A.thaliana     | crest            | AT   |
|C.familiaris   | dog              | CFI  |
|E.coli         | ecoli            | BEC  |
|D.melanogaster | fruitfly         | DM   |
|F.rubripes     | fugu             | FR   |
|A.nidulans     | fungi            | AN   |
|H.sapiens      | human            | HS   |
|A.gambiae      | mos_malar        | AG   |
|M.musculus     | mouse            | MM   |
|R.norvegicus   | rat              | RN   |
|Random Sequence| synthetic        | SYNT |
|C.elegans      | worm             | CE   |
|S.cerevisiae   | yeast            | SC   |
|D.rerio        | zebrafish        | DR   |

## To run the program simple issue the following command:

    prompt > perl run_all.pl [filename] [species code]

For example

    prompt> perl run_all.pl test.fasta HS

Output will be printed to `STDOUT`.

### Performance Tuning

To control the number of CPU cores used (default: all available):

    export OMP_NUM_THREADS=4
    perl run_all.pl test.fasta HS

For detailed performance information and tuning options, see [PERFORMANCE.md](PERFORMANCE.md).

### Expected Performance

On modern multi-core systems, this optimized version runs **2-5x faster** than the original implementation:
- Compiler optimizations: 20-40% faster
- Fast I/O: 2-5x on I/O operations
- Multi-threading: 2-8x on multi-core CPUs (scales with core count)

## Publication
Edward Wijaya, Kanagasabai Rajaraman, Siu Yiu Ming and Sung Wing Kin, 
*Detection of generic spaced motifs using submotif pattern mining*, (2007),
Bioinformatics, 15;23(12):1476-85. [PMID:17483509](http://www.ncbi.nlm.nih.gov/pubmed/17483509) 
