::

                        ||||||            ||        ||
                      ||    ||                      ||
                      ||        ||    ||  ||    ||||||    ||||
                      ||  ||||  ||    ||  ||  ||    ||  ||    ||
                      ||    ||  ||    ||  ||  ||    ||  ||    ||
                        ||||||    ||||||  ||    ||||||    ||||

Installation - option 1 (Conda)
-------------------------------

Make sure you have have Conda installed on your computer.

Create conda environment
~~~~~~~~~~~~~~~~~~~~~~~~

::

    conda env create -f environment.yml

Environment with Guido
~~~~~~~~~~~~~~~~~~~~~~

After creating the environment you can activate it by using:

::

    conda activate guido

To get out of the environment use:

::

    conda deactivate

Installation - option 2 (pip)
-----------------------------

Install dependencies with pip
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

::

    pip install -r requirements.txt

Download *Anopheles gambiae* data
---------------------------------

Download ZIP from:
https://imperialcollegelondon.box.com/v/guido-agam-data and extract the
archive into ``guido/`` folder.

Usage
-----

::

    usage: guido [-h] [--sequence-file SEQUENCE] [--region REGION] [--gene GENE]
                 [--variants VARIATION_STORE] [--conservation CONSERVATION_STORE]
                 [--pam PAM] [--threads N_THREADS]
                 [--max-flanking MAX_FLANKING_LENGTH]
                 [--min-flanking MIN_FLANKING_LENGTH]
                 [--length-weight LENGTH_WEIGHT] [--max-offtargets MAX_OFFTARGETS]
                 [--n-patterns N_PATTERNS] [--disable-mmej [DISABLE_MMEJ]]
                 [--disable-off-targets [DISABLE_OFFTARGETS]]
                 [--output-folder OUTPUT_FOLDER] [--feature-type FEATURE]
                 [--dump DUMP]

    Guido is a fashionable gRNA designer, specialised in gene drives.

    optional arguments:
      -h, --help            show this help message and exit
      --region REGION, -r REGION
                            Region in AgamP4 genome [2L:1530-1590].
      --gene GENE, -G GENE  Genome of interest (AgamP4.7 geneset).
      --variants VARIATION_STORE, -v VARIATION_STORE
                            VCF file with variants.
      --conservation CONSERVATION_STORE, -c CONSERVATION_STORE
                            Path to Zarr store with conservation data.
      --pam PAM, -P PAM     Protospacer adjacent motif (IUPAC format)
      --threads N_THREADS, -t N_THREADS
                            Number of threads used.
      --max-flanking MAX_FLANKING_LENGTH, -M MAX_FLANKING_LENGTH
                            Max length of flanking region.
      --min-flanking MIN_FLANKING_LENGTH, -m MIN_FLANKING_LENGTH
                            Min length of flanking region.
      --length-weight LENGTH_WEIGHT, -w LENGTH_WEIGHT
                            Length weight - used in scoring.
      --max-offtargets MAX_OFFTARGETS
                            Max number of reported offtargets
      --n-patterns N_PATTERNS, -p N_PATTERNS
                            Number of MH patterns used in guide evaluation.
      --disable-mmej [DISABLE_MMEJ]
                            Disable MMEJ prediction.
      --disable-off-targets [DISABLE_OFFTARGETS]
                            Disable off-targets search.
      --output-folder OUTPUT_FOLDER, -o OUTPUT_FOLDER
                            Output folder.
      --feature-type FEATURE, -f FEATURE
                            Type of genomic feature to focus guide search on.
      --dump DUMP           Dump pickled cut_sites object to the output folder.

Output
------

List of guides
~~~~~~~~~~~~~~

::

    guide_sequence  genomic_location    exon_name   strand  off_target_analysis MMEJ_score  MMEJ_sum_score  MMEJ_top_score  SNP_count   wt_prob SNP_info    MMEJ_out_of_frame_del
    GGAGCAGTTCAGCAGCGCGGCGG X:1300558-1300580   AGAP000080-RD-E7    -   [0, 0, 0, 1]    100.0   2082.9  634.7   1   0.9762  1300563:G/[T]([0.0238]) +   +   +   +   +
    GCTGCCGCTGCCCGGAGAGACGG X:1300426-1300449   AGAP000080-RD-E7    +   [0, 0, 0, 0]    100.0   1649.2  603.0   0   1.0     +   +   +   +   +
    GGTTGCCGCCTCCGTCTCTCCGG X:1300438-1300460   AGAP000080-RD-E7    -   [0, 0, 0, 0]    100.0   1649.2  603.0   0   1.0     +   +   +   +   +
    GTTGCCGCCTCCGTCTCTCCGGG X:1300437-1300459   AGAP000080-RD-E7    -   [0, 0, 0, 0]    100.0   1515.2  469.0   0   1.0     +   +   +   +   +
    GTTCGGAATGGGAGGGAAGTCGG X:1300498-1300521   AGAP000080-RD-E7    +   [0, 0, 0, 1]    100.0   1484.3  352.5   0   1.0     +   +   +   +   +
    CTGACCTTTGAACGGCTCTCAGG X:1300905-1300927   AGAP000080-RD-E7    -   [0, 0, 0, 0]    100.0   1436.9  402.0   0   1.0     +   +   +   +   +
    GCCTCCGTCTCTCCGGGCAGCGG X:1300431-1300453   AGAP000080-RD-E7    -   [0, 0, 0, 0]    100.0   1409.4  423.0   0   1.0     +   +   +   +   +
    GGGCTTGAAGCTGCTGTTTGCGG X:1300621-1300644   AGAP000080-RD-E7    +   [0, 0, 0, 0]    100.0   1383.5  352.5   0   1.0     +   +   +   +   +
    ACCTTTGAACGGCTCTCAGGCGG X:1300902-1300924   AGAP000080-RD-E7    -   [0, 0, 0, 0]    100.0   1355.6  402.0   0   1.0     +   +   +   +   +
    CGCTTCTGGTGTAGGTCGCCCGG X:1300779-1300802   AGAP000080-RD-E7    +   [0, 0, 0, 0]    100.0   1019.5  268.0   0   1.0     +   +   +   +   +

Detailed list
~~~~~~~~~~~~~

::

    Guide: GGAGCAGTTCAGCAGCGCGGCGG  Location: X:1300558-1300580 Strand: -
    MMEJ score: 100.0   MMEJ sum score: 2082.9  MMEJ top score: 634.7

    Top MMEJ patterns
    Pattern Score   Deletion size   Produces out-of-frame deletion  MH seq  Deletion seq
    GCGGTCGCGAACAGCTTCAAGTCGGAGCAGTTCA------+++++GCAGCGGTTGCCAATTACGCGCTCGGGACCTTCAA    634.7   11  +   GCAGCG  GCAGCGCGGCG
    GCGGTCGCGAACAGCTTCAAGTCGGAGCAGTTCAGCA---++GCGGCAGCGGTTGCCAATTACGCGCTCGGGACCTTCAA    467.4   5   +   GCG GCGCG
    GCGGTCGCGAACAGCTTCAAGTCGGAGCAGTTCAGCAG--CGGCGGCAGCGGTTGCCAATTACGCGCTCGGGACCTTCAA    362.0   2   +   CG  CG
    GCGGTCGCGAACAGCTTCAAGTCGGAGCAGTT--------++++++CAGCGGTTGCCAATTACGCGCTCGGGACCTTCAA    347.9   14  +   CAGC    CAGCAGCGCGGCGG
    GCGGTCGCGAACAGCTTCAAGTCGGA--------------+++++GCAGCGGTTGCCAATTACGCGCTCGGGACCTTCAA    270.9   19  +   GCAG    GCAGTTCAGCAGCGCGGCG

    Variants
    Position    Ref/Alt Frequency
    1300563 G/[T]   0.0238

    WT allele probability: 0.9762

    Off-targets
    Number off-targets with
    [0, 1, 2, 3] mismatches
    [0, 0, 0, 1]

    Chromosome  Start position  Strand  Mismatches
    2L:33877342 +   1: T>G, 2: T>A, 15: G>C
