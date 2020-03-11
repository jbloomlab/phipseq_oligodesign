# phipseq_oligodesign

Phage immuprecipitation sequencing (PhIP-Seq) is an experimental method to characterize viral exposure history from serum samples. 
This system has previously been termed "VirScan" and is described in [this](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4844011/) 2015 paper by Xu, et al.
In PhIP-Seq, peptides that fully cover the proteome of the target viruses are displayed on the surface of phage.
Those phage are then immunoprecipitated with serum and the phage pulled-down by the serum are sequenced. 
This method identifies what peptides the serum reacts with, which can be used to identify what viruses the serum contains antibodies against, and, thus, what viruses the individual whose serum was tested was exposed to.

The initial step in any PhIP-Seq experiment is designing the DNA oligos to encode the peptides to be displayed on the surface of the phage. 
These oligos encode overlapping peptides that tile along the complete proteome of the virsues of interest.
These oligos also include 5' and 3' adaptor sequences for cloning into the T7 phage for library creation.
This repository contains a script (`phipseq_oligodesign.py`) for desiging oligos for phage immunoprecipitation and sequencing (PhIP-Seq) experiments. 
This script was written by Kate Crawford in [the Bloom lab](https://research.fhcrc.org/bloom/en.html).

### Using `phipseq_oligodesign.py`

Note: This repository uses a conda environment to install necessary packages. 
To install this environment, you must have conda installed. 
If you do not, please look into installing `miniconda` [here](https://docs.conda.io/en/latest/miniconda.html).

First, clone this repository to your local computer using:

    > git clone https://github.com/jbloomlab/phipseq_oligodesign.git

Then, from within the `phipseq_oligodesign` directory, create the conda environment and activate it by running:

    > conda env create -f environment.yml

    > conda activate phipseq

Now, everything should be set-up to run the script and design oligos.

In order to design oligos, all protein or dna sequences for design oligos from must be separate `.fasta` files within a single directory. No other files should be present in this input directory and all files must be of one sequence type (either protein or dna). 

To run the `phipseq_oligodesign.py` script run:
    
    > python phipseq_oligoodesign.py input_dir output_dir seq_type

Where `seq_type` is either 'protein' or 'dna'.

Additional settings (such as oligo length, overlap length, etc.) can be specified using optional arguments. For additional information on these arguments and to see their default values run:

    > python phipseq_oligodesign.py -h


### Output

The `phipseq_oligodesign.py` script creates 4 output files. 

The first three files are all `csv`s with the following columns:
* `Index`
* `Name`: Name of the sequence file the oligo comes from
* `Oligo`: DNA sequence for the oligo, including adaptor sequences
* `Prot`: Protein sequence for the oligo (not including the adaptor sequences)
* `Prot_Start`: Site in the full-length amino acid sequence where this oligo starts

The first three files are:
* `oligos_all.txt`: a csv containing all designed oligos
* `oligos_nodups.txt`: a csv containing designed oligos with no duplicate protein sequences. If there are duplicate sequences, the first oligo is kept.
* `oligos_dups.txt`: a csv containing all duplicated oligos.

The fourth file is:
* `oligos_dup_counts.txt`: a tsv containing information about how many times each duplicated sequence was duplicated and which sequences it came from.
