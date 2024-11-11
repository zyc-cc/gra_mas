
If you however want to compile GraphAligner yourself, run these:

- Install miniconda https://conda.io/projects/conda/en/latest/user-guide/install/index.html
- `git clone https://github.com/zyc-cc/gra_mas.git`
- `cd GraphAligner`
- `git submodule update --init --recursive`
- `conda env create -f CondaEnvironment_linux.yml` or `conda env create -f CondaEnvironment_osx.yml`
- `source activate GraphAligner`
- `make bin/GraphAligner`

gatb-core and minimap2 need to be compiled manually.


GraphAligner -g short_reads.gfa -f long_reads.fa --corrected-out corrected.fa -x dbg -r 0.8 --input-reads `

### Parameters

- `-g` input graph. Format .gfa / .vg/.fa
- `-f` input reads. Format .fasta / .fastq / .fasta.gz / .fastq.gz. You can input multiple files with `-f file1 -f file2 ...` or `-f file1 file2 ...`
- `-t` number of aligner threads. The program also uses two IO threads in addition to these.
- `-a` output file name. Format .gam or .json
- `-x` parameter preset. Use `-x vg` for aligning to variation graphs and other simple graphs, and `-x dbg` for aligning to de Bruijn graphs.

All parameters below are optional.

- `--precise-clipping` use arg as the identity threshold for a valid alignment. Recommended to be less than the accuracy of the reads, for example 0.75 for ONT, 0.9 for HiFi, 0.95 for assembly-to-assembly.
- `--min-alignment-score` discard alignments whose score is less than this.
- `--multimap-score-fraction` alignment score fraction for including secondary alignments. Alignments whose alignment score is less than arg as a fraction of the best scoring overlapping alignment per read are discarded. Lower values include more poor secondary alignments and higher values less.

Seeding:

- `--seeds-minimizer-density` For a read of length `n`, use the `arg * n` most unique seeds
- `--seeds-minimizer-length` k-mer size for minimizer seeds
- `--seeds-minimizer-windowsize` Window size for minimizer seeds
- `--seeds-mum-count` MUM seeds. Use the n longest maximal unique matches. -1 for all MUMs
- `--seeds-mem-count` MEM seeds. Use the n longest maximal exact matches. -1 for all MEMs
- `--seeds-mxm-length` MUM/MEM minimum length. Don't use MUMs/MEMs shorter than n
- `--seeds-mxm-cache-prefix` MUM/MEM file cache prefix. Store the MUM/MEM index into disk for reuse. Recommended unless you are sure you won't align to the same graph multiple times

Extension:

- `-b` alignment bandwidth. Unlike in linear alignment, this is the score difference between the minimum score in a row and the score where a cell falls out of the band. Values recommended to be between 1-35.
- `-C` tangle effort. Determines how much effort GraphAligner spends on tangled areas. Higher values use more CPU and memory and have a higher chance of aligning through tangles. Lower values are faster but might return an inoptimal or a partial alignment. Use for complex graphs (eg. de Bruijn graphs of mammalian genomes) to limit the runtime in difficult areas. Values recommended to be between 1'000 - 500'000.
- `-k` k-mer size
- `--input-reads` -g the input file is a short read file
- `-r` msa filter threshold(default:0.08)，for large datasets the value should be larger or the effect of msa will be very weak
