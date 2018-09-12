# SQANTI2

Last Updated: 2018/09/12

Private repo for Liz's modified SQANTI. The original [SQANTI](https://bitbucket.org/ConesaLab/sqanti) is by Ana Conesa lab.

Until the SQANTI authors have finalized their agreement with Liz on how to integrate Liz's changes, the script names are modified to reflect this as a temporary working version.

For example, `sqanti_qc.py` is named currently `sqanti_qc2.py`.

## Prerequisite

* Perl
* Minimap2 

### Python-related libraries

* Python (2.7)
* pysam
* psutil
* bx-python
* BioPython
* BCBioGFF
* [cDNA_Cupcake](https://github.com/Magdoll/cDNA_Cupcake/wiki/Cupcake-ToFU:-supporting-scripts-for-Iso-Seq-after-clustering-step#install)

### R-related libraries

* R (>= 3.4.0)
* R packages (for `sqanti_qc2.py`): ggplot2, scales, reshape, gridExtra, grid, dplyr

## Installing Python dependencies

I recommend using Anaconda which makes installing all the Python packages much easier. If you already have Anaconda installed because you use [Iso-Seq3](https://github.com/PacificBiosciences/IsoSeq_SA3nUP/wiki/Tutorial:-Installing-and-Running-Iso-Seq-3-using-Conda) or [cDNA_Cupcake](https://github.com/Magdoll/cDNA_Cupcake/wiki/Cupcake-ToFU:-supporting-scripts-for-Iso-Seq-after-clustering-step), you can activate your current environment and directly go to step (4).

(1)  Here's the generic Anaconda installation for [Linux environment](http://docs.continuum.io/anaconda/install/#linux-install). Currently only Linux environment supported.

```
bash ~/Downloads/Anaconda2-5.2.0-Linux-x86_64.sh
export PATH=$HOME/anaconda5.2/bin:$PATH
conda -V
conda update conda
```

(2) Create a virutal environment. I will call it `anaCogent5.2`. Type `y` to agree to the interactive questions.

```
conda create -n anaCogent5.2 python=2.7 anaconda
source activate anaCogent5.2
```

(3) Once you have activated the virtualenv, you should see your prompt changing to something like this:

```
(anaCogent5.2)-bash-4.1$
```

(4) Install additional required libraries:

```
conda install -n anaCogent5.2 pysam
conda install -n anaCogent5.2 psutil
conda install -n anaCogent5.2 biopython
conda install -n anaCogent5.2 -c http://conda.anaconda.org/cgat bx-python
conda install -n anaCogent5.2 -c bioconda bcbiogff
```

If you don't already have [cDNA_Cupcake](https://github.com/Magdoll/cDNA_Cupcake/wiki/Cupcake-ToFU:-supporting-scripts-for-Iso-Seq-after-clustering-step#install) installed, you can do that now:

```
$ git clone https://github.com/Magdoll/cDNA_Cupcake.git
$ cd cDNA_Cupcake
$ python setup.py build
$ python setup.py install
```

## Running SQANTI2

Activate the Anaconda environment. Make sure minimap2 works. Add `cDNA_Cupcake/sequence` to `$PYTHONPATH`.

```
$ source activate anaCogent5.2
(anaCogent5.2)-bash-4.1$ export PYTHONPATH=$PYTHONPATH:<path_to>/cDNA_Cupcake/sequence/
(anaCogent5.2)-bash-4.1$ minimap2 --version
2.11-r797
```

### Input to SQANTI QC

* Iso-Seq output. Preferably already mapped to the genome and [collapsed to unique transcripts](https://github.com/Magdoll/cDNA_Cupcake/wiki/Cupcake-ToFU:-supporting-scripts-for-Iso-Seq-after-clustering-step#collapse). (FASTA/FASTQ)
* Reference annotation in GTF format. For example [GENCODE](https://www.gencodegenes.org/releases/current.html) or [CHESS](http://ccb.jhu.edu/chess/).
* Reference genome, in FASTA format. For example hg38. *Make sure your annotation GTF is based on the correct ref genome version!*

### Running SQANTI QC

The script usage is:

```
python sqanti_qc2.py [-t cpus] [--skipORF] <input_fasta> <annotation_gtf> <genome_fasta>
```

If you don't feel like running the ORF prediction part, use `--skipORF`. Just know that all your transcripts will be annotated as non-coding.


For example:

```
python sqanti_qc2.py -t 30 touse.rep.fasta ~/share/gencode/gencode.v28.annotation.gtf /pbi/dept/bifx/etseng/genomes/hg38/hg38.fa
```

### SQANTI QC output


