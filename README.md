# SQANTI2

Last Updated: 07/24/2019

## What is SQANTI2

SQANTI2 is an updated version of SQANTI ([publication](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5848618/), [code repository](https://bitbucket.org/ConesaLab/sqanti)). 


New features implemented in SQANTI2 not available in SQANTI:

* NMD detection -- new field `predicted_NMD` in classification output.
* Intron retention --- marked with `intron_retention` in `subcategory` field in classification output.
* CAGE peak --- new fields `dist_peak` and `within_peak` in classification output. Must provide CAGE peak data.
* polyA motif --- new field `polyA_motif` in classification output. Must provide polyA motif list.
* CDS-annotated GFF --- SQANTI2 outputs a `xxxx.cds.gff` GFF file that annotates CDS regions.


![sqanti2workflow](https://github.com/Magdoll/images_public/blob/master/SQANTI2_figures/sqanti2_workflow.png)

## Updates

2019.07.25 updated to version 3.5. Checks for Cupcake version (8.1+)

2019.07.24 updated to versoin 3.4. Now `--gtf` input option works with collapsed GFF format.

2019.07.23 updated to version 3.3. Added CDS for GFF support and IR fix.

2019.07.19 updated to version 3.2. `sqanti_qc2.py` fusion mode surpressed ORF prediction for the meantime. Minor mod to gmst to work on small input.

2019.07.17 updated to version 3.1. `sqanti_qc2.py` now supports GTF input `--gtf` again. Minor change to NIC subtype categorization naming.

2019.07.16 updated to version 3.0. now use Bioconda install of `gtfToGenePred` and `gffread`.

2019.07.12 updated to version 2.9. `sqanti_qc2.py` now annotates NMD prediction.

2019.06.19 updated to version 2.8. `sqanti_qc2.py` now works with fusion transcripts (must have ID `PBfusion.X`) using the `--is_fusion` option

2019.06.02 updated to version 2.7. Fixed GMAP option bug + added distance to closest annotated start/end for the gene (not ref isoform) and filtering afterwards.

2019.05.02 updated to version 2.6. Added deSALT aligner support using `--aligner=deSALT`. 

2019.03.18 minor typo fixed for version 2.5. updated doc for `sqanti_filter2.py`

2019.03.10 updated to version 2.5. Fixed `sqanti_filter2.py` missing fusion category also using polyA_motif as part of filtering.

2019.02.27 updated to version 2.4.  Added polyA motif finding.

2019.02.27 updated to version 2.3. `junction_category` fixed to check for (ss5,ss3) pairs in provided GTF.

2019.02.26 updated to version 2.2. added support for CAGE peak (FANTOM5) and Intropolis junction BED. 

2018.10.15 updated to version 1.1. modified use of SAM to GFF with added `source` parameter.

Private repo for Liz's modified SQANTI. The original [SQANTI](https://bitbucket.org/ConesaLab/sqanti) is by Ana Conesa lab.

![](https://github.com/Magdoll/images_public/blob/master/github_isoseq3_wiki_figures/wiki_SQANTI_sample_output.png)

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

(2) Create a virutal environment. I will call it `anaCogent3`. Type `y` to agree to the interactive questions.

```
conda create -n anaCogent3 python=2.7 anaconda
source activate anaCogent3
```

(3) Once you have activated the virtualenv, you should see your prompt changing to something like this:

```
(anaCogent3)-bash-4.1$
```

(4) Install additional required libraries:

```
conda install -n anaCogent3 -c bioconda pysam
conda install -n anaCogent3 psutil
conda install -n anaCogent3 biopython
conda install -n anaCogent3 -c http://conda.anaconda.org/cgat bx-python
conda install -n anaCogent3 -c bioconda bcbiogff
```

We also need to install [gtfToGenePred](https://bioconda.github.io/recipes/ucsc-gtftogenepred/README.html) and [gffread](https://bioconda.github.io/recipes/gffread/README.html).

```
conda install -n anaCogent3 -c bioconda ucsc-gtftogenepred openssl=1.0
conda install -n anaCogent3 -c bioconda gffread
```

If you don't already have [cDNA_Cupcake](https://github.com/Magdoll/cDNA_Cupcake/wiki/Cupcake-ToFU:-supporting-scripts-for-Iso-Seq-after-clustering-step#install) installed, you can do that now:

```
$ git clone https://github.com/Magdoll/cDNA_Cupcake.git
$ cd cDNA_Cupcake
$ python setup.py build
$ python setup.py install
```

No installation for SQANTI2 itself is required. The scripts can be run directly.


## Running SQANTI2

Activate the Anaconda environment. Make sure minimap2 works. Add `cDNA_Cupcake/sequence` to `$PYTHONPATH`.

```
$ source activate anaCogent3
(anaCogent3)-bash-4.1$ export PYTHONPATH=$PYTHONPATH:<path_to>/cDNA_Cupcake/sequence/
(anaCogent3)-bash-4.1$ minimap2 --version
2.15-r905
```

#### Input to SQANTI QC

* *Iso-Seq output*. Preferably already mapped to the genome and [collapsed to unique transcripts](https://github.com/Magdoll/cDNA_Cupcake/wiki/Cupcake-ToFU:-supporting-scripts-for-Iso-Seq-after-clustering-step#collapse). (FASTA/FASTQ)
* *Reference annotation* in GTF format. For example [GENCODE](https://www.gencodegenes.org/releases/current.html) or [CHESS](http://ccb.jhu.edu/chess/).
* *Reference genome*, in FASTA format. For example hg38. *Make sure your annotation GTF is based on the correct ref genome version!*

Optionally:

* CAGE Peak data (from FANTOM5). I've provided a version of [CAGE Peak for hg38 genome](https://github.com/Magdoll/images_public/blob/master/SQANTI2_support_data/hg38.cage_peak_phase1and2combined_coord.bed.gz) which was originally from [FANTOM5](http://fantom.gsc.riken.jp/5/datafiles/latest/extra/CAGE_peaks/). 

* [Intropolis](https://github.com/nellore/intropolis/blob/master/README.md) Junction BED file. I've provided a version of [Intropolis for hg38 genome, modified into STAR junction format](https://github.com/Magdoll/images_public/tree/master/SQANTI2_support_data).

* polyA motif list. A ranked  list of polyA motifs to find upstream of the 3' end site. See [human polyA list](https://raw.githubusercontent.com/Magdoll/images_public/master/SQANTI2_support_data/human.polyA.list.txt) for an example.

### Running SQANTI QC

The script usage is:

```
python sqanti_qc2.py [-t cpus] [--skipORF] [-c shortread_STAR_junction_out] 
     [--cage_peak CAGE_PEAK_BED]
     [--polyA_motif_list POLYA_LIST]
     [--aligner_choice=minimap2,deSALT]
     [--is_fusion]
     <input_fasta> <annotation_gtf> <genome_fasta>
```

If you don't feel like running the ORF prediction part, use `--skipORF`. Just know that all your transcripts will be annotated as non-coding.
If you have short read data, you can run STAR to get the junction file (usually called `SJ.out.tab`, see [STAR manual](https://github.com/alexdobin/STAR/blob/master/doc/STARmanual.pdf)) and supply it to SQANTI2.

If `--aligner_choice=minimap2`, the minimap2 parameter used currently is: `minimap2 -ax splice --secondary=no -C5 -O6,24 -B4 -uf`
If `--aligner_choice=deSALT`, the deSALT parameter used currently is: `deSALT aln -x ccs`. 

You can look at the [`MINIMAP2_CMD` and `DESALT_CMD` in `sqanti_qc2.py` for the full command format](https://github.com/Magdoll/SQANTI2/blob/master/sqanti_qc2.py#L61).


For example:

```
python sqanti_qc2.py -t 30 example/touse.rep.fasta gencode.v29.annotation.gtf hg38_noalt.fa \
      --cage_peak hg38.cage_peak_phase1and2combined_coord.bed \
      --coverage Public_Intronpolis/intropolis.v1.hg19_with_liftover_to_hg38.tsv.modified
```

If you have multiple bed files, you can use file patterns:

```
python sqanti_qc2.py -t 30 example/touse.rep.fasta gencode.v29.annotation.gtf hg38_noalt.fa \
      --cage_peak hg38.cage_peak_phase1and2combined_coord.bed \
      --coverage "JunctionBeds/samples.*.junctions.bed"
```


For fusion transcripts, you must use the `--is_fusion` option for `sqanti_qc2.py` to work properly. Furthermore, the IDs in the input FASTA/FASTQ *must* have the format `PBfusion.X`, as is output by [`fusion_finder.py` in Cupcake](https://github.com/Magdoll/cDNA_Cupcake/wiki/Cupcake-ToFU:-supporting-scripts-for-Iso-Seq-after-clustering-step#fusion).


### SQANTI QC output

You can look at the [example](https://github.com/Magdoll/SQANTI2/tree/master/example) subfolder for a sample output. The PDF file shows all the figures drawn using R script [SQANTI_report2.R](https://github.com/Magdoll/SQANTI2/blob/master/utilities/SQANTI_report2.R), taking the `_classification.txt` and `_junctions.txt` as the two input. If you know R well, you are free to modify the R script to add new figures! I will be constantly adding new figures as well.

Detailed explanation of `_classification.txt` and `_junctions.txt` <a href="#explain">below</a>.


### Filtering Isoforms using SQANTI


I've made a lightweight filtering script based on SQANTI2 output that filters for two things: (a) intra-priming and (b) short read junction support.  

The script usage is:

```
python sqanti_filter2.py <classification_txt> <input_fasta> <input_sam>
         [-a INTRAPRIMING] [-c MIN_COV] [-m MAX_DIST_TO_KNOWN_END]
```

where `-a` determines the fraction of genomic 'A's above which the isoform will be filtered. The default is `-a 0.8`. 

`-m` sets the maximum distance to an annotated 3' end (the `diff_to_gene_TTS` field in classification output) to offset the intrapriming rule.

`-c` is the filter for the minimum short read junction support (looking at the `min_cov` field in `.classification.txt`), and can only be used if you have short read data.


For example:

```
python sqanti_filter2.py touse.rep_classification.txt \
                         touse.rep.renamed_corrected.fasta \
                         touse.rep.renamed_corrected.sam
```

The current filtering rules are as follow:

* If a transcript is FSM, ISM, or NIC, then it is kept unless the 3' end is unreliable.
* If a transcript is NNC, genic, genomic, then it is kept only if all of below are true: (1) 3' end is reliable (2) does not have a junction that is labeled as RTSwitching (3) all junctions are either canonical or has short read coverage above `-c` threshold

A transcript 3' end is *unreliable* if both of the following are true: the genomic 'A' content is above the `-a` threshold *and* the distance to the known 3' end exceeeds the `-m` threshold.


<a name="explain"/>

### SQANTI2 Output Explanation


SQANTI/SQANTI2 categorizes each isoform by finding the best matching reference transcript in the following order:

* FSM (*Full Splice Match*): meaning the reference and query isoform have the same number of exons and each internal junction agree. The exact 5' start and 3' end can differ by any amount.

* ISM (*Incomplete Splice Match*): the query isoform has 5' exons than the reference, but each internal junction agree. The exact 5' start and 3' end can differ by any amount.

* NIC (*Novel In Catalog*): the query isoform does not have a FSM or ISM match, but is using a combination of known donor/acceptor sites.

* NNC (*Novel Not in Catalog*): the query isoform does not have a FSM or ISM match, and has at least one donor or acceptor site that is not annotated.

* *Antisense*: the query isoform does not have overlap a same-strand reference gene but is anti-sense to an annotated gene. 

* *Genic Intron*: the query isoform is completely contained within an annotated intron.

* *Genic Genomic*: the query isoform overlaps with introns and exons.

* *Intergenic*: the query isoform is in the intergenic region.

![sqanti_explain](https://github.com/Magdoll/images_public/blob/master/SQANTI2_figures/sqanti2_classification.png)


#### Classification Output Explanation

The output `.classification.txt` has the following fields:

1. `isoform`: the isoform ID. Usually in `PB.X.Y` format.
2. `chrom`: chromosome.
3. `strand`: strand.
4. `length`: isoform length.
5. `exons`: number of exons.
6. `structural_category`: one of the categories ["full-splice_match", "incomplete-splice_match", "novel_in_catalog", "novel_not_in_catalog", "genic", "antisense", "fusion", "intergenic", "genic_intron"]
7. `associated_gene`: the reference gene name.
8. `associated_transcript`: the reference transcript name.
9. `ref_length`: reference transcript length.
10. `ref_exons`: reference transcript number of exons.
11. `diff_to_TSS`: distance of query isoform 5' start to reference transcript start end. Negative value means query starts downstream of reference.
12. `diff_to_TTS`: distance of query isoform 3' end to reference annotated end site. Negative value means query ends upstream of reference.
13. `diff_to_gene_TSS`: distance of query isoform 5' start to the closest start end of *any* transcripts of the matching gene. This field is different from `diff_to_TSS` since it's looking at all annotated starts of a gene. Negative value means query starts downstream of reference.
14. `diff_to_gene_TTS`: distance of query isoform 3' end to the closest end of *any* transcripts of the matching gene. Negative value means query ends upstream of reference.
13. `subcategory`: additional splicing categorization, separated by semi-colons. Categories include: `mono-exon`, `multi-exon`. Intron rentention is marked with `intron_retention`. 
14. `RTS_stage`: TRUE if one of the junctions could be a RT switching artifact.
15. `all_canonical`: TRUE if all junctions have canonical splice sites.
16. `min_sample_cov`: sample with minimum coverage.
17. `min_cov`: minimum junction coverage based on short read STAR junction output file. NA if no short read given.
18. `min_cov_pos`: the junction that had the fewest coverage. NA if no short read data given.
19. `sd_cov`: standard deviation of junction coverage counts from short read data. NA if no short read data given.
20. `FL`: currently always NA. I will add back this information later.
21. `n_indels`: total number of indels based on alignment.
22. `n_indels_junc`: number of junctions in this isoform that have alignment indels near the junction site (indicating potentially unreliable junctions).
23. `bite`: TRUE if all junctions match reference junctions completely.
24. `iso_exp`: currently always NA. I will add back this information later.
25. `gene_exp`: currently always NA. I will add back this information later.
26. `ratio_exp`: currently always NA. I will add back this information later.
27. `FSM_class`: ignore this field for now.
28. `ORF_length`: predicted ORF length.
29. `CDS_length`: predicted CDS length. 
30. `CDS_start`: CDS start.
31. `CDS_end`: CDS end.
32. `CDS_genomic_start`: genomic coordinate of the CDS start. If on - strand, this coord will be greater than the end.
33. `CDS_genomic_end`: genomic coordinate of the CDS end. If on - strand, this coord will be smaller than the start.
34. `predicted_NMD`: TRUE if there's a predicted ORF and CDS ends before the last junction; FALSE if otherwise. NA if non-coding.
35. `perc_A_downstreamTTS`: percent of genomic "A"s in the downstream 20 bp window. If this number if high (say > 0.8), the 3' end site of this isoform is probably not reliable.
36. `dist_peak`: distance to closest TSS based on CAGE Peak data. Negative means upstream of TSS and positive means downstream of TSS. Strand-specific. SQANTI2 only searches for nearby CAGE Peaks within 10000 bp of the PacBio transcript start site. Will be `NA` if none are found within 10000 bp.
37. `within_peak`: TRUE if the PacBio transcript start site is within a CAGE Peak. 
38. `polyA_motif`: if `--polyA_motif_list` is given, shows the top ranking polyA motif found within 50 bp upstream of end.
39. `polyA_dist`: if `--polyA_motif_list` is given, shows the location of the  last base of the hexamer. Position 0 is the putative poly(A) site. This distance is hence always negative because it is upstream. 

### Junction Output Explanation


THe `.junctions.txt` file shows every junction for every PB isoform. NOTE because of this the *same* junction might appear multiple times if they are shared by multiple PB isoforms. 

1. `isoform`: Isoform ID
2. `junction_number`: The i-th junction of the isoform
3. `chrom`: Chromosome 
4. `strand`: Strand
5. `genomic_start_coord`: Start of the junction (1-based), note that if on - strand, this would be the junction acceptor site instead.
6. `genomic_end_coord`: End of the junction (1-based), note that if on - strand, this would be the junction donor site instead.
7. `transcript_coord`: Currently not implemented. Ignore.
8. `junction_category`: `known` if the (donor-acceptor) combination is annotated in the GTF file, `novel` otherwise. Note that it is possible to have a `novel` junction even though both the donor and acceptor site are known, since the combination might be novel.
9. `start_site_category`: `known` if the junction start site is annotated. If on - strand, this is actually the donor site.
10. `end_site_category`: `known` if the junction end site is annotated. If on - strand, this is actually the acceptor site.
11. `diff_to_Ref_start_site`: distance to closest annotated junction start site. If on - strand, this is actually the donor site.
12. `diff_to_Ref_end_site`: distance to closest annotated junction end site. If on - strand, this is actually the acceptor site.
13. `bite_junction`: TRUE if either or both the junction start/end site matches annotation.
14. `splice_site`: Splice motif.
15. `RTS_junction`: TRUE if junction is predicted to a template switching artifact.
16. `indel_near_junct`: TRUE if there is alignment indel error near the junction site, indicating potential junction incorrectness.
17. `sample_with_cov`: If `--coverage` (short read junction coverage info) is provided, shows the number of samples (cov files) that have short read that support this junction.
18. `total_coverage`: Total number of short read support from all samples that cover this junction.

