#!/usr/bin/env python
# SQANTI: Structural and Quality Annotation of Novel Transcript Isoforms
# Authors: Lorena de la Fuente, Hector del Risco, Cecile Pereira and Manuel Tardaguila
# Modified by Liz (etseng@pacb.com) currently as SQANTI2 working version

__author__  = "etseng@pacb.com"
__version__ = '2.5'

import os, re, sys, subprocess, timeit, glob
import itertools
import bisect
import argparse
import math
from collections import defaultdict, Counter
from csv import DictWriter, DictReader

utilitiesPath =  os.path.dirname(os.path.realpath(__file__))+"/utilities/" 
sys.path.insert(0, utilitiesPath)
from rt_switching import rts
from indels_annot import calc_indels_from_sam
import distutils.spawn


try:
    from Bio.Seq import Seq
    from Bio import SeqIO
    from Bio.SeqRecord import SeqRecord
except ImportError:
    print >> sys.stderr, "Unable to import Biopython! Please make sure Biopython is installed."
    sys.exit(-1)

try:
    from bx.intervals import Interval, IntervalTree
except ImportError:
    print >> sys.stderr, "Unable to import bx-python! Please make sure bx-python is installed."
    sys.exit(-1)

try:
    from BCBio import GFF as BCBio_GFF
except ImportError:
    print >> sys.stderr, "Unable to import BCBio! Please make sure bcbiogff is installed."
    sys.exit(-1)

try:
    from err_correct_w_genome import err_correct
    from sam_to_gff3 import convert_sam_to_gff3
    from STAR import STARJunctionReader
    from BED import LazyBEDPointReader
except ImportError:
    print >> sys.stderr, "Unable to import err_correct_w_genome or sam_to_gff3.py! Please make sure cDNA_Cupcake/sequence/ is in $PYTHONPATH."
    sys.exit(-1)

try:
    from cupcake.tofu.compare_junctions import compare_junctions
except ImportError:
    print >> sys.stderr, "Unable to import cupcake.tofu! Please make sure you install cupcake."
    sys.exit(-1)


GMAP_CMD = "gmap -D --cross-species -n 1 --max-intronlength-middle=2000000 --max-intronlength-ends=2000000 -L 3000000 -f samse -t {cpus} {dir} -d {name} -z {sense} {i} > {o}"
MINIMAP2_CMD = "minimap2 -ax splice --secondary=no -C5 -u{sense} -t {cpus} {g} {i} > {o}"

GMSP_PROG = os.path.join(utilitiesPath, "gmst", "gmst.pl")
GMST_CMD = "perl " + GMSP_PROG + " -faa --strand direct --fnn --output {o} {i}"



FIELDS_JUNC = ['isoform', 'chrom', 'strand', 'junction_number', 'genomic_start_coord',
                   'genomic_end_coord', 'transcript_coord', 'junction_category',
                   'start_site_category', 'end_site_category', 'diff_to_Ref_start_site',
                   'diff_to_Ref_end_site', 'bite_junction', 'splice_site', 'canonical',
                   'RTS_junction', 'indel_near_junct',
                   'phyloP_start', 'phyloP_end', 'sample_with_cov', "total_coverage"] #+coverage_header

FIELDS_CLASS = ['isoform', 'chrom', 'strand', 'length',  'exons',  'structural_category',
                'associated_gene', 'associated_transcript',  'ref_length', 'ref_exons',
                'diff_to_TSS', 'diff_to_TTS', 'subcategory', 'RTS_stage', 'all_canonical',
                'min_sample_cov', 'min_cov', 'min_cov_pos',  'sd_cov', 'FL', 'n_indels',
                'n_indels_junc',  'bite',  'iso_exp', 'gene_exp',  'ratio_exp',
                'FSM_class',   'coding', 'ORF_length', 'CDS_length', 'CDS_start',
                'CDS_end', 'perc_A_downstream_TTS', 'dist_to_cage_peak', 'within_cage_peak',
                'polyA_motif', 'polyA_dist']

RSCRIPTPATH = distutils.spawn.find_executable('Rscript')
RSCRIPT_REPORT = 'SQANTI_report2.R'

if os.system(RSCRIPTPATH + " --version")!=0:
    print >> sys.stderr, "Rscript executable not found! Abort!"
    sys.exit(-1)

class genePredReader(object):
    def __init__(self, filename):
        self.f = open(filename)

    def __iter__(self):
        return self

    def next(self):
        line = self.f.readline().strip()
        if len(line) == 0:
            raise StopIteration
        return genePredRecord.from_line(line)


class genePredRecord(object):
    def __init__(self, id, chrom, strand, txStart, txEnd, cdsStart, cdsEnd, exonCount, exonStarts, exonEnds, gene=None):
        self.id = id
        self.chrom = chrom
        self.strand = strand
        self.txStart = txStart
        self.txEnd = txEnd
        self.cdsStart = cdsStart
        self.cdsEnd = cdsEnd
        self.exonCount = exonCount
        self.exonStarts = exonStarts   # 0-based starts
        self.exonEnds = exonEnds       # 1-based ends
        self.gene = gene

        self.length = 0
        self.exons = []

        for s,e in zip(exonStarts, exonEnds):
            self.length += e-s
            self.exons.append(Interval(s, e))

        self.junctions = [(self.exonEnds[i],self.exonStarts[i+1]) for i in xrange(self.exonCount-1)]

    @property
    def segments(self):
        return self.exons


    @classmethod
    def from_line(cls, line):
        raw = line.strip().split('\t')
        return cls(id=raw[0],
                              chrom=raw[1],
                              strand=raw[2],
                              txStart=int(raw[3]),
                              txEnd=int(raw[4]),
                              cdsStart=int(raw[5]),
                              cdsEnd=int(raw[6]),
                              exonCount=int(raw[7]),
                              exonStarts=map(int, raw[8][:-1].split(',')),  #exonStarts string has extra , at end
                              exonEnds=map(int, raw[9][:-1].split(',')),     #exonEnds string has extra , at end
                              gene=raw[11] if len(raw)>=12 else None,
                              )

    def get_splice_site(self, genome_dict, i):
        """
        Return the donor-acceptor site (ex: GTAG) for the i-th junction
        :param i: 0-based junction index
        :param genome_dict: dict of chrom --> SeqRecord
        :return: splice site pattern, ex: "GTAG", "GCAG" etc
        """
        assert 0 <= i < self.exonCount-1

        d = self.exonEnds[i]
        a = self.exonStarts[i+1]

        seq_d = genome_dict[self.chrom].seq[d:d+2]
        seq_a = genome_dict[self.chrom].seq[a-2:a]

        if self.strand == '+':
            return (str(seq_d)+str(seq_a)).upper()
        else:
            return (str(seq_a.reverse_complement())+str(seq_d.reverse_complement())).upper()



class myQueryTranscripts:
    def __init__(self, id, tss_diff, tts_diff, num_exons, length, str_class, subtype=None,
                 genes=None, transcripts=None, chrom=None, strand=None, bite ="NA",
                 RT_switching ="????", canonical="NA", min_cov ="NA",
                 min_cov_pos ="NA", min_samp_cov="NA", sd ="NA", FL ="NA",
                 nIndels ="NA", nIndelsJunc ="NA", proteinID=None,
                 ORFlen="NA", CDS_start="NA", CDS_end="NA",
                 isoExp ="NA", geneExp ="NA", coding ="non_coding",
                 refLen ="NA", refExons ="NA",
                 FSM_class = None, percAdownTTS = None,
                 dist_cage='NA', within_cage='NA',
                 polyA_motif='NA', polyA_dist='NA'):

        self.id  = id
        self.tss_diff    = tss_diff
        self.tts_diff    = tts_diff
        self.genes 		 = genes if genes is not None else []
        self.AS_genes    = set()   # ref genes that are hit on the opposite strand
        self.transcripts = transcripts if transcripts is not None else []
        self.num_exons = num_exons
        self.length      = length
        self.str_class   = str_class  	# structural classification of the isoform.
        self.chrom       = chrom
        self.strand 	 = strand
        self.subtype 	 = subtype
        self.RT_switching= RT_switching
        self.canonical   = canonical
        self.min_samp_cov = min_samp_cov
        self.min_cov     = min_cov
        self.min_cov_pos = min_cov_pos
        self.sd 	     = sd
        self.proteinID   = proteinID
        self.ORFlen      = ORFlen
        self.CDS_start   = CDS_start
        self.CDS_end     = CDS_end
        self.coding      = coding
        self.FL          = FL
        self.nIndels     = nIndels
        self.nIndelsJunc = nIndelsJunc
        self.isoExp      = isoExp
        self.geneExp     = geneExp
        self.refLen      = refLen
        self.refExons    = refExons
        self.FSM_class   = FSM_class
        self.bite        = bite
        self.percAdownTTS = percAdownTTS
        self.dist_cage   = dist_cage
        self.within_cage = within_cage
        self.polyA_motif = polyA_motif
        self.polyA_dist  = polyA_dist

    def get_total_diff(self):
        return abs(self.tss_diff)+abs(self.tts_diff)

    def modify(self, ref_transcript, ref_gene, tss_diff, tts_diff, refLen, refExons):
        self.transcripts = [ref_transcript]
        self.genes = [ref_gene]
        self.tss_diff = tss_diff
        self.tts_diff = tts_diff
        self.refLen = refLen
        self.refExons = refExons

    def geneName(self):
        geneName = "_".join(set(self.genes))
        return geneName

    def ratioExp(self):
        if self.geneExp == 0 or self.geneExp == "NA":
            return "NA"
        else:
            ratio = float(self.isoExp)/float(self.geneExp)
        return(ratio)

    def CDSlen(self):
        if self.coding == "coding":
            return(str(int(self.CDS_end) - int(self.CDS_start) + 1))
        else:
            return("NA")

    def __str__(self):
        return "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s" % (self.chrom, self.strand,
                                                                                                                                                           str(self.length), str(self.num_exons),
                                                                                                                                                           str(self.str_class), "_".join(set(self.genes)),
                                                                                                                                                           self.id, str(self.refLen), str(self.refExons),
                                                                                                                                                           str(self.tss_diff), str(self.tts_diff),
                                                                                                                                                           self.subtype, self.RT_switching,
                                                                                                                                                           self.canonical, str(self.min_samp_cov),
                                                                                                                                                           str(self.min_cov), str(self.min_cov_pos),
                                                                                                                                                           str(self.sd), str(self.FL), str(self.nIndels),
                                                                                                                                                           str(self.nIndelsJunc), self.bite, str(self.isoExp),
                                                                                                                                                           str(self.geneExp), str(self.ratioExp()),
                                                                                                                                                           self.FSM_class, self.coding, str(self.ORFlen),
                                                                                                                                                           str(self.CDSlen()), str(self.CDS_start), str(self.CDS_end),
                                                                                                                                                           str(self.percAdownTTS),
                                                                                                                                                           str(self.dist_cage),
                                                                                                                                                           str(self.within_cage),
                                                                                                                                                           str(self.polyA_motif),
                                                                                                                                                           str(self.polyA_dist))


    def as_dict(self):
        return {'isoform': self.id,
         'chrom': self.chrom,
         'strand': self.strand,
         'length': self.length,
         'exons': self.num_exons,
         'structural_category': self.str_class,
         'associated_gene': "_".join(set(self.genes)),
         'associated_transcript': "_".join(set(self.transcripts)),
         'ref_length': self.refLen,
         'ref_exons': self.refExons,
         'diff_to_TSS': self.tss_diff,
         'diff_to_TTS': self.tts_diff,
         'subcategory': self.subtype,
         'RTS_stage': self.RT_switching,
         'all_canonical': self.canonical,
         'min_sample_cov': self.min_samp_cov,
         'min_cov': self.min_cov,
         'min_cov_pos': self.min_cov_pos,
         'sd_cov': self.sd,
         'FL': self.FL,
         'n_indels': self.nIndels,
         'n_indels_junc': self.nIndelsJunc,
         'bite': self.bite,
         'iso_exp': self.isoExp,
         'gene_exp': self.geneExp,
         'ratio_exp': self.ratioExp(),
         'FSM_class': self.FSM_class,
         'coding': self.coding,
         'ORF_length': self.ORFlen,
         'CDS_length': self.CDSlen(),
         'CDS_start': self.CDS_start,
         'CDS_end': self.CDS_end,
         'perc_A_downstream_TTS': self.percAdownTTS,
         'dist_to_cage_peak': self.dist_cage,
         'within_cage_peak': self.within_cage,
         'polyA_motif': self.polyA_motif,
         'polyA_dist': self.polyA_dist
         }

class myQueryProteins:

    def __init__(self, cds_start, cds_end, orf_length, proteinID="NA"):
        self.orf_length  = orf_length
        self.cds_start   = cds_start
        self.cds_end     = cds_end
        self.proteinID   = proteinID

def correctionPlusORFpred(args, genome_dict):
    """
    Use the reference genome to correct the sequences (unless a pre-corrected GTF is given)
    """
    global corrORF
    global corrGTF
    global corrSAM
    global corrFASTA

    corrPathPrefix = args.dir+"/"+os.path.splitext(os.path.basename(args.isoforms))[0]
    corrGTF = corrPathPrefix +"_corrected.gtf"
    corrSAM = corrPathPrefix +"_corrected.sam"
    corrFASTA = corrPathPrefix +"_corrected.fasta"
    corrORF =  corrPathPrefix +"_corrected.faa"


    # Step 1. IF GFF or GTF is provided, make it into a genome-based fasta
    #         IF sequence is provided, align as SAM then correct with genome
    if os.path.exists(corrFASTA):
        print >> sys.stderr, "Error corrected FASTA {0} already exists. Using it...".format(corrFASTA)
    else:
        if not args.gtf:
            if os.path.exists(corrSAM):
                print >> sys.stderr, "Aligned SAM {0} already exists. Using it...".format(corrSAM)
            else:
                if args.aligner_choice == "gmap":
                    print >> sys.stdout, "****Aligning reads with GMAP..."
                    cmd = GMAP_CMD.format(cpus=args.gmap_threads,
                                          dir=os.path.dirname(args.gmap_index),
                                          name=os.path.basename(args.gmap_index),
                                          sense=args.sense,
                                          i=args.isoforms,
                                          o=corrSAM)
                elif args.aligner_choice == "minimap2":
                    print >> sys.stdout, "****Aligning reads with Minimap2..."
                    cmd = MINIMAP2_CMD.format(cpus=args.gmap_threads,
                                              sense=args.sense,
                                              g=args.genome,
                                              i=args.isoforms,
                                              o=corrSAM)
                if subprocess.check_call(cmd, shell=True)!=0:
                    print >> sys.stderr, "ERROR running alignment cmd: {0}".format(cmd)
                    sys.exit(-1)
            # error correct the genome (input: corrSAM, output: corrFASTA)
            err_correct(args.genome, corrSAM, corrFASTA, genome_dict=genome_dict)
            # convert SAM to GFF --> GTF
            convert_sam_to_gff3(corrSAM, corrGTF+'.tmp', source=os.path.basename(args.genome).split('.')[0])  # convert SAM to GFF3
            cmd = "{u}/gffread {o}.tmp -T -o {o}".format(u=utilitiesPath, o=corrGTF)
            if subprocess.check_call(cmd, shell=True)!=0:
                print >> sys.stderr, "ERROR running cmd: {0}".format(cmd)
                sys.exit(-1)
        else:
            print >> sys.stdout, "Skipping aligning of sequences because gtf file was provided."

            ind = 0
            with open(args.isoforms, 'r') as isoforms_gtf:
                for line in isoforms_gtf:
                    if line[0] != "#" and len(line.split("\t"))!=9:
                        sys.stderr.write("\nERROR: input isoforms file with not gtf format.\n")
                        sys.exit()
                    elif len(line.split("\t"))==9:
                        ind += 1
                if ind==0:
                        sys.stderr.write("\nERROR: gtf has not annotation lines.\n")
                        sys.exit()


            # GFF to GTF (in case the user provides gff instead of gtf)
            corrGTF_tpm = corrGTF+".tmp"
            try:
                subprocess.call([utilitiesPath+"gffread", args.isoforms , '-T', '-o', corrGTF_tpm])
            except (RuntimeError, TypeError, NameError):
                sys.stderr.write('ERROR: File %s without gtf/gff format.\n' % args.isoforms)
                raise SystemExit(1)


            # check if gtf chromosomes inside genome file
            with open(corrGTF, 'w') as corrGTF_out:
                with open(corrGTF_tpm, 'r') as isoforms_gtf:
                    for line in isoforms_gtf:
                        if line[0] != "#":
                            chrom = line.split("\t")[0]
                            type = line.split("\t")[2]
                            if chrom not in genome_dict.keys():
                                sys.stderr.write("\nERROR: gtf \"%s\" chromosome not found in genome reference file.\n" % (chrom))
                                sys.exit()
                            elif type=="exon":
                                corrGTF_out.write(line)
            os.remove(corrGTF_tpm)

            if not os.path.exists(corrSAM):
                sys.stdout.write("\nIndels will be not calculated since you ran SQANTI without alignment step (SQANTI with gtf format as transcriptome input).\n")

            # GTF to FASTA
            subprocess.call([utilitiesPath+"gffread", corrGTF, '-g', args.genome, '-w', corrFASTA])

    # ORF generation
    print >> sys.stdout, "**** Predicting ORF sequences..."

    gmst_dir = os.path.join(args.dir, "GMST")
    gmst_pre = os.path.join(gmst_dir, "GMST_tmp")
    if not os.path.exists(gmst_dir):
        os.makedirs(gmst_dir)


    # sequence ID example: PB.2.1 gene_4|GeneMark.hmm|264_aa|+|888|1682
    gmst_rex = re.compile('(\S+\t\S+\|GeneMark.hmm)\|(\d+)_aa\|(\S)\|(\d+)\|(\d+)')
    orfDict = {}  # GMST seq id --> myQueryProteins object
    if args.skipORF:
        print >> sys.stderr, "WARNING: Skipping ORF prediction because user requested it. All isoforms will be non-coding!"
    elif os.path.exists(corrORF):
        print >> sys.stderr, "ORF file {0} already exists. Using it....".format(corrORF)
        for r in SeqIO.parse(open(corrORF), 'fasta'):
            # now process ORFs into myQueryProtein objects
            m = gmst_rex.match(r.description)
            if m is None:
                print >> sys.stderr, "Expected GMST output IDs to be of format '<pbid> gene_4|GeneMark.hmm|<orf>_aa|<strand>|<cds_start>|<cds_end>' but instead saw: {0}! Abort!".format(r.description)
                sys.exit(-1)
            orf_length = int(m.group(2))
            cds_start = int(m.group(4))
            cds_end = int(m.group(5))
            orfDict[r.id] = myQueryProteins(cds_start, cds_end, orf_length, proteinID=r.id)
    else:
        cmd = GMST_CMD.format(i=corrFASTA, o=gmst_pre)
        if subprocess.check_call(cmd, shell=True, cwd=gmst_dir)!=0:
            print >> sys.stderr, "ERROR running GMST cmd: {0}".format(cmd)
            sys.exit(-1)
        # Modifying ORF sequences by removing sequence before ATG
        with open(corrORF, "w") as f:
            for r in SeqIO.parse(open(gmst_pre+'.faa'), 'fasta'):
                m = gmst_rex.match(r.description)
                if m is None:
                    print >> sys.stderr, "Expected GMST output IDs to be of format '<pbid> gene_4|GeneMark.hmm|<orf>_aa|<strand>|<cds_start>|<cds_end>' but instead saw: {0}! Abort!".format(r.description)
                    sys.exit(-1)
                id_pre = m.group(1)
                orf_length = int(m.group(2))
                orf_strand = m.group(3)
                cds_start = int(m.group(4))
                cds_end = int(m.group(5))
                pos = r.seq.find('M')
                if pos!=-1:
                    # must modify both the sequence ID and the sequence
                    orf_length -= pos
                    cds_start += pos*3
                    newid = "{0}|{1}_aa|{2}|{3}|{4}".format(id_pre, orf_length, orf_strand, cds_start, cds_end)
                    newseq = r.seq.tostring()[pos:]
                    orfDict[r.id] = myQueryProteins(cds_start, cds_end, orf_length, proteinID=newid)
                    f.write(">{0}\n{1}\n".format(newid, newseq))
                else:
                    new_rec = r
                    orfDict[r.id] = myQueryProteins(cds_start, cds_end, orf_length, proteinID=r.id)
                    f.write(">{0}\n{1}\n".format(new_rec.description, new_rec.seq))

    if len(orfDict) == 0:
        print >> sys.stderr, "WARNING: All input isoforms were predicted as non-coding"


    return(orfDict)


def reference_parser(args, genome_chroms):
    """
    Read the reference GTF file
    :param args:
    :param genome_chroms: list of chromosome names from the genome fasta, used for sanity checking
    :return: (refs_1exon_by_chr, refs_exons_by_chr, junctions_by_chr, junctions_by_gene)
    """
    global referenceFiles

    referenceFiles = os.path.join(args.dir, "refAnnotation_"+args.output+".genePred")
    print >> sys.stdout, "**** Parsing Reference Transcriptome...."

    if os.path.exists(referenceFiles):
        print >> sys.stdout, "{0} already exists. Using it.".format(referenceFiles)
    else:
        ## gtf to genePred
        if not args.geneid:
            subprocess.call([utilitiesPath+"gtfToGenePred", args.annotation, referenceFiles, '-genePredExt', '-allErrors', '-ignoreGroupsWithoutExons', '-geneNameAsName2'])
        else:
            subprocess.call([utilitiesPath+"gtfToGenePred", args.annotation, referenceFiles, '-genePredExt', '-allErrors', '-ignoreGroupsWithoutExons'])

    ## parse reference annotation
    # 1. ignore all miRNAs (< 200 bp)
    # 2. separately store single exon and multi-exon references
    refs_1exon_by_chr = defaultdict(lambda: IntervalTree()) #
    refs_exons_by_chr = defaultdict(lambda: IntervalTree())
    # store donors as the exon end (1-based) and acceptor as the exon start (0-based)
    # will convert the sets to sorted list later
    junctions_by_chr = defaultdict(lambda: {'donors': set(), 'acceptors': set(), 'da_pairs': set()})
    # dict of gene name --> set of junctions (don't need to record chromosome)
    junctions_by_gene = defaultdict(lambda: set())

    for r in genePredReader(referenceFiles):
        if r.length < 200: continue # ignore miRNAs
        if r.exonCount == 1:
            refs_1exon_by_chr[r.chrom].insert(r.txStart, r.txEnd, r)
        else:
            refs_exons_by_chr[r.chrom].insert(r.txStart, r.txEnd, r)
            # only store junctions for multi-exon transcripts
            for d, a in r.junctions:
                junctions_by_chr[r.chrom]['donors'].add(d)
                junctions_by_chr[r.chrom]['acceptors'].add(a)
                junctions_by_chr[r.chrom]['da_pairs'].add((d,a))
                junctions_by_gene[r.gene].add((d,a))

    # check that all genes' chromosomes are in the genome file
    ref_chroms = set(refs_1exon_by_chr.keys()).union(refs_exons_by_chr.keys())
    diff = ref_chroms.difference(genome_chroms)
    if len(diff) > 0:
        print >> sys.stderr, "WARNING: ref annotation contains chromosomes not in genome: {0}\n".format(",".join(diff))

    # convert the content of junctions_by_chr to sorted list
    for k in junctions_by_chr:
        junctions_by_chr[k]['donors'] = list(junctions_by_chr[k]['donors'])
        junctions_by_chr[k]['donors'].sort()
        junctions_by_chr[k]['acceptors'] = list(junctions_by_chr[k]['acceptors'])
        junctions_by_chr[k]['acceptors'].sort()
        junctions_by_chr[k]['da_pairs'] = list(junctions_by_chr[k]['da_pairs'])
        junctions_by_chr[k]['da_pairs'].sort()

    return dict(refs_1exon_by_chr), dict(refs_exons_by_chr), dict(junctions_by_chr), dict(junctions_by_gene)


def isoforms_parser(args):
    """
    Parse input isoforms (GTF) to dict (chr --> sorted list)
    """
    global queryFile
    queryFile = os.path.splitext(corrGTF)[0] +".genePred"

    print >> sys.stderr, "**** Parsing Isoforms...."

    # gtf to genePred
    cmd = utilitiesPath+"gtfToGenePred {0} {1} -genePredExt -allErrors -ignoreGroupsWithoutExons".format(\
        corrGTF, queryFile)
    if subprocess.check_call(cmd, shell=True)!=0:
        print >> sys.stderr, "ERROR running cmd: {0}".format(cmd)
        sys.exit(-1)


    isoforms_list = defaultdict(lambda: []) # chr --> list to be sorted later

    for r in genePredReader(queryFile):
        isoforms_list[r.chrom].append(r)

    for k in isoforms_list:
        isoforms_list[k].sort(key=lambda r: r.txStart)

    return isoforms_list


def STARcov_parser(coverageFiles): # just valid with unstrand-specific RNA-seq protocols.
    """
    :param coverageFiles: comma-separated list of STAR junction output files or a directory containing junction files
    :return: list of samples, dict of (chrom,strand) --> (0-based start, 1-based end) --> {dict of sample -> unique reads supporting this junction}
    """
    cov_files = glob.glob(coverageFiles)

    print >> sys.stderr, "Input pattern: {0}. The following files found and to be read as junctions:\n{1}".format(\
        coverageFiles, "\n".join(cov_files) )

    cov_by_chrom_strand = defaultdict(lambda: defaultdict(lambda: defaultdict(lambda: 0)))
    undefined_strand_count = 0
    all_read = 0
    samples = []
    for file in cov_files:
        prefix = os.path.basename(file[:file.rfind('.')]) # use this as sample name
        samples.append(prefix)
        for r in STARJunctionReader(file):
            if r.strand == 'NA':
                # undefined strand, so we put them in BOTH strands otherwise we'll lose all non-canonical junctions from STAR
                cov_by_chrom_strand[(r.chrom, '+')][(r.start, r.end)][prefix] = r.unique_count + r.multi_count
                cov_by_chrom_strand[(r.chrom, '-')][(r.start, r.end)][prefix] = r.unique_count + r.multi_count
                undefined_strand_count += 1
            else:
                cov_by_chrom_strand[(r.chrom, r.strand)][(r.start, r.end)][prefix] = r.unique_count + r.multi_count
            all_read += 1
    print >> sys.stderr, "{0} junctions read. {1} junctions added to both strands because no strand information from STAR.".format(all_read, undefined_strand_count)

    return samples, cov_by_chrom_strand


def expression_parser(expressionFile):

    try:
        p = open(expressionFile, "r")
    except IOError:
        sys.stderr.write('ERROR: Unable to read %s expression file\n' % expressionFile)
        raise SystemExit(1)
    try:
        header = p.readline()
        exp_dicc = {}

        for line in p:
            pbid = line.split()[0]
            if len(pbid.split("|"))>2 and pbid[0:2]=="PB": # PacBio fasta header (including chained format)
                pbid_mod = pbid.split("|")[0].split(" ")[0]
            elif len(pbid.split("|"))>4: # Refseq fasta header
                pbid_mod = pbid.split("|")[3]
            else:
                pbid_mod = pbid.split()[0] # Ensembl fasta header
            mean = sum([float(i) for i in line.rstrip().split()[1:]])/len(line.rstrip().split()[1:])
            exp_dicc[pbid_mod] = mean
        p.close()

    except IOError:
        sys.stderr.write('File %s without expression matrix format' % expressionFile)
        raise SystemExit(1)

    return(exp_dicc)



def transcriptsKnownSpliceSites(refs_1exon_by_chr, refs_exons_by_chr, trec, genome_dict, nPolyA):
    """
    :param refs_1exon_by_chr: dict of single exon references (chr -> IntervalTree)
    :param refs_exons_by_chr: dict of multi exon references (chr -> IntervalTree)
    :param trec: id record (genePredRecord) to be compared against reference
    :param genome_dict: dict of genome (chrom --> SeqRecord)
    :param nPolyA: window size to look for polyA
    :return: myQueryTranscripts object that indicates the best reference hit
    """
    def get_diff_tss_tts(trec, ref):
        if trec.strand == '+':
            diff_tss = trec.txStart - ref.txStart
            diff_tts = ref.txEnd - trec.txEnd
        else:
            diff_tts = trec.txStart - ref.txStart
            diff_tss = ref.txEnd - trec.txEnd
        return diff_tss, diff_tts

    def categorize_incomplete_matches(trec, ref):
        """
        intron_retention --- at least one trec exon covers at least two adjacent ref exons
        complete --- all junctions agree and is not IR
        5prime_fragment --- all junctions agree but trec has less 5' exons
        3prime_fragment --- all junctions agree but trec has less 3' exons
        internal_fragment --- all junctions agree but trec has less 5' and 3' exons
        """
        # check intron retention
        ref_exon_tree = IntervalTree()
        for i,e in enumerate(ref.exons): ref_exon_tree.insert(e.start, e.end, i)
        for e in trec.exons:
            if len(ref_exon_tree.find(e.start, e.end)) > 1: # multiple ref exons covered
                return "intron_retention"

        agree_front = trec.junctions[0]==ref.junctions[0]
        agree_end   = trec.junctions[-1]==ref.junctions[-1]
        if agree_front:
            if agree_end:
                return "complete"
            else: # front agrees, end does not
                return ("3prime_fragment" if trec.strand=='+' else '5prime_fragment')
        else:
            if agree_end: # front does not agree, end agrees
                return ("5prime_fragment" if trec.strand=='+' else '3prime_fragment')
            else:
                return "internal_fragment"

    # Transcript information for a single query id and comparison with reference.

    # Intra-priming: calculate percentage of "A"s right after the end
    if trec.strand == "+":
        pos_TTS = trec.exonEnds[-1]
        seq_downTTS = str(genome_dict[trec.chrom].seq[pos_TTS:pos_TTS+nPolyA]).upper()
    else: # id on - strand
        pos_TTS = trec.exonStarts[0]
        seq_downTTS = str(genome_dict[trec.chrom].seq[pos_TTS-nPolyA:pos_TTS].reverse_complement()).upper()

    percA = float(seq_downTTS.count('A'))/nPolyA*100


    isoform_hit = myQueryTranscripts(id=trec.id, tts_diff="NA", tss_diff="NA",\
                                    num_exons=trec.exonCount,
                                    length=trec.length,
                                    str_class="", \
                                    chrom=trec.chrom,
                                    strand=trec.strand, \
                                    subtype="no_subcategory",\
                                    percAdownTTS=str(percA))



    ##***************************************##
    ########### SPLICED TRANSCRIPTS ###########
    ##***************************************##

    if trec.exonCount >= 2:
        if trec.chrom not in refs_exons_by_chr: return isoform_hit  # return blank hit result
        for ref in refs_exons_by_chr[trec.chrom].find(trec.txStart, trec.txEnd):
            if trec.strand != ref.strand:
                # opposite strand, just record it in AS_genes
                isoform_hit.AS_genes.add(ref.gene)
                continue
            match_type = compare_junctions(trec, ref, internal_fuzzy_max_dist=0, max_5_diff=999999, max_3_diff=999999)

            if match_type not in ('exact', 'subset', 'partial', 'concordant', 'super', 'nomatch'):
                raise Exception, "Unknown match category {0}!".format(match_type)

            diff_tss, diff_tts = get_diff_tss_tts(trec, ref)

            # #############################
            # SQANTI's full-splice_match
            # #############################
            if match_type == "exact":
                subtype = "multi-exon"
                if isoform_hit.str_class != 'full-splice_match': # prev hit is not as good as this one, replace it!
                    isoform_hit = myQueryTranscripts(trec.id, diff_tss, diff_tts, trec.exonCount, trec.length,
                                                     str_class="full-splice_match",
                                                     subtype=subtype,
                                                     chrom=trec.chrom,
                                                     strand=trec.strand,
                                                     genes=[ref.gene],
                                                     transcripts=[ref.id],
                                                     refLen = ref.length,
                                                     refExons= ref.exonCount,
                                                     percAdownTTS=str(percA))
                elif abs(diff_tss)+abs(diff_tts) < isoform_hit.get_total_diff(): # prev hit is FSM however this one is better
                    isoform_hit.modify(ref.id, ref.gene, diff_tss, diff_tts, ref.length, ref.exonCount)

            # #######################################################
            # SQANTI's incomplete-splice_match
            # (only check if don't already have a FSM match)
            # #######################################################
            elif isoform_hit.str_class!='full-splice_match' and match_type == "subset":
                subtype = categorize_incomplete_matches(trec, ref)

                if isoform_hit.str_class != 'incomplete-splice_match': # prev hit is not as good as this one, replace it!
                    isoform_hit = myQueryTranscripts(trec.id, diff_tss, diff_tts, trec.exonCount, trec.length,
                                                     str_class="incomplete-splice_match",
                                                     subtype=subtype,
                                                     chrom=trec.chrom,
                                                     strand=trec.strand,
                                                     genes=[ref.gene],
                                                     transcripts=[ref.id],
                                                     refLen = ref.length,
                                                     refExons= ref.exonCount,
                                                     percAdownTTS=str(percA))
                elif abs(diff_tss)+abs(diff_tts) < isoform_hit.get_total_diff():
                    isoform_hit.modify(ref.id, ref.gene, diff_tss, diff_tts, ref.length, ref.exonCount)
                    isoform_hit.subtype = subtype
            # #######################################################
            # Some kind of junctio nmatch that isn't ISM/FSM
            # #######################################################
            elif match_type in ('partial', 'concordant', 'super', 'nomatch') and isoform_hit.str_class not in ('full-splice_match', 'incomplete-splice_match'):
                if isoform_hit.str_class=="":
                    isoform_hit = myQueryTranscripts(trec.id, "NA", "NA", trec.exonCount, trec.length,
                                                     str_class="anyKnownSpliceSite",
                                                     subtype="no_subcategory",
                                                     chrom=trec.chrom,
                                                     strand=trec.strand,
                                                     genes=[ref.gene],
                                                     transcripts=["novel"],
                                                     refLen=ref.length,
                                                     refExons=ref.exonCount,
                                                     percAdownTTS=str(percA))
    ##***************************************####
    ########### UNSPLICED TRANSCRIPTS ###########
    ##***************************************####
    else: # single exon id
        if trec.chrom in refs_1exon_by_chr:
            for ref in refs_1exon_by_chr[trec.chrom].find(trec.txStart, trec.txEnd):
                if ref.strand != trec.strand:
                    # opposite strand, just record it in AS_genes
                    isoform_hit.AS_genes.add(ref.gene)
                    continue
                diff_tss, diff_tts = get_diff_tss_tts(trec, ref)

                # see if there's already an existing match AND if so, if this one is better
                if isoform_hit.str_class == "": # no match so far
                    isoform_hit = myQueryTranscripts(trec.id, diff_tss, diff_tts, 1, trec.length, "full-splice_match",
                                                            subtype="mono-exon",
                                                            chrom=trec.chrom,
                                                            strand=trec.strand,
                                                            genes=[ref.gene],
                                                            transcripts=[ref.id],
                                                            refLen=ref.length,
                                                            refExons = ref.exonCount,
                                                            percAdownTTS=str(percA))
                elif abs(diff_tss)+abs(diff_tts) < isoform_hit.get_total_diff():
                    isoform_hit.modify(ref.id, ref.gene, diff_tss, diff_tts, ref.length, ref.exonCount)

        if isoform_hit.str_class == "" and trec.chrom in refs_exons_by_chr:
            # no hits to single exon genes, let's see if it hits multi-exon genes
            # (1) if it overlaps with a ref exon and is contained in an exon, we call it ISM
            # (2) else, if it is completely within a ref gene start-end region, we call it NIC by intron retention
            for ref in refs_exons_by_chr[trec.chrom].find(trec.txStart, trec.txEnd):
                if ref.strand != trec.strand:
                    # opposite strand, just record it in AS_genes
                    isoform_hit.AS_genes.add(ref.gene)
                    continue
                diff_tss, diff_tts = get_diff_tss_tts(trec, ref)

                # ToDo: debug why these two are worse categorized in the new sqanti fix
                if trec.id in ('PB.2009.1', 'PB.1276.3'):
                    print "DEBUG", trec.id
                    print "query:", trec.txStart, trec.txEnd
                    print "ref:", ref.id, "exons:", ref.exons
                    for e in ref.exons:
                        print e.start <= trec.txStart < trec.txEnd <= e.end
                    #raw_input()

                for e in ref.exons:
                    if e.start <= trec.txStart < trec.txEnd <= e.end:
                        isoform_hit.str_class = "incomplete-splice_match"
                        isoform_hit.subtype = "mono-exon"
                        isoform_hit.modify(ref.id, ref.gene, diff_tss, diff_tts, ref.length, ref.exonCount)
                        # this is as good a match as it gets, we can stop the search here
                        return isoform_hit

                # if we haven't exited here, then ISM hit is not found yet
                # instead check if it's NIC by intron retention
                # but we don't exit here since the next gene could be a ISM hit
                if ref.txStart <= trec.txStart < trec.txEnd <= ref.txEnd:
                    isoform_hit.str_class = "novel_in_catalog"
                    isoform_hit.subtype = "mono-exon_by_intron_retention"
                    isoform_hit.modify("novel", ref.gene, diff_tss, diff_tts, ref.length, ref.exonCount)
                    return isoform_hit

                # if we get to here, means neither ISM nor NIC, so just add a ref gene and categorize further later
                isoform_hit.genes.append(ref.gene)

    return isoform_hit


def novelIsoformsKnownGenes(isoforms_hit, trec, junctions_by_chr, junctions_by_gene):
    """
    At this point: definitely not FSM or ISM, see if it is NIC, NNC, or fusion
    :return isoforms_hit: updated isoforms hit (myQueryTranscripts object)
    """
    ref_genes = list(set(isoforms_hit.genes))

    #
    # at this point, we have already found matching genes/transcripts
    # hence we do not need to update refLen or refExon
    # or tss_diff and tts_diff (always set to "NA" for non-FSM/ISM matches)
    #
    isoforms_hit.transcripts = ["novel"]
    if len(ref_genes) == 1:
        # hits exactly one gene, must be either NIC or NNC
        ref_gene_junctions = junctions_by_gene[ref_genes[0]]
        # 1. check if all donors/acceptor sites are known (regardless of which ref gene it came from)
        # 2. check if this query isoform uses a subset of the junctions from the single ref hit
        all_junctions_known = True
        all_junctions_in_hit_ref = True
        for d,a in trec.junctions:
            all_junctions_known = all_junctions_known and (d in junctions_by_chr[trec.chrom]['donors']) and (a in junctions_by_chr[trec.chrom]['acceptors'])
            all_junctions_in_hit_ref = all_junctions_in_hit_ref and ((d,a) in ref_gene_junctions)
        if all_junctions_known:
            isoforms_hit.str_class="novel_in_catalog"
            if all_junctions_in_hit_ref:
                isoforms_hit.subtype = "combination_of_known_junctions"
            else:
                isoforms_hit.subtype = "no_combination_of_known_junctions"
        else:
            isoforms_hit.str_class="novel_not_in_catalog"
            isoforms_hit.subtype = "any annotated donor/acceptor"
    else: # see if it is fusion
        # list of a ref junctions from all genes, including potential shared junctions
        all_ref_junctions = list(itertools.chain(junctions_by_gene[ref_gene] for ref_gene in ref_genes))

        # (junction index) --> number of refs that have this junction
        junction_ref_hit = dict((i, all_ref_junctions.count(junc)) for i,junc in enumerate(trec.junctions))

        # if the same query junction appears in more than one of the hit references, it is not a fusion
        if max(junction_ref_hit.itervalues()) > 1:
            isoforms_hit.str_class = "moreJunctions"
        else:
            isoforms_hit.str_class = "fusion"
            isoforms_hit.subtype = "mono-exon" if trec.exonCount==1 else "multi-exon"

    return isoforms_hit

def associationOverlapping(isoforms_hit, trec, junctions_by_chr):

    # at this point: definitely not FSM or ISM or NIC
    # possibly (in order of preference assignment):
    #  - NNC  (multi-exon and overlaps some ref on same strand, dun care if junctions are known)
    #  - antisense  (on opp strand of a known gene)
    #  - genic (overlaps a combination of exons and introns), ignore strand
    #  - genic_intron  (completely within an intron), ignore strand
    #  - intergenic (does not overlap a gene at all), ignore strand

    isoforms_hit.str_class = "intergenic"
    isoforms_hit.transcripts = ["novel"]
    isoforms_hit.subtype = "mono-exon" if trec.exonCount==1 else "multi-exon"

    if len(isoforms_hit.genes) == 0:
        # completely no overlap with any genes on the same strand
        # check if it is anti-sense to a known gene, otherwise it's genic_intron or intergenic
        if len(isoforms_hit.AS_genes) == 0 and trec.chrom in junctions_by_chr:
            # no hit even on opp strand
            # see if it is completely contained within a junction
            da_pairs = junctions_by_chr[trec.chrom]['da_pairs']
            i = bisect.bisect_left((trec.txStart, trec.txEnd), da_pairs)
            while i < len(da_pairs) and da_pairs[i][0] <= trec.txStart:
                if da_pairs[i][0] <= trec.txStart <= trec.txStart <= da_pairs[i][1]:
                    isoforms_hit.str_class = "genic_intron"
                    break
                i += 1
        else:
            # hits one or more genes on the opposite strand
            isoforms_hit.str_class = "antisense"
            isoforms_hit.genes = ["novelGene_{g}_AS".format(g=g) for g in isoforms_hit.AS_genes]
    else:
        # overlaps with one or more genes on the same strand
        if trec.exonCount >= 2:
            # multi-exon and has a same strand gene hit, must be NNC
            isoforms_hit.str_class = "novel_not_in_catalog"
            isoforms_hit.subtype = "not any annotated donor/acceptor"
        else:
            # single exon, must be genic
            isoforms_hit.str_class = "genic"

    return isoforms_hit


def write_junctionInfo(trec, junctions_by_chr, accepted_canonical_sites, indelInfo, genome_dict, fout, covInf=None, covNames=None, phyloP_reader=None):
    """
    :param trec: query isoform genePredRecord
    :param junctions_by_chr: dict of chr -> {'donors': <sorted list of donors>, 'acceptors': <sorted list of acceptors>, 'da_pairs': <sorted list of junctions>}
    :param accepted_canonical_sites: list of accepted canonical splice sites
    :param indelInfo: indels near junction information, dict of pbid --> list of junctions near indel (in Interval format)
    :param genome_dict: genome fasta dict
    :param fout: DictWriter handle
    :param covInf: (optional) junction coverage information, dict of (chrom,strand) -> (0-based start,1-based end) -> dict of {sample -> unique read count}
    :param covNames: (optional) list of sample names for the junction coverage information
    :param phyloP_reader: (optional) dict of (chrom,0-based coord) --> phyloP score

    Write a record for each junction in query isoform
    """
    def find_closest_in_list(lst, pos):
        i = bisect.bisect_left(lst, pos)
        if i == 0:
            return lst[0]-pos
        elif i == len(lst):
            return lst[-1]-pos
        else:
            a, b = lst[i-1]-pos, lst[i]-pos
            if abs(a) < abs(b): return a
            else: return b

    if trec.chrom not in junctions_by_chr:
        # nothing to do
        return

    # go through each trec junction
    for junction_index, (d, a) in enumerate(trec.junctions):
        # NOTE: donor just means the start, not adjusted for strand
        # find the closest junction start site
        min_diff_s = -find_closest_in_list(junctions_by_chr[trec.chrom]['donors'], d)
        # find the closest junction end site
        min_diff_e = find_closest_in_list(junctions_by_chr[trec.chrom]['acceptors'], a)

        splice_site = trec.get_splice_site(genome_dict, junction_index)

        indel_near_junction = "NA"
        if indelInfo is not None:
            indel_near_junction = "TRUE" if (trec.id in indelInfo and Interval(d,a) in indelInfo[trec.id]) else "FALSE"

        sample_cov = defaultdict(lambda: 0)  # sample -> unique count for this junction
        if covInf is not None:
            sample_cov = covInf[(trec.chrom, trec.strand)][(d,a)]

        # if phyloP score dict exists, give the triplet score of (last base in donor exon), donor site -- similarly for acceptor
        phyloP_start, phyloP_end = 'NA', 'NA'
        if phyloP_reader is not None:
            phyloP_start = ",".join(map(str, [phyloP_reader.get_pos(trec.chrom, d-1), phyloP_reader.get_pos(trec.chrom, d), phyloP_reader.get_pos(trec.chrom, d+1)]))
            phyloP_end = ",".join(map(str, [phyloP_reader.get_pos(trec.chrom, a-1), phyloP_reader.get_pos(trec.chrom, a),
                                              phyloP_reader.get_pos(trec.chrom, a+1)]))

        qj = {'isoform': trec.id,
              'junction_number': "junction_"+str(junction_index+1),
              "chrom": trec.chrom,
              "strand": trec.strand,
              "genomic_start_coord": d+1,  # write out as 1-based start
              "genomic_end_coord": a,      # already is 1-based end
              "transcript_coord": "?????",  # this is where the exon ends w.r.t to id sequence, ToDo: implement later
              "junction_category": "known" if ((d,a) in junctions_by_chr[trec.chrom]['da_pairs']) else "novel",
              "start_site_category": "known" if min_diff_s==0 else "novel",
              "end_site_category": "known" if min_diff_e==0 else "novel",
              "diff_to_Ref_start_site": min_diff_s,
              "diff_to_Ref_end_site": min_diff_e,
              "bite_junction": "TRUE" if (min_diff_s==0 or min_diff_e==0) else "FALSE",
              "splice_site": splice_site,
              "canonical": "canonical" if splice_site in accepted_canonical_sites else "non_canonical",
              "RTS_junction": "????", # First write ???? in _tmp, later is TRUE/FALSE
              "indel_near_junct": indel_near_junction,
              "phyloP_start": phyloP_start,
              "phyloP_end": phyloP_end,
              "sample_with_cov": sum(cov!=0 for cov in sample_cov.itervalues()) if covInf is not None else "NA",
              "total_coverage": sum(sample_cov.itervalues()) if covInf is not None else "NA"}

        if covInf is not None:
            for sample in covNames:
                qj[sample] = sample_cov[sample]

        fout.writerow(qj)


def isoformClassification(args, isoforms_by_chr, refs_1exon_by_chr, refs_exons_by_chr, junctions_by_chr, junctions_by_gene, genome_dict, indelsJunc):

    ## read coverage files if provided

    if args.coverage is not None:
        print >> sys.stdout, "**** Reading Splice Junctions coverage files."
        SJcovNames, SJcovInfo = STARcov_parser(args.coverage)
        fields_junc_cur = FIELDS_JUNC + SJcovNames # add the samples to the header
    else:
        SJcovNames, SJcovInfo = None, None
        print >> sys.stdout, "Splice Junction Coverage files not provided."
        fields_junc_cur = FIELDS_JUNC

    if args.cage_peak is not None:
        print >> sys.stdout, "**** Reading CAGE Peak data."
        cage_peak_obj = CAGEPeak(args.cage_peak)
    else:
        cage_peak_obj = None


    if args.polyA_motif_list is not None:
        print >> sys.stdout, "**** Reading PolyA motif list."
        polyA_motif_list = []
        for line in open(args.polyA_motif_list):
            x = line.strip().upper().replace('U', 'A')
            if any(s not in ('A','T','C','G') for s in x):
                print >> sys.stderr, "PolyA motif must be A/T/C/G only! Saw: {0}. Abort!".format(x)
                sys.exit(-1)
            polyA_motif_list.append(x)
    else:
        polyA_motif_list = None


    if args.phyloP_bed is not None:
        print >> sys.stdout, "**** Reading PhyloP BED file."
        phyloP_reader = LazyBEDPointReader(args.phyloP_bed)
    else:
        phyloP_reader = None

    # running classification
    print >> sys.stdout, "**** Performing Classification of Isoforms...."


    accepted_canonical_sites = list(args.sites.split(","))

    outputPathPrefix = args.dir+"/"+args.output

    outputClassPath = outputPathPrefix+"_classification.txt"
    fout_class = DictWriter(open(outputClassPath+"_tmp", "w"), fieldnames=FIELDS_CLASS, delimiter='\t')
    fout_class.writeheader()

    outputJuncPath = outputPathPrefix+"_junctions.txt"
    fout_junc = DictWriter(open(outputJuncPath+"_tmp", "w"), fieldnames=fields_junc_cur, delimiter='\t')
    fout_junc.writeheader()

    isoforms_info = {}
    novel_gene_index = 1

    for chrom,records in isoforms_by_chr.iteritems():
        for rec in records:
            # Find best reference hit
            isoform_hit = transcriptsKnownSpliceSites(refs_1exon_by_chr, refs_exons_by_chr, rec, genome_dict, nPolyA=args.window)


            if isoform_hit.str_class == "anyKnownSpliceSite":
                # not FSM or ISM --> see if it is NIC, NNC, or fusion
                isoform_hit = novelIsoformsKnownGenes(isoform_hit, rec, junctions_by_chr, junctions_by_gene)
            elif isoform_hit.str_class == "":
                # possibly NNC, genic, genic intron, anti-sense, or intergenic
                isoform_hit = associationOverlapping(isoform_hit, rec, junctions_by_chr)

            # write out junction information
            write_junctionInfo(rec, junctions_by_chr, accepted_canonical_sites, indelsJunc, genome_dict, fout_junc, covInf=SJcovInfo, covNames=SJcovNames, phyloP_reader=phyloP_reader)

            if isoform_hit.str_class in ("intergenic", "genic_intron"):
                # Liz: I don't find it necessary to cluster these novel genes. They should already be always non-overlapping.
                isoform_hit.genes = ['novelGene_' + str(novel_gene_index)]
                isoform_hit.transcripts = ['novel']
                novel_gene_index += 1

            isoforms_info[rec.id] = isoform_hit

            # look at Cage Peak info (if available)
            if cage_peak_obj is not None:
                if rec.strand == '+':
                    within_cage, dist_cage = cage_peak_obj.find(rec.chrom, rec.strand, rec.txStart)
                else:
                    within_cage, dist_cage = cage_peak_obj.find(rec.chrom, rec.strand, rec.txEnd)
                isoform_hit.within_cage = within_cage
                isoform_hit.dist_cage = dist_cage

            # polyA motif finding: look within 50 bp upstream of 3' end for the highest ranking polyA motif signal (user provided)
            if polyA_motif_list is not None:
                if rec.strand == '+':
                    polyA_motif, polyA_dist = find_polyA_motif(str(genome_dict[rec.chrom][rec.txEnd-50:rec.txEnd].seq), polyA_motif_list)
                else:
                    polyA_motif, polyA_dist = find_polyA_motif(str(genome_dict[rec.chrom][rec.txStart:rec.txStart+50].reverse_complement().seq), polyA_motif_list)
                isoform_hit.polyA_motif = polyA_motif
                isoform_hit.polyA_dist = polyA_dist

            fout_class.writerow(isoform_hit.as_dict())

    return isoforms_info

def FLcount_parser (fl_files):

    # may be one file, files separated by comma or a directory where file are.
    if os.path.isdir(fl_files)==True:
        cov_paths = [os.path.join(fl_files,fn) for fn in next(os.walk(fl_files))[2]]
    else:
        cov_paths = fl_files.split(",")

    fl_dicc = {}
    for path in cov_paths:
        if not os.path.exists(path):
            sys.stdout.write("\nERROR: Unable to read %s FL count files.\n" % path)
        else:
            p = open(path.strip(), "r")
            for line in p:
                if line[0] != "#" and line[0:4] != "pbid":
                    pbid = line.split("\t")[0]
                    fl = int(line.split("\t")[1])
                    if pbid in fl_dicc:
                        fl_dicc[pbid] = fl_dicc[pbid] + fl
                    else:
                        fl_dicc[pbid] = fl
            p.close()

    return(fl_dicc)


def pstdev(data):
    """Calculates the population standard deviation."""
    n = len(data)
    mean = sum(data)*1. / n  # mean
    var = sum(pow(x - mean, 2) for x in data) / n  # variance
    return math.sqrt(var)  # standard deviation


def find_polyA_motif(genome_seq, polyA_motif_list):
    """
    :param genome_seq: genomic sequence to search polyA motifs from, must already be oriented
    :param polyA_motif_list: ranked list of motifs to find, report the top one found
    :return: polyA_motif, polyA_dist (how many bases upstream is this found)
    """
    for motif in polyA_motif_list:
        i = genome_seq.find(motif)
        if i >= 0:
            return motif, -(len(genome_seq)-i-len(motif)+1)
    return 'NA', 'NA'

def run(args):

    start3 = timeit.default_timer()

    print >> sys.stdout, "**** Parsing provided files...."

    print >> sys.stdout, "Reading genome fasta {0}....".format(args.genome)
    # NOTE: can't use LazyFastaReader because inefficient. Bring the whole genome in!
    genome_dict = dict((r.name, r) for r in SeqIO.parse(open(args.genome), 'fasta'))

    ## correction of sequences and ORF prediction (if gtf provided instead of fasta file, correction of sequences will be skipped)
    orfDict = correctionPlusORFpred(args, genome_dict)

    ## parse reference id (GTF) to dicts
    refs_1exon_by_chr, refs_exons_by_chr, junctions_by_chr, junctions_by_gene = reference_parser(args, genome_dict.keys())

    ## parse query isoforms
    isoforms_by_chr = isoforms_parser(args)

    ## Run indel computation if sam exists
    # indelsJunc: dict of pbid --> list of junctions near indel (in Interval format)
    # indelsTotal: dict of pbid --> total indels count
    if os.path.exists(corrSAM):
        (indelsJunc, indelsTotal) = calc_indels_from_sam(corrSAM)
    else:
        indelsJunc = None
        indelsTotal = None

    # isoform classification + intra-priming + id and junction characterization
    isoforms_info = isoformClassification(args, isoforms_by_chr, refs_1exon_by_chr, refs_exons_by_chr, junctions_by_chr, junctions_by_gene, genome_dict, indelsJunc)

    print >> sys.stdout, "Number of classified isoforms: {0}".format(len(isoforms_info))

    outputPathPrefix = os.path.join(args.dir, args.output)
    outputClassPath = outputPathPrefix + "_classification.txt"
    outputJuncPath = outputPathPrefix + "_junctions.txt"

    ## RT-switching computation
    print >> sys.stderr, "**** RT-switching computation...."

    # RTS_info: dict of (pbid) -> list of RT junction. if RTS_info[pbid] == [], means all junctions are non-RT.
    RTS_info = rts([outputJuncPath+"_tmp", args.genome, "-a"], genome_dict)
    for pbid in isoforms_info:
        if pbid in RTS_info and len(RTS_info[pbid]) > 0:
            isoforms_info[pbid].RT_switching = "TRUE"
        else:
            isoforms_info[pbid].RT_switching = "FALSE"

    # ORF information
    for pbid in orfDict:
        if pbid in isoforms_info:
            isoforms_info[pbid].coding = "coding"
            isoforms_info[pbid].ORFlen = orfDict[pbid].orf_length
            isoforms_info[pbid].CDS_start = orfDict[pbid].cds_start
            isoforms_info[pbid].CDS_end = orfDict[pbid].cds_end

    ## FSM classification
    geneFSM_dicc = defaultdict(lambda: [])
    for iso in isoforms_info:
        gene = isoforms_info[iso].geneName()  # if multi-gene, returns "geneA_geneB_geneC..."
        geneFSM_dicc[gene].append(isoforms_info[iso].str_class)

    ## FL count files
    # ToDo: to look over this section later. Ignore now.
    if args.fl_count:
        sys.stdout.write("\n**** Reading Full-length read abundance files.\n")
        fl = FLcount_parser(args.fl_count)
        for i in fl:
            if i not in isoforms_info:
                sys.stderr.write("\nWARNING: %s ID found in FL abundance files but not in transcriptome.\n" %i)
        if len(fl)>0:
            for x in isoforms_info:
                if x in fl:
                    isoforms_info[x].FL = fl[x]
                else:
                    isoforms_info[x].FL = 0
    else:
        sys.stdout.write("\nFull-length read abundance files not provided.\n")


    ## Isoform expression information
    # ToDo: to look over this section later. Ignore now.
    if args.expression:
        sys.stdout.write("\n**** Reading Isoform Expression Information.\n")
        exp_dicc = expression_parser(args.expression)
        geneExp_dicc = {}
        for iso in isoforms_info:
            if iso not in exp_dicc:
                exp_dicc[iso] = 0
                sys.stdout.write("\nIsoform %s not found in expression matrix. Nule expression associated.\n" %(iso))
            gene = isoforms_info[iso].geneName()
            if gene not in geneExp_dicc:
                geneExp_dicc[gene] = exp_dicc[iso]
            else:
                geneExp_dicc[gene] = geneExp_dicc[gene]+exp_dicc[iso]
    else:
        exp_dicc = None
        geneExp_dicc = None
        sys.stdout.write("\nIsoforms expression files not provided.\n")


    # ToDo: to look over this section later. Ignore now.
    ## Adding indel, FSM class and expression information
    for iso in isoforms_info:
        gene = isoforms_info[iso].geneName()
        if exp_dicc!=None and geneExp_dicc!=None:
            isoforms_info[iso].geneExp = geneExp_dicc[gene]
            isoforms_info[iso].isoExp = exp_dicc[iso]
        if len(geneFSM_dicc[gene])==1:
            isoforms_info[iso].FSM_class = "A"
        elif "full-splice_match" in geneFSM_dicc[gene]:
            isoforms_info[iso].FSM_class = "C"
        else:
            isoforms_info[iso].FSM_class = "B"

    if indelsTotal is not None:
        for iso in isoforms_info:
            if iso in indelsTotal:
                isoforms_info[iso].nIndels = indelsTotal[iso]
            else:
                isoforms_info[iso].nIndels = 0


    ## Read junction files and create attributes per id
    # Read the junction information to fill in several remaining unfilled fields in classification
    # (1) "canonical": is "canonical" if all junctions are canonical, otherwise "non_canonical"
    # (2) "bite": is TRUE if any of the junction "bite_junction" field is TRUE
    # ToDO: assign other fields later, right now I only did the .canonical and .bite field

    reader = DictReader(open(outputJuncPath+"_tmp"), delimiter='\t')
    fields_junc_cur = reader.fieldnames
    sj_covs_by_isoform = defaultdict(lambda: [])  # pbid --> list of total_cov for each junction so we can calculate SD later
    for r in reader:
        # only need to do assignment if:
        # (1) the .canonical field is still "NA"
        # (2) the junction is non-canonical
        assert r['canonical'] in ('canonical', 'non_canonical')
        if (isoforms_info[r['isoform']].canonical == 'NA') or \
            (r['canonical'] == 'non_canonical'):
            isoforms_info[r['isoform']].canonical = r['canonical']

        if (isoforms_info[r['isoform']].bite == 'NA') or (r['bite_junction'] == 'TRUE'):
            isoforms_info[r['isoform']].bite = r['bite_junction']

        if r['indel_near_junct'] == 'TRUE':
            if isoforms_info[r['isoform']].nIndelsJunc == 'NA':
                isoforms_info[r['isoform']].nIndelsJunc = 0
            isoforms_info[r['isoform']].nIndelsJunc += 1

        # min_cov: min( total_cov[j] for each junction j in this isoform )
        # min_cov_pos: the junction [j] that attributed to argmin(total_cov[j])
        # min_sample_cov: min( sample_cov[j] for each junction in this isoform )
        # sd_cov: sd( total_cov[j] for each junction j in this isoform )
        if r['sample_with_cov'] != 'NA':
            sample_with_cov = int(r['sample_with_cov'])
            if (isoforms_info[r['isoform']].min_samp_cov == 'NA') or (isoforms_info[r['isoform']].min_samp_cov > sample_with_cov):
                isoforms_info[r['isoform']].min_samp_cov = sample_with_cov

        if r['total_coverage'] != 'NA':
            total_cov = int(r['total_coverage'])
            sj_covs_by_isoform[r['isoform']].append(total_cov)
            if (isoforms_info[r['isoform']].min_cov == 'NA') or (isoforms_info[r['isoform']].min_cov > total_cov):
                isoforms_info[r['isoform']].min_cov = total_cov
                isoforms_info[r['isoform']].min_cov_pos = r['junction_number']


    for pbid, covs in sj_covs_by_isoform.iteritems():
        isoforms_info[pbid].sd = pstdev(covs)

    #### Printing output file:

    print >> sys.stderr, "**** Writing output files...."

    with open(outputClassPath, 'w') as h:
        fout_class = DictWriter(h, fieldnames=FIELDS_CLASS, delimiter='\t')
        fout_class.writeheader()
        for r in isoforms_info.itervalues():
            fout_class.writerow(r.as_dict())

    # Now that RTS info is obtained, we can write the final junctions.txt
    with open(outputJuncPath, 'w') as h:
        fout_junc = DictWriter(h, fieldnames=fields_junc_cur, delimiter='\t')
        fout_junc.writeheader()
        for r in DictReader(open(outputJuncPath+"_tmp"), delimiter='\t'):
            if r['isoform'] in RTS_info:
                if r['junction_number'] in RTS_info[r['isoform']]:
                    r['RTS_junction'] = 'TRUE'
                else:
                    r['RTS_junction'] = 'FALSE'
            fout_junc.writerow(r)

    ## Generating report

    print >> sys.stderr, "**** Generating SQANTI report...."
    cmd = RSCRIPTPATH + " {d}/{f} {c} {j}".format(d=utilitiesPath, f=RSCRIPT_REPORT, c=outputClassPath, j=outputJuncPath)
    if subprocess.check_call(cmd, shell=True)!=0:
        print >> sys.stderr, "ERROR running command: {0}".format(cmd)
        sys.exit(-1)
    stop3 = timeit.default_timer()

    print >> sys.stderr, "Removing temporary files...."
    os.remove(outputClassPath+"_tmp")
    os.remove(outputJuncPath+"_tmp")


    print >> sys.stderr, "SQANTI complete in {0} sec.".format(stop3 - start3)


def rename_isoform_seqids(input_fasta):
    """
    Rename input isoform fasta/fastq, which is usually mapped, collapsed Iso-Seq data with IDs like:

    PB.1.1|chr1:10-100|xxxxxx

    to just being "PB.1.1"

    :param input_fasta: Could be either fasta or fastq, autodetect.
    :return: output fasta with the cleaned up sequence ID
    """
    type = 'fasta'
    with open(input_fasta) as h:
        if h.readline().startswith('@'): type = 'fastq'
    f = open(input_fasta[:input_fasta.rfind('.')]+'.renamed.fasta', 'w')
    for r in SeqIO.parse(open(input_fasta), type):
        if r.id.startswith('PB.'):  # PacBio fasta header
            newid = r.id.split('|')[0]
        else:
            raw = r.id.split('|')
            if len(raw) > 4:  # RefSeq fasta header
                newid = raw[3]
            else:
                newid = r.id.split()[0]  # Ensembl fasta header
        f.write(">{0}\n{1}\n".format(newid, r.seq))
    f.close()
    return f.name


class CAGEPeak:
    def __init__(self, cage_bed_filename):
        self.cage_bed_filename = cage_bed_filename
        self.cage_peaks = defaultdict(lambda: IntervalTree()) # (chrom,strand) --> intervals of peaks

        self.read_bed()

    def read_bed(self):
        for line in open(self.cage_bed_filename):
            raw = line.strip().split()
            chrom = raw[0]
            start0 = int(raw[1])
            end1 = int(raw[2])
            strand = raw[5]
            tss0 = int(raw[6])
            self.cage_peaks[(chrom,strand)].insert(start0, end1, (tss0, start0, end1))

    def find(self, chrom, strand, query, search_window=10000):
        """
        :param start0: 0-based start of the 5' end to query
        :return: <True/False falls within a cage peak>, <nearest dist to TSS>
        dist to TSS is 0 if right on spot
        dist to TSS is + if downstream, - if upstream (watch for strand!!!)
        """
        within_peak, dist_peak = False, 'NA'
        for (tss0,start0,end1) in self.cage_peaks[(chrom,strand)].find(query-search_window, query+search_window):
            if not within_peak:
                within_peak, dist_peak = (start0<=query<end1), (query - tss0) * (-1 if strand=='-' else +1)
            else:
                d = (query - tss0) * (-1 if strand=='-' else +1)
                if abs(d) < abs(dist_peak):
                    within_peak, dist_peak = (start0<=query<end1), d
        return within_peak, dist_peak



def main():

    global utilitiesPath

    #arguments
    parser = argparse.ArgumentParser(description="Structural and Quality Annotation of Novel Transcript Isoforms")
    parser.add_argument('isoforms', help='\tIsoforms (Fasta/fastq or gtf format; By default "fasta/fastq". GTF if specified -g option)')
    parser.add_argument('annotation', help='\t\tReference annotation file (GTF format)')
    parser.add_argument('genome', help='\t\tReference genome (Fasta format)')
    parser.add_argument('--cage_peak', help='\t\tFANTOM5 Cage Peak (BED format, optional)')
    parser.add_argument("--polyA_motif_list", help="\t\tRanked list of polyA motifs (text, optional)")
    parser.add_argument("--phyloP_bed", help="\t\tPhyloP BED for conservation score (BED, optional)")
    parser.add_argument("--skipORF", default=False, action="store_true", help="\t\tSkip ORF prediction (to save time)")
    parser.add_argument('-g', '--gtf', help='\t\tUse when running SQANTI by using as input a gtf of isoforms', action='store_true')
    parser.add_argument('-e','--expression', help='\t\tExpression matrix', required=False)
    parser.add_argument('-x','--gmap_index', help='\t\tPath and prefix of the reference index created by gmap_build. Mandatory unless -g option is specified.')
    parser.add_argument('-t', '--gmap_threads', help='\t\tNumber of threads used during alignment by GMAP or Minimap2.', required=False, default="1", type=int)
    parser.add_argument('-z', '--sense', help='\t\tOption that helps GMAP/Minimap2 know that the exons in you cDNA sequences are in the correct sense. Applicable just when you have a high quality set of cDNA sequences', required=False, action='store_true')
    parser.add_argument('-o','--output', help='\t\tPrefix for output files.', required=False)
    parser.add_argument('-d','--dir', help='\t\tDirectory for output files. Default: Directory where the script was run.', required=False)
    parser.add_argument('-c','--coverage', help='\t\tJunction coverage files (provide a single file or a file pattern, ex: "mydir/*.junctions").', required=False)
    parser.add_argument('-s','--sites', default="ATAC,GCAG,GTAG", help='\t\tSet of splice sites to be considered as canonical (comma-separated list of splice sites). Default: GTAG,GCAG,ATAC.', required=False)
    parser.add_argument('-w','--window', default="20", help='\t\tSize of the window in the genomic DNA screened for Adenine content downstream of TTS', required=False, type=int)
    parser.add_argument('--geneid', help='\t\tUse gene_id tag from GTF to define genes. Default: gene_name used to define genes', default=False, action='store_true')
    parser.add_argument('-fl', '--fl_count', help='\t\tFull-length PacBio abundance files (comma-separated list of PacBio abundance files generated by PacBio or directory where there are in).', required=False)
    parser.add_argument("-v", "--version", help="Display program version number.", action='version', version='SQANTI2 '+str(__version__))

    args = parser.parse_args()

    if args.gtf:
        print >> sys.stderr, "--gtf option currently not supported."
        sys.exit(-1)

    # path and prefix for output files
    if args.output is None:
        args.output = os.path.splitext(os.path.basename(args.isoforms))[0]

    if args.dir is None:
        args.dir = os.getcwd()
    else:
        if not os.path.isdir(os.path.abspath(args.dir)):
            print >> sys.stderr, "ERROR: {0} directory doesn't exist. Abort!".format(args.dir)
            sys.exit()
        else:
            args.dir = os.path.abspath(args.dir)

    args.genome = os.path.abspath(args.genome)
    if not os.path.isfile(args.genome):
        print >> sys.stderr, "ERROR: genome fasta {0} doesn't exist. Abort!".format(args.genome)
        sys.exit()

    if not args.gtf:
        args.aligner_choice = "minimap2"
        if args.gmap_index is not None:
            args.aligner_choice = "gmap"
            if not os.path.isdir(os.path.abspath(args.gmap_index)):
                print >> sys.stderr, "GMAP index {0} doesn't exist! Abort.".format(args.gmap_index)
                sys.exit()
            else:
                print("Aligner choice: GMAP.")
        else:
            print("Aligner choice: Minimap2.")


    args.isoforms = os.path.abspath(args.isoforms)
    if not os.path.isfile(args.isoforms):
        print >> sys.stderr, "ERROR: Input isoforms {0} doesn't exist. Abort!".format(args.isoforms)
        sys.exit()

    print >> sys.stderr, "Cleaning up isoform IDs..."
    args.isoforms = rename_isoform_seqids(args.isoforms)
    print >> sys.stderr, "Cleaned up isoform fasta file written to: {0}".format(args.isoforms)


    args.annotation = os.path.abspath(args.annotation)
    if not os.path.isfile(args.annotation):
        print >> sys.stderr, "ERROR: Annotation doesn't exist. Abort!".format(args.annotation)
        sys.exit()

    if args.sense:
        if args.aligner_choice == "gmap":
            args.sense = "sense_force"
        else:
            args.sense = "auto"
    else:
        if args.aligner_choice == "minimap2":
            args.sense = "f"
        else:
            args.sense = "b"

    # Running functionality
    print >> sys.stdout, "**** Running SQANTI..."
    run(args)


if __name__ == "__main__":
    main()
