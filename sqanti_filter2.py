__author__  = "etseng@pacb.com"
__version__ = '3.8'

"""
Lightweight filtering of SQANTI by using .classification.txt output

Only keep Iso-Seq isoforms if:
The isoform is FSM, ISM, or NIC and (does not have intrapriming or has polyA_motif)
The isoform is NNC, does not have intrapriming/or polyA motif, not RT-switching, and all junctions are either all canonical or short-read-supported
The isoform is antisense, intergenic, genic, does not have intrapriming/or polyA motif, not RT-switching, and all junctions are either all canonical or short-read-supported

"""

import os, sys, argparse, subprocess
import distutils.spawn
from csv import DictReader, DictWriter
from Bio import SeqIO
from cupcake.io.BioReaders import GMAPSAMReader

utilitiesPath =  os.path.dirname(os.path.realpath(__file__))+"/utilities/"
RSCRIPTPATH = distutils.spawn.find_executable('Rscript')
RSCRIPT_REPORT = 'SQANTI_report2.R'

if os.system(RSCRIPTPATH + " --version")!=0:
    print >> sys.stderr, "Rscript executable not found! Abort!"
    sys.exit(-1)


CATEGORY_DICT = {'full-splice_match': 'FSM',
                 'incomplete-splice_match': 'ISM',
                 'novel_in_catalog': 'NIC',
                 'novel_not_in_catalog': 'NNC',
                 'antisense': 'AS',
                 'intergenic': 'intergenic',
                 'genic_intron': 'intron',
                 'genic': 'genic',
                 'fusion': 'fusion'}

def sqanti_filter_lite(args):

    fafq_type = 'fasta'
    with open(args.isoforms) as h:
        if h.readline().startswith('@'): fafq_type = 'fastq'

    prefix = args.sqanti_class[:args.sqanti_class.rfind('.')]

    fcsv = open(prefix + '.filtered_lite_reasons.txt', 'w')
    fcsv.write("# classification: {0}\n".format(args.sqanti_class))
    fcsv.write("# isoform: {0}\n".format(args.isoforms))
    fcsv.write("# intrapriming cutoff: {0}\n".format(args.intrapriming))
    fcsv.write("# min_cov cutoff: {0}\n".format(args.min_cov))
    fcsv.write("filtered_isoform,reason\n")

    fout = open(prefix + '.filtered_lite.' + fafq_type, 'w')

    seqids_to_keep = []
    total_count = 0
    for r in DictReader(open(args.sqanti_class), delimiter='\t'):
        total_count += 1
        filter_flag, filter_msg = False, ""
        percA = float(r['perc_A_downstream_TTS']) / 100
        assert 0 <= percA <= 1
        min_cov = float(r['min_cov']) if r['min_cov']!='NA' else None
        num_exon = int(r['exons'])
        is_RTS = r['RTS_stage'] == 'TRUE'
        is_canonical = r['all_canonical']=='canonical'

        cat = CATEGORY_DICT[r['structural_category']]

        if cat in ['FSM', 'ISM', 'NIC']:
            if (percA >= args.intrapriming and r['polyA_motif']=='NA' and
                    (r['diff_to_gene_TSS']=='NA' or abs(int(r['diff_to_gene_TTS'])) > args.max_dist_to_known_end)):
                filter_flag, filter_msg = True, "IntraPriming"
        else:
            if (percA >= args.intrapriming and r['polyA_motif']=='NA' and
                    (r['diff_to_gene_TSS']=='NA' or abs(int(r['diff_to_gene_TTS'])) > args.max_dist_to_known_end)):
                filter_flag, filter_msg = True, "IntraPriming"
            elif is_RTS:
                filter_flag, filter_msg = True, "RTSwitching"
            elif (not is_canonical) and (min_cov is None or (min_cov is not None and min_cov < args.min_cov)):
                filter_flag, filter_msg = True, "LowCoverage/Non-Canonical"

        if not filter_flag:
            seqids_to_keep.append(r['isoform'])
        else:
            fcsv.write("{0},{1}\n".format(r['isoform'], filter_msg))

    print >> sys.stdout, "{0} isoforms read from {1}. {2} to be kept.".format(total_count, args.sqanti_class, len(seqids_to_keep))

    for r in SeqIO.parse(open(args.isoforms), fafq_type):
        if r.id in seqids_to_keep:
            SeqIO.write(r, fout, fafq_type)
    fout.close()
    print >> sys.stdout, "Output written to: {0}".format(fout.name)


    # write out a new .classification.txt, .junctions.txt
    outputClassPath = prefix + '.filtered_lite_classification.txt'
    with open(outputClassPath, 'w') as f:
        reader = DictReader(open(args.sqanti_class), delimiter='\t')
        writer = DictWriter(f, reader.fieldnames, delimiter='\t')
        writer.writeheader()
        for r in reader:
            if r['isoform'] in seqids_to_keep:
                writer.writerow(r)
        print >> sys.stdout, "Output written to: {0}".format(f.name)

    outputJuncPath = prefix + '.filtered_lite_junctions.txt'
    with open(outputJuncPath, 'w') as f:
        reader = DictReader(open(args.sqanti_class.replace('_classification', '_junctions')), delimiter='\t')
        writer = DictWriter(f, reader.fieldnames, delimiter='\t')
        writer.writeheader()
        for r in reader:
            if r['isoform'] in seqids_to_keep:
                writer.writerow(r)
        print >> sys.stdout, "Output written to: {0}".format(f.name)

    outputSam = prefix + '.filtered_lite.sam'
    with open(outputSam, 'w') as f:
        reader = GMAPSAMReader(args.sam_file, True)
        f.write(reader.header)
        for r in reader:
            if r.qID in seqids_to_keep:
                f.write(r.record_line + '\n')
        print >> sys.stdout, "Output written to: {0}".format(f.name)


    print >> sys.stderr, "**** Generating SQANTI report...."
    cmd = RSCRIPTPATH + " {d}/{f} {c} {j}".format(d=utilitiesPath, f=RSCRIPT_REPORT, c=outputClassPath, j=outputJuncPath)
    if subprocess.check_call(cmd, shell=True)!=0:
        print >> sys.stderr, "ERROR running command: {0}".format(cmd)
        sys.exit(-1)



def main():

    parser = argparse.ArgumentParser(description="Filtering of Isoforms based on SQANTI attributes")
    parser.add_argument('sqanti_class', help='\t\tSQANTI classification output file.')
    parser.add_argument('isoforms', help='\t\tfasta/fastq isoform file to be filtered by SQANTI')
    parser.add_argument('sam_file', help='\t\tSAM alignment of the input fasta/fastq')
    parser.add_argument('-a',"--intrapriming", type=float, default=0.8, help='\t\tAdenine percentage at genomic 3\' end to flag an isoform as intra-priming (default: 0.8)')
    parser.add_argument('-m',"--max_dist_to_known_end", type=int, default=50, help="\t\tMaximum distance to an annotated 3' end to preserve as a valid 3' end and not filter out (default: 50bp)")
    parser.add_argument("-c", "--min_cov", type=int, default=3, help="\t\tMinimum junction coverage for each isoform (only used if min_cov field is not 'NA'), default: 3")
    #parser.add_argument("--always_keep_canonical", default=False, action="store_true", help="Always keep isoforms with all canonical junctions, regardless of other criteria. (default: False)")

    args = parser.parse_args()

    assert 0 <= args.intrapriming <= 1

    args.sqanti_class = os.path.abspath(args.sqanti_class)
    if not os.path.isfile(args.sqanti_class):
        print >> sys.stderr, "ERROR: {0} doesn't exist. Abort!".format(args.sqanti_class)
        sys.exit(-1)

    if args.isoforms is not None and not os.path.isfile(args.isoforms):
        print >> sys.stderr, "ERROR: {0} doesn't exist. Abort!".format(args.isoforms)
        sys.exit(-1)

    # if args.dir is None:
    #     args.dir = os.getcwd() + "/Filter_out"
    #     if not os.path.exists(args.dir):
    #         os.makedirs(args.dir)
    # else:
    #     args.dir = os.path.abspath(args.dir)
    #     if not os.path.exists(args.dir):
    #         os.makedirs(args.dir)

    print >> sys.stdout, "\nRunning SQANTI filtering...\n"

    sqanti_filter_lite(args)



if __name__ == "__main__":
    main()
