#!/usr/bin/env python

import logging
import os
import pysam
import preprocess_and_align
import breakseq_core
import breakseq_post
import compute_zygosity
import gen_vcf
from _version import __version__


def add_options(main_parser):
    preprocess_and_align.add_options(main_parser)
    breakseq_core.add_options(main_parser)
    breakseq_post.add_options(main_parser)
    compute_zygosity.add_options(main_parser)
    gen_vcf.add_options(main_parser)

    main_parser.add_argument("--nthreads", help="Number of processes to use for parallelism", type=int, default=1)
    main_parser.add_argument("--bams", help="Alignment BAMs", nargs="+", required=True, default=[])
    main_parser.add_argument("--work", help="Working directory", default="work")
    main_parser.add_argument("--chromosomes", nargs="+", help="List of chromosomes to process", default=[])
    main_parser.add_argument("--reference", help="Reference FASTA", required=True)
    main_parser.add_argument("--sample", help="Sample name. Leave unspecified to infer sample name from BAMs.")
    main_parser.add_argument("--keep_temp", help="Keep temporary files", action="store_true")
    main_parser.add_argument('--version', action='version', version='%(prog)s ' + __version__)


def infer_sample(bam):
    samfile = pysam.Samfile(bam, "rb")
    if "RG" not in samfile.header:
        raise Exception("Unable to infer sample name from %s since RG is missing" % bam)
    samples = list(set([item["SM"] for item in samfile.header["RG"]]))
    if len(samples) > 1:
        raise Exception("Multiple samples found: %s" % (", ".join(samples)))
    samfile.close()
    return samples[0]


def breakseq2_workflow(sample=None, bplib=None, bwa=None, samtools=None, bams=[], work="work", chromosomes=[],
                       nthreads=1, min_span=breakseq_core.DEFAULT_MIN_SPAN,
                       min_overlap=compute_zygosity.DEFAULT_MIN_OVERLAP, reference=None, keep_temp=False, window=compute_zygosity.DEFAULT_WINDOW):
    func_logger = logging.getLogger(breakseq2_workflow.__name__)

    bams = [os.path.abspath(bam) for bam in args.bams]

    if not bams:
        func_logger.error("No BAMs specified so nothing to do")
        return

    if not sample:
        sample = infer_sample(bams[0])

    if not os.path.isdir(work):
        func_logger.info("Created working directory %s" % work)
        os.makedirs(work)

    aligned_bams = preprocess_and_align.parallel_preprocess_and_align(bplib, bwa, samtools, bams, work, chromosomes,
                                                                      nthreads, keep_temp)

    if not aligned_bams:
        func_logger.warn("Read-extraction and alignment generated nothing")
        return

    breakseq_core.breakseq_core(aligned_bams, "%s/breakseq.out" % work, min_span=min_span)
    breakseq_post.generate_final_gff(["%s/breakseq.out" % work], "%s/breakseq.gff" % work)
    compute_zygosity.compute_zygosity(bams, window, "%s/breakseq.gff" % work, "%s/breakseq_genotyped.gff" % work,
                                      min_overlap)
    gen_vcf.gff_to_vcf(reference, "%s/breakseq_genotyped.gff" % work, sample, "%s/breakseq.vcf" % work)
