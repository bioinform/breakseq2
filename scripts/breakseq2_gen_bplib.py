#!/usr/bin/env python

import argparse
from breakseq2 import breakseq_index, _version

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generate breakpoint library FASTA from breakpoint GFF",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    breakseq_index.add_options(parser)
    parser.add_argument("--reference", help="Reference FASTA", required=True)
    parser.add_argument("--output", help="Output FASTA to generate. Leave unspecified for stdout")
    parser.add_argument('--version', action='version', version='%(prog)s ' + _version.__version__)
    args = parser.parse_args()

    breakseq_index.generate_bplib(args.bplib_gff, args.reference, args.output, args.junction_length, args.format_version)
