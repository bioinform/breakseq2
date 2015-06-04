#!/usr/bin/env python

import argparse
import sys
from biopy.io import Fasta
from biopy.io import SV

DEFAULT_JUNCTION_LENGTH = 200


def add_options(main_parser):
    local_parser = main_parser.add_argument_group("Breakpoint library FASTA generation options")
    local_parser.add_argument("--bplib_gff", help="Breakpoint GFF input", required=False)
    local_parser.add_argument("--junction_length", help="Junction length", type=int, default=DEFAULT_JUNCTION_LENGTH)


def get_seq(sv, sv_type, seq):
    return ">%s:%s:%s\n%s" % (sv.id, sv.size(), sv_type, seq)


def generate_bplib(gff, reference, output, junction_length=DEFAULT_JUNCTION_LENGTH):
    outfd = open(output, "w") if output else sys.stdout
    insertion_sequence_file = gff.replace(".gff", "") + ".ins";
    if not os.path.isfile(insertion_sequence_file):
        raise Exception("Insertion sequence file %s missing" % insertion_sequence_file)
    for sv in SV.parse(gff, Fasta.Seqs(reference, junction_length)):
        flanks = sv.get_flanks()
        if sv.is_insertion():
            if flanks[0] is None or flanks[1] is None:
                raise Exception("No inserted sequence found for insertion %s" % sv.id)
            outfd.write("%s\n" % (get_seq(sv, "A", flanks[0])))
            outfd.write("%s\n" % (get_seq(sv, "B", flanks[1])))
        if sv.is_deletion():
            outfd.write("%s\n" % (get_seq(sv, "C", flanks[2])))
    outfd.close()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generate breakpoint library FASTA from breakpoint GFF",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    add_options(parser)
    parser.add_argument("--reference", help="Reference FASTA", required=True)
    parser.add_argument("--output", help="Output FASTA to generate. Leave unspecified for stdout")
    args = parser.parse_args()

    generate_bplib(args.bplib_gff, args.reference, args.output, args.junction_length)

