#!/usr/bin/env python

import argparse
from biopy.io import Fasta
from biopy.io import SV

DEFAULT_JUNCTION_LENGTH = 200


def add_options(main_parser):
    local_parser = main_parser.add_argument_group("Breakpoint library FASTA generation options")
    local_parser.add_argument("--bplib_gff", help="Breakpoint GFF input", required=True)
    local_parser.add_argument("--junction_length", help="Junction length", type=int, default=DEFAULT_JUNCTION_LENGTH)


def print_seq(sv, sv_type, seq):
    print ">%s:%s:%s\n%s" % (sv.id, sv.size(), sv_type, seq)


def generate_bplib(gff, reference, junction_length=DEFAULT_JUNCTION_LENGTH):
    for sv in SV.parse(gff, Fasta.Seqs(reference, junction_length)):
        flanks = sv.get_flanks()
        if sv.is_insertion():
            if flanks[0] is None or flanks[1] is None:
                raise Exception("No inserted sequence found for insertion %s" % sv.id)
            print_seq(sv, "A", flanks[0])
            print_seq(sv, "B", flanks[1])
        if sv.is_deletion():
            print_seq(sv, "C", flanks[2])


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generate breakpoint library FASTA from breakpoint GFF",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    add_options(parser)
    parser.add_argument("--reference", help="Reference FASTA", required=True)
    args = parser.parse_args()

    generate_bplib(args.bplib_gff, args.reference, args.junction_length)

