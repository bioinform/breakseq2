import argparse
import sys
import os
from biopy.io import Fasta
from biopy.io import SV

DEFAULT_JUNCTION_LENGTH = 200
SUPPORTED_FORMAT_VERSIONS = ["1", "2"]
DEFAULT_FORMAT_VERSION = "2"

def add_options(main_parser):
    local_parser = main_parser.add_argument_group("Breakpoint library FASTA generation options")
    local_parser.add_argument("--bplib_gff", help="Breakpoint GFF input", required=False)
    local_parser.add_argument("--junction_length", help="Junction length", type=int, default=DEFAULT_JUNCTION_LENGTH)
    local_parser.add_argument("--format_version", help="Version of breakpoint library format to use", options=breakseq_index.SUPPORTED_FORMAT_VERSIONS, default=breakseq_index.DEFAULT_FORMAT_VERSION) 


def get_seq(sv, jn_type, seq, format_version):
    if format_version == "1":
        return ">%s:%s\n%s" % (sv.id, jn_type, seq)
    return ">%s:%s:%s\n%s" % (sv.id, sv.size(), jn_type, seq)


def generate_bplib(gff, reference, output, junction_length=DEFAULT_JUNCTION_LENGTH, format_version=DEFAULT_FORMAT_VERSION):
    outfd = open(output, "w") if output else sys.stdout
    insertion_sequence_file = gff.replace(".gff", "") + ".ins";
    if not os.path.isfile(insertion_sequence_file):
        raise Exception("Insertion sequence file %s missing" % insertion_sequence_file)
    for sv in SV.parse(gff, Fasta.Seqs(reference, junction_length)):
        flanks = sv.get_flanks()
        if sv.is_insertion():
            if flanks[0] is None or flanks[1] is None:
                raise Exception("No inserted sequence found for insertion %s" % sv.id)
            outfd.write("%s\n" % (get_seq(sv, "A", flanks[0], format_version)))
            outfd.write("%s\n" % (get_seq(sv, "B", flanks[1], format_version)))
        if sv.is_deletion():
            outfd.write("%s\n" % (get_seq(sv, "C", flanks[2], format_version)))
    outfd.close()
