#!/usr/bin/env python3

import sys
import argparse
import csv
import re

def get_args():
    help_usage="A script to convert a tab-delimited file generated with 'vcf-to-tab' in VCFtools into a fasta file.\nv1.1 (June 1, 2020) by Satoshi Ansai."
    parser = argparse.ArgumentParser(usage=help_usage)

    # If there is no std-in, an input file is required
    if sys.stdin.isatty():
        parser.add_argument("input", help="an input tab file", type=str)

    parser.add_argument("-m", "--miss", help="Set the substitute letter(s) for missing bases. Default: \"-\".", default='-', type=str)
    parser.add_argument("-i", "--iupac", help="Heterozygous SNPs are output as IUPAC code.", action="store_true")
    parser.add_argument("-j", "--homo", help="Either of two variants in each heterozygous site is output.", action="store_true")

    args = parser.parse_args()

    return(args)

def main():
    args = get_args()

    # read an input file or stdin and store as "l"
    if hasattr(args, 'input'):
        with open(args.input) as f:
            reader = csv.reader(f, delimiter="\t")
            l = [row for row in reader]
    else:
        reader = csv.reader(sys.stdin, delimiter="\t")
        l = [row for row in reader]

    l_T = [list(x) for x in zip(*l)] # transpose

    def multiple_replace(text, adict):
        rx = re.compile('|'.join(map(re.escape,adict)))
        def one_xlat(match):
            return adict[match.group(0)]
        return rx.sub(one_xlat, text)

    if args.iupac:
        adict = {
        "G/G"   :"G",
        "C/C"   :"C",
        "T/T"   :"T",
        "A/A"   :"A",
        "G/T"   :"K",
        "T/G"   :"K",
        "A/C"   :"M",
        "C/A"   :"M",
        "C/G"   :"S",
        "G/C"   :"S",
        "A/G"   :"R",
        "G/A"   :"R",
        "A/T"   :"W",
        "T/A"   :"W",
        "C/T"   :"Y",
        "T/C"   :"Y",
        "./."   :args.miss,
        }
    elif args.homo:
        adict = {
        "G/G"   :"G",
        "C/C"   :"C",
        "T/T"   :"T",
        "A/A"   :"A",
        "G/T"   :"G",
        "T/G"   :"T",
        "A/C"   :"A",
        "C/A"   :"C",
        "C/G"   :"C",
        "G/C"   :"G",
        "A/G"   :"A",
        "G/A"   :"G",
        "A/T"   :"A",
        "T/A"   :"T",
        "C/T"   :"C",
        "T/C"   :"T",
        "./."   :args.miss,
        }
    else:
        adict = {
        "G/G"   :"G",
        "C/C"   :"C",
        "T/T"   :"T",
        "A/A"   :"A",
        "G/T"   :"N",
        "T/G"   :"N",
        "A/C"   :"N",
        "C/A"   :"N",
        "C/G"   :"N",
        "G/C"   :"N",
        "A/G"   :"N",
        "G/A"   :"N",
        "A/T"   :"N",
        "T/A"   :"N",
        "C/T"   :"N",
        "T/C"   :"N",
        "./."   :args.miss,
        }

    # output as a multi-sample fasta
    for line in l_T[3:]:
        print(">", line[0], sep="")
        line_rep = [multiple_replace(s, adict) for s in line[1:]]
        print(''.join(line_rep))

if __name__ == '__main__':
    main()
