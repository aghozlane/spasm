#!/usr/bin/env python
# -*- coding: utf-8 -*-
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#    A copy of the GNU General Public License is available at
#    http://www.gnu.org/licenses/gpl-3.0.html


"""Extract Metagenemark produced sequence."""

from __future__ import print_function
import argparse
import os
import sys
import re


__author__ = "Amine Ghozlane"
__copyright__ = "Copyright 2014, INRA"
__license__ = "GPL"
__version__ = "1.0.0"
__maintainer__ = "Amine Ghozlane"
__email__ = "amine.ghozlane@jouy.inra.fr"
__status__ = "Developpement"


def isfile(path):
    """Check if path is an existing file.
      :Parameters:
          path: Path to the file
    """
    if not os.path.isfile(path):
        if os.path.isdir(path):
            msg = "{0} is a directory".format(path)
        else:
            msg = "{0} does not exist.".format(path)
        raise argparse.ArgumentTypeError(msg)
    return path


def isdir(path):
    """Check if path is an existing file.
      :Parameters:
          path: Path to the file
    """
    if not os.path.isdir(path):
        if os.path.isfile(path):
            msg = "{0} is a file.".format(path)
        else:
            msg = "{0} does not exist.".format(path)
        raise argparse.ArgumentTypeError(msg)
    return path


def get_arguments():
    """Retrieves the arguments of the program.
      Returns: An object that contains the arguments
    """
    # Parsing arguments
    parser = argparse.ArgumentParser(description=__doc__, usage="{0} -h"
                                     .format(sys.argv[0]))
    parser.add_argument('-i', dest='metagenemark_file', type=isfile,
                        required=True,
                        help="Metagenemark lst file (with gene and protein)")
    parser.add_argument('-a', dest='protein_file', type=str,
                        help="Protein output file")
    parser.add_argument('-d', dest='dna_file', type=str,
                        help="Gene output file")
    parser.add_argument('-o', dest='output_name', type=str,
                        default=os.curdir + os.sep + 'sample',
                        help="Output file name")
    return parser.parse_args()


def parse_metagenemark(metagenemark_file):
    """Parse MetaGeneMark file output to generate prot and nt file
    """
    gene_dict = {}
    protein_dict = {}
    gene_seq = False
    protein_seq = False
    info = {}
    node_reg = re.compile(r"FASTA\s+definition\s+line:\s+(\S+)")
    property_reg = re.compile(r"\s*([0-9]+)\s+[+-]\s+(\S+)\s+(\S+)\s+[0-9]+\s+\S+")
    property_dict = {}
    try:
        with open(metagenemark_file, 'rt') as mgm_file:
            for line in mgm_file:
                property_match = property_reg.match(line)
                node_match = node_reg.match(line)
                if node_match: # Get the node
                    del(property_dict)
                    property_dict = {}
                    new_node = node_match.group(1)
                if property_match: # Get completeness
                    partial = ""
                    if property_match.group(2).startswith("<"):
                        partial += "1"
                    else:
                        partial += "0"
                    if property_match.group(3).startswith(">"):
                        partial += "1"
                    else:
                        partial += "0"
                    property_dict[new_node + "_" + property_match.group(1)] = partial
                if line.startswith('>'): #if header detected
                    title = line[1:].replace("\n", "").replace("\r", "").split("\t")
                    gene = title[0].split('|')
                    gene_number = gene[0].split('_')[1]
                    # node + gene_number (like prodigal)
                    geneid = title[1][1:] + "_" + gene_number
                    # Hmm algo, length, strand, begin, end, partial
                    info[geneid] = gene[1:] + [property_dict[geneid]]
                    if gene[2].endswith('_nt'): # DNA sequence
                        gene_dict[geneid] = ""
                        gene_seq = True
                        protein_seq = False
                    else: # Protein sequence
                        protein_dict[geneid] = ""
                        protein_seq = True
                        gene_seq = False
                elif gene_seq and len(line) > 1:
                    gene_dict[geneid] += line.replace("\n", "").replace("\r", "")
                elif protein_seq and len(line) > 1:
                    protein_dict[geneid] += line.replace("\n", "").replace("\r", "")
                else:
                    gene_seq = False
                    protein_seq = False
            assert(len(gene_dict)  == len(protein_dict)
                   and len(gene_dict) == len(info) and len(gene_dict) > 0)
    except IOError:
        sys.exit("Error cannot open {0}".format(metagenemark_file))
    except AssertionError:
        sys.exit("There is something wrong with "
                 "the parsing of {0}".format(metagenemark_file))
    return gene_dict, protein_dict, info


def fill(text, width=80):
    """Split text"""
    return os.linesep.join(text[i:i+width] for i in xrange(0, len(text), width))


def write_data(data_dict, info, output_filename):
    """Write gene and protein file
    """
    try:
        with open(output_filename, "wt") as output:
            for item in data_dict:
                output.write(">{1} # {2[3]} # {2[4]} # {2[2]} # {2[1]};"
                             "partial={2[5]};{2[0]}{0}{3}{0}".format(
                                 os.linesep, item, info[item],
                                 fill(data_dict[item])))
    except IOError:
        sys.exit("Error cannot open {0}".format(output_filename))


#==============================================================
# Main program
#==============================================================
def main():
    """
    Main program function
    """
    # Get arguments
    args = get_arguments()
    # Parse Metagenemark file
    (gene_dict, protein_dict, info) = parse_metagenemark(args.metagenemark_file)
    # Extract data
    if args.dna_file:
        write_data(gene_dict, info, args.dna_file)
    else:
        write_data(gene_dict, info, args.output_name + "_gene.fna")
    if args.protein_file:
        write_data(protein_dict, info, args.protein_file)
    else:
        write_data(protein_dict, info, args.output_name + "_protein.faa")


if __name__ == '__main__':
    main()
