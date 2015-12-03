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

"""Filter the contigs depending on their length."""

from __future__ import print_function
import os
import sys
from argparse import ArgumentParser


__author__ = "Mathieu Almeida, Amine Ghozlane"
__copyright__ = "Copyright 2014, INRA"
__license__ = "GPL"
__version__ = "1.0.0"
__maintainer__ = "Mathieu Almeida, Amine Ghozlane"
__email__ = "mathieu.almeida77@jouy.inra.fr, amine.ghozlane@jouy.inra.fr"
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


#===================
#parameters
#===================
def config_parameters():
    """Extract program options
    """
    parser = ArgumentParser(description=__doc__,
                            usage="{0} -h [options] [arg]".format(sys.argv[0]))
    parser.add_argument('-f', '--readFasta', dest='fasta_filename', type=isfile,
                        required=True, help='Input fasta file.')
    parser.add_argument('-t', '--trim', dest='trimsize', type=int, default=300,
                        help='Keep the contigs with trimSize or '
                        'more nucleotides.')
    parser.add_argument('-s', '--sample', dest='samplename', type=str,
                        default='', help='Define a sample name')
    parser.add_argument('-o', '--outfile', dest='out_filename', type=str,
                        default=None, help='Define ouput filename')
    return parser.parse_args()


#===========================================
#remove contigs or scaffolds shorter than trimSize
#===========================================
def parse_fasta_file(fasta_filename):
    """Parse the fasta file and remove contigs or scaffolds shorter than
    trimsize
    """
    nb_contig = 0
    fragment =''
    header=''
    seq_dict = {}
    try:
        with open(fasta_filename, "rt") as fastafile:
            for line in fastafile:
                if line[0] == '>':
                    nb_contig = nb_contig + 1
                    header = line[1:].strip()
                    seq_dict[header] = ""
                elif len(line) > 0 and len(header) > 0:
                    seq_dict[header] += line.strip().replace(os.linesep, "")
                elif len(line) > 0 and len(header) == 0:
                    sys.exit("Strange sequence before any header")
            assert(len(seq_dict) > 0)
    except IOError:
        sys.exit("Error cannot open {0}".format(fasta_filename))
    except AssertionError:
        sys.exit("Error no sequence read from {0}".format(fasta_filename))
    return seq_dict

def write_sequence(seq_dict, trimsize, out_filename):
    """
    """
    nbcontig_output = 0
    try:
        with open(out_filename, "wt") as outfile:
            for header in seq_dict:
                if len(seq_dict[header]) >= trimsize:
                    outfile.write(">{1}{0}{2}{0}".format(
                                  os.linesep, header, seq_dict[header]))
                    nbcontig_output = nbcontig_output + 1
    except IOError:
        sys.exit("Error cannot open {0}".format(out_filename))
    print('Number of sequence in the output: {0}'.format(nbcontig_output))


#===================
#MAIN
#===================
def main():
    """Main function
    """
    args = config_parameters()
    seq_dict = parse_fasta_file(args.fasta_filename)
    print('Number of sequence in the input: {0}'.format(len(seq_dict)))
    if not args.out_filename:
        args.out_filename = ('.'.join(args.fasta_filename.split('.')[0:-1]) +
                             '.' + str(args.trimsize) + '.fasta')
    write_sequence(seq_dict, args.trimsize, args.out_filename)


if __name__ == '__main__':
    main()
