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


"""Extract Cazy annotation."""


from __future__ import print_function
import argparse
import os
import sys
import csv


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
    parser = argparse.ArgumentParser(description=__doc__, usage=
                                     "{0} -h".format(sys.argv[0]))
    parser.add_argument('-i', dest='hmmer_result', type=isfile,
                        required=True, help='Hmmer result against cazy hmm '
                        'datablase (type --tblout with hmmer)')
    parser.add_argument('-f', dest='faminfo', type=isfile,
                        required=True, help='Family information file')
    parser.add_argument('-o', dest='output_file', type=str, help='Output file')
    parser.add_argument('-r', dest='results', type=isdir,
                        default=os.curdir + os.sep,
                        help='Path to result directory.')
    return parser.parse_args()


def get_faminfo(faminfo_file):
    """Parse family information
    """
    faminfo_dict = {}
    try:
        with open(faminfo_file, "rt") as faminfo:
            faminfo_reader = csv.reader(faminfo, delimiter="\t")
            # Pass header
            faminfo_reader.next()
            for line in faminfo_reader:
                nb_item = len(line)
                if nb_item < 5:
                    faminfo_dict[line[0]] = (
                        line[1:] + ["NA" for i in xrange(5 - nb_item)])
                else:
                    faminfo_dict[line[0]] = ["NA" if not e else e
                                             for e in line[1:]]
            assert(len(faminfo_dict) > 0)
    except IOError:
        sys.exit("Error cannot open {0}".format(faminfo_file))
    except AssertionError:
        sys.exit("Error nothing read from {0}".format(faminfo_file))
    return faminfo_dict


def get_hmmer_result(hmmer_result_file, faminfo_desc):
    """Parse hmmer result and set the annotation
    """
    annotation_result = []
    try:
        with open(hmmer_result_file, "rt") as hmmer_result:
            for line in hmmer_result:
                if line[0] != "#":
                    splited_line = line.split()
                    family = splited_line[0].split(".")[0]
                    if family in faminfo_desc:
                        annotation_result += [[splited_line[2], family] +
                                              faminfo_desc[family]]
                    else:
                        annotation_result += [[splited_line[2], family, "NA",
                                               "NA", "NA", "NA"]]
            assert(len(annotation_result) > 0)
    except IOError:
        sys.exit("Error cannot open {0}".format(hmmer_result_file))
    except AssertionError:
        sys.exit("Error nothing read from {0}".format(hmmer_result_file))
    return annotation_result


def write_result(annotation_result, results, output_file):
    """Write the result table
    """
    if not output_file:
        output_file = results + os.sep + "egg_nog_annotation.txt"
    try:
        with open(output_file, "wt") as output:
            output_writer = csv.writer(output, delimiter='\t')
            output_writer.writerow(["SequenceName", "Family", "ncbi-cdd",
                                    "cazy-class", "cazy-note",
                                    "cazy-activities"])
            output_writer.writerows(annotation_result)
    except IOError:
        sys.exit("Error cannot open {0}".format(output_file))


#==============================================================
# Main program
#==============================================================
def main():
    """
    Main program function
    """
    # Get arguments
    args = get_arguments()
    # Parse faminfo
    faminfo_desc = get_faminfo(args.faminfo)
    # Parse hmmer result
    annotation_result = get_hmmer_result(args.hmmer_result, faminfo_desc)
    # Write result
    write_result(annotation_result, args.results, args.output_file)


if __name__ == '__main__':
    main()
