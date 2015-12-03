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

"""Extract a list of sequence from the catalogue."""


from __future__ import print_function
import argparse
import sys
import os
import csv
import bisect


__author__ = "Amine Ghozlane"
__copyright__ = "Copyright 2014, INRA"
__credits__ = ["Amine Ghozlane"]
__license__ = "GPL"
__version__ = "1.0.0"
__maintainer__ = "Amine Ghozlane"
__email__ = "amine.ghozlane@jouy.inra.fr"
__status__ = "Developpement"


def isfile(path):
    """Check if path is an existing file.
      Arguments:
          path: Path to the file
    """
    # from Jonathan Barnoud
    if not os.path.isfile(path):
        if os.path.isdir(path):
            msg = "{0} is a directory".format(path)
        else:
            msg = "{0} does not exist.".format(path)
        raise argparse.ArgumentTypeError(msg)
    return path


def isdir(path):
    """Check if path is an existing file.
      Arguments:
          path: Path to the file
    """
    if not os.path.isdir(path):
        if os.path.isfile(path):
            msg = "{0} is a file".format(path)
        else:
            msg = "{0} does not exist.".format(path)
        raise argparse.ArgumentTypeError(msg)
    return path


def getArguments():
    """Retrieves the arguments of the program.
      Returns: An object that contains the arguments
    """
    # Parsing arguments
    parser = argparse.ArgumentParser(description=__doc__, usage=
                                     "{0} -h".format(sys.argv[0]))
    parser.set_defaults(results=".{0}".format(os.sep))
    parser.add_argument('-i', dest='list_sequences_file', type=isfile,
                        required=True,
                        help='List of sequence to extract.')
    parser.add_argument('-d', dest='catalogue_file', type=isfile,
                        required=True,
                        help='Database query.')
    parser.add_argument('-o', dest='output_file', type=str,
                        help='Ouput file.')
    parser.add_argument('-r', dest='results', type=isdir,
                        help='Path to result directory.')
    return parser.parse_args()


def extract_interest_elements(list_sequences_file):
    """Get a list of the element of interest
    """
    list_seq = []
    try:
        with open(list_sequences_file, "rt") as list_sequences:
            list_sequences_reader = csv.reader(list_sequences)
            for line in list_sequences_reader:
                 list_seq.append(line[0])
            assert(list_sequences > 0)
            list_seq.sort()
    except IOError:
        sys.exit("Error cannot the file : {0}".format(list_sequences_file))
    except AssertionError:
        sys.exit("Error no element detected in the file : {0}"
                 .format(list_sequences_file))
    return list_seq


# def is_selected(header, list_sequences):
#     """
#     """
#     for element in list_sequences:
#         if element in header:
#             return element
#     return None


def get_element(name, input_list):
    """Search name in input list
      Arguments:
        input_list: List
        name: Search criteria
    """
    # Searching the node with its name
    i = bisect.bisect_left(input_list, name)
    # Object has been found
    if(i != len(input_list) and input_list[i] == name):
        return input_list[i]
    return None


def extract_catalogue_sequence(list_sequences, catalogue_file):
    """
    """
    grab_sequence = False
    interest_sequence = {}
    try:
        with open(catalogue_file, "rt") as catalogue:
            for line in catalogue:
                if line[0] == ">":
                    grab_sequence = False
		    selection = get_element(line[1:].split(" ")[0], list_sequences)
                    if selection:
                        interest_sequence[selection] = ""
                        grab_sequence = True
                elif grab_sequence and len(line):
                    interest_sequence[selection] += line.replace("\n", "").replace("\r", "")
            assert(len(interest_sequence) > 0)
    except IOError:
        sys.exit("Error cannot the file : {0}".format(catalogue_file))
    except AssertionError:
        sys.exit("Error no element detected in the file : {0}"
                 .format(catalogue_file))
    return interest_sequence


def fill(text, width=80):
    """Split text"""
    return os.linesep.join(text[i:i+width] for i in xrange(0, len(text), width))


def write_interest_sequence(interest_sequence, output_file):
    """
    """
    try:
        with open(output_file, "wt") as output:
            for key in interest_sequence:
                output.write(">{1}{0}{2}{0}".format(
                                os.linesep, key, fill(interest_sequence[key])))
    except IOError:
        sys.exit("Error cannot open {0}".format(output_file))


#==============================================================
# Main program
#==============================================================
def main():
    """
    Main program function
    """
    # Get the arguments
    args = getArguments()
    # Get List of sequence of interest
    print("Get List of sequence of interest ...")
    list_sequences = extract_interest_elements(args.list_sequences_file)
    # Extract catalogue sequence
    print("Extract sequences from the catalogue...")
    interest_sequence = extract_catalogue_sequence(list_sequences,
                                                   args.catalogue_file)
    # Write sequences
    print("Write sequences to {0}".format(args.output_file))
    write_interest_sequence(interest_sequence, args.output_file)
    print("Done.")


if __name__ == '__main__':
    main()
