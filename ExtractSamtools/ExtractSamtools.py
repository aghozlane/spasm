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

"""Extract the annotation in case of blast on imomi database."""

from __future__ import print_function
import os
import sys
import argparse
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


#===================
# parameters
#===================
def get_arguments():
    """Extract program options
    """
    parser = argparse.ArgumentParser(description=__doc__,
                                     usage="{0} -h [options] [arg]"
                                     .format(sys.argv[0]))
    parser.add_argument('-i', dest='input_sam', type=isfile,
                        required=True, help='Input samtools count.')
    parser.add_argument('-a', dest='annotation_file', type=isfile,
                        required=True,
                        help='Input blast result file from ExtractNCBIDB.')
    parser.add_argument('-o', '--output_file', dest='detailed_file', type=str,
                        default=(os.curdir + os.sep + "abundance_samtools.txt"),
                        help='Detailed samtools abundance')
    parser.add_argument('-k', '--output_krona', dest='output_krona', type=str,
                        default=(os.curdir + os.sep + "krona_samtools.txt"),
                        help='Output file for krona')
    return parser.parse_args()


def parse_sam_count(input_sam):
    """Parse sam result
    """
    sam_dict = {}
    try:
        with open(input_sam, "rt") as sam:
            sam_reader = csv.reader(sam, delimiter='\t')
            for line in sam_reader:
                if line[0] == "*":
                    sam_dict["Unmapped_reads"] = [line[3]]
                else:
                    sam_dict[line[0]] = [line[2]]
    except IOError:
        sys.exit("Error cannot open {0}".format(input_sam))
    return sam_dict


def parse_annotation(sam_dict, annotation_file):
    """Match the annotation with count table
    """
    try:
        with open(annotation_file, "rt") as annotation:
            annotation_reader = csv.reader(annotation, delimiter='\t')
            # Pass header
            annotation_reader.next()
            for line in annotation_reader:
                if line[0] in sam_dict:
                    sam_dict[line[0]] += line[1:8]
                else:
                    sys.exit("New gene in NCBI annotation and not in sam")
    except IOError:
        sys.exit("Error cannot open {0}".format(annotation_file))
    return sam_dict


def write_krona(sam_dict, output_file, detailed=False):
    """Write krona result
    """
    try:
        with open(output_file, "wt") as output:
            output_writer = csv.writer(output, delimiter="\t")
            for key in sam_dict:
                if detailed:
                    output_writer.writerow([key] + sam_dict[key])
                else:
                    # Add ummapped reads count
                    if key == "Unmapped_reads":
                        output_writer.writerow(sam_dict[key] +
                                               ["Unmapped reads"])
                    else:
                        output_writer.writerow(sam_dict[key])
    except IOError:
        sys.exit("Error cannot open {0}".format(output_file))


#===================
# MAIN
#===================
def main():
    """Main program
    """
    args = get_arguments()
    sam_dict = parse_sam_count(args.input_sam)
    sam_dict = parse_annotation(sam_dict, args.annotation_file)
    write_krona(sam_dict, args.output_krona)
    write_krona(sam_dict, args.detailed_file, True)

if __name__ == "__main__":
    main()
# END
