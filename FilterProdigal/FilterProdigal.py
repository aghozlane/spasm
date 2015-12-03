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

"""Remove incomplete set of fragment in a prodigual file."""


from __future__ import print_function
import argparse
import sys
import os


__author__ = "Mathieu Almeida, Amine Ghozlane"
__copyright__ = "Copyright 2014, INRA"
__license__ = "GPL"
__version__ = "1.0.0"
__maintainer__ = "Mathieu Almeida,  Amine Ghozlane"
__email__ = "mathieu.almeida@jouy.inra.fr, amine.ghozlane@jouy.inra.fr"
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


def config_parameters():
    """Extract program options
    """
    parser = argparse.ArgumentParser(description=__doc__,
                            usage="{0} -h [options] [arg]".format(sys.argv[0]))
    parser.add_argument('-i', '--fragment_file', dest='fragment_file',
                        type=isfile, required=True,
                        help='Input fragment file. Warning : only *.fasta')
    parser.add_argument('-c', '--choices', dest='choices',
                        type=str, nargs="+", default=["11"],
                        choices=["00", "01","10","11"],
                        help="""Remove fragment incomplete in :
                                '5 term: 10 ; 3 term: 01 ; both term: 11""")
    parser.add_argument('-o', '--output', dest='output_file', type=str,
                        help='Remove fragment depending on the choices.')
    parser.add_argument('-r', dest='results', type=isdir,
                        default=os.curdir + os.sep, help='Path to result '
                        'directory.')
    return parser.parse_args()


def select_fragment(fragment_file_name, choices):
    """Remove fragment that match cut choices condition :
       missing lack 3 term, lac 5 term, lack both term
    """
    print ('fragment type to remove : {0}'.format(choices))
    fragment_list = []
    prodigal_header = ""
    fragment_seq = ""
    fragment_status = ""
    try:
        with open(fragment_file_name, 'rt') as fragment_file:
            for line in fragment_file:
                if line.startswith('>'):
                    if prodigal_header == "":
                        prodigal_header = line.strip()
                        fragment_status = (
                            prodigal_header.split(";")[1].split("=")[1])
                    else:
                        #Test fragment and save it
                        #print (fragment_status)
                        if fragment_status not in choices:
                            fragment_list += [[prodigal_header , fragment_seq]]
                            #print (prodigal_header, fragment_seq)
                        #Reset fragment
                        fragment_seq = ""
                        prodigal_header = line.strip()
                        fragment_status = (
                             prodigal_header.split(";")[1].split("=")[1])
                else:
                    fragment_seq = fragment_seq + line.strip()
            if fragment_status not in choices:
                fragment_list += [[prodigal_header , fragment_seq]]
            assert(len(fragment_list) > 0)
    except IOError:
        sys.exit("Error cannot open {0}".format(fragment_file))
    except AssertionError:
        sys.exit("Error nothing selected for {0}".format(fragment_file))
    return fragment_list


def fill(text, width=80):
    """Split text"""
    return os.linesep.join(text[i:i+width] for i in xrange(0, len(text), width))


def write_result(fragment_list, results, output_file):
    """Write the result table
    """
    if not output_file:
        output_file = results + os.sep + "prodigal_filtered.txt"
    try:
        with open(output_file, "wt") as output:
            for fragment in fragment_list:
                output.write("{1}{0}{2}{0}".format(
                    os.linesep, fragment[0], fill(fragment[1])))
    except IOError:
        sys.exit("Error cannot open {0}".format(output_file))


def main():
    """Main program function
    """
    #Extract arguments
    args = config_parameters()
    #Start fragment selection
    fragment_list = select_fragment(args.fragment_file, args.choices)
    #Write Result
    write_result(fragment_list, args.results, args.output_file)


if __name__ == '__main__':
    main()
