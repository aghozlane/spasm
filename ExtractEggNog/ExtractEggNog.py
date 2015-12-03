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


"""Extract EggNogg annotation."""


from __future__ import print_function
import argparse
import os
import sys
import csv
import re
import glob
import multiprocessing as mp
import bisect
import subprocess


__author__ = "Amine Ghozlane"
__copyright__ = "Copyright 2014, INRA"
__license__ = "GPL"
__version__ = "1.0.0"
__maintainer__ = "Amine Ghozlane"
__email__ = "amine.ghozlane@jouy.inra.fr"
__status__ = "Developpement"


class FullPaths(argparse.Action):
    """Expand user- and relative-paths"""
    def __call__(self, parser, namespace, values, option_string=None):
        # If false then its a list
        if not isinstance(values, basestring):
            out = []
            for val in values:
                if os.path.isfile(val):
                    out += [os.path.abspath(os.path.expanduser(val))]
                elif os.path.isdir(val):
                    out += [os.path.abspath(os.path.expanduser(val)) + os.sep]
                else:
                    out += [val]
            setattr(namespace, self.dest, out)
        # Value is a string
        else:
            if os.path.isfile(values):
                setattr(namespace, self.dest,
                        os.path.abspath(os.path.expanduser(values)))
            elif os.path.isdir(values):
                setattr(namespace, self.dest,
                        os.path.abspath(os.path.expanduser(values)) + os.sep)


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
                                     "{0} -h (see also "
                                     "ftp://eggnog.embl.de/eggNOG/3.0/)"
                                     .format(sys.argv[0]))
    parser.add_argument('-b', dest='blast_result_file', type=isfile, 
                        required=True, action=FullPaths,
                        help='Blast result against eggnog '
                        ' sequence datablase (type outfmt 6 with blast+)')
    parser.add_argument('-il', dest='identity_limit', type=float, default=0.0,
                        help='Identity percent threshold on the reference')
    parser.add_argument('-cl', dest='coverage_limit', type=float, default=0.0,
                        help='Coverage percent threshold on the reference')
    parser.add_argument('-nb', dest='nbest', type=int, default=0,
                        help='Number of best selected (default:None, '
                        'based on the best value : identity + coverage)')
    parser.add_argument('-m', dest='nog_members', type=isdir, required=True,
                        action=FullPaths,
                        help='Directory with all NOG members')
    parser.add_argument('-n', dest='nog_description', type=isdir,
                        required=True, action=FullPaths,
                        help='Directory with all NOG description')
    parser.add_argument('-a', dest='nog_funccat_desc', type=isfile,
                        required=True, action=FullPaths,
                        help='Funccat description file '
                        '(eggnogv3.funccats.txt)')
    parser.add_argument('-f', dest='nog_funccat', type=isdir, action=FullPaths,
                        required=True, help='Directory with all NOG funccat')
    parser.add_argument('-o', dest='output_file', type=str,
                        help='Output file')
    parser.add_argument('-r', dest='results', type=isdir, action=FullPaths,
                        default=os.curdir + os.sep, help='Path to result '
                        'directory.')
    return parser.parse_args()


def getfiles(directory, type_file):
    """Get all files of one type in a given directory
    """
    return glob.glob('{0}{1}*.{2}.txt'.format(directory, os.sep, type_file))


def get_description(nog_description_files):
    """Get the annotation for a COG/NOG/KOG id
    """
    nog_desc = {}
    try:
        for nog_description_file in nog_description_files:
            with open(nog_description_file, "rt") as nog_description:
                nog_description_reader = csv.reader(nog_description,
                                                    delimiter='\t')
                # Pass header
                nog_description_reader.next()
                for line in nog_description_reader:
                    if len(line) == 2:
                        nog_desc[line[0]] = line[1]
                    else:
                        print("Strange, more than two lines :{0}{1}"
                              .format(os.linesep, line), file=sys.stderr)
                assert(len(nog_desc) > 0)
    except IOError:
        sys.exit("Error cannot open {0}".format(nog_description_file))
    except AssertionError:
        sys.exit("Error nothing read from {0}".format(nog_description_file))
    return nog_desc


def get_members(nog_members_files):
    """Get all the sequences id belonging to a COG/NOG/KOG id
    """
    nog_mem = {}
    try:
        for nog_members_file in nog_members_files:
            with open(nog_members_file, "rt") as nog_members:
                nog_members_reader = csv.reader(nog_members, delimiter='\t')
                # Pass header
                nog_members_reader.next()
                for line in nog_members_reader:
                    if len(line) == 4:
                        if line[1] in nog_mem:
                            nog_mem[line[1]] += [line[0]]
                        else:
                            nog_mem[line[1]] = [line[0]]
                assert(len(nog_mem) > 0)
    except IOError:
        sys.exit("Error cannot open {0}".format(nog_members_file))
    except AssertionError:
        sys.exit("Error nothing read from {0}".format(nog_members_file))
    return nog_mem


def get_funccat_description(nog_funccat_desc_file):
    """Load Functionnal category description
    """
    nog_fun_desc = {}
    desc_regex = re.compile(r"\s+\[([A-Z])\]\s+(.+)")
    try:
        with open(nog_funccat_desc_file, "rt") as nog_funccat_desc:
            for line in nog_funccat_desc:
                desc_match = desc_regex.match(line)
                if desc_match:
                    nog_fun_desc[desc_match.group(1)] = desc_match.group(2)
            assert(len(nog_fun_desc) > 0)
    except IOError:
        sys.exit("Error cannot open {0}".format(nog_funccat_desc_file))
    except AssertionError:
        sys.exit("Error nothing read from {0}".format(nog_funccat_desc_file))
    return nog_fun_desc


def get_funccat(nog_funccat_files):
    """Get the funcid for each COG/NOG/KOG id
    """
    nog_func = {}
    try:
        for nog_funccat_file in nog_funccat_files:
            with open(nog_funccat_file, "rt") as nog_funccat:
                nog_funccat_reader = csv.reader(nog_funccat, delimiter='\t')
                # Pass header
                nog_funccat_reader.next()
                for line in nog_funccat_reader:
                    if len(line) == 2:
                        nog_func[line[0]] = line[1]
                assert(len(nog_func) > 0)
    except IOError:
        sys.exit("Error cannot open {0}".format(nog_funccat_file))
    except AssertionError:
        sys.exit("Error nothing read from {0}".format(nog_funccat_file))
    return nog_func


def get_blast_result(blast_result_file):
    """Parse and directly annotate the sequence
    """
    blast_dict = {}
    try:
        with open(blast_result_file, "rt") as blast_result:
            blast_result_reader = csv.reader(blast_result, delimiter='\t')
            for line in blast_result_reader:
                # Get target in database, identity, coverage and evalue
                
                ##TO CHECK
                if line[0] in blast_dict:
                    blast_dict[line[0]].append([line[1], float(line[10]),
                                                float(line[11]), float(line[12])])
                else:
                    blast_dict[line[0]] = [[line[1], float(line[10]),
                                            float(line[11]), float(line[12])]]
            assert(len(blast_dict) > 0)
    except IOError:
        sys.exit("Error cannot open {0}".format(blast_result_file))
    except AssertionError:
        sys.exit("Error nothing read from {0}".format(blast_result_file))
    return blast_dict



def set_annotation(blast_dict, nog_desc, nog_mem, nog_funccat,
                   nog_funccat_desc):
    """Get annotation
    """
    annotation_result = {}
    try:
        for seq in blast_dict:
            annotation_result[seq] = []
            for blast in blast_dict[seq]:
                ref = blast[0]
                # UGLY
                identity = blast[1]
                coverage = blast[2]
                evalue = blast[3]
                if ref in nog_mem:
                    cog_list = nog_mem[ref]
                    for cog_id in cog_list:
                        if cog_id in nog_funccat and cog_id in nog_desc:
                            annotation_result[seq] += [[ref, cog_id,
                                                "{0}".format(", ".join(
                                                [nog_funccat_desc[ident]
                                                 if ident in nog_funccat_desc
                                                 else "Function unknown"
                                                for ident in nog_funccat[cog_id]]
                                                )), nog_desc[cog_id], evalue,
                                                identity, coverage]]
                        elif cog_id in nog_desc:
                            annotation_result[seq] += [[ref , cog_id, "NA",
                                                    nog_desc[cog_id], evalue]]
                        elif cog_id in nog_funccat:
                            annotation_result[seq] += [[ref, cog_id,
                                                "{0}".format(", ".join(
                                                [nog_funccat_desc[ident]
                                                 if ident in nog_funccat_desc
                                                 else "Function unknown"
                                                for ident in nog_funccat[cog_id]]
                                                )), "NA", evalue, identity,
                                                coverage]]
                        else:
                            annotation_result[seq] += [[ref , cog_id, "NA", "NA",
                                                        evalue, identity,
                                                        coverage]]
                else:
                    annotation_result[seq] += [[ref, "NA", "NA", "NA", evalue,
                                                identity, coverage]]
        assert(len(annotation_result) > 0)
    except AssertionError:
        sys.exit("No annotation found at all !")
    return annotation_result


def write_result(annotation_result, identity_limit, coverage_limit,
                 nbest, results, output_file):
    """Write the result table
    """
    data = []
    if not output_file:
        output_file = results + os.sep + "egg_nog_annotation.txt"
    try:
        with open(output_file, "wt") as output:
            output_writer = csv.writer(output, delimiter='\t')
            output_writer.writerow(["SequenceName", "ReferenceGene",
                                    "COG/NOG/KOG", "Functional category",
                                    "Annotation", "evalue", "Identity",
                                    "Coverage"])
            #output_writer.writerows(annotation_result)
            for seq in annotation_result:
                for ref_set in annotation_result[seq]:
                    if ref_set[-1] >= coverage_limit:
                        data += [ref_set]
                # Sort depending on the identity and coverage
                data.sort(key=lambda x: x[0][-2] + x[0][-1])
                if nbest > 0:
                    short_set = data[0: nbest]
                else:
                    short_set = data
                for line in short_set:
                    output_writer.writerow([seq] + line)
                del(data)
                data = []
    except IOError:
        sys.exit("Error cannot open {0}".format(output_file))


def write_result_check(annotation_result, result_alignment, identity_limit,
                       coverage_limit, nbest, results, output_file):
    """Write the result table with alignment data
    """
    data = []
    if not output_file:
        output_file = results + os.sep + "egg_nog_annotation.txt"
    try:
        with open(output_file, "wt") as output:
            output_writer = csv.writer(output, delimiter='\t')
            output_writer.writerow(["SequenceName", "ReferenceProtein",
                                    "COG/NOG/KOG", "Functional category",
                                    "Annotation", "evalue", "Identity",
                                    "Similarity", "Coverage"])
            #for annot in annotation_result:
            for seq in annotation_result:
                for ref_set in annotation_result[seq]:
                    ref = ref_set[0]
                    ali = result_alignment[tuple([seq, ref])]
                    if ali[1] >= identity_limit and ali[2] >= coverage_limit:
                        data += [ref_set[0:6] + ali]
                # Sort depending on the identity and coverage
                data.sort(key=lambda x: x[0][-3] + x[0][-1], reverse=True)
                if nbest > 0:
                    short_set = data[0: nbest]
                else:
                    short_set = data
                for line in short_set:
                    output_writer.writerow([seq] + line)
                del(data)
                data = []
    except IOError:
        sys.exit("Error cannot open {0}".format(output_file))


#==============================================================
# Main program
#==============================================================
def main():
    """
    Main program function
    """
    result_alignment = None
    # Get arguments
    args = get_arguments()
    # Parse description
    print("Load description")
    nog_desc = get_description(getfiles(args.nog_description, "description"))
    # Parse members
    print("Load members")
    nog_mem = get_members(getfiles(args.nog_members, "members"))
    # Parse funccats description
    print("Load desc func")
    nog_funccat_desc = get_funccat_description(args.nog_funccat_desc)
    # Parse funccat
    print("Load func")
    nog_funccat = get_funccat(getfiles(args.nog_funccat, "funccat"))
    # Parse blast result
    print("Load blast")
    blast_dict = get_blast_result(args.blast_result_file)
    # Set annotation
    print("Set annotation")
    annotation_result = set_annotation(blast_dict, nog_desc, nog_mem,
                                       nog_funccat, nog_funccat_desc)
    # Write result
    if result_alignment:
        write_result_check(annotation_result, result_alignment, args.identity_limit,
                           args.coverage_limit, args.nbest, args.results,
                           args.output_file)
    else:
        write_result(annotation_result, args.identity_limit, args.coverage_limit,
                     args.nbest, args.results, args.output_file)


if __name__ == '__main__':
    main()
