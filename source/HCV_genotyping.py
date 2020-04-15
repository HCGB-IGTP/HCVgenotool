#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
HCV genotyping Tool.
"""
#------------------------------------------------------------------------------
# Script name: HCV_genotyping.py
# Instructions: python3 HCV_genotyping.py [OPTIONS] \
#               -i  <path/to/identity_df.csv> \
#               -s <path/to/align.sam> \
#               -o <path/to/output.txt>
#   Required arguments:
#       -i, --identity_df <path/to/identity_df.csv>: CSV file with pairwise
#           alignment identities in the form of a DataFrame with sequences as
#           rows and columns, computed for all the HCV genomes used as
#           references in the BWA alignment.
#       -s, --sam <path/to/align.sam>: sam alignemnt file or comma sepparated
#           list of sam files (one sample per file). Tested only for
#           single-read alignment comming from a pair-end overlapping sequencing
#           in which both pairs were joined using
#           "ea-utils/fastq-join v.1.01.759" (see
#           https://expressionanalysis.github.io/ea-utils/).
#       -o, --outfile <path/to/output.txt>: tab separated file with header, one
#           sample per column.
#   [OPTIONS]:
#       -h --help: prints help message.
#
# Description: starting from one or more BWA sam alignment file, coming from an
#   Illumina NGS HCV amplicon sequencing, this tool reports all genotypes, as
#   well as their read numbers, percentages and mean reference divergenece, for
#   all the processed samples. Only one sample for .sam file.
#
# HCVgenotool
# Copyright (C) 2018  David Piñeyro
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.
#
# Author: David Piñeyro <igtp.genbioinfo@gmail.com>
# License: GPL-3.0-or-later
# Version: 0.0.1
# Date: February 21st, 2018
#------------------------------------------------------------------------------
# IMPORTS
import argparse
import operator
import os
import pandas
import math
#------------------------------------------------------------------------------
# GLOBAL VARIABLES
SCRIPT_DESCRIPTION = """
HCV_genotyping.py
Version: 0.1
Author: David Piñeyro
"""
#------------------------------------------------------------------------------
# FUNCTIONS
def arguments():
    """This function uses argparse functionality to collect arguments."""
    parser = argparse.ArgumentParser(
                formatter_class=argparse.RawDescriptionHelpFormatter,
                description= SCRIPT_DESCRIPTION)
    parser.add_argument('-i', '--identity_df',
                        metavar='<path/to/identity_df.csv>', type=str,
                        required=True,
                        help="""CSV file with pairwise alignment identities in
                             the form of a DataFrame with sequences as rows and
                             columns, computed for all the HCV genomes used as a
                             references in the BWA alignment.""")
    # E.g.: sam1.sam,sam2.sam,sam3.sam  # Three samples.
    # E.g.: align.sam                   # One sample.
    parser.add_argument('-s', '--sam',
                        metavar='<path/to/align.sam>', type=str, required=True,
                        help="""sam alignemnt file or comma sepparated
                                list of sam files (one sample per file).
                                Tested only for single-read alignment comming
                                from a pair-end overlapping sequencing in which
                                both pairs were joined using
                                "ea-utils/fastq-join v.1.01.759" (see
                                https://expressionanalysis.github.io/ea-utils/).
                                """)
    parser.add_argument('-o', '--outfile',
                        metavar='<path/to/output.txt>', type=str, required=True,
                        help="""tab separated file with header, one sample per
                                column.""")

    args = parser.parse_args()
    print('Input sam file/s:', args.sam)
    print('Reference genomes identity Data Frame:', args.identity_df)
    return args
def aligned(l):
    """This function takes a sam alignment line as a list and checks whether
    the alignment line is suitable to create an Alignment object. Requeriments:
        - The sequence should be aligned.
        - The alignment should be a primary alignment.
        - Read should not fails platform/vendor quality checks.
        - Read is not PCR or optical duplicate.
        - Alignment is not supplementary alignment.
        - If more than one optimal alignment exists, all of the optimally
          aligned references should belong to the same HCV subtype.

    Arguments:
        l[list]: a list, comming from a sam alignment line, splitted by '\t'.
    Output:
        [True/False]: True if alignment is valid, False otherwise.
    """
    bit_flag = format(int(l[1]), '012b')  # Returns string of 12 0s and 1s.
    read_unmapped = bit_flag[-3] == '1'
    not_primary_alignment = bit_flag[-9] == '1'
    quality_failed = bit_flag[-10] == '1'
    duplicate = bit_flag[-11] == '1'
    suppl_align = bit_flag[-12] == '1'
    # Discarding reads.
    if read_unmapped or \
       not_primary_alignment or \
       quality_failed or \
       duplicate or \
       suppl_align:
        return False
    else:
        # Obtaining optimal and suboptimal alignment scores.
        optimal = ''
        suboptimal = ''
        # This loop assumes that AS and XS fields are present. At least AS
        # should be present.
        for f in l:
            if f[:5] == 'AS:i:':
                optimal = f.split(':')[-1]
            if f[:5] == 'XS:i:':
                suboptimal = f.split(':')[-1]
        if optimal == suboptimal:
            # Checking if the aligned HCV sequence is from the same subtype.
            # HCV nomenclature is NX_aaa, where N is type, X is subtype and
            # aaa is the clone name. If NX is the same for optimal and
            # suboptimal alignment, alignment is considered valid.
            optimal_subtype = l[2].split('_')[0]
            suboptimal_subtype = ''
            for f in l:
                if f[:5] == 'XA:Z:':
                    subtype = f.split(':')[-1]
                    subtype = subtype.split(',')[0]
                    suboptimal_subtype = subtype.split('_')[0]
                    break
            if optimal_subtype == suboptimal_subtype:
                return True
            else:
                return False
        else:
            return True
def collect_aligns(sam):
    """This function collects alignments from a sam file and creates
    Alignment objects for each of them. Sam flag is read to discard non-primary
    alignments, if any.

    Arguments:
        sam[string]: sam alignment file name.
    Output:
        aligns[list]: list of Alignment objects.
    """
    sam_file = open(sam, 'r')
    aligns = []  # List to return with Alignment objects.
    discarded = 0  # Number of discarded lines of sam file (non-aligned, etc).
    for l in sam_file:
        if l[0] != '@':  # If is not a comment line.
            l = l.rstrip().split('\t')
            # Discarding not aligned and not primary alignment.
            if aligned(l):
                align = Alignment(l)
                aligns.append(align)
            else:
                discarded += 1
    sam_file.close()
    # Printing some statistics.
    print('Collecting alignments from file {}'.format(sam))
    print('\tAlignments collected:', len(aligns))
    print('\tDiscarded alignments:', discarded)
    return aligns
def cigar_to_end(start, cigar):
    """This function takes a start position and a cigar string and calculates
    the real end position of the alignment.

    Arguments:
        start [int]: alignment start (1-based, leftmost) position.
        cigar [str]: CIGAR string.
    Output:
        end [int]: alignment end position.
    """
    to_add = 0  # variable to add to start in order to get end.
    # List of CIGAR operations that consumes reference.
    consumes_ref = ['M', 'D', 'N', '=', 'X']
    prev_int = ''  # String to store numeric values.
    for c in cigar:
        if c.isnumeric():
            prev_int += c
        else:
            if c in consumes_ref:
                to_add += int(prev_int)
            prev_int = ''
    end = start + to_add
    return end
def calculate_genotypes(aligns):
    """This function takes all Alignments from a sample and creates all the
    possible Genotype objects for this sample.

    Arguments:
        aligns [list]: list of Alignment objects, coming from a single sample.
    Output:
        genotypes [list]: list of Genotype objects.
    """
    genotypes = []
    for a in aligns:
        # Alignment relevant information.
        reference = a.ref
        start_pos = a.pos
        end_pos = cigar_to_end(start_pos, a.cigar)
        nm = a.edit_dist
        # Collect genotypes stored.
        stored_genotypes = [g.reference for g in genotypes]
        # Genotype already present in genotypes list. Just update.
        if reference in stored_genotypes:
            stored = False  # To control alignment was stored.
            for g in genotypes:
                if g.reference == reference:
                    # Check alignment positions. If it's not identical, create
                    # a new genotype object with the same reference.
                    if g.start_pos == start_pos and g.end_pos == end_pos:
                        # Update an existing genotype.
                        g.count += 1
                        if nm in g.nm_distribution:
                            g.nm_distribution[nm] += 1
                        else:
                            g.nm_distribution[nm] = 1
                        stored = True
                        break
            if not stored:
                # Reference was already present but with different positions.
                # Create a new Genotype object.
                genotype = Genotype(reference, 1, start_pos, end_pos, {nm: 1})
                genotypes.append(genotype)
        else:  # Reference was not present in stored_genotypes.
            genotype = Genotype(reference, 1, start_pos, end_pos, {nm: 1})
            genotypes.append(genotype)
    # Printing some statistics.
    # Calculate the number of seqs and reference of the main genotype.
    n_seqs, main_g = max([[g.count, g.reference] for g in genotypes])
    # Percentage of sequences.
    percent = round(n_seqs/len(aligns) * 100, 2)
    print('\tMain genotype: {}, with {} seqs ({}% of total).'.format(main_g,
                                                                     n_seqs,
                                                                     percent))
    print()
    # Computing identity distribution for each genotype.
    for g in genotypes:
        g.compute_identity_dist()
    return genotypes
def group_genotypes(g_list):
    """Group genotypes by their subtype.

    Arguments:
        g_list [list]: a list of Genotype objects.
    Output:
        g_groups [list]: a list of Genotype_group objects.
    """
    g_groups = []
    for g in g_list:
        # Create the first Genotype_group object.
        if len(g_groups) == 0:
            g_group = Genotype_group(g)
            g_groups.append(g_group)
        else:
            # Update existing Genotype_group.
            if g.reference.split('_')[0] in [gr.group for gr in g_groups]:
                i = 0
                while i < len(g_groups):
                    if g.reference.split('_')[0] == g_groups[i].group:
                        g_groups[i].update(g)
                        break
                    i += 1
            else:
                # Creating new Genotype_group.
                g_group = Genotype_group(g)
                g_groups.append(g_group)
    return g_groups
def output_genotypes(s, s_files, identity_df, outfile):
    """Writes the table with genotypes for all the samples.

    Arguments:
        s [list]: list of Genotype_group objects, one for each sample.
        s_files [tuple]: a tuple of strings with samples filenames.
        identity_df [string]: filename of the identity dataframe.
        outfile [string]: filename for the output file.
    Output:
        out [file]: a table containing all the samples genotyping information in
            .csv format.
    """
    # Reading identity_df and loading it into a DataFrame.
    id_df = pandas.read_csv(identity_df, index_col=0)
    # Creating the output file.
    out = open(outfile, 'w')
    # Writing legend:
    out.write('Sample: sample id.\n')
    out.write('Reads: number of filtered reads.\n')
    out.write('Percentage: percentage of total sample filtered reads.\n')
    out.write('HCV_subtype: closest subtype selected by sequence alignment.\n')
    out.write('Min_identity: minimal reference identity percentage.\n')
    out.write('Max_identity: maximum reference identity percentage.\n')
    out.write('Avg_identity: average reference identity percentage.\n')
    out.write('SD_identity: reference identity percentage standard deviation.\n')
    out.write('Mode_identity: most frequent reference identity percentage.\n')
    out.write('Min_subtype_identity: minimal reference subtype identity percentage.\n')
    out.write('Refs: reference sequences to which reads were found aligned.\n')
    # Writing column names.
    out.write('Sample,Reads,Percentage,HCV_subtype,Min_identity,Max_identity,Avg_identity,SD_identity,Mode_identity,Min_subtype_identity,HCV_refs\n')
    # Extracting sample names (only works on linux).
    sample_names = [os.path.basename(f).split('.')[0]
                    for f in s_files.split(',')]
    # s and s_files are in the same order.
    i = 0
    while i < len(sample_names):
        name = sample_names[i]
        genotypes = s[i]
        n_reads = sum([g.count for g in genotypes])
        # Writing one line per genotype found.
        for g in genotypes:
            out.write('{},'.format(name))
            out.write('{},'.format(g.count))
            out.write('{},'.format(round(g.count/n_reads*100, 2)))
            out.write('{},'.format(g.group))
            out.write('{},'.format(g.min_identity()))
            out.write('{},'.format(g.max_identity()))
            out.write('{},'.format(round(g.mean_identity(), 2)))
            out.write('{},'.format(round(g.sd_identity(), 2)))
            out.write('{},'.format(round(g.mode_identity(), 2)))
            out.write('{},'.format(id_df.loc[g.group, g.group]))
            out.write('{}'.format('; '.join(g.reference_set)))
            out.write('\n')
        i += 1
    out.close()
#------------------------------------------------------------------------------
# CLASSES
class Alignment(object):
    """This class is designed to store alignment lines from sam file.

    Attributes:
        seq_name [str]: aligned sequence name.
        flag [int]: sam flag (as it appears in sam file).
        ref [str]: reference sequence mapped.
        pos [int]: 1-based leftmost POSition/coordinate of clipped sequence.
        mapq [int]: MAPping Quality (Phred-scaled).
        cigar [str]: extended CIGAR string.
        seq [str]: query SEQuence on the same strand as the reference.
        qual [str]: query QUALity (ASCII-33 gives the Phred base quality) of seq
        edit_dist [int]: number of changes compared to reference.
        align_score [int]: alignment score.
    Methods:
        __init__(self, l):
            l[list]: a list, comming from a sam alignment line,splitted by '\t'.
    """
    def __init__(self, l):
        self.seq_name = l[0]
        self.flag = int(l[1])
        self.ref = l[2]
        self.pos = int(l[3])
        self.mapq = int(l[4])
        self.cigar = l[5]
        self.seq = l[9]
        self.qual = l[10]
        self.edit_dist = [int(f.split(':')[-1]) for f in l if f[:5] == 'NM:i:'][0]
        self.align_score = [int(f.split(':')[-1]) for f in l if f[:5] == 'AS:i:'][0]
class Genotype(object):
    """This class stores genotypes information.

    Attributes:
        reference [str]: name of the reference sequence.
        count [int]: number of sequences (reads) with this genotype.
        start_pos [int]: 1-based leftmost POSition/coordinate of aligned seq.
        end_pos [int]: 1-based rightmost POSition/coordinate of aligned seq.
        alignment_length [int]: end_pos - start_pos.
        nm_distribution [dict]:
            key [int]: edit distance, it is: number of transformations
                (insertions, deletions or substitutions) required to transform
                the aligned read into the reference.
            value [int]: number of sequences with this edit distance.
            E.g: {0: 10, 1: 100, 2: 345, 5: 445}
        identity_distribution [dict]: attribute computed after
            Genotype.compute_identity_dist() method execution. It is the
            reference-read pairwise identity percentage distribution. Calculated
            from nm_distribution, being:
                identity = (1 - (edit distance)/alignment_length)*100
            key [float]: identity percentage, it is: percentage of matched
                positions in the alignment.
            value [int]: number of sequences with this identity percentage.
            E.g: {100.00: 10, 99.00: 100, 98.00: 345, 95.00: 445}
    Methods:
        compute_identity_dist(self): computes identity distribution from
            nm_distribution (edit distances distribution).
    """
    def __init__(self, reference, count, start_pos, end_pos, nm_distribution):
        self.reference = reference
        self.count = count
        self.start_pos = start_pos
        self.end_pos = end_pos
        self.alignment_length = end_pos - start_pos
        self.nm_distribution = nm_distribution
    def compute_identity_dist(self):
        """Computes identity distribution from nm_distribution (edit distances
        distribution).
        """
        self.identity_distribution = {
            round((1-(k/(self.end_pos - self.start_pos)))*100, 2):v
            for k, v in self.nm_distribution.items()}
class Genotype_group(object):
    """Creates a genotype group object to store all the genotypes of the same
    HCV subtype found in a sample. The object is initialized with a single
    Genotype object. Then, each Genotype from the same HCV subtype (group)
    updates the Genotype_group object, adding its information.

    Attributes:
        group [string]: each HCV reference sequence is named following this
            rule: Nx_ZZZZZZ being N the type, x the subtype and ZZZZZZ the
            id of each particular reference sequence. Thus, selecting the
            part before '_' char gives the type and subtype of HCV sequence.
            The attribute group stores the type and subtype of a given
            Genotype object.
        reference_set [set]: a set (unordered set of unique elements) with the
            refernce sequence ids of the Genotypes grouped by Genotype_group.
        count [int]: number of reads found for this Genotype_group.
        al_length_set [set]: set with the different alignment lengths found
            for this Genotype_group.
        nm_distribution [dict]:
            key [int]: edit distance, it is: number of transformations
                (insertions, deletions or substitutions) required to transform
                the aligned read into the reference.
            value [int]: number of sequences with this edit distance.
            E.g: {0: 10, 1: 100, 2: 345, 5: 445}
        identity_distribution [dict]: attribute computed after
            Genotype.compute_identity_dist() method execution. It is the
            reference-read pairwise identity percentage distribution. Calculated
            from nm_distribution, being:
                identity = (1 - (edit distance)/alignment_length)*100
            key [float]: identity percentage, it is: percentage of matched
                positions in the alignment.
            value [int]: number of sequences with this identity percentage.
            E.g: {100.00: 10, 99.00: 100, 98.00: 345, 95.00: 445}
    Methods:
        update(self, Genotype): adds information from another Genotype object
            to the existing Genotype_group object.
        min_nm(self): returns the smallest key [int] from self.nm_distribution.
        min_identity(self): returns the smallest key [float] from
            self.identity_distribution.
        max_nm(self): returns the largest key [int] from self.nm_distribution.
        max_identity(self): returns the largest key [float] from
            self.identity_distribution.
        mode_nm(self): returns the self.nm_distribution key [int] with the
            largest value.
        mode_identity(self): returns the self.identity_distribution key [float]
            with the largest value.
        mean_nm(self): returns the mean of all the self.nm_distribution keys.
        mean_identity(self): returns the mean of all the
            self.identity_distribution keys.
        sd_nm(self): returns the corrected sample standard deviation (using
            n - 1 as denominator) of all the self.nm_distribution keys.
        sd_identity(self): returns the corrected sample standard deviation
            (using n - 1 as denominator) of all the self.identity_distribution
            keys.
    """
    def __init__(self, genotype):
        self.group = genotype.reference.split('_')[0]
        self.reference_set = set([genotype.reference])
        self.count = genotype.count
        self.al_length_set = set([genotype.alignment_length])
        self.nm_distribution = genotype.nm_distribution
        self.identity_distribution = genotype.identity_distribution
    def update(self, genotype):
        """Adds information of an additional Genotype object to the existing
        Genotype_group. It updates the following Genotype_group attributes:
            self.reference_set
            self.count
            self.al_length_set
            self.nm_distribution
            self.identity_distribution
        Arguments:
            genotype [Genotype]: the Genotype object with the information to add.
        """
        self.reference_set.add(genotype.reference)
        self.count += genotype.count
        self.al_length_set.add(genotype.alignment_length)
        # Updating nm_distribution:
        for k, v in genotype.nm_distribution.items():
            if k in self.nm_distribution:
                self.nm_distribution[k] += v
            else:
                self.nm_distribution[k] = v
        # Updating identity_distribution:
        for k, v in genotype.identity_distribution.items():
            if k in self.identity_distribution:
                self.identity_distribution[k] += v
            else:
                self.identity_distribution[k] = v
    def min_nm(self):
        """Returns the smallest key [int] from self.nm_distribution."""
        return min(self.nm_distribution)
    def min_identity(self):
        """Returns the smallest key [float] from self.identity_distribution."""
        return min(self.identity_distribution)
    def max_nm(self):
        """Returns the largest key [int] from self.nm_distribution."""
        return max(self.nm_distribution)
    def max_identity(self):
        """Returns the largest key [float] from self.identity_distribution."""
        return max(self.identity_distribution)
    def mode_nm(self):
        """Returns the self.nm_distribution key [int] with the largest value,
        i.e. the most frequent key."""
        # Using operator.itemgetter to use dict value as the key for max().
        # Note that if there are more than one key with the largest value, this
        # method only returns one of these keys.
        return max(self.nm_distribution.items(), key=operator.itemgetter(1))[0]
    def mode_identity(self):
        """Returns the self.identity_distribution key [float] with the largest
        value, i.e. the most frequent key.
        """
        return max(self.identity_distribution.items(),
                   key=operator.itemgetter(1))[0]
    def mean_nm(self):
        """Returns the mean of all the self.nm_distribution keys."""
        c = 0
        for k, v in self.nm_distribution.items():
            c += k*v
        return c / sum(self.nm_distribution.values())
    def mean_identity(self):
        """Returns the mean of all the self.identity_distribution keys."""
        c = 0
        for k, v in self.identity_distribution.items():
            c += k*v
        return c / sum(self.identity_distribution.values())
    def sd_nm(self):
        """Returns the corrected sample standard deviation (using n - 1 as
        denominator) of all the self.nm_distribution keys."""
        # List to store all the distribution values.
        l = []
        for k, v in self.nm_distribution.items():
            l_temp = [k for i in range(v)]
            l += l_temp
        # Return 'NA' when len(l) == 1.
        if len(l) == 1:
            return 0
        mean = sum(l) / len(l)
        sd = math.sqrt((1/(len(l) - 1)) * sum([(x - mean)**2 for x in l]))
        return sd
    def sd_identity(self):
        """Returns the corrected sample standard deviation (using n - 1 as
        denominator) of all the self.identity_distribution keys."""
        # List to store all the distribution values.
        l = []
        for k, v in self.identity_distribution.items():
            l_temp = [k for i in range(v)]
            l += l_temp
        # Return 'NA' when len(l) == 1.
        if len(l) == 1:
            return 0
        mean = sum(l) / len(l)
        sd = math.sqrt((1/(len(l) - 1)) * sum([(x - mean)**2 for x in l]))
        return sd
#------------------------------------------------------------------------------
# MAIN PROGRAM
def main_program():
    """
    This is the main program.
    """
    # args.reference, args.sam, args.outfile
    args = arguments()
    # Process each sample separately. Store results as a list of
    # Genotype_group objects.
    samples = []
    for sam in args.sam.split(','):
        # Create a list Alignment objects.
        aligns = collect_aligns(sam)
        # Create a list of Genotype objects from aligns.
        genotypes = calculate_genotypes(aligns)
        # Group genotypes by subtype and create a list of Genotype_group objects
        # for this sample.
        g_groups = group_genotypes(genotypes)
        samples.append(g_groups)
    # Write the output file.
    output_genotypes(samples, args.sam, args.identity_df, args.outfile)
#------------------------------------------------------------------------------
# Conditional to run the script
if __name__ == '__main__':
    main_program()
