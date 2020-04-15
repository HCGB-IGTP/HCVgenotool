#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Pairwaise alignment tool for HCV genome fragments.
"""
#------------------------------------------------------------------------------
# Script name: HCV_amplicon_identity.py
# Instructions: python3 HCV_amplicon_identity.py [OPTIONS] \
#               -p <path/to/primers.fasta> \
#               -r <path/to/reference.fasta> \
#               -o <path/to/basename_>
#   Required arguments:
#       -p, --primers, <path/to/primers.fasta>: fasta file with the two primers
#           used to amplify genome fragments (make amplicons) before sequencing.
#       -r, --reference <path/to/reference.fasta>: fasta file with all the HCV
#           reference genomes.
#       -o, --outfile <path/to/basename_>: output file base name for all the
#           output files:
#               basename_amplicons.fasta,
#               basename_alignments.txt,
#               basename_identity_df_all.csv,
#               basename_identity_df_per_group.csv.
#   [OPTIONS]:
#       -h --help: prints help message.
#
# Description: taking a pair of primer sequences and all HCV genomes, pcr
#   amplicons are obtained and aligned. As a result, several files are created:
#       basename_amplicons.fasta: amplicons determined as if a pcr was performed
#           with the given primers for each HCV genome.
#       basename_alignments.txt: all pairwise global alignments of the
#           amplicons.
#       basename_identity_df_all.csv: CSV table with all the amplicons pairwise
#           identities.
#       basename_identity_df_per_group.csv: CSV table with per group MINIMAL
#           pairwise identity.
# note: "identity" is defined as number of matches / alignment length * 100.
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
import time
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import pairwise2
#------------------------------------------------------------------------------
# GLOBAL VARIABLES
SCRIPT_DESCRIPTION = """
HCV_amplicon_aligner.py
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
    parser.add_argument('-p', '--primers',
                        metavar='<path/to/primers.fasta>', type=str,
                        required=True,
                        help="""fasta file with the two primers used to
                             amplify genome fragments (make amplicons) before
                             sequencing.""")
    # E.g.: sam1.sam,sam2.sam,sam3.sam  # Three samples.
    # E.g.: align.sam                   # One sample.
    parser.add_argument('-r', '--reference',
                        metavar='<path/to/reference.fasta>', type=str,
                        required=True,
                        help="""fasta file with all the HCV reference genomes.
                                """)
    parser.add_argument('-o', '--outfile',
                        metavar='<path/to/basename_>', type=str, required=True,
                        help="""output file base name for all the output files:
                             basename_amplicons.fasta,
                             basename_alignments.txt,
                             basename_identity_df_all.csv,
                             basename_identity_df_per_group.csv.""")

    args = parser.parse_args()
    print('Reference genomes from file:', args.reference)
    print('Primers used for amplicons:', args.primers)
    return args
def EDNAfull_matrix():
    """EDNAfull matrix used by EMBOSS water and needle algorithms for nucleotide
    alignments. For protein alignments use blosum62 from biopython.

    This matrix was created by Todd Lowe   12/10/92

    Uses ambiguous nucleotide codes, probabilities rounded to nearest integer

    Lowest score = -4, Highest score = 5

    Available at: <ftp://ftp.ncbi.nih.gov/blast/matrices/NUC.4.4>
    """
    # Matrix nucleotides.
    bases = ["A", "T", "G", "C", "S", "W", "R", "Y",
             "K", "M", "B", "V", "H", "D", "N"]
    # Matrix scores (it's a symetric matrix).
    scores = [5, -4, -4, -4, -4,  1,  1, -4, -4,  1, -4, -1, -1, -1, -2,
             -4,  5, -4, -4, -4,  1, -4,  1,  1, -4, -1, -4, -1, -1, -2,
             -4, -4,  5, -4,  1, -4,  1, -4,  1, -4, -1, -1, -4, -1, -2,
             -4, -4, -4,  5,  1, -4, -4,  1, -4,  1, -1, -1, -1, -4, -2,
             -4, -4,  1,  1, -1, -4, -2, -2, -2, -2, -1, -1, -3, -3, -1,
              1,  1, -4, -4, -4, -1, -2, -2, -2, -2, -3, -3, -1, -1, -1,
              1, -4,  1, -4, -2, -2, -1, -4, -2, -2, -3, -1, -3, -1, -1,
             -4,  1, -4,  1, -2, -2, -4, -1, -2, -2, -1, -3, -1, -3, -1,
             -4,  1,  1, -4, -2, -2, -2, -2, -1, -4, -1, -3, -3, -1, -1,
              1, -4, -4,  1, -2, -2, -2, -2, -4, -1, -3, -1, -1, -3, -1,
             -4, -1, -1, -1, -1, -3, -3, -1, -1, -3, -1, -2, -2, -2, -1,
             -1, -4, -1, -1, -1, -3, -1, -3, -3, -1, -2, -1, -2, -2, -1,
             -1, -1, -4, -1, -3, -1, -3, -1, -3, -1, -2, -2, -1, -2, -1,
             -1, -1, -1, -4, -3, -1, -1, -3, -1, -3, -2, -2, -2, -1, -1,
             -2, -2, -2, -2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1]
    matrix = {}
    counter = 0
    for a in bases:
        for b in bases:
            matrix[(a, b)] = scores[counter]
            counter += 1
    return matrix
def virtual_pcr(genomes, primers):
    """Outputs the amplicon resulted using the given primers for each genome.

    Arguments:
        genomes [list]: list of Bio.SeqRecord objects with reference DNA
            sequences collected from a fasta file using SeqIO.parse().
        primers [list]: list of Bio.SeqRecord objects with DNA primer sequences
            collected from a fasta file using SeqIO.parse(). It must contain
            only 2 sequences.
    Output:
        amplicons [list]: list of Bio.SeqRecord objects with DNA amplicons
            resulted from a virtual pcr using primers with reference genomes.
    """
    # Check for only two primers.
    if len(primers) != 2:
        raise Exception('{} primer sequences supplied, only 2 allowed.'.format(
                        len(primers)))
    # Computing all the two primers alignments in both, forward and reverse
    # complement.
    primers_al = []
    for p in primers:
        # Smith-Waterman local alignment.
        # Primers are supposed to be in a 5' to 3' orientation. As 3' is the
        # most important part for a PCR and taking into account that many
        # primers come with flags at 5', to avoid noisy alignment, only the
        # last 30 bases are taken (if len(primer) <= 30, all primer is taken).
        # Primer in forward orientation
        fwd = [pairwise2.align.localds(p.seq[-30:], g.seq, EDNAfull_matrix(),
               -10, -0.5, one_alignment_only=True) for g in genomes]
        # Primer in reverse complement orientation.
        p_rev = p.reverse_complement(id=p.id+'_rc')
        rev = [pairwise2.align.localds(p_rev.seq[:30], g.seq, EDNAfull_matrix(),
               -10, -0.5, one_alignment_only=True) for g in genomes]
        # Convert each alignment result in a Align_formatted object.
        fwd_formatted = [Align_formatted(al) for al in fwd]
        rev_formatted = [Align_formatted(al) for al in rev]
        # For each genome, select the best scoring align (fwd or rev).
        best = []
        for g in range(len(genomes)):
            if fwd_formatted[g].score >= rev_formatted[g].score:
                best.append(fwd_formatted[g])
            else:
                best.append(rev_formatted[g])
        # For each primer, append the best alignments.
        primers_al.append(best)
    # Compute slicing positions for genomes.
    slicings = []
    for g in range(len(genomes)):
        gaps_1 = primers_al[0][g].gaps('b')
        gaps_2 = primers_al[1][g].gaps('b')
        p1_start, p1_end = primers_al[0][g].al_start, primers_al[0][g].al_end
        p2_start, p2_end = primers_al[1][g].al_start, primers_al[1][g].al_end
        # Gaps up to relevant position should be subtracted.
        p1_start -= len([gap for gap in gaps_1 if gap < p1_start])
        p1_end -= len([gap for gap in gaps_1 if gap < p1_end])
        p2_start -= len([gap for gap in gaps_2 if gap < p2_start])
        p2_end -= len([gap for gap in gaps_2 if gap < p2_end])
        # The smallest end will be the start position of the amplicon. The
        # bigest start will be the end position of the amplicon.
        al_start = min(p1_end, p2_end)
        al_end = max(p1_start, p2_start)
        s = al_start, al_end
        slicings.append(s)
    # Slicing seqs and creating SeqRecord objects
    amplicons = []
    for i, s in enumerate(slicings):
        amplicon = SeqRecord(seq = genomes[i].seq[s[0]:s[1]],
                             id = genomes[i].id,
                             description = "amplicon")
        amplicons.append(amplicon)
    return amplicons
def nw_global_alignment(seqs):
    """Performs pairwise global alignments for each sequence in seqs against
    each other. It uses Bio.pairwise2.align.globalds() function.

    Arguments:
        seqs [list]: list of Bio.SeqRecord objects with DNA sequences.
    Output:
        alignments [list]: list of Align_formatted objects with global
            alignments results. It is a 2D list. First dimension are the
            sequences and the second dimension are the alignments of this
            sequence.
    """
    alignments = []
    i = 0  # index counter
    while i < len(seqs):
        seq_a = seqs[i].seq
        seq_a_id = seqs[i].id
        # Generate a list with the rest of the indexes, to iterate over.
        rest_i = [idx for idx in range(len(seqs)) if idx != i]
        als = [] # list to store alignments with this seq_a.
        for idx in rest_i:
            seq_b = seqs[idx].seq
            seq_b_id = seqs[idx].id
            al = pairwise2.align.globalds(seq_a, seq_b, EDNAfull_matrix(), -10,
                                          -0.5, one_alignment_only=True)
            # Formatting alignment.
            al_formatted = Align_formatted(al, seq_a_id = seq_a_id,
                                           seq_b_id = seq_b_id)
            als.append(al_formatted)
        alignments.append(als)
        i += 1
    return alignments
def write_pairwise_al(al, output_file):
    """Writes all pairwise alignments, with some alignment statistics and info.

    Parameters:
        al [list]: a 2D list (list of list) with all pairwise alignments
            (Align_formatted objects) for each sequence.
        output_file [string]: output file name.
    Output:
        output_f [file]: a file with all pairwise alignments and alignment
            information.
    """
    # Open output file to write.
    output_f = open(output_file, 'w')
    # Writing heading information.
    output_f.write('#######################################\n')
    output_f.write('# Needleman-Wunsch global alignments. #\n')
    output_f.write('#######################################\n')
    # Loop to iterate over each sequence.
    i = 0
    while i < len(al):
        seq_a_id = al[i][0].seq_a_id
        output_f.write('-' * 80)
        output_f.write('\n')
        output_f.write('Pairwise global alignments for sequence {}.\n'.format(
                       seq_a_id))
        # Loop to iterate over each alignment.
        for align in al[i]:
            seq_b_id = align.seq_b_id
            output_f.write('\n{} vs {} global alignment.\n'.format(seq_a_id,
                                                                   seq_b_id))
            # Writing out alignment statistics.
            output_f.write(
                'Alignment score: {}.\n'.format(align.score))
            output_f.write(
                'Alignment length: {}.\n'.format(align.alignment_length()))
            output_f.write(
                'Alignment start position: {}.\n'.format(align.al_start))
            output_f.write(
                'Alignment end position: {}.\n'.format(align.al_end))
            output_f.write(
                'Numbers of mismatches: {}.\n'.format(len(align.mismatches())))
            output_f.write(
                'Gaps in {}: {}.\n'.format(seq_a_id, len(align.gaps('a'))))
            output_f.write(
                'Gaps in {}: {}.\n'.format(seq_b_id, len(align.gaps('b'))))
            output_f.write(
                'Identity (matches/alignment_length): {}%.\n'.format(
               round(len(align.matches()) / align.alignment_length() * 100, 4)))
            # Calculate match string that will show match '|' or mismatch '*' or
            # ' ' for a gap.
            match = ''
            for a, b in zip(align.seq_a, align.seq_b):
                    if a == '-' or b == '-':
                        match += ' '
                    elif a == b:
                        match += '|'
                    else:
                        match += '*'
            # First, separate lines by each 80 char.
            if len(match) < 80:  # In case length is less than 80 char.
                initial_idx = [0]
                final_idx = [len(match)]
            else:  # Indexes to cut.
                initial_idx = range(0, len(match), 80)
                final_idx = [idx for idx in initial_idx[1:]]
                # Correct the last idx.
                final_idx.append(len(match))
            # Write alignment in clustal alike format.
            for s, e in zip(initial_idx, final_idx):
                output_f.write(align.seq_a[s:e])
                output_f.write('\n')
                output_f.write(match[s:e])
                output_f.write('\n')
                output_f.write(align.seq_b[s:e])
                output_f.write('\n')
        output_f.write('-' * 80)
        output_f.write('\n')
        i += 1
    # Close file.
    output_f.close()
def write_identity_matrix(al, all_file, groups_file):
    """Writes identity matrices for a list of alignments.

    Arguments:
        al [list]: a 2D list (list of list) with all pairwise alignments
            (Align_formatted objects) for each sequence.
        all_file [string]: filename for all sequences identity matrix.
        groups_file [string]: filename for group identity matrix.
    Output:
        a_file [file]: a file with the pairwise identities for all the
            sequences. It is a csv table.
        g_file [file]: a file with the pairwise MINIMAL identities for each
            group of sequences. Sequence names are in the form of XY_ZZZZZ,
            where X is the group and Y is the subtype. ZZZZZ is the seq
            identifier. So, two sequences are of the same group if XY = XY for
            both. It is a csv table.
    """
    # Collecting all sequence names.
    id_list = [a[0].seq_a_id for a in al]
    # Collecting all group names. Caution! comming from an unordered set!.
    group_list = list({i.split('_')[0] for i in id_list})
    ############################################################################
    # Writing out a_file.
    a_file = open(all_file, 'w')
    # Writing column names.
    a_file.write(',' + ','.join(id_list) + '\n')
    # Writing each line.
    i = 0
    while i < len(al):
        # Seq id of the line to write.
        seq_id = al[i][0].seq_a_id
        a_file.write(seq_id)
        j = 0
        while j < len(al[i]):
            identity = str(round(len(al[i][j].matches()) /
                                 al[i][j].alignment_length() * 100, 2))
            # Writing the case of self-alignment, which is not in al[i] list.
            if j == len(al[i])-1 and i == len(al)-1:  # the very last alignment.
                a_file.write(',' + identity + ',100.00')
            elif j == i:  # Not the last alignment.
                a_file.write(',100.00,' + identity)
            else:
                a_file.write(',' + identity)
            j += 1
        a_file.write('\n')
        i += 1
    a_file.close()
    ############################################################################
    # Writing out g_file.
    g_file = open(groups_file, 'w')
    # Writing column names.
    g_file.write(',' + ','.join(group_list) + '\n')
    # A 2D list, each file and each column is a group. The values are the
    # MINIMUM idenity found for each pair of groups.
    all_groups_identities = []
    for g1 in group_list:
        group_identities = []  # Each group against the others.
        for g2 in group_list:
            identities = [round(len(a.matches())/a.alignment_length()*100, 2)
                          for seq in al for a in seq
                          if a.seq_a_id.split('_')[0] == g1 and
                          a.seq_b_id.split('_')[0] == g2]
            # Only self identity found. Only one seq in this group.
            if len(identities) == 0:
                identities = 'NA'
            else: # Select only the min value. Convert to a string.
                identities = str(min(identities))
            group_identities.append(identities)
        all_groups_identities.append(group_identities)
    # Writing each line.
    i = 0
    while i < len(all_groups_identities):
        g_file.write(group_list[i] + ',' + ','.join(all_groups_identities[i]))
        g_file.write('\n')
        i += 1
    g_file.close()
#------------------------------------------------------------------------------
# CLASSES
class Align_formatted(object):
    """This class provides a convenience object to store and use
    Bio.pairwise2.align output.

    Attributes:
        seq_a [string]: first sequence aligned.
        seq_b [string]: second sequence aligned.
        score [float]: alignment score.
        al_start [int]: alignment start position in the alignment, that is,
            in both sequences counting gaps.
        al_end [int]: alignment end position.
        seq_a_id [string]: id for seq a. Optional. Default = None.
        seq_b_id [string]: id for seq b. Optional. Default = None.
    Methods:
        gaps(self, sequence, gap_symbol): returns a list with each gap position
            in a sequence [list].
        matches(self): returns a list with the matched positions [list].
        mismatches(self): returns a list with the mismatched positions [int].
        alignment_length(self): returns al_end - al_start [int].
    """
    def __init__(self, alignment, seq_a_id = None, seq_b_id = None):
        # Alignment is the return of Bio.pairwise.align.globalxx or .localxx.
        # It is always a list of tuples in which each tuple is an alignment. The
        # first tuple is the best scoring alignment. Therefore, the first
        # tuple is the one that is taken.
        self.seq_a = alignment[0][0]
        self.seq_b = alignment[0][1]
        self.score = alignment[0][2]
        self.al_start = alignment[0][3]
        self.al_end = alignment[0][4]
        self.seq_a_id = seq_a_id
        self.seq_b_id = seq_b_id
    def gaps(self, sequence='a', gap_symbol = '-'):
        """Locates all gap positions in an aligned sequence.

        Arguments:
            sequence [char]: can be either 'a' or 'b' meaning sequence a or b.
                Default = 'a'.
            gap_symbol [char]: the symbol used to represent gaps in the
                alignment. Default = '-'.
        Output:
            gap symbol positions in sequence [list].
        """
        if sequence == 'a':
            seq = self.seq_a
        else:
            seq = self.seq_b
        gaps = [i for i, b in enumerate(seq) if b == gap_symbol]
        return gaps
    def matches(self):
        """Locates all matched positions in an alignment.

        Output:
            Matched positions in an alignment [list].
        """
        match = [i for i, s in enumerate(zip(self.seq_a, self.seq_b))
                       if s[0] == s[1]]
        return match
    def mismatches(self, gap_symbol = '-'):
        """Locates all mismatches positions in an alignment. A gap is not
        considered as a mismatch.

        Arguments:
            gap_symbol [char]: the symbol used to represent gaps in the
                alignment. Default = '-'.
        Output:
            mismatched positions in an alignment [list.]
        """
        mismatch = [i for i, s in enumerate(zip(self.seq_a, self.seq_b))
                       if s[0] != s[1] and gap_symbol not in s]
        return mismatch
    def alignment_length(self):
        """Computes alignment length."""
        return self.al_end - self.al_start

#------------------------------------------------------------------------------
# MAIN PROGRAM
def main_program():
    print("Starting amplicon alignments...")
    time0 = time.time()  # Start time.
    # Collecting argument values.
    args = arguments()
    # Parsing reference sequences and primers as a list of SeqRecord objects.
    genomes = list(SeqIO.parse(args.reference, "fasta"))
    primers = list(SeqIO.parse(args.primers, "fasta"))
    # Defining output filenames.
    basename = args.outfile
    amplicons_file = basename + '_amplicons.fasta'
    alignments_file = basename + '_alignments.txt'
    identity_df_all_file = basename + '_identity_df_all.csv'
    identity_df_groups_file = basename + '_identity_df_per_group.csv'
    # Align primers against references to create amplicons.
    amplicons = virtual_pcr(genomes, primers)
    # Output amplicons file.
    SeqIO.write(amplicons, amplicons_file, "fasta")
    # Optimal global pairwise alignment (Needleman–Wunsch) for all amplicons.
    amplicon_aligns = nw_global_alignment(amplicons)
    # Write alignments file.
    write_pairwise_al(amplicon_aligns, alignments_file)
    # writing identity matrix for each sequence and for subtype.
    write_identity_matrix(amplicon_aligns, identity_df_all_file,
                          identity_df_groups_file)
    print("\t...amplicon alignments performed in {} seconds".format(time.time()-
                                                                    time0))
#------------------------------------------------------------------------------
# Conditional to run the script
if __name__ == '__main__':
    main_program()
