#!/usr/bin/env python
# DADA2 PIPELINE: script 1


import gzip
import argparse
import numpy as np


def get_args():
    """ get file names for read 1, read 2, quality score threshhold, and tolerance for mismatch """
    parser = argparse.ArgumentParser(description=' File name to be used')
    parser.add_argument('-R', '--Read', type=str,help='File Name for Read 1 as a gzip file', required=True)
    parser.add_argument('-r', '--read', type=str, help='File Name for Read 2 as a gzip file', required=True)
    parser.add_argument('-I', '--Index', type=str, help='File Name for index 1 as a gzip file', required=True)
    parser.add_argument('-i', '--index', type=str, help='File Name for index 2 as a gzip file', required=True)
    parser.add_argument('-b', '--map', type=str, help='Mapping File as .csv ', required=True)
    parser.add_argument('-P', '--forwardPrimer', type=str, help='Forward primer sequence', required=True)
    parser.add_argument('-p', '--reversePrimer', type=str, help='Reverse Primer Sequence', required=True)
    parser.add_argument('-t', '--tolerance', type=str, required=True \
                        , help='Amount of mismatched nucleotides allowed in forwardPrimer and reversePrimer, must be'\
                               'less than 5')
    parser.add_argument('-q', '--quality', type=str, help='quality score threshold', required=True)
    return parser.parse_args()


args= get_args()
r1 = args.Read
r2 = args.read
i1 = args.Index
i2 = args.index
b = args.map
q = args.quality
fw = args.forwardPrimer
rv = args.reversePrimer
t = args.tolerance


def barcodes(b):
    """takes barcodes and sample IDs from a csv separated file and returns an array of known combinations and samples"""
    with open(b,'r') as bar:
        barcodes = []
        map = []
        for line in bar:
            line = line.strip()
            x = line.split(',')[0]
            y = line.split(',')[1]
            z = line.split(',')[2]
            barcodes.append([y, z])
            map.append(x)
        barcodes.pop(0)
        map.pop(0)
        barcodes[:] = [x for x in barcodes if x != ['', '']]
        #removes empty indexes
        map = tuple(map)
        barcodes = tuple(barcodes)

        return barcodes, map


bar, sample = barcodes(b)

def rev(score):
    """Reverses a given string, used to reverse quality score line in this program"""
    return score[::-1]


def revcomp(seq):
    # revcomp creates a reverse compliment of a given string of nucleotides
    rev = seq[::-1]
    new = ''
    for n in rev:
        if n == 'A':
            new += 'T'
        elif n == 'T':
            new += 'A'
        elif n == 'C':
            new += 'G'
        elif n == 'G':
            new += 'C'
        elif n == 'N':
            new += 'N'

    return new


def qual(q,seq):
    #takes a quality score and checks if the mean phred score of a given score line is greater and returns true if it is
    #and false if not
    sum=0
    for score in seq:
        # sums all quality scores
        n = (((ord(score) - 33)))
        sum += n
    if (int(sum)/len(seq))>=int(q):
        return True
    else:
        return False


def tolerance(t,r,fw):
    if int(t) >= sum(a != b for a,b in zip(fw, r)):
        return 0


def despace(t,l1,s1,l2,s2,fw,rv):
    """removes spacers an returns properly oriented read and q score lines
     t is tolerance, l is read line, s is the quality score line
     returns """
    win = len(fw + rv)
    for i in range(win):
        # defines the size of the serach window
        seq1 = l1[i:i+len(fw)]
        # defines the search sequence matching fw length
        seq2 = l1[i:i+len(rv)]
        # defines the search sequence matching rv length
        if tolerance(t, seq1, fw) == 0:
            # checks if the search window matches with the given tolerances
            for x in range(win):
                # searches for the rv sequence in the second read
                seq3 = l2[x:x+len(rv)]
                if tolerance(t, seq3, rv) == 0:
                    return l1[i:], s1[i:], l2[x:], s2[x:], True
        elif tolerance(t, seq2, rv) == 0:
            # checks if the search window matches with the given tolerances
            for x in range(win):
                # searches for the rv sequence in the second read
                seq3 = l2[x:x+len(fw)]
                if tolerance(t, seq3, fw) == 0:
                    return l1[i:], s1[i:], l2[x:], s2[x:], False
    return l1, s1, l2, s2, None
    #if nothing mathces, returns None

#despace(0,'AAAAACAGCCGCGGTAATTCCAGCTCCCGACCGACCCAGCCCTTAGAGCCAATCCTTATCCCGAAGTTACGGATCCGGCTTGCCGACTTCCCTTACCTACATTGTTCCAACATGCCAGAGGCTGTT','AAAAAAAAAAAAAAAAAAAAAAAAAAAFFFJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJFJJJJJJJ','GAATCCGACAGTTTATACAAGAGTTTACATTCCGTGGCCTTTTGCCACACATAGAGAAAACAATTCGTCAGTTAAATGATCAGCTAATATCAAGGAAAGGTGAACCCAAACACTTTGGTTTCCAAAA','AAAFFJJFJFJJFFJJJJJJFJJFJ7JJJJJJJJJFJJAAJJJFJJJJJJJJJJJFJJJJJJJJFJJJJ7FA7FJJJA<JJJFAJJJJJJJFJJ<<<FFJFAAAAAAAAAAAAAAAAAAAAAAAAAA','CAGCCGCGGTAATTCCAGCT','GAACCCAAACACTTTGGTTTCC')

def demult(r1,r2,i1,i2,bar,q,t,fw,rv,sample):
    #takes read and index files record by record, checks if indexes contain an N, and if it does places it in low quality,
    #checks if average read quality is above or equal to user input threshold, if not places record in low quality files,
    #checks if index 1 is the reverse compliment of index two, if not places in unmatched index files, if true places
    #in that barcode file. All headers are appended to include the barcodes at the start. Also includes counters  All
    # headers are appended to include the barcodes at the start. Returns counts for each index, unknown indexes,
    # and index hopped indices as well as creating files

    with gzip.open(r1, 'rt') as r1, gzip.open(r2, 'rt') as r2, gzip.open(i1, 'rt') as i1, gzip.open(i2, 'rt') as i2, \
    open('unknown_R1.fastq', 'a') as uf, open('unknown_R2.fastq', 'a') as ur:
        #open read files and unknown files
        codes = []
        for x in sample:
            # generates forward and reverse files for each sample
            codes.append(x+'_R1.fastq')
            codes.append(x+'_R2.fastq')
        files = {code: open(code,'a') for code in codes}
        #opens all sample files to write to
        counts = {'unknown': 0}
        ln = 0
        r1temp = []
        r2temp = []
        i1temp = []
        i2temp = []
        i2_rev = []
        #set temp variables
        for r1_line, r2_line, i1_line, i2_line in zip(r1,r2,i1,i2):
            ln += 1
            r1temp.append(r1_line.strip())
            r2temp.append(r2_line.strip())
            i1temp.append(i1_line.strip())
            i2temp.append(i2_line.strip())
            # add record lines to temp files
            if ln % 4 == 2:
                #if the record is the read line, revcomps the r3 barcode and stores it
                i1_rev = revcomp(i1_line.strip())
            if ln % 4 == 0:
                #checks if a full record is stored
                r1temp[1], r1temp[3], r2temp[1], r2temp[3], ord = despace(t, r1temp[1], r1temp[3], r2temp[1], r2temp[3], fw, rv)
                # removes adapter and checks orientation
                if ord == None:
                    # if the fw and rv seqs were not found, places read in unknown
                    uf.write(r1temp[0] + '\n' + r1temp[1] + '\n' + r1temp[2] + '\n' + r1temp[3] + '\n')
                    ur.write(r2temp[0]+ '\n' + r2temp[1] + '\n' + r2temp[2] + '\n' + r2temp[3] + '\n')
                    counts['unknown'] += 1
                    r1temp = []
                    r2temp = []
                    i1temp = []
                    i2temp = []
                elif qual(q, r1temp[3]) and qual(q, r2temp[3]) and qual(q, i1temp[3]) and qual(q, i2temp[3]) is True:
                    # tests for quality score passing
                    if [i2temp[1], i1_rev] in bar:
                        # checks if barcodes are known in correct orientation
                        pair = [i2temp[1],i1_rev]
                        name = sample[bar.index(pair)]
                        # pulls the correct sample name out of the list
                        if ord is True:
                            # if r1 is fw and r2 is rv places in correct files
                            files[(name+"_R1.fastq")].write(r1temp[0] + '\n' + r1temp[1] + '\n' + r1temp[2] + '\n' + r1temp[3] + '\n')
                            files[(name+"_R2.fastq")].write(r2temp[0]+ '\n' + r2temp[1] + '\n' + r2temp[2] + '\n' + r2temp[3] + '\n')
                            if name in counts:
                                counts[name] += 1
                                r1temp = []
                                r2temp = []
                                i1temp = []
                                i2temp = []
                            else:
                                counts[name] = 1
                                r1temp = []
                                r2temp = []
                                i1temp = []
                                i2temp = []
                        else:
                            # if r1 is rv and r2 is fw, places in appropriate files
                            files[(name+"_R2.fastq")].write(r1temp[0] + '\n' + r1temp[1] + '\n' + r1temp[2] + '\n' + r1temp[3] + '\n')
                            files[(name+"_R1.fastq")].write(r2temp[0] + '\n' + r2temp[1] + '\n' + r2temp[2] + '\n' + r2temp[3] + '\n')
                            if name in counts:
                                counts[name] += 1
                                r1temp = []
                                r2temp = []
                                i1temp = []
                                i2temp = []
                            else:
                                counts[name] = 1
                                r1temp = []
                                r2temp = []
                                i1temp = []
                                i2temp = []
                    elif[i1_rev, i2temp[1]] in bar:
                        # if barcode is in opposite orientation
                        pair = [i1_rev, i2temp[1]]
                        name = sample[bar.index(pair)]
                        # pulls the correct sample name out
                        if ord is True:
                            # if r1 is fw and r2 is rv places in correct files
                            files[(name + "_R1.fastq")].write(r1temp[0] + '\n' + r1temp[1] + '\n' + r1temp[2] + '\n' + r1temp[3] + '\n')
                            files[(name + "_R2.fastq")].write(r2temp[0] + '\n' + r2temp[1] + '\n' + r2temp[2] + '\n' + r2temp[3] + '\n')
                            if name in counts:
                                counts[name] += 1
                                r1temp = []
                                r2temp = []
                                i1temp = []
                                i2temp = []
                            else:
                                counts[name] = 1
                                r1temp = []
                                r2temp = []
                                i1temp = []
                                i2temp = []
                        else:
                            # if r1 is rv and r2 is fw, places in appropriate files
                            files[(name + "_R2.fastq")].write(r1temp[0] + '\n' + r1temp[1] + '\n' + r1temp[2] + '\n' + r1temp[3] + '\n')
                            files[(name + "_R1.fastq")].write(r2temp[0] + '\n' + r2temp[1] + '\n' + r2temp[2] + '\n' + r2temp[3] + '\n')
                            if name in counts:
                                counts[name] += 1
                                r1temp = []
                                r2temp = []
                                i1temp = []
                                i2temp = []
                            else:
                                counts[name] = 1
                                r1temp = []
                                r2temp = []
                                i1temp = []
                                i2temp = []
                    else:
                        uf.write(r1temp[0] + '\n' + r1temp[1] + '\n' + r1temp[2] + '\n' + r1temp[3] + '\n')
                        ur.write(r2temp[0] + '\n' + r2temp[1] + '\n' + r2temp[2] + '\n' + r2temp[3] + '\n')
                        counts['unknown'] += 1
                        r1temp = []
                        r2temp = []
                        i1temp = []
                        i2temp = []
                else:
                    uf.write(r1temp[0] + '\n' + r1temp[1] + '\n' + r1temp[2] + '\n' + r1temp[3] + '\n')
                    ur.write(r2temp[0] + '\n' + r2temp[1] + '\n' + r2temp[2] + '\n' + r2temp[3] + '\n')
                    counts['unknown'] += 1
                    r1temp = []
                    r2temp = []
                    i1temp = []
                    i2temp = []
    for file in files.values():
        file.close()
    print("Barcode\tCount")
    for item in counts:
        print(str(item)+'\t'+str(counts[item]))
    return counts


counts = demult(r1, r2, i1, i2, bar, q, t, fw, rv, sample)

