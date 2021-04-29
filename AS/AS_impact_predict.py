#!/usr/bin/env python3
# coding = utf-8
# -*- coding: utf-8 -*-
#@Author: 'Lvmj'

# from Bio import Seq
import re, os, collections, argparse
from Bio import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
##This software is for prediction of putative impact of AS which parse by ASTALAVISTA

parser = argparse.ArgumentParser(description='This software is for prediction of putative impact of AS which parse by ASTALAVISTA.')
parser.add_argument('-g','--gtf', dest='gtf', metavar='', help='reference gtf', required=True)
parser.add_argument('-a','--ag', dest='ag', metavar='', help='assembled gtf', required=True)
parser.add_argument('-as','--asg', dest='asg', metavar='', help='astalavista gtf', required=True)
parser.add_argument('-c','--comp', dest='comp', metavar='', help='assembled to reference genes text file', required=True)
parser.add_argument('-r','--reference', dest='reference', metavar='', help='reference fasta', required=True)
parser.add_argument('-o','--output', dest='output', metavar='', help='output file', required=True)
args = parser.parse_args()


def findall(s, ch):
    return [i for i, ltr in enumerate(s) if ltr == ch]


def DNA_complement(seq):
    complement = {"A":"T", "G":"C", "T":"A", "C":"G", "N":"N"}
    seq = seq.upper()
    seq = list(seq)
    rev_seq = ''
    for i in seq:
        rev_seq += complement[i]
    return rev_seq[::-1]


fa = args.reference  # reference fasta
seq = {}
for Seq_record in SeqIO.parse(fa,"fasta"):
    seq[Seq_record.id] = Seq_record.seq


Gid2Euid = {}  #assembled ID to reference ID
fin1 = open(args.comp, 'r')
fin1.readline()
for line in fin1:
    line = line.strip("\n")
    _line = line.split("\t")
    Gid = _line[0]
    Euid = _line[3]    #EuID in fourth column
    if len(Euid.split(',')) == 1 and Euid.startswith('Eu'):  #only the one to one match
        Gid2Euid[Gid] = Euid[0:-5]
        # print(Euid[0:-5])
# print(len(Gid2Euid))

cds_infor = collections.defaultdict(list)  #CDS infomation of reference
fin2 = open(args.gtf, 'r')
for line in fin2:
    line = line.strip("\n")
    _line = line.split("\t")
    if _line[2] == "CDS":
        chro = _line[0]
        start = int(_line[3])
        end = int(_line[4])
        strand = _line[6]
        tuid = re.match(r'.*transcript_id \"(.+?)\"', _line[8]).group(1)[0:-5]
        # gid = re.match(r'.*gene_id \"(.+?)\"', _line[8]).group(1)[0:-5]
        cds_infor[tuid].extend([start, end])


fin3 = open(args.ag, 'r')
cds_seq = {}  #cds sequences
direction = {}
tid2gid = {}  # transcript to gene id
for line in fin3:
    line = line.strip("\n")
    _line = line.split("\t")
    if _line[2] == "exon":
        chro = _line[0]
        start = int(_line[3])
        end = int(_line[4])
        strand = _line[6]
        tuid = re.match(r'.*transcript_id \"(.+?)\"', _line[8]).group(1)
        # print(tuid)
        # print(tuid)
        direction[tuid] = strand
        gid = re.match(r'.*gene_id \"(.+?)\"', _line[8]).group(1)
        tid2gid[tuid] = gid
        if gid in Gid2Euid:
            # print(gid, tuid)
            if strand == '-':  # strand -
                atg = int(cds_infor[Gid2Euid[gid]+".1"][-1])  #atg at the end
                # print(Gid2Euid[gid],gid)
                if end <= atg:
                    if tuid not in cds_seq:
                        cds_seq[tuid] = seq[chro][(start-1):end]
                        # print(tuid, cds_seq[tuid])
                    else:
                        cds_seq[tuid] = cds_seq[tuid] + seq[chro][(start-1):end]
                else:
                    if tuid not in cds_seq:
                        if start <= atg:
                            cds_seq[tuid] = seq[chro][(start-1):atg]
                    else:  # 接前面外显子序列
                        if start <= atg:
                            cds_seq[tuid] = cds_seq[tuid] + seq[chro][(start-1):atg]

            else:  # strand +
                atg = int(cds_infor[Gid2Euid[gid]+".1"][0])
                # print(Gid2Euid[gid],gid)
                if end >= atg:
                    if tuid not in cds_seq:
                        cds_seq[tuid] = seq[chro][(atg-1):end]
                    else:
                        cds_seq[tuid] = cds_seq[tuid] + seq[chro][(start-1):end]

def protein(id):  #transcrition sequence
    if direction[id] == "-":
        cds = DNA_complement(cds_seq[id])
    else:
        cds = cds_seq[id]

    dna = Seq.Seq(str(cds))
    mrna = dna.transcribe()
    protein = mrna.translate()
    return(protein)

All_id = []
fin4 = open(args.asg, 'r')
fot = open(args.output, 'w')
fot.write("AS Events\tGid\tEuid\tAS type\tPutative impact\n")
for line in fin4:
    line = line.strip("\n")
    _line = line.split("\t")
    chro = _line[0]
    start = int(_line[3])
    end = int(_line[4])
    strand = _line[6]
    geneid = re.match(r'.*gene_id \"(.+?)\"', _line[8]).group(1)
    splice_chain = re.match(r'.*splice_chain \"(.+?)\"', _line[8]).group(1)
    as_event = geneid+";"+splice_chain
    tuid1, tuid2 = re.match(r'.*transcript_id \"(.+?)\"', _line[8]).group(1).split(",")[0].split("/")[0], \
                   re.match(r'.*transcript_id \"(.+?)\"', _line[8]).group(1).split(",")[1].split("/")[0]
    structure = re.match(r'.*structure \"(.+?)\"', _line[8]).group(1)
    Gid = tid2gid[tuid1]
############RI#############
    if structure == "0,1^2-":  #RI 内含子保留
        if strand == "-":
            splice_end, splice_start = int(re.match(r'.*splice_chain \"(.+?)\"', _line[8]).group(1).split(",")[1].split("^")[0]), \
                                       int(re.match(r'.*splice_chain \"(.+?)\"', _line[8]).group(1).split(",")[1].split("^")[1][0:-1])
        else:
            splice_start, splice_end = int(re.match(r'.*splice_chain \"(.+?)\"', _line[8]).group(1).split(",")[1].split("^")[0]), \
                                       int(re.match(r'.*splice_chain \"(.+?)\"', _line[8]).group(1).split(",")[1].split("^")[1][0:-1])
        # print(tid2gid[tuid1], cds_infor[tid2gid[tuid1]])
        if tid2gid[tuid1] in Gid2Euid:  # if annotated genes
            if Gid2Euid[tid2gid[tuid1]] not in All_id:
                All_id.append(Gid2Euid[tid2gid[tuid1]])
            Euid = Gid2Euid[tid2gid[tuid1]]+".1"
            if splice_start > cds_infor[Euid][0] and splice_end < cds_infor[Euid][-1]:  # if the CDS in UTR
                splice_dis = abs(splice_end - splice_start) - 1
                if splice_dis % 3 == 0:   #frame shift detect
                    if tuid1 in cds_seq and tuid2 in cds_seq:
                        try:
                            pro_dis = protein(tuid1).find("*") - protein(tuid2).find("*")  #first stop codon
                            # splice_dis = splice_end - splice_start - 1
                            if pro_dis == splice_dis/3:
                                fot.write(as_event+"\t"+Gid+"\t"+Euid+"\tRI\tInsertion/Deletion\n")
                            else:
                                fot.write(as_event+"\t"+Gid+"\t"+Euid+"\tRI\tPTC\n")
                        except:
                            fot.write(as_event+"\t"+Gid+"\t"+Euid + "\tRI\tNA\n")
                    else:
                        fot.write(as_event+"\t"+Gid+"\t"+Euid + "\tRI\tstart or stop codon lost\n")
                else:
                    fot.write(as_event+"\t"+Gid+"\t"+Euid + "\tRI\tframe shift\n")

            elif (splice_start < cds_infor[Euid][-1] and splice_end > cds_infor[Euid][-1]) or \
                (splice_start < cds_infor[Euid][0] and splice_end > cds_infor[Euid][0]):
                fot.write(as_event+"\t"+Gid+"\t"+Euid + "\tRI\tstart or stop codon lost\n")
            else:
                fot.write(as_event+"\t"+Gid+"\t"+Euid + "\tRI\tUTR\n")
        else:
            fot.write(as_event+"\t"+Gid+"\tNA\tRI\tNA\n")

##########A3########
    elif structure == "1-,2-":  #A3(AA)
        if strand == "-":
            splice_end, strand_start = int(re.match(r'.*splice_chain \"(.+?)\"', _line[8]).group(1).split(",")[0][0:-1]), \
                                       int(re.match(r'.*splice_chain \"(.+?)\"', _line[8]).group(1).split(",")[1][0:-1])
        else:
            splice_start, splice_end = int(re.match(r'.*splice_chain \"(.+?)\"', _line[8]).group(1).split(",")[0][0:-1]), \
                                       int(re.match(r'.*splice_chain \"(.+?)\"', _line[8]).group(1).split(",")[1][0:-1])
        # print(tid2gid[tuid1], cds_infor[tid2gid[tuid1]])
        if tid2gid[tuid1] in Gid2Euid:  #判断annotated gene
            if Gid2Euid[tid2gid[tuid1]] not in All_id:
                All_id.append(Gid2Euid[tid2gid[tuid1]])
            Euid = Gid2Euid[tid2gid[tuid1]]+".1"
            if splice_start > cds_infor[Euid][0] and splice_end < cds_infor[Euid][-1]:
                splice_dis = abs(splice_end - splice_start)
                if splice_dis % 3 == 0:
                    if tuid1 in cds_seq and tuid2 in cds_seq:
                        try:
                            pro_dis = protein(tuid1).find("*") - protein(tuid2).find("*")
                            # splice_dis = splice_end - splice_start - 1
                            if pro_dis == splice_dis/3:
                                fot.write(as_event+"\t"+Gid+"\t"+Euid+"\tA3\tInsertion/Deletion\n")
                            else:
                                fot.write(as_event+"\t"+Gid+"\t"+Euid+"\tA3\tPTC\n")
                        except:
                            fot.write(as_event+"\t"+Gid+"\t"+Euid + "\tA3\tNA\n")
                    else:
                        fot.write(as_event+"\t"+Gid+"\t"+Euid + "\tA3\tstart or stop codon lost\n")
                else:
                    fot.write(as_event+"\t"+Gid+"\t"+Euid + "\tA3\tframe shift\n")
            elif (splice_start < cds_infor[Euid][-1] and splice_end > cds_infor[Euid][-1]) or \
                (splice_start < cds_infor[Euid][0] and splice_end > cds_infor[Euid][0]):
                fot.write(as_event+"\t"+Gid+"\t"+Euid + "\tA3\tstart or stop codon lost\n")
            else:
                fot.write(as_event+"\t"+Gid+"\t"+Euid + "\tA3\tUTR\n")
        else:
            fot.write(as_event+"\t"+Gid+"\tNA\tA3\tNA\n")

##########A5#####
    elif structure == "1^,2^":  #A5(AD)
        if strand == "-":
            splice_end, strand_start = int(re.match(r'.*splice_chain \"(.+?)\"', _line[8]).group(1).split(",")[0][0:-1]), \
                                       int(re.match(r'.*splice_chain \"(.+?)\"', _line[8]).group(1).split(",")[1][0:-1])
        else:
            splice_start, splice_end = int(re.match(r'.*splice_chain \"(.+?)\"', _line[8]).group(1).split(",")[0][0:-1]), \
                                       int(re.match(r'.*splice_chain \"(.+?)\"', _line[8]).group(1).split(",")[1][0:-1])
        # print(tid2gid[tuid1], cds_infor[tid2gid[tuid1]])
        if tid2gid[tuid1] in Gid2Euid:  # annotated gene
            if Gid2Euid[tid2gid[tuid1]] not in All_id:
                All_id.append(Gid2Euid[tid2gid[tuid1]])
            Euid = Gid2Euid[tid2gid[tuid1]]+".1"
            if splice_start > cds_infor[Euid][0] and splice_end < cds_infor[Euid][-1]:
                splice_dis = abs(splice_end - splice_start)
                if splice_dis % 3 == 0:
                    if tuid1 in cds_seq and tuid2 in cds_seq:
                        try:
                            pro_dis = protein(tuid2).find("*") - protein(tuid1).find("*")
                            # splice_dis = splice_end - splice_start - 1
                            if pro_dis == splice_dis/3:
                                fot.write(as_event+"\t"+Gid+"\t"+Euid+"\tA5\tInsertion/Deletion\n")
                            else:
                                fot.write(as_event+"\t"+Gid+"\t"+Euid+"\tA5\tPTC\n")
                        except:
                            fot.write(as_event+"\t"+Gid+"\t"+Euid + "\tA5\tNA\n")
                    else:
                        fot.write(as_event+"\t"+Gid+"\t"+Euid + "\tA5\tstart or stop codon lost\n")
                else:
                    fot.write(as_event+"\t"+Gid+"\t"+Euid + "\tA5\tframe shift\n")
            elif (splice_start < cds_infor[Euid][-1] and splice_end > cds_infor[Euid][-1]) or \
                (splice_start < cds_infor[Euid][0] and splice_end > cds_infor[Euid][0]):
                fot.write(as_event+"\t"+Gid+"\t"+Euid + "\tA5\tstart or stop codon lost\n")
            else:
                fot.write(as_event+"\t"+Gid+"\t"+Euid + "\tA5\tUTR\n")
        else:
            fot.write(as_event+"\t"+Gid+"\tNA\tA5\tNA\n")

##########SE#####
    elif structure == "0,1-2^":  #A5
        if strand == "-":
            splice_end, splice_start = int(re.match(r'.*splice_chain \"(.+?)\"', _line[8]).group(1).split(",")[1].split("-")[0]), \
                                       int(re.match(r'.*splice_chain \"(.+?)\"', _line[8]).group(1).split(",")[1].split("-")[1][0:-1])
        else:
            splice_start, splice_end = int(re.match(r'.*splice_chain \"(.+?)\"', _line[8]).group(1).split(",")[1].split("-")[0]), \
                                       int(re.match(r'.*splice_chain \"(.+?)\"', _line[8]).group(1).split(",")[1].split("-")[1][0:-1])
        # print(tid2gid[tuid1], cds_infor[tid2gid[tuid1]])
        if tid2gid[tuid1] in Gid2Euid:
            if Gid2Euid[tid2gid[tuid1]] not in All_id:
                All_id.append(Gid2Euid[tid2gid[tuid1]])
            Euid = Gid2Euid[tid2gid[tuid1]]+".1"
            if splice_start > cds_infor[Euid][0] and splice_end < cds_infor[Euid][-1]:
                splice_dis = abs(splice_end - splice_start)+1
                if splice_dis % 3 == 0:
                    if tuid1 in cds_seq and tuid2 in cds_seq:
                        try:
                            pro_dis = protein(tuid2).find("*") - protein(tuid1).find("*")
                            # splice_dis = splice_end - splice_start - 1
                            if pro_dis == splice_dis/3:
                                fot.write(as_event+"\t"+Gid+"\t"+Euid+"\tSE\tInsertion/Deletion\n")
                            else:
                                fot.write(as_event+"\t"+Gid+"\t"+Euid+"\tSE\tPTC\n")
                        except:
                            fot.write(as_event+"\t"+Gid+"\t"+Euid + "\tSE\tNA\n")
                    else:
                        fot.write(as_event+"\t"+Gid+"\t"+Euid + "\tSE\tstart or stop codon lost\n")
                else:
                    fot.write(as_event+"\t"+Gid+"\t"+Euid + "\tSE\tframe shift\n")
            elif (splice_start < cds_infor[Euid][-1] and splice_end > cds_infor[Euid][-1]) or \
                (splice_start < cds_infor[Euid][0] and splice_end > cds_infor[Euid][0]):
                fot.write(as_event+"\t"+Gid+"\t"+Euid + "\tSE\tstart or stop codon lost\n")
            else:
                fot.write(as_event+"\t"+Gid+"\t"+Euid + "\tSE\tUTR\n")
        else:
            fot.write(as_event+"\t"+Gid+"\tNA\tSE\tNA\n")


##########MXE#####
    elif structure == "1-2^,3-4^":  #MXE
        if strand == "-":
            splice_end2, splice_start2, splice_end1, splice_start1 = \
                int(re.match(r'.*splice_chain \"(.+?)\"', _line[8]).group(1).split(",")[0].split("-")[0]), \
                                                                     int(re.match(r'.*splice_chain \"(.+?)\"', _line[8]).
                                                                        group(1).split(",")[0].split("-")[1][0:-1]),\
                                                                     int(re.match(r'.*splice_chain \"(.+?)\"', _line[8]).
                                                                        group(1).split(",")[1].split("-")[0]), \
                                                                     int(re.match(r'.*splice_chain \"(.+?)\"', _line[8]).
                                                                        group(1).split(",")[1].split("-")[1][0:-1]),

        else:
            splice_start1, splice_end1, splice_start2, splice_end2 = \
                int(re.match(r'.*splice_chain \"(.+?)\"', _line[8]).group(1).split(",")[0].split("-")[0]), \
                                                                     int(re.match(r'.*splice_chain \"(.+?)\"', _line[8]).
                                                                        group(1).split(",")[0].split("-")[1][0:-1]),\
                                                                     int(re.match(r'.*splice_chain \"(.+?)\"', _line[8]).
                                                                        group(1).split(",")[1].split("-")[0]), \
                                                                     int(re.match(r'.*splice_chain \"(.+?)\"', _line[8]).
                                                                        group(1).split(",")[1].split("-")[1][0:-1]),
        # print(tid2gid[tuid1], cds_infor[tid2gid[tuid1]])
        if tid2gid[tuid1] in Gid2Euid:
            if Gid2Euid[tid2gid[tuid1]] not in All_id:
                All_id.append(Gid2Euid[tid2gid[tuid1]])
            Euid = Gid2Euid[tid2gid[tuid1]]+".1"
            if splice_start1 > cds_infor[Euid][0] and splice_end2 < cds_infor[Euid][-1]:
                splice_dis = abs(splice_end1 - splice_start1)-abs(splice_end2 - splice_start2)
                if splice_dis % 3 == 0:
                    if tuid1 in cds_seq and tuid2 in cds_seq:
                        try:
                            pro_dis = protein(tuid2).find("*") - protein(tuid1).find("*")
                            # splice_dis = splice_end - splice_start - 1
                            if pro_dis == splice_dis/3:
                                fot.write(as_event+"\t"+Gid+"\t"+Euid+"\tMXE\tInsertion/Deletion\n")
                            else:
                                fot.write(as_event+"\t"+Gid+"\t"+Euid+"\tMXE\tPTC\n")
                        except:
                            fot.write(as_event+"\t"+Gid+"\t"+Euid + "\tMXE\tNA\n")
                    else:
                        fot.write(as_event+"\t"+Gid+"\t"+Euid + "\tMXE\tstart or stop codon lost\n")
                else:
                    fot.write(as_event+"\t"+Gid+"\t"+Euid + "\tMXE\tframe shift\n")
            elif (splice_start1 < cds_infor[Euid][-1] and splice_end2 > cds_infor[Euid][-1]) or \
                (splice_start1 < cds_infor[Euid][0] and splice_end2 > cds_infor[Euid][0]):
                fot.write(as_event+"\t"+Gid+"\t"+Euid + "\tMXE\tstart or stop codon lost\n")
            else:
                fot.write(as_event+"\t"+Gid+"\t"+Euid + "\tMXE\tUTR\n")
        else:
            fot.write(as_event+"\t"+Gid+"\tNA\tMXE\tNA\n")
# print(len(All_id))



