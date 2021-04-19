#!/usr/bin/env python3
# coding = utf-8
# -*- coding: utf-8 -*-
#@Author: 'Lvmj'

import re, argparse
from collections import defaultdict
from Bio import SeqIO
from sys import argv


parser = argparse.ArgumentParser(description='Taco result gtf filter.')
parser.add_argument('-i','--input', dest='input', metavar='', help='input gtf', required=True)
parser.add_argument('-r','--reference', dest='reference', metavar='', help='reference fasta', required=True)
parser.add_argument('-o','--output', dest='output', metavar='', help='output gtf', required=True)
args = parser.parse_args()

fot = open(args.output, 'w')
fin = args.input
Seq = {}
fa = args.reference ##### 参考fasta
for Seq_record in SeqIO.parse(fa,"fasta"):
    Seq[Seq_record.id] = Seq_record.seq


def exonseq(tuple,Chr):
    tuple = sorted(tuple)
    exseq = []
    for i in range(0, len(tuple), 2):
        exseqs = str(Seq[Chr][int(tuple[i])-1:int(tuple[i+1])]).upper()
        exseq.append(exseqs)
    return exseq


def intronseq(tuple, Chr):
    tuple = sorted(tuple)
    inseq = []
    if len(tuple) == 2:
        inseq = 0
    else:
        for i in range(0, len(tuple)-2, 2):
            inseqs = str(Seq[Chr][int(tuple[i+1])-1:int(tuple[i+2])]).upper()
            inseq.append(inseqs)
    return inseq


def exonlen(tuple):
    tuple=sorted(tuple)
    exlen=[]
    for i in range(0, len(tuple), 2):
        exlens=int(tuple[i+1])-int(tuple[i])+1
        exlen.append(exlens)
    return exlen


def intronlen(tuple):
    tuple = sorted(tuple)
    inlen = []
    if len(tuple)==2:
        inlen=[0]
    else:
        for i in range(0, len(tuple)-2, 2):
            inlens=int(tuple[i+2])-int(tuple[i+1])+1
            inlen.append(inlens)
    return inlen


with open(fin, 'r') as f:
    startends = defaultdict(list)
    tiChr = {}
    tigi = {}
    allinfor = {}
    transline = {}
    for line in f:
        line = line.strip()
        _line = line.split('\t')
        if re.search(r'Chr.+|chr.+|scaffold|Pp', _line[0]):
            if _line[2] == 'transcript':
                infor = _line[8]
                matchtrid = re.match(r'.+transcript_id\ \"(TU\d+).+', infor)
                if matchtrid:
                    trid = matchtrid.group(1)
                transline[trid] = line

            elif _line[2] == 'exon':
                infor = _line[8]
                matchtrid = re.match(r'.+transcript_id\ \"(TU\d+).+', infor)
                if matchtrid:
                    trid = matchtrid.group(1)
                matchgid = re.match(r'.+gene_id\ \"(G\d+).+',infor)
                if matchgid:
                    gid = matchgid.group(1)
                tigi[trid] = gid
                Chr = _line[0]
                tiChr[trid] = Chr
                sta = int(_line[3])
                end = int(_line[4])
                stainfo = Chr+str(sta)+trid
                allinfor[stainfo] = line.strip()
                strand = _line[6]
                if strand == ".":
                    strand = "+"
                startends[trid].append(sta)
                startends[trid].append(end)

# 过滤
startends1 = defaultdict(list)
Tuinfor = defaultdict(list)
for j in startends:
    exons = exonlen(startends[j])
    introns = intronlen(startends[j])
    exon_seq = exonseq(startends[j], tiChr[j])
    intron_seq = intronseq(startends[j], tiChr[j])
    if max(introns) >150000:    # over 150kb
        if j in transline:
            transline.pop(j)
    else:
        for k in range(0, len(exons)):
            stainfo1 = tiChr[j]+str(startends[j][2*k])+j
            atcg = {}
            for m in str(exon_seq[k]):
                atcg[m] = str(exon_seq[k]).count(m)
            if 'T' in atcg:
                percentT = atcg['T']/len(exon_seq[k])

            if 'A' in atcg:
                percentA = atcg['A']/len(exon_seq[k])

            if exons[k] < 10:      # filter exon < 10 bp
                if j in transline:
                    transline.pop(j)
                continue
            elif percentT > 0.7 and exons[k] < 50:      # filter ployT
                if j in transline:
                    transline.pop(j)
                continue

            elif percentA > 0.7 and exons[k] < 50:      # filter polyA
                if j in transline:
                    transline.pop(j)
                continue

            elif re.search(r'TCTCTCTCTCTCT|GAGAGAGAGAGAG',str(exon_seq[k])) and exons[k] < 50:  # AG TC repeat
                if j in transline:
                    transline.pop(j)
                continue
            else:
                _line1 = allinfor[stainfo1].split('\t')
                if re.search(r'Chr.+|chr.+|scaffold|Pp', _line1[0]):
                    if _line1[2] == 'exon':
                        infor1 = _line1[8]
                        trid1 = re.match(r'.+transcript_id\ \"(TU\d+).+', infor1).group(1)
                        Tuinfor[trid1].append(allinfor[stainfo1])


for l in transline:
    fot.write(transline[l]+"\n")
    for m in Tuinfor[l]:
        fot.write(str(m)+"\n")


