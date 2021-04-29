#!/usr/bin/python3
# coding = utf-8
# -*- coding: utf-8 -*-


path = "work_dir"  

file1 = open(path + "Include_reads.txt", 'w')
file2 = open(path + "skipped_reads.txt", 'w')
file3 = open(path + "Include_level.txt", 'w')
file4 = open(path + "rMATS_dpsi.txt", 'w')

samples = ['sample1',      ####sample names，ordered
           'sample2',
           '...'
           ]

### print titles
titles1 = []
titles2 = []
titles3 = []
titles4 = []

for sample in samples:
    Iread = "IncludeReads_"+str(sample)
    Sread = "SkippedReads_"+str(sample)
    Ilevel = "IncludeLevel_"+str(sample)
    titles1.append(Iread)
    titles2.append(Sread)
    titles3.append(Ilevel)

for i in range(len(samples)):
    for j in range(len(samples)):
        if i < j:
            dpsi = "dpsi_%s_vs_%s" % (samples[i], samples[j])
            pvalue = "pvalue_%s_vs_%s" % (samples[i], samples[j])
            fdr = "FDR_%s_vs_%s" % (samples[i], samples[j])
            dpsis = "\t".join([dpsi, pvalue, fdr])
            titles4.append(dpsis)

file1.write("EventID\t"+"\t".join(titles1)+"\n")
file2.write("EventID\t"+"\t".join(titles2)+"\n")
file3.write("EventID\t"+"\t".join(titles3)+"\n")
file4.write("EventID\t"+"\t".join(titles4)+"\n")

eventids = {}
names = []

##### intersect of AS events
for i in range(len(samples)):
    for j in range(len(samples)):
        if i < j:
            name = "%s_vs_%s" %(samples[i], samples[j])
            names.append(name)
            file = "%s%s_vs_%s/Combine.MATS.JC.txt" %(path,samples[i],samples[j])
            eventids[name] = []
            f = open(file, 'r')
            f.readline()
            for line in f:
                line = line.strip()
                _line = line.split("\t")
                eid = ";".join([str(_line[23]), str(_line[1]).strip("\""), ";".join(_line[5:11])])
                if eid not in eventids[name]:
                    eventids[name].append(eid)
# intersects ids
list_ids = []
for name in eventids:
    list_ids.append(set(eventids[name]))

set = set.intersection(*list_ids)
Eids = list(set)  #获得所有基因id

################################## output
dic_Ir = {}
dic_Sr = {}
dic_Il = {}
for i in range(len(samples)-1):
    if i == len(samples)-2:
        file = "%s%s_vs_%s/Combine.MATS.JC.txt" % (path, samples[i], samples[i+1])
        f = open(file, 'r')
        f.readline()
        for line in f:
            line = line.strip()
            _line = line.split("\t")
            eid = ";".join([str(_line[23]), str(_line[1]).strip("\""), ";".join(_line[5:11])])
            Ir1 = str(sum(map(int, _line[12].split(','))))
            Ir2 = str(sum(map(int, _line[14].split(','))))
            Sr1 = str(sum(map(int, _line[13].split(','))))
            Sr2 = str(sum(map(int, _line[15].split(','))))
            try:
                Il1 = str(sum(map(float, _line[20].split(',')))/3)
                Il2 = str(sum(map(float, _line[21].split(',')))/3)
            except ValueError:
                Il1, Il2 = "NA", "NA"
            if eid in Eids:
                if eid in dic_Ir:
                    dic_Ir[eid].append(Ir1)
                    dic_Ir[eid].append(Ir2)
                    dic_Sr[eid].append(Sr1)
                    dic_Sr[eid].append(Sr2)
                    dic_Il[eid].append(Il1)
                    dic_Il[eid].append(Il2)
                else:
                    dic_Ir[eid] = [Ir1, Ir2]
                    dic_Sr[eid] = [Sr1, Sr2]
                    dic_Il[eid] = [Il1, Il2]
    else:
        file = "%s%s_vs_%s/Combine.MATS.JC.txt" % (path, samples[i], samples[i+1])
        f = open(file, 'r')
        f.readline()
        for line in f:
            line = line.strip()
            _line = line.split("\t")
            eid = ";".join([str(_line[23]), str(_line[1]).strip("\""), ";".join(_line[5:11])])
            Ir = str(sum(map(int, _line[12].split(','))))
            Sr = str(sum(map(int, _line[13].split(','))))
            try:
                Il = str(sum(map(float, _line[20].split(',')))/3)
            except ValueError:
                Il = "NA"
            if eid in Eids:
                if eid in dic_Ir:
                    dic_Ir[eid].append(Ir)
                    dic_Sr[eid].append(Sr)
                    dic_Il[eid].append(Il)
                else:
                    dic_Ir[eid] = [Ir]
                    dic_Sr[eid] = [Sr]
                    dic_Il[eid] = [Il]
    f.close()

dic_dpsi = {}
dic_Pv = {}
dic_Fdr = {}
dic_dpsis = {}
for i in range(len(samples)):
    for j in range(len(samples)):
        if i < j:
            file = "%s%s_vs_%s/Combine.MATS.JC.txt" % (path, samples[i], samples[j])
            f = open(file, "r")
            f.readline()
            for line in f:
                line = line.strip()
                _line = line.split()
                eid = ";".join([str(_line[23]), str(_line[1]).strip("\""), ";".join(_line[5:11])])
                dpsi = _line[22]
                Pv = _line[18]
                Fdr = _line[19]
                dpsis = "\t".join(map(str, [dpsi, Pv, Fdr]))
                # print(dpsis)
                if eid in dic_dpsis:
                    dic_dpsi[eid].append(dpsi)
                    dic_Pv[eid].append(Pv)
                    dic_Fdr[eid].append(Fdr)
                    dic_dpsis[eid].append(dpsis)
                else:
                    dic_dpsi[eid] = [dpsi]
                    dic_Pv[eid] = [Pv]
                    dic_Fdr[eid] = [Fdr]
                    dic_dpsis[eid] = [dpsis]

for eid in dic_Ir:
    file1.write(eid+"\t"+"\t".join(dic_Ir[eid])+"\n")
    file2.write(eid+"\t"+"\t".join(dic_Sr[eid])+"\n")
    file3.write(eid+"\t"+"\t".join(dic_Il[eid])+"\n")
    file4.write(eid+"\t"+"\t".join(dic_dpsis[eid])+"\n")
