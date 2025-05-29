#!/usr/bin/env python
# Description: processing a pileup-like sequence file and calculate base frequencies.
import argparse
import re

parser = argparse.ArgumentParser()
parser.add_argument("infile", help="input a sequence file", type=str)
parser.add_argument("outfile", help="output a file", type=str)
args = parser.parse_args()

frin = args.infile
fwout = args.outfile

pileup = open(frin, 'r', encoding='utf-8')
results = open(fwout, "w")

header="\t".join(["chrom","position","ref","totalnoN","A","T","C","G","N","A_fre","T_fre","C_fre","G_fre")
results.write(header + "\n")

for line in pileup.readlines():
    linesplit = line.split()
    chrom = linesplit[0]
    position = linesplit[1]
    ref = linesplit[2]
    total = linesplit[3]
    insert = re.findall(r'\+[0-9]+\w+',linesplit[4],re.IGNORECASE)
    deletion = re.findall(r'-[0-9]+\w+',linesplit[4],re.IGNORECASE)
    skip = re.findall('r<>',linesplit[4],re.IGNORECASE)
    
    list1=[]
    if insert:
        for i1 in insert:
            j1 = "".join(re.findall(r'\d+', i1))
            k1 = i1[0:int(j1) + len(j1) + 1]
            list1.append(k1)

        a1 = [];t1 = [];c1 = [];g1 = []
        for h1 in list1:
            h1 = h1.lower()
            f1 = h1.count("a")
            a1.append(f1)

            f2 = h1.count("t")
            t1.append(f2)

            f3 = h1.count("g")
            g1.append(f3)

            f4 = h1.count("c")
            c1.append(f4)
    else:
        a1 = [0, 0];t1 = [0, 0];c1 = [0, 0];g1 = [0, 0]

    list2 = []
    if deletion:
        for i2 in deletion:
            j2 ="".join(re.findall(r'\d+',i2))
            k2 = i2[0:int(j2) + len(j2)+1]
            list2.append(k2)

        a2 = [];t2 = [];c2 = [];g2 = []
        for h2 in list2:
            h2 = h2.lower()
            a = h2.count("a")
            a2.append(a)

            t = h2.count("t")
            t2.append(t)

            g = h2.count("g")
            g2.append(g)

            c = h2.count("c")
            c2.append(c)
    else:
        a2=[0,0];t2=[0,0];c2=[0,0];g2=[0,0]
    
    suma = sum(a1)+sum(a2)
    sumt = sum(t1)+sum(t2)
    sumg = sum(g1)+sum(g2)
    sumc = sum(c1)+sum(c2)
    
    A = str(linesplit[4].count("A") + linesplit[4].count("a")-len(re.findall(r'\^A',linesplit[4],re.IGNORECASE))-suma)
    T = str(linesplit[4].count("T") + linesplit[4].count("t")-len(re.findall(r'\^T',linesplit[4],re.IGNORECASE))-sumt)
    C = str(linesplit[4].count("C") + linesplit[4].count("c")-len(re.findall(r'\^C',linesplit[4],re.IGNORECASE))-sumc)
    G = str(linesplit[4].count("G") + linesplit[4].count("g")-len(re.findall(r'\^G',linesplit[4],re.IGNORECASE))-sumg)
    N = str(linesplit[4].count("*")+linesplit[4].count("<")+linesplit[4].count(">")-len(re.findall(r'\^[<>]',linesplit[4],re.IGNORECASE))-len(re.findall(r'\^\*',linesplit[4],re.IGNORECASE)))
    
    if ref == "a" or ref == "A":
        A = str(linesplit[4].count(".") + linesplit[4].count(",")-len(re.findall(r'\^[,.]',linesplit[4],re.IGNORECASE)))
    elif ref == "t" or ref == "T":
        T = str(linesplit[4].count(".") + linesplit[4].count(",")-len(re.findall(r'\^[,.]',linesplit[4],re.IGNORECASE)))
    elif ref == "c" or ref == "C":
        C = str(linesplit[4].count(".") + linesplit[4].count(",")-len(re.findall(r'\^[,.]',linesplit[4],re.IGNORECASE)))
    elif ref == "g" or ref == "G":
        G = str(linesplit[4].count(".") + linesplit[4].count(",")-len(re.findall(r'\^[,.]',linesplit[4],re.IGNORECASE)))
    else:
        A = "0";T = "0";C = "0";G = "0";N = "0";total= "0"

    if total == str(0) or int(A)+int(T)+int(C)+int(G)==0:
        A_fre = str('%.2f%%' % 0)
        T_fre = str('%.2f%%' % 0)
        C_fre = str('%.2f%%' % 0)
        G_fre = str('%.2f%%' % 0)
    else:
        total = str(int(A)+int(T)+int(C)+int(G))
        A_fre = str('%.2f%%' % (int(A) / int(total) * 100)) 
        T_fre = str('%.2f%%' % (int(T) / int(total) * 100))
        C_fre = str('%.2f%%' % (int(C) / int(total) * 100))
        G_fre = str('%.2f%%' % (int(G) / int(total) * 100))
        
    output = "\t".join([chrom, position, ref, total, A, T, C, G, N, A_fre, T_fre, C_fre, G_fre])
    results.write(output + "\n")
print("------Work Done------")
results.close()