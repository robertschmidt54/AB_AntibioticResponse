# -*- coding: utf-8 -*-
"""
Created on Thu Sep  9 09:54:38 2021

@author: zenak
"""

gffFile = "../RefGenome/GCF_001593425.2_ASM159342v2_genomic.gff"
regionOfInterest = "gene"

with open(gffFile, 'r') as infile:
    with open("../gffID_to_GeneID.txt", 'w+') as outfile:
        outfile.write("GeneID\tGenedbRef\tGeneName\tLocus_tag\tStart\tEnd\n")
        gffLines = infile.readlines()
        for line in gffLines:
            if line.startswith("#"):
                continue
            
            linarray = line.split('\t')
            if linarray[2] != regionOfInterest:
                continue
            
            start = linarray[3]
            end = linarray[4]            
            infoString = linarray[8]
                        
            infoArray = infoString.split(';')
            #print(infoArray)
            infoDict = dict(elem.strip().split('=', 1) for elem in infoArray)
            print(infoDict)
            geneID = infoDict['ID']
            GeneRefNum = infoDict['gbkey']
            GeneName = infoDict['Name']
            Locus_tag = infoDict['locus_tag']
            #old_locus_tag = infoArray[6].split("=")[1]
            
            outrow = geneID+"\t"+GeneRefNum+"\t"+GeneName+"\t"+Locus_tag+"\t"+start+"\t"+end+"\n"
            outfile.write(outrow)
    outfile.close()
infile.close()