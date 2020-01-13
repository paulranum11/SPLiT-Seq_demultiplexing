#!/usr/bin/python
#Written by Dipankar / Dumaatravaie on 28th August, 2019 (version python 3 )
#Collapse OligoDT and RandomHexamer Barcodes from the same well
import sys
import os 
import os.path

#Container declaration
RanHex_BARCODES={}
Odt_BARCODES={}

index_RanHex={}
index_Odt={}

#input and output files and directory
RanHex_file="./RanHex.txt"
Odt_file="./OligoDT.txt"
OUTPUT_DIR="results"

#read all the files and stores the file names in a container
files=os.listdir(OUTPUT_DIR)

# Read in the RandomHexamer and OligoDT RT primer barcode sequences and store them
# in the relative containers 
inf=open(RanHex_file,"r")
counter=0
for l in inf:
	RanHex_BARCODES[l.strip()]=[]
	index_RanHex[counter]=l.strip()
	counter+=1
inf.close()

counter=0
inf=open(Odt_file,"r")
for l in inf:
	Odt_BARCODES[l.strip()]=[]
	index_Odt[counter]=l.strip()
	counter+=1
inf.close()

# Walk through all RanHex and Odt barcodes in pairs
# These pairs reflect the barcode pairs added to the first 48 wells of the ROUND1 plate

for file in files:
	p=file.split(".")[0].split("-")
	if p[0] in RanHex_BARCODES:
		RanHex_BARCODES[p[0]].append(p[1]+"-"+p[2])
	else:
		if p[0] in Odt_BARCODES:
			Odt_BARCODES[p[0]].append(p[1]+"-"+p[2])
counter=0
for i in range(0,len(index_RanHex)):
	Odt_hit={}
	RanHex_hit={}

	p=RanHex_BARCODES[index_RanHex[i]]
	for k in range(0,len(p)):
		RanHex_hit[p[k]]=1
		
	p=Odt_BARCODES[index_Odt[i]]
	for k in range(0,len(p)):
		Odt_hit[p[k]]=1
	
	for k,v in RanHex_hit.items():
		if k in Odt_hit:
			R=index_RanHex[i]+"-"+k+".fastq"
			O=index_Odt[i]+"-"+k+".fastq"
			print ("collapsing \n"+O+" and \n"+R)
			os.system("cat ./%s/%s >> ./%s/%s"%(OUTPUT_DIR,R,OUTPUT_DIR,O))
			os.system("rm ./%s/%s"%(OUTPUT_DIR,R))
			counter+=1

print ("\nTotal Pairs of Files Collapsed = "+str(counter))
print ("\nCleaning up \n" )		
	
	
