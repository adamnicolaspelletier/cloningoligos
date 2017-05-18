#!/usr/bin/env python

import primer3
import pandas as pd
from Bio.Restriction import *
from Bio import SeqIO
from Bio.Alphabet.IUPAC import IUPACAmbiguousDNA
amb = IUPACAmbiguousDNA()
from Bio.Align.Applications import MuscleCommandline
from Bio import AlignIO
from StringIO import StringIO

###############################################################################################################################
###############################################################################################################################
## User Input

cdnafile = "cdna_seq_test.txt"
ccdsfile = "ccds_seq_test.txt"

anneal = 21 		#write the number of nucleotides in the annealing portion of the oligo. 21+ recommended
deltaTM = 5			#acceptable temperature diffference between F and R oligos
minimumTM = 62		# Minimum acceptable Tm for any oligo. 55 is optimal. 
maximumTM = 65		# Minimum acceptable Tm for any oligo. 75 is a big max, under 70 would be optimal.

fivepr = ""  # write the sequence you want to add in the 5' portion of the restriction enzyme site. 6 nucleotides recommended, and mostly A and T recommended
threepr = "" # write the sequence you want to add in the 3' portion of the restriction enzyme site. 6 nucleotides recommended, and mostly A and T recommended
kozak = "ACC"





###############################################################################################################################
###############################################################################################################################

def fastaconvert(fastalist):  #convert a conventional fasta file into a list of IDs and whole sequences (merges the 50 characters per line)
	a = ""
	x = ""
	z = []
	for i in fastalist:
		if ">" in i:
			a = a.replace("\n","")
			z.append(a)
			x += i
			x = x.replace("\n","")
			z.append(x)
			a = ""
			x = ""
		else:
			a += i
	del z[0]
	return z

# def idsplitter(ID):
# 	ID.pop()
# 	realID = []
# 	for i in ID:
# 		name = str(i)
# 		idsplit = name.split("|")
# 		realID.append(idsplit[2])
# 	return realID

# sitelist = []
# enzlistname = ""


# for i in enzlist:
# 	sitelist.append(i.site)
# 	enzlistname += str(i)



# seq = ""
# seqlist = []


# ID = []
# idsplit = []
# name = ""


# fastafinal = fastaconvert(fastalist)

# a = ""
# x = []    #list of gene ID
# v = ""
# z = []	#lsit of fasta sequence one-liners

# for i in fastafinal:
# 	if ">" in i:
# 		a += i
# 		a = a.replace("\n","")
# 		x.append(a) #creates a list of ID
# 		a = ""

# 	else:
# 		v += i
# 		v = v.replace("\n","")
# 		z.append(v) #creates a list of sequences
# 		v = ""
# y = idsplitter(x)
# print y

def reversecomp(rprimsequence): ## make a complement version of the sequence, and reverse it so it has the proper orientation
	a = ""
	tempzrev = rprimsequence
	tempzrev = tempzrev.replace("T","X")
	tempzrev = tempzrev.replace("A","T")
	tempzrev = tempzrev.replace("X","A")
	tempzrev = tempzrev.replace("C","Y")
	tempzrev = tempzrev.replace("G","C")
	tempzrev = tempzrev.replace("Y","G")
	templist = list(tempzrev)
	templist.reverse()
	for i in templist:
		a += i
	return a

def fseqgen(seq, anneal):  #Generates a Forward primer based on a complete fasta string, towards the beggining for the length of the anneal parameter
	ftemp = ""
	seqlist = list(seq)
	flist = seqlist[0:anneal]
	for item in flist:
		ftemp += item
	return ftemp

def rseqgen(seq, anneal):  #Generates a Reverse primer based on a complete fasta string, towards the end for the length of the anneal parameter
	rtemp = ""
	seqlist = list(seq)
	rlist = seqlist[len(seqlist)-anneal:len(seqlist)]
	for item in rlist:
		rtemp += item
	rprim = reversecomp(rtemp)
	return rprim

def fprimergen(seq, anneal, sitelist, fivepr, kozak):  #Generates a Forward primer based on a complete fasta string, towards the beggining for the length of the anneal parameter
	ftemp = ""
	seqlist = list(seq)
	flist = seqlist[0:anneal]
	for item in flist:
		ftemp += item
	fprim = fivepr + sitelist[0] + kozak + ftemp
	return fprim

def rprimergen(seq, anneal, sitelist, threepr):  #Generates a Reverse primer based on a complete fasta string, towards the end for the length of the anneal parameter
	rtemp = ""
	seqlist = list(seq)
	rlist = seqlist[len(seqlist)-anneal:len(seqlist)]
	for item in rlist:
		rtemp += item
	rprimint = rtemp + sitelist[1] + threepr
	rprim = reversecomp(rprimint)
	return rprim


def primer3cloning(ccdsfile, cdnafile):
	with open("internal_files/fastaoutput.txt", "w") as tempoutput:	
		tempoutput.write("")

	for record in SeqIO.parse(ccdsfile, "fasta"):
		with open("internal_files/fastaoutput.txt", "a") as tempoutput:
				ccds_seq = record.seq		
				tempoutput.write('>' + record.id + '\n'+ str(ccds_seq.upper()) + '\n')
		ccdsid = record.id

	for recordc in SeqIO.parse(cdnafile, "fasta"):
		with open("internal_files/fastaoutput.txt", "a") as tempoutput:
				cdna_seq = recordc.seq		
				tempoutput.write('>' + recordc.id + '\n'+ str(cdna_seq.upper()))
		cdnaid = recordc.id

	muscle_cline = MuscleCommandline(input="internal_files/fastaoutput.txt") 
	stdout, stderr = muscle_cline()
	muscle_cline.gapopen = -10.0
	align = AlignIO.read(StringIO(stdout), "fasta")
	AlignIO.write([align], "internal_files/muscleoutput.txt", "fasta")

	for recordf in SeqIO.parse("internal_files/muscleoutput.txt", "fasta"):
		if recordf.id == ccdsid:
			ccds = list(str(recordf.seq))

		count = 0
		gene = 0
		for i in ccds:
			if i == "-" and gene == 0:
				count += 1
			elif i == "-":
				pass
			else:
				gene += 1
	bleh = count + gene
	print ccds[bleh-3:bleh]
	print count, gene
		



primer3cloning(ccdsfile,cdnafile)