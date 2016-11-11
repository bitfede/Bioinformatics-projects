#!/usr/bin/python

from Bio import SeqIO

#
# EXERCISE 1 - write code to parse info of human.fa
#

record = SeqIO.parse(open("human.fa"), "fasta")

#the next two lines print each sequence 
# for x in record:
# 	print x.id
# 	print x.seq

print "-----"

#
# EXERCISE 2 - parse infomrration of synonyms and location by parsing the protein-coding_gene.txt file
#

inputfile = open('protein-coding_gene.txt')
outputfile = open('results.fa', 'w')

#store lines in a list
fileInfo = inputfile.readlines()

for i in range(0, len(fileInfo)):
	fileInfo[i] = fileInfo[i].split('\t')



#this helps me to see how table is organized
# print "######"
# for x in range(0, len(fileInfo[0])):
# 	print("->" + fileInfo[0][x]) #these are the labels (ID, Synonims etc etc)
# 	print fileInfo[100][x]			 #in this case we look at the hundereth's line of the file
# print "#####"



#
# EXERCISE 3 - write a separate fasta file. In particular, the information
#              of each entry should be formatted by writing:
#							 > gene name|synonyms (separated by commas)|gene location.
#

#create a dictionary with sequence name as key and all the info as value
#let's look for matches from our human.fa
record = SeqIO.parse(open("human.fa"), "fasta")
for gene in record:
	print gene.id
	for info in fileInfo:
		# print gene.id + "MATCHING? " + info[1]
		if (gene.id == info[1]):
			outputfile.write(">" + gene.id + " | " + info[6] + " | " + info[7] + "\n")
			sequence = str(gene.seq)
			for x in range(0, len(sequence)):
				if (x % 60 == 0 and x != 0):
					outputfile.write("\n")
				outputfile.write(sequence[x])
			outputfile.write("\n")

print "-----------------"

