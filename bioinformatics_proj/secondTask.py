#!/usr/local/bin/python

import Bio
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
from Bio import SeqIO

#
# EXERCISE 1 - set up a Python code that parses the human sequences
#

record = SeqIO.parse(open("human.fa"), "fasta")

#opening output file
outputfile = open('BLASTresults.txt', 'w')

# decided to use xxxxxxxx beacuse ...

#
# EXERCISE 2, 3 , 4 -  prints the human sequence ID, the mouse ID of the most
#											 similar homolog and the alignment, E-value and bitscore
#	 									 	 to a file of your choice.
#



for gene in record:
	result_handle = NCBIWWW.qblast("blastp", "refseq_protein", gene.seq, matrix_name="BLOSUM80", expect=0.001, format_type='XML')
	result = Bio.Blast.NCBIXML.read(result_handle)
	for alignment in result.alignments:
		for hsp in alignment.hsps:
			if "Mus musculus" in alignment.title:  #we are looking for mouse
				print "HUMAN GENE:     " +  gene.name + "\n"
				outputfile.write("HUMAN GENE:     " +  gene.name + "\n")
				print "MOUSE ID:			 " +  alignment.title + "\n"
				outputfile.write("MOUSE ID:			 " +  alignment.title + "\n")

				print "E-value:        " +  str(hsp.expect) + "\n"
				outputfile.write("E-value:        " +  str(hsp.expect) + "\n")
				print "BITSCORE:       " +  str(hsp.bits)
				outputfile.write("BITSCORE:       " +  str(hsp.bits) + "\n")
				outputfile.write("----------------------------------\n")
				print "-------------------------------------"
