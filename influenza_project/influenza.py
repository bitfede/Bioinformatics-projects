import tstatistics as stat
import pprint

######################################
# 	1) peform t-test                 #
######################################

# set filenames
filename_con = "Influenza_con_expr.txt"
filename_inf = "Influenza_inf_expr.txt"

# open files to read
f_con = open(filename_con)
f_inf = open(filename_inf)

# create two dictionaries that will contain the data for each gene
con_expr = {}
inf_expr = {}

#reading the first line on both files but we are not saving it
f_con.readline() 
f_inf.readline() 

#the next two for loops will parse the data and store it in their respective dictionaries

#control data
for line in f_con:
	line = line.rstrip("\n")
	line = line.split(" ")
	con_expr[line[0]] = map(float, line[1:len(line)]) # length of data is 25

#influenza data
for line in f_inf:
	line = line.rstrip("\n")
	line = line.split(" ")
	inf_expr[line[0]] = map(float, line[1:len(line)])

#create a dictionary that will hold the ttest results
ttest = {}

#create a dictionary that will hold just the p-values
pvalues = {}

for k in con_expr:
	#perform ttest
	ttest[k] = stat.tutest(con_expr[k] , inf_expr[k])
	#extract p-value
	pvalues[k] = ttest[k][1]

#pretty printer, used to write to my .txt files
pp = pprint.PrettyPrinter(indent=4)


######################################
# 	2) get p-value and correct it    #
######################################

# pp.pprint(pvalues) #  <--- uncomment this to print to stdout  (I piped this output to get my two txt files)
# pp.pprint(pvalues) #  same as the line before

# NOW we need to CORRECT the p values

# let's perform a Bonferroni correction
bonf_alpha = 0.05													#our alpha value
bonf_N     = 60													#(number of statistical tests)

bonf_val 	 = bonf_alpha / bonf_N

print bonf_val		#DEBUG checking corrected value

# now we will save in this new dictionary ONLY the genes that we are interested in
filtered_list = {}

#looping thru the genes and filtering them by P-value
for k in pvalues:
	if pvalues[k] < bonf_val:
		filtered_list[k] = pvalues[k]


######################################
# 	3) list the output               #
######################################

#create formatted output
outputline = "ID P-VALUE\n"
for k in filtered_list:
	outputline = outputline + k + " " + str(filtered_list[k]) + "\n"

print outputline  #piped it straight to gene_over-expressed.txt