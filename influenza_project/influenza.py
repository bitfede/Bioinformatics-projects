import tstatistics as stat
import pprint

filename_con = "Influenza_con_expr.txt"
filename_inf = "Influenza_inf_expr.txt"

f_con = open(filename_con)
f_inf = open(filename_inf)

con_expr = {}
inf_expr = {}

f_con.readline() #we do not care about the first lines
f_inf.readline() 

for line in f_con:
	line = line.rstrip("\n")
	line = line.split(" ")
	con_expr[line[0]] = map(float, line[1:len(line)]) # length of data is 25

for line in f_inf:
	line = line.rstrip("\n")
	line = line.split(" ")
	inf_expr[line[0]] = map(float, line[1:len(line)])


ttest = {}
pvalues = {}

for k in con_expr:
	ttest[k] = stat.tutest(con_expr[k] , inf_expr[k])
	pvalues[k] = ttest[k][1]

pp = pprint.PrettyPrinter(indent=4)

pp.pprint(pvalues)