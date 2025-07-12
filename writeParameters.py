#require python 3.6 or above
# usage
# python writeParameters.py -i <paramExampleFile> -p <paraList> -n <run number>
import csv
from optparse import OptionParser

parser = OptionParser()
parser.add_option("-p", "--parameters", dest="paraList",\
					help="csv file for all the parameter combinations", metavar="FILE")
parser.add_option("-i", "--input", dest="paramExampleFile",\
					help="input template file", metavar="FILE")
parser.add_option("-n", "--number", dest="outNumber",type = "int",\
					help="the run number")
parser.add_option("-x", "--prefix", dest="prefix",type = "string",\
					help="prefix for output filenames")
                  
(options, args) = parser.parse_args()


if __name__ == "__main__":
	#step 1 read in parameter combination file
	allParamFile = open(options.paraList,newline = '')
	allParam = csv.DictReader(allParamFile)
	
	#step 2 read in input file template
	prototype = open(options.paramExampleFile,"r").read()
	
	#iterate through the list and find the correct run number
	for row in allParam:
		if row['num'] == str(options.outNumber):
			if row['LOAD_FROM_CHECKPOINT'] == 'FALSE':
				row['SAMPLE_DB_FILENAME'] = options.prefix + "_" + row['num'] + "_beforeIRS_sd.sqlite"
				row['CHECKPOINT_SAVE_FILENAME'] = options.prefix + "_" + row['num'] + "_beforeIRS_cp.sqlite"
				row['CHECKPOINT_LOAD_FILENAME'] = ''
			else:
				row['SAMPLE_DB_FILENAME'] = options.prefix + "_" + row['num'] + "_afterIRS_sd.sqlite"
				row['CHECKPOINT_SAVE_FILENAME'] = options.prefix + "_" + row['num'] + "_afterIRS_cp.sqlite"
				row['CHECKPOINT_LOAD_FILENAME'] = options.prefix + "_" + row['num'] + "_beforeIRS_cp.sqlite"
			
			out = open(options.prefix + "_" + row['num'] + "_input.py", "w")
			out.write(prototype.format(**row))
			out.close()
			break
	

