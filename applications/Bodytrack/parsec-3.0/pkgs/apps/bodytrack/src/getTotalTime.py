import sys

def getTime(starterString, fileName):
	fh = open(fileName)
	fileLines = fh.readlines()
	totalTime = 0

	for line in fileLines:
		if starterString in line:
			totalTime += float(line.split(starterString)[1].replace("s", ""))

	return totalTime


if __name__ == "__main__":
	fileName = "log5.txt"
	
	if len(sys.argv) > 1:
		fileName = sys.argv[1]

	print("Total Time taken by loop inner error inside = {0}".format(getTime("Time taken by loop inner error inside = ",  fileName)))
	print("Total Time taken by loop outer error inside = {0}".format(getTime("Time taken by loop outer error inside = ",  fileName)))
	print("Total Time taken by loop inner error edge = {0}".format(getTime("Time taken by loop inner error edge = ",  fileName)))
	print("Total Time taken by loop outer error edge = {0}".format(getTime("Time taken by loop outer error edge = ",  fileName)))
	print("Total Time taken by loop inner inside error  = {0}".format(getTime("Time taken by loop inner inside error  = ",  fileName)))
	print("Total Time taken by loop pf outer inside error = {0}".format(getTime("Time taken by loop pf outer inside error = ",  fileName)))
	print("Total Time taken by edge error loop 1 = {0}".format(getTime("Time taken by edge error loop 1 = ",  fileName)))
	print("Total Time taken by loop edge error loop 2 = {0}".format(getTime("Time taken by loop edge error loop 2 = ",  fileName)))
	print("Total Time taken by loop pf update = {0}".format(getTime("Time taken by loop pf update = ",  fileName)))
	