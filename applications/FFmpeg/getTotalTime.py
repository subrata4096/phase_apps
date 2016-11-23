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
	fileName = "timeFinal.txt"
	
	if len(sys.argv) > 1:
		fileName = sys.argv[1]

	print("Total deflate time = {0}".format(getTime("Time taken by deflate = ",  fileName)))
	print("Total vblur time = {0}".format(getTime("Time taken by vblur = ",  fileName)))
	print("Total embross time = {0}".format(getTime("Time taken by embross = ",  fileName)))
	print("bench: utime = {0}".format(getTime("bench: utime=",  fileName)))
	# print("Total h815scale time = {0}".format(getTime("Time taken by h815scale = ",  fileName)))
	# print("Total swscale time = {0}".format(getTime("Time taken by swscale = ",  fileName)))
	# print("Total hscale time = {0}".format(getTime("Time taken by hscale = ",  fileName)))