import sys

def minFunct(starterString, fileName):
	fh = open(fileName)
	fileLines = fh.readlines()
	totalTime = 0
	minParticles = 1000

	for line in fileLines:
		if starterString in line:
			indexParticles = line.index(starterString)
			endStrIndex = line.index(",", indexParticles)
			numParticles = int(line[indexParticles:endStrIndex].replace(starterString, "").strip())
			if numParticles < minParticles:
				minParticles = numParticles

	return minParticles

def invalidParticles(starterString, fileName):
	fh = open(fileName)
	fileLines = fh.readlines()
	totalTime = 0
	numTimesInvalid = 0

	for line in fileLines:
		if starterString in line:
			numTimesInvalid += 1
			
	return numTimesInvalid


if __name__ == "__main__":
	fileName = "log5.txt"
	
	if len(sys.argv) > 1:
		fileName = sys.argv[1]

	print("Min particles at any point of time: {0}".format(minFunct("particles = ",  fileName)))
	print("Number of reiterations done: {0}".format(invalidParticles("Not enough valid particles - Resampling!!!",  fileName)))

	