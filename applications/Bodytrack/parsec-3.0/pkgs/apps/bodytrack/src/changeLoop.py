#Manish's script
import re
import os
import shutil
import subprocess
import time

#Initialize variables, regex
MAXLOOPSKIP = 7
OUTPUTDIR = "./approx_vid_out"

errorDict = {}
timeDict = {}
errorPattern = re.compile(".*\*\*\*\sError:\s(\d+\.\d+).*")
timePattern = re.compile("real\s*(\d*)m(\d*.\d*)s")

# APPROXTYPESARRAY = ["LOOPCHANGE_IEINSIDE_OUTER", "LOOPCHANGE_IEINSIDE_INNER", "LOOPCHANGE_IEEDGE_OUTER", "LOOPCHANGE_IEEDGE_INNER","LOOPCHANGE_INSIDEERROR_OUTER","LOOPCHANGE_INSIDEERROR_INNER", "LOOPCHANGE_EDGEERROR_1","LOOPCHANGE_EDGEERROR_2"]   #Removed update one
APPROXTYPESARRAY = ["LOOPCHANGE_INSIDEERROR_OUTER", "LOOPCHANGE_INSIDEERROR_INNER", "LOOPCHANGE_EDGEERROR_1","LOOPCHANGE_EDGEERROR_2"]   #Removed update one
# PHASETYPESARRAY = ["PHASE_INSIDEERROR_OUTER", "PHASE_INSIDEERROR_INNER", "PHASE_EDGEERROR_1","PHASE_EDGEERROR_2"]   #Removed update one
LOGFILENAME = "./approx_log.txt"
PHASETOAPPROXIMATE = 0

# TIMESTRINGS = ["Time taken by loop inner error inside = ", "Time taken by loop outer error inside = ", "Time taken by loop inner error edge = ", "Time taken by loop outer error edge = ", "Time taken by loop inner inside error  = ", "Time taken by loop pf outer inside error = ", "Time taken by edge error loop 1 = ", "Time taken by loop edge error loop 2 = ", "Time taken by loop pf update = ",]     
TIMESTRINGS = ["Time taken by loop pf outer inside error = ", "Time taken by loop inner inside error  = " , "Time taken by edge error loop 1 = ", "Time taken by loop edge error loop 2 = "]     
LOOPSTRINGS = ["while loop minvalid, iter = ", "for loop annealing layers, iter = ", "Frames "]
ANNEALINGLAYERS = [5, 100, 25, 75, 3]
MINPARTICLES    = [4000, 5, 1000, 500, 8000]
NUMFRAMES = 5
PHASEOPTIONS = [1, 2, 3, 4, -1]
#PHASEOPTIONS = [1, 2, 3, 4, 5, 6, 7, 8, -1]
#PHASEOPTIONS = [1, 2, -1]
NUMPHASES             = 4
MAX_ANNEALING_LAYERS  = 63
MIN_ANNEALING_LAYERS  = 3
STEP_ANNEALING_LAYERS = -10

MAX_MIN_PARTICLES     = 8000
MIN_MIN_PARTICLES     = 500
STEP_MIN_PARTICLES    = -500

# returns the time/iterations as a float using the starterString
def getTotalVal(starterString, fileStr, typeInt=0, timeVal=0):
    # fh = open(fileName)
    # fileLines = fh.readlines()
    fileLines = fileStr.split("\n")
    totalVal = 0
    numSampled = 0

    for line in fileLines:
        if starterString in line:
            numSampled += 1
            if typeInt == 1:
                totalVal += int(line.split(starterString)[1].replace("s", ""))
                continue

            totalVal += float(line.split(starterString)[1].replace("s", ""))

    if timeVal:
        # print("time for {0}, time={1}, numSampled={2}, finalVal={3}".format(starterString, totalVal, numSampled, totalVal/numSampled))
        totalVal = totalVal / numSampled
    return totalVal


# Runs bodytrack and returns dictionary containing time take for each string and different iteration values
def runBodytrack(annealingLayers, minParticles):
    outputValDict = {} 
    outputValDict["time"] = [] # First key in dictionary is an array of all the time values gotten from the output
    outputValDict["iter"] = [] # Second key in dictionary is an array of all the iteration values it will need

    overallTimeTaken = -1
    #  ../obj/amd64-linux.gcc-serial/TrackingBenchmark/bodytrack ../inputs/sequenceB_261 4 261 1000 5 0 1
    #args -> input, numberCams, frames, min particles, annealing layers, ..,..
    # print(minParticles)
    # print(annealingLayers)
    runProc = subprocess.Popen(['../obj/amd64-linux.gcc-serial/TrackingBenchmark/bodytrack',  '../inputs/sequenceB_261', '4', str(NUMFRAMES) , str(minParticles), str(annealingLayers), '0', '1'],stdout=subprocess.PIPE,stderr=subprocess.PIPE)
    output, error = runProc.communicate()
    # print(output)
    overallTimeTaken = getTotalVal("Number of HW instructions Used=", output, typeInt=1)

    for timeIndex, timeString in enumerate(TIMESTRINGS):
        outputValDict["time"].append(getTotalVal(timeString, output, timeVal=1))

    for loopStrIndex, loopString in enumerate(LOOPSTRINGS):
        outputValDict["iter"].append(getTotalVal(loopString, output, typeInt=1))

    outputValDict["time"].append(overallTimeTaken)

    # timeMatchObj = timePattern.search(output)

    # if(timeMatchObj):
    #     timeTaken = float(timeMatchObj.group(1)) * 60 + float(timeMatchObj.group(2))
    # else:
    #     print("Failed to get time amount for loopskip " + str(i))

    return outputValDict


# Runs qos script and returns error
def getError(inputFileName, referenceFileName):
    error = -1

    #Comparison script
    p = subprocess.Popen(['python', 'qos.py', inputFileName, referenceFileName],stdout=subprocess.PIPE,stderr=subprocess.PIPE)
    output, error = p.communicate()
    #Get error
    errorMatchObj = errorPattern.match(output)
    if(errorMatchObj):
        error = float(errorMatchObj.group(1))
    else:
        print("Failed to get error amount")

    return error


def runBodytrackAndGetError(annealingLayers, minParticles):
    outputValDict = runBodytrack(annealingLayers, minParticles)
    error = getError('../inputs/sequenceB_261/poses.txt', OUTPUTDIR + '/correctoutput.txt')

    return outputValDict, error


# initializes all necessary environment variables
def initEnvironmentVariables():    
    os.environ["LOOPITERSKIP_PFUPDATE"]              ='1'
    os.environ["LOOPCHANGE_PFUPDATE"]                ='1'

    os.environ["LOOPITERSKIP_IEINSIDE_OUTER"]        ='1'
    os.environ["LOOPCHANGE_IEINSIDE_OUTER"]          ='1'

    os.environ["LOOPITERSKIP_IEINSIDE_INNER"]        ='1'
    os.environ["LOOPCHANGE_IEINSIDE_INNER"]          ='1'

    os.environ["LOOPITERSKIP_IEEDGE_OUTER"]          ='1'
    os.environ["LOOPCHANGE_IEEDGE_OUTER"]            ='1'

    os.environ["LOOPITERSKIP_IEEDGE_INNER"]          ='1'
    os.environ["LOOPCHANGE_IEEDGE_INNER"]            ='1'

    os.environ["LOOPCHANGE_INSIDEERROR_OUTER"]       ='1'
    os.environ["LOOPITERSKIP_INSIDEERROR_OUTER"]     ='1'

    os.environ["LOOPCHANGE_INSIDEERROR_INNER"]       ='1'
    os.environ["LOOPITERSKIP_INSIDEERROR_INNER"]     ='1'

    os.environ["LOOPCHANGE_EDGEERROR_1"]             ='1'
    os.environ["LOOPITERSKIP_EDGEERROR_1"]           ='1'

    os.environ["LOOPCHANGE_EDGEERROR_2"]             ='1'
    os.environ["LOOPITERSKIP_EDGEERROR_2"]           ='1'

    os.environ["PHASETOAPPROXIMATE"]                 ='1'
    os.environ["TOTALPHASES"]                        =str(NUMPHASES)
    os.environ["TOTALITERATIONS"]                    ='1'
    os.environ["PARTICLE_FACTOR_CHANGE"]             ='1'

# write to log file
def writeToLog(text):
    print(text.strip("\n"))
    fh = open(LOGFILENAME, "a")
    fh.write(text)
    fh.close()
    return


# Runs base case, copies it to reference file location and returns time taken
def getBaseValue(annealingLayers, minParticles):
    outputValDict = runBodytrack(annealingLayers, minParticles)
    os.system("cp ../inputs/sequenceB_261/poses.txt " + OUTPUTDIR + "/correctoutput.txt")

    return outputValDict

# Initialize everything
def init():
    #Make tempdirectoy
    if not os.path.isdir(OUTPUTDIR):
        os.mkdir(OUTPUTDIR)
    initEnvironmentVariables()


def runIndividualApproximations():
    init() 

    for inputIndex in range(len(ANNEALINGLAYERS) - 1,  len(ANNEALINGLAYERS)):       
        annealingLayers = ANNEALINGLAYERS[inputIndex]
        minParticles    = MINPARTICLES[inputIndex]
        initValBaseDict = getBaseValue(annealingLayers, minParticles)
        writeToLog("###!!!####INIT {0}, {1}, {2}, {3}, {4}\n".format(annealingLayers, minParticles, initValBaseDict["iter"][0], initValBaseDict["iter"][1], initValBaseDict["iter"][2])) 
        for approxIndex, approxType in enumerate(APPROXTYPESARRAY):
            initEnvironmentVariables()
            for approxLevel in range(2, MAXLOOPSKIP):
                os.environ[approxType] = str(approxLevel)
                outputValDict, error = runBodytrackAndGetError(annealingLayers, minParticles)
                localSpeedUp = initValBaseDict["time"][approxIndex] / outputValDict["time"][approxIndex]
                globalSpeedUp = initValBaseDict["time"][-1] / outputValDict["time"][-1]
                approxLevelToPrint = 1 / float(approxLevel)
                if "LOOPCHANGE" in approxType:
                    approxLevelToPrint = 1 - approxLevelToPrint

                writeToLog("{0}, {1}, {2}, {3}, {4}, {5}, {6}, {7}, {8}, {9}\n".format(annealingLayers, minParticles, approxType, approxLevelToPrint, localSpeedUp, globalSpeedUp, error, outputValDict["iter"][0], outputValDict["iter"][1], outputValDict["iter"][2]))


def runDifferentParams():
    global LOGFILENAME

    LOGFILENAME = "./combinedApproximationLog_diffParams.txt"

    initValBaseDict = getBaseValue(MAX_ANNEALING_LAYERS, MAX_MIN_PARTICLES )
    writeToLog("###!!!####INIT {0}, {1}, {2}, {3}, {4}\n".format(MAX_ANNEALING_LAYERS, MAX_MIN_PARTICLES, initValBaseDict["iter"][0], initValBaseDict["iter"][1], initValBaseDict["iter"][2])) 

    for annealingLayer in range(MAX_ANNEALING_LAYERS, MIN_ANNEALING_LAYERS, STEP_ANNEALING_LAYERS):    
        for minParticle in range(MAX_MIN_PARTICLES, MIN_MIN_PARTICLES, STEP_MIN_PARTICLES):       
            annealingLayers = str(annealingLayer)
            minParticles    = str(minParticle)
            # for approxIndex, approxType in enumerate(APPROXTYPESARRAY):
            #     initEnvironmentVariables()
                # for approxLevel in range(2, MAXLOOPSKIP):
            # os.environ[approxType] = str(approxLevel)
            outputValDict, error = runBodytrackAndGetError(annealingLayers, minParticles)
            # localSpeedUp = initValBaseDict["time"][approxIndex] / outputValDict["time"][approxIndex]
            globalSpeedUp = initValBaseDict["time"][-1] / outputValDict["time"][-1]
            # approxLevelToPrint = 1 / float(approxLevel)
            # if "LOOPCHANGE" in approxType:
            #     approxLevelToPrint = 1 - approxLevelToPrint

            writeToLog("{0}, {1}, {2}, {3}, {4}, {5}, {6}\n".format(annealingLayers, minParticles, globalSpeedUp, error, outputValDict["iter"][0], outputValDict["iter"][1], outputValDict["iter"][2]))


def combineApproximations():
    global OUTPUTDIR
    global LOGFILENAME

    OUTPUTDIR = "./approx_vid_out_combined"
    LOGFILENAME = "./combinedApproximationLog_phases.txt"
    init()
    # initValBaseDict = getBaseValue()

    # approxCounterMax = (MAXLOOPSKIP - 2) ** (len(APPROXTYPESARRAY)) 

    # for inputIndex in range(0, len(ANNEALINGLAYERS)):

    inputIndex = 3

    os.environ["PHASETOAPPROXIMATE"] = str(0)
    os.environ[APPROXTYPESARRAY[0]] = str(0)
    #annealingLayers = ANNEALINGLAYERS[inputIndex]
    #minParticles    = MINPARTICLES[inputIndex]
    annealingLayers = 23
    minParticles    = 2000
    initValBaseDict = getBaseValue(annealingLayers, minParticles)

    writeToLog("###!!!####INIT {0}, {1}, {2}, {3}, {4}, {5}\n".format(annealingLayers, minParticles, initValBaseDict["iter"][0], initValBaseDict["iter"][1], initValBaseDict["iter"][2], initValBaseDict["time"][-1])) 

    os.environ["TOTALITERATIONS"] = str(initValBaseDict["iter"][0])

    for approx0Value in range(1, MAXLOOPSKIP):
        for phase in PHASEOPTIONS:   
            for approx1Value in range(1, MAXLOOPSKIP):
                for approx2Value in range(1, MAXLOOPSKIP):
                    for approx3Value in range(1, MAXLOOPSKIP):
                        # for approx4Value in range(20, 40, 30):
                        os.environ[APPROXTYPESARRAY[0]] = str(approx0Value)
                        os.environ[APPROXTYPESARRAY[1]] = str(approx1Value)
                        os.environ[APPROXTYPESARRAY[2]] = str(approx2Value)
                        os.environ[APPROXTYPESARRAY[3]] = str(approx3Value)
                        # os.environ["PARTICLE_FACTOR_CHANGE"] = str(approx4Value)
                        os.environ["PHASETOAPPROXIMATE"] = str(phase)
                        outputValDict, error = runBodytrackAndGetError(annealingLayers, minParticles)
                        # localSpeedUp = initValBaseDict["time"][approxIndex] / outputValDict["time"][approxIndex]
                        # globalSpeedUp = initValBaseDict["time"][-1] / outputValDict["time"][-1]
                        approx0ValueToPrint = 1 / float(approx0Value)
                        if "LOOPCHANGE" in APPROXTYPESARRAY[0]:
                            approx0ValueToPrint = 1 - approx0ValueToPrint

                        approx1ValueToPrint = 1 / float(approx1Value)
                        if "LOOPCHANGE" in APPROXTYPESARRAY[1]:
                            approx1ValueToPrint = 1 - approx1ValueToPrint

                        approx2ValueToPrint = 1 / float(approx2Value)
                        if "LOOPCHANGE" in APPROXTYPESARRAY[2]:
                            approx2ValueToPrint = 1 - approx2ValueToPrint

                        approx3ValueToPrint = 1 / float(approx3Value)
                        if "LOOPCHANGE" in APPROXTYPESARRAY[3]:
                            approx3ValueToPrint = 1 - approx3ValueToPrint

                        writeToLog("{0}, {1}, {2}, {3}, {4}, {5}, {6}, {7}\n".format(phase, approx0ValueToPrint, approx1ValueToPrint, approx2ValueToPrint, approx3ValueToPrint, error, outputValDict["time"][-1], outputValDict["iter"][0]))



    # for approxCounter in range(0, approxCounterMax):
    #     tempApproxCounter = approxCounter
    #     for approxIndex,approxType in reversed(enumerate(APPROXTYPESARRAY)):
    #         envVariableValue = (approxCounter / ((MAXLOOPSKIP - 2) ** (approxIndex))) % (MAXLOOPSKIP - 2)
    #         os.environ[approxType] = str(envVariableValue)
    #         if(envVariableValue > 0):
    #             tempApproxCounter -= ((MAXLOOPSKIP - 2) ** (approxIndex) + envVariableValue)





def main():
    # global MAXLOOPSKIP
    #Clean the bin and compile the code
    # os.system("bash testChangeLoop  > /dev/null 2>&1")
    # for i in range(5):
    # runIndividualApproximations()
    # runDifferentParams()
    # MAXLOOPSKIP = 7
    combineApproximations()

if __name__ == "__main__":
    main()
