#Manish's script
import re
import os
import shutil
import subprocess
import time
import random

BACKUP_INITIAL_LOG = "./approximationLog_bup.txt"
INITIAL_LOG = "./approximationLog.txt"
OUTPUT_LOG  = "./fixed_approximationLog.txt"
BITRATE_PATTERN_COMPILED = re.compile('[.\s]*Duration.*start.* bitrate:\s(\d*)\skb')
BASEDIRECTORY_VIDEOS = "./approx_vid_out/"
APPROXTYPESARRAY  = ['DEFLATE', 'EMBOSS', 'BOXBLUR']

# puts the command in the bash file
def setCommand(command, fileName="test.sh"):
    fh = open(fileName, "w")
    fh.write(command)
    fh.close()
    return


def writeToNewLog(command):
    OutputFile = open(OUTPUT_LOG, "a")
    OutputFile.write(command)
    OutputFile.close()


# Runs a given command and returns output and error
def runFFPROBECommand(fileToProbe, fileName="test.sh"):
    setCommand("./build/bin/ffprobe " + fileToProbe, fileName)
    runProc = subprocess.Popen(['bash', 'test.sh'],stdout=subprocess.PIPE,stderr=subprocess.PIPE)
    output, error = runProc.communicate()
    return output,error


def getFPS(videoPath, fileName="test.sh"):
    setCommand("./build/bin/ffprobe -v error -count_frames -select_streams v:0   -show_entries stream=nb_read_frames -of default=nokey=1:noprint_wrappers=1 " + videoPath)
    runProc = subprocess.Popen(['bash', 'test.sh'],stdout=subprocess.PIPE,stderr=subprocess.PIPE)
    output, error = runProc.communicate()
    return output.strip('\n')

# Reads in input line, gets video and updates fps
def addFPS():
    InputFile = open(INITIAL_LOG)
    inputLines = InputFile.readlines()
    newBitRate = "-1"
    InputFile.close()
    os.rename(INITIAL_LOG, BACKUP_INITIAL_LOG)
    fpsDict = {}
    
    for inputLine in inputLines:
        inputLineArray = inputLine.split(',')
        videoPath = os.path.join(BASEDIRECTORY_VIDEOS, inputLineArray[5].strip(), inputLineArray[2].strip(), inputLineArray[4].strip().replace(":", "_"), "input.mp4")
        if videoPath not in fpsDict.keys():
            fpsDict[videoPath] = getFPS(videoPath)
        inputLineArray.insert(6,' ' + fpsDict[videoPath])
        writeToNewLog(','.join(inputLineArray))
    os.rename(OUTPUT_LOG, INITIAL_LOG)
    print("DONE WITH INPUT FPS FIXING")
    return


 # DEFLATE, 0.5, 10, 1011, 1280:720, static_adele, 
  # approx_vid_out/fast_nypace/10/1024_640/input.mp4  
# Reads in input line, gets video and updates bitrate
def fixInputBitrate():
    InputFile = open(INITIAL_LOG)
    inputLines = InputFile.readlines()
    newBitRate = "-1"
    InputFile.close()
    os.rename(INITIAL_LOG, BACKUP_INITIAL_LOG)
    
    for inputLine in inputLines:
        inputLineArray = inputLine.split(',')
        videoPath = os.path.join(BASEDIRECTORY_VIDEOS, inputLineArray[5].strip(), inputLineArray[2].strip(), inputLineArray[4].strip().replace(":", "_"), "input.mp4")
        output, error = runFFPROBECommand(videoPath)
        bitrateMatchObj = BITRATE_PATTERN_COMPILED.search(error)
        if(bitrateMatchObj):
            newBitRate = bitrateMatchObj.group(1)
        inputLineArray[3] = ' ' + newBitRate
        writeToNewLog(','.join(inputLineArray))

    os.rename(OUTPUT_LOG, INITIAL_LOG)
    print("DONE WITH INPUT BITRATE FIXING")
    return

# adds output bitrate to log
 # DEFLATE, 0.5, 10, 1011, 1280:720, static_adele, 
# approx_vid_out/fast_nypace/10/1024_640/0/globalApprox_2.mp4
# add global followed by local
def fixOutputBitrate():
    InputFile = open(INITIAL_LOG)
    inputLines = InputFile.readlines()
    newBitRate = "-1"
    InputFile.close()

    for inputLine in inputLines:

        # Get the global approximated video's bitrate
        inputLineArray = inputLine.split(',')
        globablVideoPath = os.path.join(BASEDIRECTORY_VIDEOS, inputLineArray[5].strip(), inputLineArray[2].strip(), inputLineArray[4].strip().replace(":", "_"),  str(APPROXTYPESARRAY.index(inputLineArray[0].replace("!", "").replace("#", "").strip())), "globalApprox_" + str(int(round(1 / float(inputLineArray[1].strip())))) + ".mp4" )
        output, error = runFFPROBECommand(globablVideoPath)
        bitrateMatchObj = BITRATE_PATTERN_COMPILED.search(error)
        if(bitrateMatchObj):
            newBitRate = bitrateMatchObj.group(1)
        inputLineArray[-1] = inputLineArray[-1].rstrip('\n')
        inputLineArray.append(' ' + newBitRate)

        # Get the local approximated video's bitrate
        localVideoPath = os.path.join(BASEDIRECTORY_VIDEOS, inputLineArray[5].strip(), inputLineArray[2].strip(), inputLineArray[4].strip().replace(":", "_"),  str(APPROXTYPESARRAY.index(inputLineArray[0].replace("!", "").replace("#", "").strip())), "localApprox_" + str(int(round(1 / float(inputLineArray[1].strip())))) + ".mp4" )
        output, error = runFFPROBECommand(localVideoPath)
        bitrateMatchObj = BITRATE_PATTERN_COMPILED.search(error)
        if(bitrateMatchObj):
            newBitRate = bitrateMatchObj.group(1)
        inputLineArray.append(' ' + newBitRate)

        # Get the base global video's bitrate
        globablVideoPath = os.path.join(BASEDIRECTORY_VIDEOS, inputLineArray[5].strip(), inputLineArray[2].strip(), inputLineArray[4].strip().replace(":", "_"),  str(APPROXTYPESARRAY.index(inputLineArray[0].replace("!", "").replace("#", "").strip())), "globalBase.mp4" )
        output, error = runFFPROBECommand(globablVideoPath)
        bitrateMatchObj = BITRATE_PATTERN_COMPILED.search(error)
        if(bitrateMatchObj):
            newBitRate = bitrateMatchObj.group(1)
        inputLineArray.append(' ' + newBitRate)

        # Get the base local video's bitrate
        localVideoPath = os.path.join(BASEDIRECTORY_VIDEOS, inputLineArray[5].strip(), inputLineArray[2].strip(), inputLineArray[4].strip().replace(":", "_"),  str(APPROXTYPESARRAY.index(inputLineArray[0].replace("!", "").replace("#", "").strip())), "localBase.mp4" )
        output, error = runFFPROBECommand(localVideoPath)
        bitrateMatchObj = BITRATE_PATTERN_COMPILED.search(error)
        if(bitrateMatchObj):
            newBitRate = bitrateMatchObj.group(1)
        inputLineArray.append(' ' + newBitRate)

        inputLineArray[-1] = inputLineArray[-1] + '\n'
        writeToNewLog(','.join(inputLineArray))

    InputFile.close()
    return

def main():
    fixInputBitrate()
    addFPS()
    # fixOutputBitrate()

if __name__ == "__main__":
    main()