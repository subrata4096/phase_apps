#Manish's script
import re
import os
import shutil
import subprocess
import time
#import matplotlib
#matplotlib.use('Agg')
#import matplotlib.pyplot as plt
#import numpy as np
#from scipy import interpolate
import random

MAXAPPROX = 6
NUMRANDOM = 200
OUTPUTDIR = "./approx_vid_out"
SSIM_PATTERN_COMPILED = re.compile('[.\s]*Parsed_ssim_0.*All:(.*)\s\(')
PSNR_PATTERN_COMPILED = re.compile('[.\s]*Parsed_psnr_1.*average:(.*)\smin')
BITRATE_PATTERN_COMPILED = re.compile('[.\s]*Duration.*start.* bitrate:\s(\d*)\skb')

# FPSARRAY = ['20','30','40','50','60']
FPSARRAY = ['10','15','20','25','25']
PHASEOPTIONS = [1, 2, 3, 4, -1]
NUM_PHASES = 4

# FHD       = "1920:1080"
# WSXGAPlus = "1680:1050"
# UXGA      = "1600:1200"
# HDPlus    = "1600:900"
# HD        = "1366:768"
# RESOLUTIONARRAY = [FHD, WSXGAPlus, UXGA, HDPlus, HD]

WXGAH720p       = "1280:720"
IPHONE5         = "1136:640"
WSVGA           = "1024:640"
PSVITA          = "960:544"
VGA             = "640:480"
RESOLUTIONARRAY = [WXGAH720p, IPHONE5, WSVGA, PSVITA, VGA]

VIDEONAMES = ["./videos/adele_static_short_500.mp4", "./videos/top_static_short.mp4", "./videos/sherlock_medium_short.mp4", "./videos/himym_medium_short.mp4", "./videos/basketball_fast_short.mp4", "./videos/nypace_fast_short.mp4", "./videos/debug_vid_orig.mp4"]
VIDEOTYPENAMES = ["static_adele", "static_top", "medium_sherlock", "medium_himym", "fast_bball", "fast_nypace", "test_debug"] 

LOGFILENAME = "approximationLog.txt"
# LOGFILE = open("approximationLog.txt", "a")

APPROXTYPESARRAY  = ['DEFLATE', 'EMBOSS', 'BOXBLUR']
APPROXCOMMANDSARR = ['deflate', 'convolution="-2 -1 0 -1 1 1 0 1 2:-2 -1 0 -1 1 1 0 1 2:-2 -1 0 -1 1 1 0 1 2:-2 -1 0 -1 1 1 0 1 2"', 'boxblur=10']
APPROXCOMMANDSARR = ['convolution="-2 -1 0 -1 1 1 0 1 2:-2 -1 0 -1 1 1 0 1 2:-2 -1 0 -1 1 1 0 1 2:-2 -1 0 -1 1 1 0 1 2"'] #DEBUGDEBUG

APPROXGLOBALVARS  = ['LOOPITERSKIP_NEIGHBOUR', 'LOOPITERSKIP_CONVOLVE', 'LOOPITERSKIP_BOXBLUR']
# APPROXGLOBALVARS  = ['LOOPITERSKIP_CONVOLVE']

# The overall time i.e. bench has to be last index in array
STARTERSTRINGS    = ["Time taken by deflate = ", "Time taken by embross = ", "Time taken by vblur = ", "bench: utime="]


# returns the time as a float using the starterString
def getTime(starterString, fileStr):
    # print("Total deflate time = {0}".format(getTime("Time taken by deflate = ",  fileName)))
    # print("Total vblur time = {0}".format(getTime("Time taken by vblur = ",  fileName)))
    # print("Total embross time = {0}".format(getTime("Time taken by embross = ",  fileName)))
    # print("bench: utime = {0}".format(getTime("bench: utime=",  fileName)))
    # fh = open(fileName)
    # fileLines = fh.readlines()
    fileLines = fileStr.split("\n")
    totalTime = 0

    for line in fileLines:
        if starterString in line:
            totalTime += float(line.split(starterString)[1].replace("s", ""))

    return totalTime


def getBitrate(fileName):    
    bitrate = -1
    output, error = runFFPROBECommand(fileName)
    bitrateMatchObj = BITRATE_PATTERN_COMPILED.search(error)
    if(bitrateMatchObj):
        bitrate = bitrateMatchObj.group(1)

    return bitrate



# returns a dict with keys 'psnr', 'ssim'
def getPSNR_SSIM(referenceFileName, outputFileName):
    outputDict = {}
    outputDict['psnr'] = 'ERROR'
    outputDict['ssim'] = 'ERROR' 

    # ./build/bin/ffmpeg -i videos/output/correct_conv.mp4 -i videos/output/approx_conv.mp4 -lavfi  "ssim;[0:v][1:v]psnr" -f null - > psnr_ssimtxt.txt
    runProc = subprocess.Popen(['./build/bin/ffmpeg', '-i', referenceFileName, '-i', outputFileName, '-lavfi', 'ssim;[0:v][1:v]psnr', '-f', 'null', '-'],stdout=subprocess.PIPE,stderr=subprocess.PIPE)
    output, error = runProc.communicate()

    psnrMatchOBj = PSNR_PATTERN_COMPILED.search(error)
    if(psnrMatchOBj):
        outputDict['psnr'] = psnrMatchOBj.group(1)


    ssimMatchObj = SSIM_PATTERN_COMPILED.search(error)
    if(ssimMatchObj):
        outputDict['ssim'] = ssimMatchObj.group(1)

    return outputDict


# puts the command in the bash file
def setCommand(command, fileName="test.sh"):
    fh = open(fileName, "w")
    fh.write(command)
    fh.close()
    return

# write to log file
def writeToLog(text):
    print(text.strip("\n"))
    fh = open(LOGFILENAME, "a")
    fh.write(text)
    fh.close()
    return


# changes the video to specified dimensions BIT RATE COMMENTED OUT AS OF NOW
def changeVideo(fileName, outputFileName, bitRate, resolution, fps):
    # '-y' force overwrite
    # runProc = subprocess.Popen(['./build/bin/ffmpeg', '-y', '-i', fileName, '-r', fps, '-b:v', bitRate, '-bufsize', '64k', '-vf', 'scale='+resolution, outputFileName],stdout=subprocess.PIPE,stderr=subprocess.PIPE)
    # runProc = subprocess.Popen(['./build/bin/ffmpeg', '-y', '-i', fileName, '-r', '30', outputFileName],stdout=subprocess.PIPE,stderr=subprocess.PIPE)
    # errors out if not done through bash, not sure why? but doesnt change a lot of stuff #'-b:v', bitRate, '-bufsize', '64k', 
    # setCommand(" ".join(['./build/bin/ffmpeg', '-y', '-r', fps, '-i', fileName, '-vf', 'scale='+resolution, outputFileName])) #DEBUGDEBUGDEBUG
    setCommand(" ".join(['cp', fileName, outputFileName])) #DEBUGDEBUGDEBUG

    runProc = subprocess.Popen(['bash', 'test.sh'],stdout=subprocess.PIPE,stderr=subprocess.PIPE)
    output, error = runProc.communicate()

    return


# Runs a given command and returns output and error
def runFFMPEGCommand(command, fileName="test.sh"):
    print command
    setCommand(command, fileName)
    runProc = subprocess.Popen(['bash', 'test.sh'],stdout=subprocess.PIPE,stderr=subprocess.PIPE)
    output, error = runProc.communicate()

    return output,error



# Runs a given command and returns output and error
def runFFPROBECommand(fileToProbe, fileName="test.sh"):
    setCommand("./build/bin/ffprobe " + fileToProbe, fileName)
    runProc = subprocess.Popen(['bash', 'test.sh'],stdout=subprocess.PIPE,stderr=subprocess.PIPE)
    output, error = runProc.communicate()
    return output,error


# runs the actual ffmpeg which is being measured
def runFFMPEGAndGetData(videoName, fps, resolution, bitrate, vidType, approxType, dirPath, frames):

    approxTypePath = os.path.join(dirPath, str(approxType))
    if not os.path.isdir(approxTypePath):
        os.mkdir(approxTypePath)

    # Filter used to get local error i.e. stop right after place being approximated
    localFilter  = ",".join(APPROXCOMMANDSARR[0:(approxType+1)])

    # Filter used to get global error 
    globalFilter = ",".join(APPROXCOMMANDSARR)

    globalBitrate = -1
    localBitrate = -1
    globalBaseBitrate = -1
    localBaseBitrate = -1

    # Runs to get base video used for local error calculation
    output,error = runFFMPEGCommand("./build/bin/ffmpeg -y -i " + videoName + " -vf " + localFilter + " " + os.path.join(approxTypePath, "localBase.mp4"))
    localBaseTime  = getTime(STARTERSTRINGS[approxType], output)

    output,error = runFFMPEGCommand("./build/bin/ffmpeg -y -i " + videoName + " -vf " + globalFilter + " " + os.path.join(approxTypePath, "globalBase.mp4"))
    globalBaseTime = getTime(STARTERSTRINGS[-1], output)

    for approxLevel in range(2, MAXAPPROX):
        os.environ[APPROXGLOBALVARS[approxType]] = str(approxLevel)
    
        output,error = runFFMPEGCommand("./build/bin/ffmpeg -y -i " + videoName + " -vf " + localFilter + " " + os.path.join(approxTypePath, "localApprox.mp4"))
        localTime  = getTime(STARTERSTRINGS[approxType], output)
    
        localError = getPSNR_SSIM(os.path.join(approxTypePath, "localBase.mp4"), os.path.join(approxTypePath, "localApprox.mp4"))

        output,error = runFFMPEGCommand("./build/bin/ffmpeg -y -i " + videoName + " -vf " + globalFilter + " " + os.path.join(approxTypePath, "globalApprox.mp4"))
        globalTime = getTime(STARTERSTRINGS[-1], output)

        globalError = getPSNR_SSIM(os.path.join(approxTypePath, "globalBase.mp4"), os.path.join(approxTypePath, "globalApprox.mp4"))

        os.rename(os.path.join(approxTypePath, "localApprox.mp4"), os.path.join(approxTypePath, "localApprox_" + str(approxLevel) + ".mp4"))
        os.rename(os.path.join(approxTypePath, "globalApprox.mp4"), os.path.join(approxTypePath, "globalApprox_" + str(approxLevel) + ".mp4"))

        # Speedup is old/new
        localSpeedUp = localBaseTime/float(localTime)
        globalSpeedup = globalBaseTime/float(globalTime)

        # ApproxFunct, APPROX_LEVEL, in_fps, in_bitrate, in_resolution, in_type, in_frames, out_path, out_localerror_psnr, out_localerror_ssim, out_globalerror_psnr, out_globalerror_ssim, out_localspeedup, out_globalspeedup, out_globalBitRate, out_localBitRate
        writeToLog("{0}, {1}, {2}, {3}, {4}, {5}, {6}, {7}, {8}, {9}, {10}, {11}, {12}, {13}\n".format(APPROXTYPESARRAY[approxType], 1/float(approxLevel), fps, bitrate, resolution, vidType, frames, "-".join(APPROXTYPESARRAY[approxType:]), localError['psnr'], localError['ssim'], globalError['psnr'], globalError['ssim'], localSpeedUp, globalSpeedup)) 

    # #Delete temp file and folder
    # shutil.rmtree(approxTypePath)



def getFrames(videoPath, fileName="test.sh"):
    setCommand("./build/bin/ffprobe -v error -count_frames -select_streams v:0   -show_entries stream=nb_read_frames -of default=nokey=1:noprint_wrappers=1 " + videoPath)
    runProc = subprocess.Popen(['bash', 'test.sh'],stdout=subprocess.PIPE,stderr=subprocess.PIPE)
    output, error = runProc.communicate()
    return output.strip('\n')

# runs the actual ffmpeg which is being measured
def randomRunFFMPEGAndGetData(videoName, fps, resolution, bitrate, vidType, dirPath, approx1ValRandom, approx2ValRandom, approx3ValRandom, frames, phase, basePath, globalBaseTime):

    approxTypePath = dirPath

    # Filter used to get global error 
    globalFilter = ",".join(APPROXCOMMANDSARR)

    os.environ[APPROXGLOBALVARS[0]] = str(approx1ValRandom)
    os.environ[APPROXGLOBALVARS[1]] = str(approx2ValRandom)
    os.environ[APPROXGLOBALVARS[2]] = str(approx3ValRandom)
    globalBitrate = -1
    globalBaseBitrate = -1

    output,error = runFFMPEGCommand("./build/bin/ffmpeg -y -i " + videoName + "  -vf " + globalFilter + " " + os.path.join(approxTypePath, "globalApprox.mp4"))
    # output,error = runFFMPEGCommand("./build/bin/ffmpeg -y -i " + videoName + " -vf " + globalFilter + " " + os.path.join(approxTypePath, "globalApprox.yuv"))
    globalTime = getTime(STARTERSTRINGS[-1], output)
    if globalTime == 0:
        print("ERROR")
        print(output)
        print(error)
        print(globalTime)
        return

    globalError = getPSNR_SSIM(os.path.join(basePath), os.path.join(approxTypePath, "globalApprox.mp4"))

    globalBitrate = getBitrate(os.path.join(approxTypePath, "globalApprox.mp4"))

    # Speedup is old/new
    globalSpeedup = globalBaseTime/float(globalTime)

    # ApproxFunct, APPROX_LEVEL, in_fps, in_bitrate, in_resolution, in_type, out_path, out_localerror_psnr, out_localerror_ssim, out_globalerror_psnr, out_globalerror_ssim, out_localspeedup, out_globalspeedup ///////////out_globalBitRate, out_localBitRate
    writeToLog("{0}, {1}, {2}, {3}, {4}, {5}, {6}, {7}, {8}, {9}, {10}, {11}\n".format(phase, 1/float(approx1ValRandom) if (approx1ValRandom != 1) else 0.0, 1/float(approx2ValRandom) if (approx2ValRandom != 1) else 0.0, 1/float(approx3ValRandom) if (approx3ValRandom != 1) else 0.0, fps, bitrate, resolution, vidType,frames, globalError['psnr'], 1 - float(globalError['ssim']), globalSpeedup))

    # #Delete temp file and folder
    # shutil.rmtree(approxTypePath)



# Initialize everything
def init():
    #Make tempdirectoy
    if not os.path.isdir(OUTPUTDIR):
        os.mkdir(OUTPUTDIR)
    initEnvironmentVariables()


# initializes all necessary environment variables
def initEnvironmentVariables():    
    os.environ["LOOPITERSKIP_BOXBLUR"]        ='1'
    os.environ["LOOPCHANGE_BOXBLUR"]          ='1'
    os.environ["LOOPITERSKIP_NEIGHBOUR"]      ='1'
    os.environ["LOOPCHANGE_NEIGHBOUR"]        ='1'
    os.environ["LOOPCHANGE_H815SCALE"]        ='1'
    os.environ["LOOPITERSKIP_H815SCALE"]      ='1'
    os.environ["LOOPCHANGE_CONVOLVE"]         ='1'
    os.environ["LOOPITERSKIP_CONVOLVE"]       ='1'
    os.environ["LOOPCHANGE_HSCALE"]           ='1'
    os.environ["LOOPITERSKIP_HSCALE"]         ='1'
    os.environ["LOOPCHANGE_SWSCALE"]          ='1'
    os.environ["LOOPITERSKIP_SWSCALE"]        ='1'
    os.environ["LOOPCHANGE_SWSCALE"]          ='1'
    os.environ["TOTAL_FRAMES_MANISH"]         ='125'
    os.environ["PHASETOAPPROXIMATE"]          ='0'

def runAllVideos():
    bitrate = "-1"

    init()
    for vidIndex, videoName in enumerate(VIDEONAMES):
        vidType = VIDEOTYPENAMES[vidIndex]
        videoPath = os.path.join(OUTPUTDIR, vidType)
        if not os.path.isdir(videoPath):
            os.mkdir(videoPath)
        for fps in FPSARRAY:
            fpsPath = os.path.join(videoPath, fps)
            if not os.path.isdir(fpsPath):
                os.mkdir(fpsPath)
            for resolution in RESOLUTIONARRAY:
                resolutionPath = os.path.join(fpsPath, resolution.replace(":", "_"))
                if not os.path.isdir(resolutionPath):
                    os.mkdir(resolutionPath)
                finalPath = resolutionPath

                changeVideo(fileName=videoName, outputFileName=os.path.join(finalPath, "input.mp4"), bitRate=bitrate, resolution=resolution, fps=fps)

                frames = getFrames(os.path.join(resolutionPath, "input.mp4"))
                bitrate = getBitrate(os.path.join(finalPath, "input.mp4"))

                for approxType in range(0, len(APPROXCOMMANDSARR)):
                    runFFMPEGAndGetData(os.path.join(finalPath, "input.mp4"), fps, resolution, bitrate, vidType, approxType, finalPath, frames)
                    os.environ[APPROXGLOBALVARS[approxType]] = str(1)

# Randomly runs testcases with combinations of approximations
def combineApprox():
    global OUTPUTDIR
    global LOGFILENAME

    OUTPUTDIR = "./approx_vid_out_combined_phases"
    LOGFILENAME = "./combinedApproximationLog_phases.txt"
    init()
    # DEFLATE, 0.2, 30, 772, 1280:720, medium_sherlock,

    count = 0
    bitrate = 01
    globalFilter = ",".join(APPROXCOMMANDSARR)

    videoName = VIDEONAMES[0]
    vidType = VIDEOTYPENAMES[0]
    resolution = RESOLUTIONARRAY[0]
    fps        = FPSARRAY[4]

    videoPath = os.path.join(OUTPUTDIR, vidType)
    if not os.path.isdir(videoPath):
        os.mkdir(videoPath)
    fpsPath = os.path.join(videoPath, fps)
    if not os.path.isdir(fpsPath):
        os.mkdir(fpsPath)
    resolutionPath = os.path.join(fpsPath, resolution.replace(":", "_"))
    if not os.path.isdir(resolutionPath):
        os.mkdir(resolutionPath)

    os.environ["TOTAL_FRAMES_MANISH"] = str(((NUM_PHASES - 1) / int(NUM_PHASES)) * int(NUM_PHASES)) 
    os.environ["PHASETOAPPROXIMATE"] = str(0)
    os.environ["TOTAL_PHASES_MANISH"] = str(NUM_PHASES)

    changeVideo(fileName=videoName, outputFileName=os.path.join(resolutionPath, "input.mp4"), bitRate=bitrate, resolution=resolution, fps=fps)
    frames = getFrames(os.path.join(resolutionPath, "input.mp4"))
    os.environ["TOTAL_FRAMES_MANISH"] = str(1000) #str(((int(frames) + NUM_PHASES - 1) / int(NUM_PHASES)) * int(NUM_PHASES)) 
    os.environ["PHASETOAPPROXIMATE"] = str(0)
    os.environ["TOTAL_PHASES_MANISH"] = str(NUM_PHASES)
    output, error = runFFPROBECommand(os.path.join(resolutionPath, "input.mp4"))
    bitrateMatchObj = BITRATE_PATTERN_COMPILED.search(error)
    if(bitrateMatchObj):
        bitrate = bitrateMatchObj.group(1)
                        


    os.environ[APPROXGLOBALVARS[0]] = str(1)
    os.environ[APPROXGLOBALVARS[1]] = str(1)
    os.environ[APPROXGLOBALVARS[2]] = str(1)

    # Runs to get base video used for local error calculation
    output,error = runFFMPEGCommand("./build/bin/ffmpeg -y -i " + os.path.join(resolutionPath, "input.mp4") + " -vf " + globalFilter + " " + os.path.join(resolutionPath, "globalBase.mp4"))
    globalBaseTime = getTime(STARTERSTRINGS[-1], output)

    for phase in PHASEOPTIONS:                      
        if phase > NUM_PHASES:
            continue     
        os.environ["PHASETOAPPROXIMATE"] = str(phase)
        for approx2ValRandom in range(2, MAXAPPROX):
        # for approx2ValRandom in range(1, MAXAPPROX):
            # for approx3ValRandom in range(1, MAXAPPROX): #DEBUGDEBUG
            approx3ValRandom = 1 #DEBUGDEBUG
            approx1ValRandom = 1 #DEBUGDEBUG
            finalPath = os.path.join(resolutionPath, str(approx1ValRandom) + "_" + str(approx2ValRandom) + "_" + str(approx3ValRandom))
            if not os.path.isdir(finalPath):
                os.mkdir(finalPath)


            randomRunFFMPEGAndGetData(os.path.join(resolutionPath, "input.mp4"), fps, resolution, bitrate, vidType, finalPath, approx1ValRandom, approx2ValRandom, approx3ValRandom, frames, phase, os.path.join(resolutionPath, "globalBase.mp4"), globalBaseTime)


# Randomly runs testcases with combinations of approximations
def randomRunVideo():
    global OUTPUTDIR
    global LOGFILENAME

    OUTPUTDIR = "./approx_vid_out_random"
    LOGFILENAME = "./randomApproximationLog.txt"
    init()

    count = 0
    bitrate = 01

    while(count < NUMRANDOM):
        vidRandom        = random.randrange(len(VIDEONAMES))
        fpsRandom        = random.randrange(len(FPSARRAY))
        resolutionRandom = random.randrange(len(RESOLUTIONARRAY))
        approx1ValRandom = random.randrange(1, MAXAPPROX)
        approx2ValRandom = random.randrange(1, MAXAPPROX)
        approx3ValRandom = random.randrange(1, MAXAPPROX)

        videoName = VIDEONAMES[vidRandom]
        vidType = VIDEOTYPENAMES[vidRandom]
        resolution = RESOLUTIONARRAY[resolutionRandom]
        fps        = FPSARRAY[fpsRandom]

        videoPath = os.path.join(OUTPUTDIR, vidType)
        if not os.path.isdir(videoPath):
            os.mkdir(videoPath)
        fpsPath = os.path.join(videoPath, fps)
        if not os.path.isdir(fpsPath):
            os.mkdir(fpsPath)
        resolutionPath = os.path.join(fpsPath, resolution.replace(":", "_"))
        if not os.path.isdir(resolutionPath):
            os.mkdir(resolutionPath)

        finalPath = os.path.join(resolutionPath, str(approx1ValRandom) + "_" + str(approx2ValRandom) + "_" + str(approx3ValRandom))
        if os.path.isdir(finalPath):
            continue
        else:
            os.mkdir(finalPath)

        changeVideo(fileName=videoName, outputFileName=os.path.join(finalPath, "input.mp4"), bitRate=bitrate, resolution=resolution, fps=fps)
        
        output, error = runFFPROBECommand(videoPath)
        bitrateMatchObj = BITRATE_PATTERN_COMPILED.search(error)
        if(bitrateMatchObj):
            bitrate = bitrateMatchObj.group(1)

        randomRunFFMPEGAndGetData(os.path.join(finalPath, "input.mp4"), fps, resolution, bitrate, vidType, finalPath, approx1ValRandom, approx2ValRandom, approx3ValRandom)
        count += 1


def main():
    # runAllVideos()
    # randomRunVideo()
    combineApprox()

if __name__ == "__main__":
    main()


# Notes1
# Not using scale, cause its algorithm is used by  default so it will add errors in 2 different places if used in addition to a filte rrather 
# than just one.
