import os
import subprocess

runProc = subprocess.Popen(['time', 'ls'],stdout=subprocess.PIPE,stderr=subprocess.PIPE)
output, error = runProc.communicate()

print(output)
