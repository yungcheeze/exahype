#!/bin/env python

import os

inputData = ""
outputPrefix = ""

#Generate default input data if no path specified or invalid path
if inputData == "" or not os.path.isfile(inputData):
    print("Generating and using default data")
    defaultData = "inputDataDefault.txt"
    if os.path.isfile(defaultData):
        os.system("rm "+defaultData)
    os.system("python generateDefaultInputData.py > "+defaultData)
    
#Run the executables and store the results
print("Executing test")
outputC = outputPrefix+"C.txt"
outputF = outputPrefix+"F.txt"
os.system("./C/mainC inputDataDefault.txt > Result/"+outputC)
os.system("./Fortran/mainFortran inputDataDefault.txt > Result/"+outputF)

#Comparing results
print("Comparing result")
CResult = list()
FResult = list()

with open("Result/"+outputC) as f:
    for line in f:  #Line is a string
        CResult.append(float(line))
        
with open("Result/"+outputF) as f:
    for line in f:  #Line is a string
        FResult.extend(map(float,line.split()))

if len(FResult) != len (CResult):
    print("ERROR: Result have different size !!!")
    print(CResult)
    print(FResult)
    quit()

sum_relError = 0.0
for (f,c) in zip(FResult,CResult):
    sum_relError =  sum_relError + (f-c)/f
    
print("Average relative error: "+str(sum_relError / len(FResult)))
