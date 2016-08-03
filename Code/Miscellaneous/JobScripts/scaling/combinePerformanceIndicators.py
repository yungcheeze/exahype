import os,sys

output  = "\"JobId\",\"Application\",\"Dim\",\"nProcs\",\"p\",\"hMax\",\"Compiler\",\"Mode\",\"optimisticADERDGSteps\",\"averageTimestepCostOfSingleADERDGStep\",\"numberOfSteps\",\"simulationRuntime\",\"costOfADERDGSteps\"\n"

filenames = os.listdir(".")
for filename in filenames:
  #print "Looking at", filename
  if os.path.isdir(filename):
    nestedfiles = os.listdir(filename + "/")
    for i in range(len(nestedfiles)):
       nestedfiles[i] = filename + "/" + nestedfiles[i]
    filenames.extend(nestedfiles)
  else:
    if (filename == "out.txt"):
      continue
      
    output += filename.split(".")[0] + ","
    # print filename
    tokens = filename.split("__")
    
    for token in tokens:
      output += token.split("_")[1] + ","
      
      
    # get scaling results
    file = open(filename, "r" )

    line = file.readline()
    line = file.readline()
    
    tokens = line.split(" ")
    # print tokens
    
    output += str(tokens[1]) + "," + str(tokens[4]) + "," + str(tokens[7]) + "," + str(tokens[10]) + "," + str(tokens[13]) + "\n"
      
print output
 
   
   
   
   
   
   
   
   
    
    