import re
import sys

def isSearchedStatisticsLine(data):
 # The space after the string is important as we do otherwise catch
 # the plot variant of the adapter as well
 SearchString = "ADERDGTimeStep "

 return re.search( SearchString, data ) and re.search( "logIterationStatistics",data);



def getEntryFromStatistics(filename,column):
  try:
    file = open(filename, "r" )
  except:
    print( "was not able to open file " + filename )
    return -1
  try:
    data = "   "
    
    while ( not data=="" and not isSearchedStatisticsLine(data) ):
        data = file.readline()
    if ( not data=="" ):
      return float(data.split("|")[column])
    else:
      return -1
  except:
    print "was not able to parse file " + filename
    return -1


def getOptimisticADERDGSteps(filename):
    return getEntryFromStatistics(filename,2)


def getAverageTimestepCostOfSingleADERDGStep(filename):
  return getEntryFromStatistics(filename,6)

def getCostOfADERDGSteps(filename):
  return getEntryFromStatistics(filename,5)

  
  
def getNumberOfSteps(filename):
  try:
    file = open(filename, "r" )
  except:
    print( "was not able to open file " + filename )
    return -1

  result = 0
  for line in file:
    if re.search("startNewTimeStep",line) and re.search( "step ", line):
      result = int( line.split( "step " )[1].split( " " )[0] )

  return result


def getSimulationRuntime(filename):
  try:
    file = open(filename, "r" )
  except:
    print( "was not able to open file " + filename )
    return -1

  SearchPattern          = "([0-9]\.?[0-9]*).*"

  result = 0.0

  for line in file:
    if re.search("startNewTimeStep",line) and re.search( "step ", line):
      m = re.search( SearchPattern, line )
      result = float( m.group(1) )

  return result

  
print "A = [A; [" , getOptimisticADERDGSteps(sys.argv[1]), ", ", getAverageTimestepCostOfSingleADERDGStep(sys.argv[1]), ", ", getNumberOfSteps(sys.argv[1]), ", ", getSimulationRuntime(sys.argv[1]), ", ", getCostOfADERDGSteps(sys.argv[1]), "]]"


