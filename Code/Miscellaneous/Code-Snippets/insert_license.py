import os
import glob
import shutil

license0 = "/**\n"
license1 = " * This file is part of the ExaHyPE project.\n"
license2 = " * Copyright (c) 2016  http://exahype.eu\n"
license3 = " * All rights reserved.\n"
license4 = " *\n"
license5 = " * The project has received funding from the European Union's Horizon \n"
license6 = " * 2020 research and innovation programme under grant agreement\n"
license7 = " * No 671698. For copyrights and licensing, please consult the webpage.\n"
license8 = " *\n"
license9 = " * Released under the BSD 3 Open Source License.\n"
license10= " * For the full license text, see LICENSE.txt\n"
license11= " **/\n"
license12= " \n"
 
filenames = os.listdir(".")
#filenames = [ glob("*.hpp"), glob("*.cpp")]
for filename in filenames:
   #print "Looking at", filename
   if os.path.isdir(filename):
      nestedfiles = os.listdir(filename + "/")
      for i in range(len(nestedfiles)):
         nestedfiles[i] = filename + "/" + nestedfiles[i]
      filenames.extend(nestedfiles)
   else:
      insert = filename[-4:] == ".cpp"
      insert = insert or filename[-4:] ==".hpp" 
      insert = insert or filename[-2:] == ".h" 
      insert = insert or filename[-5:] == ".cpph"
      insert = insert or filename[-2:] == ".c"
#      insert = insert or filename[-4:] == ".f90"
      insert = insert or filename[-8:] == ".exahype"
      if insert:
         file = open(filename, 'r+')
         if file.readline() != license0:
            print "Insert license in file \"" + filename + "\""
            file.seek(0)
            text = file.read()
            file.seek(0)
            file.write(license0)
            file.write(license1)
            file.write(license2)
            file.write(license3)
            file.write(license4)
            file.write(license5)
            file.write(license6)
            file.write(license7)
            file.write(license8)
            file.write(license9)
            file.write(license10)
            file.write(license11)
            file.write(license12)
            file.write(text)
         file.close()
      else:
        print "DONT Insert license in file \"" + filename + "\""
