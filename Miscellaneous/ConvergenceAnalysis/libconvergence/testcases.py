#!/usr/bin/env python

# an arena for simple unit tests

from convergence_helpers import MethodActions

class B:
   chain = MethodActions()

   @chain.add(1, "The infamous foo")
   def doFoo(self):
      pass
   @chain.add(2, "The bar")
   def doBar(self):
      pass
   def __init__(self,):
      def information(i,key):
          print("Number %d: Calling %s" % (i, self.chain.storage[key]))
      self.chain.sortedcall(self, information)

B()
