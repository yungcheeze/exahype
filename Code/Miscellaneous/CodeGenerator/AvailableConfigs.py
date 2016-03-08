#!/bin/env python
##
# Comments on the meaning of the arguments are copied from 
# https://github.com/hfp/libxsmm
# 

architectures = [ 
                    "noarch", 
                    "wsm",               # westmere
                    "snb",               # sandy bridge
                    "hsw",               # haswell
                    "knc",               # knights corner
                    "knl"                # knights landing
                ]

prefetchStrategies = [
                        "nopf",              # no prefetching at all, just 3 inputs (*A, *B, *C)
                        "pfsigonly",         # just prefetching signature, 6 inputs (*A, *B, *C, *A', *B', *C')
                        "BL2viaC",           # uses accesses to *C to prefetch *B'
                        "AL2",               # uses accesses to *A to prefetch *A'
                        "curAL2",            # prefetches current *A ahead in the kernel
                        "AL2_BL2viaC",       # combines AL2 and BL2viaC
                        "curAL2_BL2viaC",    # combines curAL2 and BL2viaC
                        "AL2jpst",           # aggressive *A' prefetch of first rows without any structure
                        "AL2jpst_BL2viaC"    # combines AL2jpst and BL2viaC
                     ]

precisions =  [
                "DP",
                "SP"
              ]

numerics   =  [
                "linear",
                "nonlinear"
              ]