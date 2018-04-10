#!/usr/bin/env python

"""An Simple Implement of Mapper in Streaming MapReduce for Friends Triangle"""

import sys

for line in sys.stdin:
    
    # split line to get ID
    
    line.strip()
    IDs = line.split()

    ## Emit (Key, Value) Pair
    ## Key: (A B C) is all the tuple such that B < C  
    ##   and A B C appear simutaneously in one file
    ## Value: 1
    ## Since value is all 1 and this is a simple streaming version
    ## We don't output the value 1 for simplicity

    for A in IDs:
        for B in IDs:
            for C in IDs:
                if A != B and A != C and B < C:
                    print "%s %s %s" % ( A, B, C )
