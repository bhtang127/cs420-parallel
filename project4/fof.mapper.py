#!/usr/bin/env python

"""An Simple Implement of Mapper in Streaming MapReduce for Friends Triangle"""

import sys

for line in sys.stdin:
    
    # split line to get ID
    
    line.strip()
    IDs = line.split()

    ## Emit (Key, Value) Pair
    ## Key: (A B C) is all the tuple such that B < C  
    ##       A B C all different and appear simutaneously in one file
    ##       and one of A B C is the name for this file 
    ## Value: 1
    ## Since value is all 1 and this is a simple streaming version
    ## We don't output the value 1 for simplicity

    L = len(IDs)
    if L > 2:
        ## ID in filename is fixed, Loop for other two IDs to form triple
        for i in range(1, L):
            for j in range(i+1, L):
                triple = sorted( [ int(IDs[0]), int(IDs[i]), int(IDs[j]) ] )
                print "%s %s %s" % ( triple[0], triple[1], triple[2] )
                print "%s %s %s" % ( triple[1], triple[0], triple[2] )
                print "%s %s %s" % ( triple[2], triple[0], triple[1] )
                
