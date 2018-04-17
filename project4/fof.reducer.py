#!/usr/bin/env python

"""An Simple Implement of Reducer in Streaming MapReduce for Friends Triangle"""

import sys
   
## Reducer is simply doing counting job
## Since mapper have already out put all
## possible legal triples of mutual friends suggest by each friends
## If a triple appear more than once, 
## then they are truly mutual friends and vice versa  
## Also, Since the output from mapper is sorted
## Therefore we can do the count serially

# Counting
## Begin with first line
line = sys.stdin.readline()
current_triple = line
current_count = 1

## Loop for serial counting
for next_triple in sys.stdin:

    if next_triple == current_triple:
        current_count += 1
        current_triple = next_triple

    else:
        if current_count > 1:
            print( current_triple.strip() )
        
        current_triple = next_triple
        current_count = 1

## End with last line
if current_count > 1:
    print( current_triple.strip() )



