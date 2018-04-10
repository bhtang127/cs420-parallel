#!/usr/bin/env python

"""An Simple Implement of Reducer in Streaming MapReduce for Friends Triangle"""

import sys

count_dict = {}

for triple in sys.stdin:
   
    ## Reducer is simply doing counting job
    ## Since mapper have already out put all
    ## possible legal triples of mutual friends suggest by each friends
    ## If a triple appear more than once, 
    ## then they are truly mutual friends and vice versa  

    # Remove the ending \n
    triple = triple.strip()

    # Counting
    if triple not in count_dict:
        count_dict[triple] = 1
    else:
        count_dict[triple] += 1

## Output result, we output them in order
for triple in sorted( count_dict.iterkeys() ):
    if count_dict[triple] > 1:
        print triple
    
