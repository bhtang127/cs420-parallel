#! /bin/bash

echo "num_threads total_fips heads elapsed" > $1;

ITER=0
while [ $ITER -lt 20 ]; do
    let ITER=ITER+1;
    NUMTHREADS=1;
    while [ $NUMTHREADS -le $3 ]; do
        java $2 $NUMTHREADS $4 | grep -Po "\d+" | awk -v RS='^$' "{print $NUMTHREADS,\$2,\$1,\$3}" >> $1;
        let NUMTHREADS=NUMTHREADS*2;
    done
done

