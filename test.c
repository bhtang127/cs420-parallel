#include <stdio.h>
#include <string.h>
#include <stdlib.h>

int main(){
    int i, DATA_LEN = 10;
    unsigned int *output_array;
    output_array = (unsigned int*) malloc ( DATA_LEN * sizeof(unsigned int));
    memset ( output_array, 0, DATA_LEN * sizeof(unsigned int));
    for(i=0;i<DATA_LEN;i++){
        printf("%d ",output_array[i]);
    }
    free(output_array);
    return 0;   
}