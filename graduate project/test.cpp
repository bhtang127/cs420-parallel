#include <iostream>
#include <vector>
#include "omp.h"

#pragma omp declare reduction(vec_insert : std::vector<int> : \
                              omp_out.insert(omp_out.end(), omp_in.begin(), omp_in.end())) \
                    // initializer(omp_priv = omp_orig)

int main(){
    std::vector<int> a(10000000);
    std::vector<int> b(10000000, 2);
    std::vector<int> mol(20000000), id, aux, tag;  
    omp_set_dynamic(0);
    omp_set_num_threads(1);
    #pragma omp parallel for 
    for(int i=0; i<a.size(); i++){
        aux = std::vector<int>( b[i], 1 );
        tag = std::vector<int>( b[i], i );
        mol[2*i]=0; mol[2*i+1] = 0;
    }  
}