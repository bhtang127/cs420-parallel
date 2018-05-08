/*  To Do:
    Use the function below to finish pure MPI(therrefore sequential) version of the project.
    Input: Length N vector with value {0,1} to indicate the original mutation state for molecules
    Procedure:
        1. Divide N molecule to M wells randomly but evenly
        2. For each well, say well i has ni molecules, then use function "sample" to calculate how many
           duplications should each molecule produce. And use function "PCR" to generate those molecules 
        3. Combine all wells, then you should get K*M molecules in total, then use function "sample" again
           to compute #duplications each molecule should produce. Then use function "PCR" to generate 
           those molecules. Now you should have S molecules, each should have tag for well it comes from 
           and if it is mutated (mutation state will be generate in function "PCR" but well id should be 
           attached by your own)
        4. Calculate the mutation rate --the ratio of mutated molecules and total molecules after step 3
    Output: Double type mutation rate
    Comments: We will do a 8 wells case, so you shouldn't parallel just at well level if num of threads
              is greater than 8 (go a little further to molecule level) 
              Also, the code below should add --std=c++11 to compile

    Waring: if a molecule is mutated before PCR, then it can only produce mutated molecule during PCR. 
            the result of PCR function actually only have meaning when the molecule is not mutated before.  
*/

#include "mpi.h"
#include <cmath>
#include <vector>
#include <set>
#include <algorithm>
#include <random>
#include <iostream>
#include <cassert>
#include <time.h>

using uint32 = unsigned int;

// Sample function
std::vector<uint32> Sample(uint32 ni, uint32 K){
    // This is a function to place K objects to ni place 
    // randomly and return a vector to indicate num of objects
    // in each place
    std::vector<uint32> retval(ni);
    std::random_device rd;
    std::mt19937 mt(rd());
    std::uniform_int_distribution<uint32> dist(0, ni-1);
    for(uint32 i=0; i < K; i++){
        retval[dist(mt)]++;
    }
    return retval;
}

///////////////////////////////////////////////////////////
// Auxiliary functions start
///////////////////////////////////////////////////////////

// Main Randomness
class Random {
public:
    Random() = default;
    Random(std::mt19937::result_type seed) : eng(seed) {}
    uint32 DrawInt(uint32 Max);
    uint32 Poisson(double mean);

private:        
    std::mt19937 eng{std::random_device{}()};
};

uint32 Random::DrawInt(uint32 Max)
{
    return std::uniform_int_distribution<uint32>{0, Max - 1}(eng);
}

uint32 Random::Poisson(double mean){
    std::poisson_distribution<uint32> dist(mean);
    return dist(eng);
}


// Aux function for sampleID
uint32 find_index(uint32 cycles, std::set<uint32> indexset, uint32 randi){
    uint32 index = 0, count = 0;
    for(; index < (1<<cycles); index++){
        if(indexset.count(index) != 0) count++;
        if(count == (randi+1)) return index; 
    }
}

// Function to sample ID without replacement from a large range
std::vector<uint32> sampleID(uint32 cycles, uint32 n_samples, Random rd){
    std::vector<uint32> retval(n_samples);
    std::set<uint32> maintain;
    uint32 randi, ind;
    if(n_samples == 1){
        retval[0] = rd.DrawInt(1<<cycles);
    }
    else if(cycles < 15){
        std::set<uint32> indexes;
        for(uint32 j=0; j < n_samples; j++){
            randi = rd.DrawInt((1<<cycles) - j);
            ind = find_index(cycles, indexes, randi);
            indexes.insert(ind);
            retval[j] = ind;
        }
        std::sort(retval.begin(), retval.end());
    }
    else{
        while(maintain.size() < n_samples){
            maintain.insert( rd.DrawInt(1<<cycles) );
        }
        retval.assign( maintain.begin(), maintain.end() );
    }
    return retval;
}

std::vector< std::set<uint32> > find_all_ancestors(std::vector<uint32> ids, 
                                                   uint32 cycles){
    std::vector< std::set<uint32> > ancestors(cycles);
    for(uint32 i=0; i<cycles; i++){
        for(uint32 j=0; j<ids.size(); j++){
            ancestors[i].insert( (uint32) (ids[j] / (1<<i)) );
        }
    }
    return ancestors;
}
////////////////////////////////////////////////////////////
// Auxiliary functions end
////////////////////////////////////////////////////////////


// Duplication function
std::vector<uint32> PCR(uint32 target, uint32 cycles, 
                      uint32 bases_per_amplicon, double error_rate){
    // When you know a molecule should produce "target" duplications
    // in the procedure, you can use this functions to generate them 
    // Return a vector of uint32 to indicate the mutation state of those
    // duplications. 1 for mutated and 0 not
    // "cycles, bases_per_amplicon, error_rate" are just parameters
    // Warning: target here should > 0

    std::vector<uint32> retval(target);

    // We do a simplified duplicate procedure
    // Actually this is not right
    // But the amount of computation is similar

    std::vector<uint32> ids;
    std::vector< std::set<uint32> > ancestors;
    uint32 ancsum;
    double poisson_mean;
    uint32 mutation;

    Random rd;

    ids = sampleID(cycles, target, rd);
    ancestors = find_all_ancestors(ids, cycles);
    ancsum = 0;
    for(uint32 k=0; k < ancestors.size(); k++){
        ancsum += ancestors[k].size();
    }
    poisson_mean = ancsum * bases_per_amplicon * error_rate;

    mutation = std::min(target, rd.Poisson(poisson_mean));
    for(uint32 i=0; i < mutation; i++){
        retval[i] = 1;
    }
    
    return retval;
}


// Main function
void process(uint32 &len, std::vector<uint32> &molecules,
             std::vector<uint32> &ids, std::vector<uint32> &ni,
             uint32 cycles, uint32 bases_per_amplicon, uint32 error_rate,
             int rank, int num_procs, uint32 n_well, uint32 K, uint32 S, int tag,
             MPI_Status &stat){
    std::vector<uint32> pcr_ret, idtag, ret_mole, ret_ids;
    uint32 accum = 0;

    if ( rank == 0 ) {
        // processing duplication
        for ( int i=0; i < len; i++ ) {
            if ( ni[i] == 0 ) continue;
            else if ( molecules[i] ) {
                pcr_ret = std::vector<uint32>( ni[i], 1 );
            }
            else {
                pcr_ret = PCR(ni[i], cycles, bases_per_amplicon, error_rate);
            }
            idtag = std::vector<uint32>( ni[i], ids[i] );
            ret_mole.insert(ret_mole.end(), pcr_ret.begin(), pcr_ret.end());
            ret_ids.insert(ret_ids.end(), idtag.begin(), idtag.end());
        } 

        // Combine work for PCR round
        molecules.clear();
        ids.clear();
        if ( tag == 1 ){
            molecules.resize(n_well * K);
            ids.resize(n_well * K);
        }
        else{
            molecules.resize(S);
            ids.resize(S);
        }
        
        for ( int i=0; i < ret_mole.size(); i++ ) {
            molecules[i] = ret_mole[i];
            ids[i] = ret_ids[i];
        }
        accum = ret_mole.size();     
        for ( int r=1; r < num_procs; r++ ) {
            MPI_Recv ( &len, 1, MPI::UNSIGNED, r, 2, MPI_COMM_WORLD, &stat );            
            MPI_Recv ( &molecules[accum], len, MPI::UNSIGNED, r, 2, MPI_COMM_WORLD, &stat );
            MPI_Recv ( &ids[accum], len, MPI::UNSIGNED, r, 2, MPI_COMM_WORLD, &stat );
            accum += len;
        }
    }
    else {
        // Reciving work
        MPI_Recv ( &len, 1, MPI::UNSIGNED, 0, 2, MPI_COMM_WORLD, &stat );
        molecules.resize(len);
        ids.resize(len);
        ni.resize(len);
        MPI_Recv ( &molecules[0], len, MPI::UNSIGNED, 0, 2, MPI_COMM_WORLD, &stat );
        MPI_Recv ( &ids[0], len, MPI::UNSIGNED, 0, 2, MPI_COMM_WORLD, &stat );
        MPI_Recv ( &ni[0], len, MPI::UNSIGNED, 0, 2, MPI_COMM_WORLD, &stat );
        
        // Processing distributed work
        for ( int i=0; i < len; i++ ) {
            if ( ni[i] == 0 ) continue;
            else if ( molecules[i] ) {
                pcr_ret = std::vector<uint32>( ni[i], 1 );
            }
            else {
                pcr_ret = PCR(ni[i], cycles, bases_per_amplicon, error_rate);
            }
            idtag = std::vector<uint32>( ni[i], ids[i] );
            ret_mole.insert(ret_mole.end(), pcr_ret.begin(), pcr_ret.end());
            ret_ids.insert(ret_ids.end(), idtag.begin(), idtag.end());
        } 

        // Send Back
        len = ret_mole.size();
        MPI_Ssend ( &len, 1, MPI::UNSIGNED, 0, 2, MPI_COMM_WORLD );        
        MPI_Ssend ( &ret_mole[0], ret_mole.size(), MPI::UNSIGNED, 0, 2, MPI_COMM_WORLD );
        MPI_Ssend ( &ret_ids[0], ret_ids.size(), MPI::UNSIGNED, 0, 2, MPI_COMM_WORLD );
    }
}

int main ( int argc, char** argv ) {

    // MPI Standard variable
    int num_procs;
    int rank, j;
    clock_t t;

    // Parameters in the design
    // we only consider this unqiue situation
    uint32 n_well = 8, n_molecule = 10000;
    uint32 cycle1 = 30, cycle2 = 15;
    uint32 K = 2e6, S = 1e8;
    uint32 bases_per_amplicon = 33;
    double error_rate = 1e-4;

    // Messaging variables
    MPI_Status stat;
    // necessary variables 
    uint32 len, accum;
    std::vector<uint32> molecules, ret_mole;
    std::vector<uint32> ids, ret_ids, ni;
    std::vector<uint32> pcr_ret;
    std::vector<uint32> idtag;
    std::vector<uint32> mut_count(n_well), total_count(n_well), aux_count(n_well);
    std::vector<double> mutation_rate(n_well);

    // MPI Setup
    if ( MPI_Init( &argc, &argv ) != MPI_SUCCESS )
        {
        printf ( "MPI_Init error\n" );
        }

    MPI_Comm_size ( MPI_COMM_WORLD, &num_procs ); // Set the num_procs
    MPI_Comm_rank ( MPI_COMM_WORLD, &rank );

    // First PCR round 

    // distributing work
    if ( rank == 0 ) {
        t = clock();
        // init status
        molecules.resize(n_molecule);
        ids.resize(n_molecule);
        ni.resize(n_molecule);
        std::vector<uint32> samples;
        for ( int i=0; i < n_well; i++ ) {
            samples = Sample(n_molecule/n_well, K);
            for ( int j=0; j < n_molecule/n_well; j++ ) {
                ni[i * n_molecule/n_well + j] = samples[j];
                ids[i * n_molecule/n_well + j] = i;
            }
        }

        len = n_molecule / num_procs;
        // distribute work to other process
        for ( int r=1; r < num_procs; r++ ) {
            MPI_Ssend ( &len, 1, MPI::UNSIGNED, r, 2, MPI_COMM_WORLD );            
            MPI_Ssend ( &molecules[r * len], len, MPI::UNSIGNED, r, 2, MPI_COMM_WORLD );
            MPI_Ssend ( &ids[r * len], len, MPI::UNSIGNED, r, 2, MPI_COMM_WORLD );
            MPI_Ssend ( &ni[r * len], len, MPI::UNSIGNED, r, 2, MPI_COMM_WORLD );                        
        }
    }

    // Processing work
    process(len, molecules, ids, ni, cycle1, bases_per_amplicon, error_rate,
            rank, num_procs, n_well, K, S, 1, stat);

    MPI_Barrier(MPI_COMM_WORLD);

    // Second PCR round

    // distributing work
    if ( rank == 0 ) {
        // sample to compute work to do
        ni = Sample(K * n_well, S);
        len = molecules.size() / num_procs;

        // distributing work to other process
        for ( int r=1; r < num_procs; r++ ) {
            MPI_Ssend ( &len, 1, MPI::UNSIGNED, r, 2, MPI_COMM_WORLD );            
            MPI_Ssend ( &molecules[r * len], len, MPI::UNSIGNED, r, 2, MPI_COMM_WORLD );
            MPI_Ssend ( &ids[r * len], len, MPI::UNSIGNED, r, 2, MPI_COMM_WORLD );
            MPI_Ssend ( &ni[r * len], len, MPI::UNSIGNED, r, 2, MPI_COMM_WORLD );                        
        }
    }

    // processing work
    process(len, molecules, ids, ni, cycle2, bases_per_amplicon, error_rate,
            rank, num_procs, n_well, K, S, 2, stat);

    MPI_Barrier(MPI_COMM_WORLD);

    // Accumulating mutation rate
    if ( rank == 0 ) {
        len = S / num_procs;
        // Distribute work
        for ( int r=1; r < num_procs; r++ ) {
            MPI_Ssend ( &len, 1, MPI::UNSIGNED, r, 2, MPI_COMM_WORLD );            
            MPI_Ssend ( &molecules[r * len], len, MPI::UNSIGNED, r, 2, MPI_COMM_WORLD );
            MPI_Ssend ( &ids[r * len], len, MPI::UNSIGNED, r, 2, MPI_COMM_WORLD );                    
        }

        // Processing distributed work
        for ( int i=0; i < len; i++ ) {
            total_count[ ids[i] ] ++;
            if ( molecules[i] ) mut_count[ ids[i] ] ++;
        }

        // Combine work
        for ( int r=1; r < num_procs; r++ ) {
            MPI_Recv ( &aux_count[0], n_well, MPI::UNSIGNED, r, 2, MPI_COMM_WORLD, &stat );
            for ( int i=0; i < n_well; i++) 
                mut_count[i] += aux_count[i];
            MPI_Recv ( &aux_count[0], n_well, MPI::UNSIGNED, r, 2, MPI_COMM_WORLD, &stat );
            for ( int i=0; i < n_well; i++) 
                total_count[i] += aux_count[i];
        }

        // Compute mutation rate each well
        for ( int i=0; i < n_well; i++){ 
            mutation_rate[i] = (double) mut_count[i] / total_count[i];
        }
        t = clock() - t;
        std::cout<<"processing time: "<<t<<std::endl;        
    }
    else {
        // Reciving work
        MPI_Recv ( &len, 1, MPI::UNSIGNED, 0, 2, MPI_COMM_WORLD, &stat );
        molecules.clear(); molecules.resize(len);
        ids.clear(); ids.resize(len);
        MPI_Recv ( &molecules[0], len, MPI::UNSIGNED, 0, 2, MPI_COMM_WORLD, &stat );
        MPI_Recv ( &ids[0], len, MPI::UNSIGNED, 0, 2, MPI_COMM_WORLD, &stat );
        
        // Processing distributed work
        for ( int i=0; i < len; i++ ) {
            total_count[ ids[i] ] ++;
            if ( molecules[i] ) mut_count[ ids[i] ] ++;
        }

        // Send Back    
        MPI_Ssend ( &mut_count[0], n_well, MPI::UNSIGNED, 0, 2, MPI_COMM_WORLD );
        MPI_Ssend ( &total_count[0], n_well, MPI::UNSIGNED, 0, 2, MPI_COMM_WORLD );
    }

    MPI_Finalize(); // finalize so I can exit
    return 0;
}
