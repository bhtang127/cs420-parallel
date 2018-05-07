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

#include <cmath>
#include <vector>
#include <set>
#include <algorithm>
#include <random>
#include <cassert>

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