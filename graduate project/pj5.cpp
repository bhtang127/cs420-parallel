#include "pre_sequencing.h"
#include "mpi.h"
#include <iostream>
#include <time.h>


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
    uint32 n_well = 8, n_molecule = 20000;
    uint32 cycle1 = 30, cycle2 = 15;
    uint32 K = 4e6, S = 2e8;
    uint32 bases_per_amplicon = 33;
    double error_rate = 1e-6;

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
