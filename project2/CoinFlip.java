////////////////////////////////////////////////////////////////////////
//
//  class: CoinFlip
//  author: Bohao Tang
//  A simple version using Java Threads to flip a coin parallelly
//
////////////////////////////////////////////////////////////////////////

import java.util.Random;

class CoinFlip implements Runnable
{
    int thread_id;
    long flips;  // the total flips number for this thread
    Random generator;
    private long heads;

    public void run ()
    {
        for ( int i=0; i<flips; i++ )
            this.heads += generator.nextInt(2);
    }

    CoinFlip ( int id, long flips )
    {
        this.thread_id = id;
        this.flips = flips;
        this.generator = new Random();
    }

    // get the head counts
    public long getHeads()
    {
        return this.heads;
    }

    public static void main ( String[] args )
    {
        if ( 2 != args.length ) 
        {
            System.out.println ("Usage: CoinFlip #threads #iterations");
            return;
        } 

        int numthreads = Integer.parseInt ( args[0] );
        long totalflips = Long.parseLong ( args[1] );
        long totalheads = 0;

        // Init
        Thread[] threads = new Thread[numthreads];
        CoinFlip[] coins = new CoinFlip[numthreads];

        long runstart = System.currentTimeMillis();

        for ( int i=0; i<numthreads; i++ )
        {
            // Allocate tasks
            coins[i] = new CoinFlip ( i, totalflips/numthreads );
            threads[i] = new Thread ( coins[i] );
            threads[i].start();
        }

        // Add up all counts
        for ( int i=0; i<numthreads; i++ )
        {
            try
            {
                threads[i].join();
                totalheads += coins[i].getHeads();
            }
            catch (InterruptedException e)
            {
                System.out.println("Thread interrupted.  Exception: " + e.toString() +
                                   " Message: " + e.getMessage()) ;
                return;
            }
        }
        
        long elapsed = System.currentTimeMillis() - runstart;
        System.out.println(totalheads + " heads in " + totalflips + " coin tosses");
        System.out.println("Elapsed time: " + elapsed + "ms");
    }

}