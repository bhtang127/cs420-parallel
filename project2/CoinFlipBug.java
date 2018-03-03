/**
 * CoinFlipBug 
 */

import java.util.Random;

class CoinFlipBug implements Runnable
{
    int thread_id;
    long flips;
    private long heads;
    static long seed = 1;

    public void run ()
    {
        for ( int i=0; i<flips; i++ )
            if ( (this.next(32) % 2) != 0 )
                ++this.heads;
    }

    CoinFlipBug ( int id, long flips )
    {
        this.thread_id = id;
        this.flips = flips;
    }

    protected synchronized int next(int bits)
    {
        seed = (seed * 0x5DEECE66DL + 0xBL) & ((1L << 48) - 1);
        return (int) (seed >>> (48 - bits));
    }

    public long getHeads()
    {
        return this.heads;
    }

    public static void main ( String[] args )
    {
        if ( 2 != args.length ) 
        {
            System.out.println ("Usage: CoinFlipBug #threads #iterations");
            return;
        } 

        int numthreads = Integer.parseInt ( args[0] );
        long totalflips = Integer.parseInt ( args[1] );
        long totalheads = 0;

        Thread[] threads = new Thread[numthreads];
        CoinFlipBug[] coins = new CoinFlipBug[numthreads];

        long runstart = System.currentTimeMillis();

        for ( int i=0; i<numthreads; i++ )
        {
            coins[i] = new CoinFlipBug ( i, totalflips/numthreads );
            threads[i] = new Thread ( coins[i] );
            threads[i].start();
        }

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