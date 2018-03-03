//For this experiment you will need to use some java functions. Specifically, System.currentTimeMillis() and Random.generate() and Random.nextlong(). The SealedDES.java program I give you in part 2 has examples. You will also need join() all your child threads from the parent thread. This is the most basic form of synchronization.
import java.util.concurrent.ThreadLocalRandom;

public class CoinFlipTest implements Runnable{
    private long headCount;
    long numThreads;
    long numFlips;
    private long tailCount;

    // Constructor for Flip
    public CoinFlipTest(long nT, long nF) {
        this.numThreads = nT;
        this.numFlips = nF;
        this.headCount = 0;
        this.tailCount = 0;
    }

    public long getHeads() {
        return this.headCount;
    }

    public long getTails() {
        return this.tailCount;
    }

    public void run() {
        for (long i = 0; i < numFlips/numThreads; i++) {
            if (ThreadLocalRandom.current().nextInt(2) == 1) {
                ++this.headCount;
            } else {
                ++this.tailCount;
            }
        }
    }

    public static void main (String[] args) {
        // Error if user doesn't put in the right arguments
        if (args.length != 2) {
            System.out.println("Usage: CoinFlipTest #threads #iterations");
            return;
        }
        long nT = Long.parseLong(args[0]);
        long nF = Long.parseLong(args[1]);

        long numHeads = 0;
        long numTails = 0;

        // Create array of threads and start timer
        Thread[] threads = new Thread[(int)nT];
        CoinFlipTest[] flippers = new CoinFlipTest[(int)nT];
        long startTime = System.currentTimeMillis();
        // Start threads
        for (int i = 0; i < nT; i++) {
            flippers[i] = new CoinFlipTest(nT, nF);
            threads[i] = new Thread(flippers[i]);
            threads[i].start();
        }

        // Join threads
        for (int i = 0; i < nT; i++) {
            try {
                threads[i].join();
                numTails += flippers[i].getTails();
                numHeads += flippers[i].getHeads();
            } catch (InterruptedException e) {
                e.printStackTrace();
            }
        }
        long runTime = System.currentTimeMillis() - startTime;

        // Print out information
        System.out.println(numHeads + " heads in " + nF + " coin tosses");
        System.out.println("Elapsed time: " + runTime);
    }
}