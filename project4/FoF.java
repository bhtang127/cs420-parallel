// Java implement for MapReduce for friends trangle 

import java.io.IOException;
import java.util.*;

import org.apache.hadoop.conf.Configuration;
import org.apache.hadoop.fs.Path;
import org.apache.hadoop.io.*;
import org.apache.hadoop.mapreduce.Job;
import org.apache.hadoop.mapreduce.Mapper;
import org.apache.hadoop.mapreduce.Reducer;
import org.apache.hadoop.mapreduce.lib.input.FileInputFormat;
import org.apache.hadoop.mapreduce.lib.output.FileOutputFormat;

public class FoF {

	public static class FoFMapper extends Mapper<LongWritable, Text, Text, IntWritable> {
		private final static IntWritable one = new IntWritable(1);
		private Text word = new Text();

		public void map(LongWritable key, Text value, Context context) throws IOException, InterruptedException 
		{
            //  Emit (Key, Value) Pair
            //  Key: (A B C) is all the tuple such that B < C  
            //        A B C all different and appear simutaneously in one file
            //        and one of A B C is the name for this file 
            //  Value: 1
			String line = value.toString();
            String [] IDs = line.split("\\s+");
            // An aux array for output, but only cost 3 int
            int[] triple = new int[3];
			if(IDs.length > 2){
                for(int i = 1; i < IDs.length; i++){
                    for(int j = i+1; j < IDs.length; j++){
                        triple[0] = Integer.parseInt( IDs[ 0 ] );
                        triple[1] = Integer.parseInt( IDs[ i ] );
                        triple[2] = Integer.parseInt( IDs[ j ] );
                        Arrays.sort(triple);
                        
                        word.set( triple[0] + " " + triple[1] + " " + triple[2] );
                        context.write(word, one);
                        
                        word.set( triple[1] + " " + triple[0] + " " + triple[2] );
                        context.write(word, one);
                        
                        word.set( triple[2] + " " + triple[0] + " " + triple[1] );
                        context.write(word, one);
                    }
                }
            }
        }
    }
	
	

	public static class FoFReducer extends Reducer<Text, IntWritable, Text, NullWritable> {
		private final static NullWritable nw = NullWritable.get();
		
        public void reduce(Text key, Iterable<IntWritable> values, Context context) throws IOException, InterruptedException
        {
            //  Reducer is simply doing counting job
            //  Since mapper have already out put all
            //  possible legal triples of mutual friends suggest by each friends
            //  If a triple appear more than once, 
            //  then they are truly mutual friends and vice versa 
			int count = 0;
            for (IntWritable val : values) {
                count += val.get();
            }
			if(count > 1){
				context.write(key, nw);
			}
		}
	}
	

	public static void main(String[] args) throws Exception {
        
        Configuration conf = new Configuration();
        Job job = Job.getInstance(conf, "Friends Triangle");

        job.setJarByClass(FoF.class);
        job.setMapperClass(FoFMapper.class);
        job.setReducerClass(FoFReducer.class);

        job.setMapOutputKeyClass(Text.class);
        job.setMapOutputValueClass(IntWritable.class);
        
        job.setOutputKeyClass(Text.class);
        job.setOutputValueClass(NullWritable.class);
        
        FileInputFormat.addInputPath(job, new Path(args[0]));
        FileOutputFormat.setOutputPath(job, new Path(args[1]));
        
        System.exit(job.waitForCompletion(true) ? 0 : 1);
	}
}

