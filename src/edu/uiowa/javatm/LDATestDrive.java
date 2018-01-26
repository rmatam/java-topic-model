package edu.uiowa.javatm;

import edu.uiowa.javatm.LDAGibbsSampler;

class LDATestDrive {
    public static void main(String[] args) {

    		String basePath, dataname;
   		//basePath = "/Users/zhiyzuo/OneDrive - University of Iowa/Data/journal-article-data/";
    		//dataname = "jasist";

		//basePath = "/Users/zhiyzuo/OneDrive - University of Iowa/Data/Stack-Exchange-Data-Dump/converted-data/";
    		//dataname = "jasist";

    		basePath = "/Users/zhiyzuo/OneDrive - University of Iowa/Data/synthetic-data/";
    		dataname = "doc1";


    		String docFilePath = basePath + "WD-" + dataname + ".csv";
    		String vocFilePath = basePath + "vocabulary-" + dataname + ".csv";
    				
    		int T = 2;	
    		int numIters = 1000;
    		double alpha = T/50.;
    		double beta = 0.1;

        long startTime = System.nanoTime();
        LDAGibbsSampler lda = new  LDAGibbsSampler(T, numIters, alpha, beta, 
												 docFilePath, vocFilePath);

		lda.fit();
        //LDAOutcome outcome = (LDAOutcome) lda.fit();
        //LDAOutcome outcome = (LDAOutcome) lda.get_outcome() ;
		LDAOutcome outcome = (LDAOutcome) lda.get_outcome();
        long endTime = System.nanoTime();
        System.out.println("Overall elapsed time: " + (endTime - startTime)/1000000000. + "s");

        lda.showTopics(12);
        outcome.showTopicDistribution();
        //System.out.println(Arrays.deepToString(lda.getPhi()));

        // save result
        //Utils.write2DArray(outcome.phi, "/Users/zhiyzuo/documents/java-topic-model/data/phi-sua-lda.csv");
        //Utils.write2DArray(outcome.theta, "/Users/zhiyzuo/documents/java-topic-model/data/theta-sua-lda.csv");
    }    
}
