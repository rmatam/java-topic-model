package edu.uiowa.javatm;

import edu.uiowa.javatm.ATGibbsSampler;
import edu.uiowa.javatm.ATOutcome;

class ATTestDrive {
    public static void main(String[] args) {

    		String docFilePath = "/Users/zhiyzuo/documents/java-topic-model/data/WD-doc1.csv";
    		String vocFilePath = "/Users/zhiyzuo/documents/java-topic-model/data/vocabulary-doc1.csv";
    		String authFilePath = "/Users/zhiyzuo/documents/java-topic-model/data/AS-doc1.csv";
    		String authArrayFilePath = "/Users/zhiyzuo/documents/java-topic-model/data/authors-doc1.csv";
    				
    		int T = 2;	
    		int numIters = 100;
    		double alpha = (double) T/50.;
    		double beta = 0.1;

        long startTime = System.nanoTime();
        ATGibbsSampler at = new ATGibbsSampler(T, numIters, alpha, beta, 
											  docFilePath, vocFilePath,
											  authFilePath,
											  authArrayFilePath);

        at.fit();

        ATOutcome outcome = (ATOutcome) at.get_outcome();
        long endTime = System.nanoTime();
        System.out.println("Overall elapsed time: " + (endTime - startTime)/1000000000. + "s");

        at.showTopics(12);
        outcome.showTopicDistribution();
        outcome.showATDistribution(2);

        //System.out.println(Arrays.deepToString(lda.getPhi()));

        // save result
        //Utils.write2DArray(outcome.phi, "/Users/zhiyzuo/git/java-topic-model/data/phi-sua-lda.csv");
        //Utils.write2DArray(outcome.theta, "/Users/zhiyzuo/git/java-topic-model/data/theta-sua-lda.csv");
    }    
}

