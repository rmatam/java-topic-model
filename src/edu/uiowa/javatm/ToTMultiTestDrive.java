package edu.uiowa.javatm;

import java.util.Arrays;

class ToTMultiTestDrive {
    public static void main(String[] args) {

    		String docFilePath = "/Users/zhiyzuo/OneDrive - University of Iowa/java-topic-model/data/WD-doc1.csv";
    		String vocFilePath = "/Users/zhiyzuo/OneDrive - University of Iowa/java-topic-model/data/vocabulary-doc1.csv";
    		String timeFilePath = "/Users/zhiyzuo/OneDrive - University of Iowa/java-topic-model/data/DK-index-doc1.csv";
    				
    		int K = 4;

    		int T = 30;	
    		int numIters = 1000;
    		double alpha = 50./T;
    		double beta = 0.1;
    		double pi = 1;

        long startTime = System.nanoTime();
        ToTMultiGibbsSampler totmulti = new ToTMultiGibbsSampler(T, K, numIters, alpha, beta, pi, 
																docFilePath, vocFilePath, timeFilePath);

		totmulti.fit();
        ToTMultiOutcome outcome = (ToTMultiOutcome) totmulti.get_outcome() ;

        long endTime = System.nanoTime();
        System.out.println("Overall elapsed time: " + (endTime - startTime)/(1e9) + "s");

        totmulti.showTopics(12);
        outcome.showTopicDistribution();
        System.out.println(Arrays.deepToString(outcome.getPsi()));
    }    
}
