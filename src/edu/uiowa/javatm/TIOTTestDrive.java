package edu.uiowa.javatm;

import java.util.Arrays;

class TIOTTestDrive {
    public static void main(String[] args) {

    		String basePath, dataname;
   		//basePath = "/Users/zhiyzuo/OneDrive - University of Iowa/Data/journal-article-data/";
    		//dataname = "jasist";

		//basePath = "/Users/zhiyzuo/OneDrive - University of Iowa/Data/Stack-Exchange-Data-Dump/converted-data/";
    		//dataname = "chess";

    		basePath = "/Users/zhiyzuo/OneDrive - University of Iowa/Data/synthetic-data/";
    		dataname = "doc1";

    		String docFilePath = basePath + "WD-" + dataname + ".csv";
    		String vocFilePath = basePath + "vocabulary-" + dataname + ".csv";
    		String timeFilePath = basePath + "DK-index-" + dataname + ".csv";
    		String citationFilePath = basePath + "citation-" + dataname + ".csv";
    				
    		int T = 2;	
    		String distType = "BETA";
    		int numIters = 1000;
    		double alpha = (double) T/50.;
    		double beta = 0.1;
    		double pi = .1;

        long startTime = System.nanoTime();
        //TIOTGibbsSampler totmulti = new TIOTGibbsSampler(T, numIters, alpha, beta, pi, 
		//												docFilePath, vocFilePath, timeFilePath, citationFilePath);
        TIOTGibbsSampler totmulti = new TIOTGibbsSampler(T, numIters, alpha, beta, pi, 
														docFilePath, vocFilePath, timeFilePath, citationFilePath,
														distType);

		totmulti.fit();
        TIOTOutcome outcome = (TIOTOutcome) totmulti.get_outcome() ;

        long endTime = System.nanoTime();
        System.out.println("Overall elapsed time: " + (endTime - startTime)/(1e-9) + "s");

        totmulti.showTopics(12);
        outcome.showTopicDistribution();
        double[][] psi = outcome.getPsi();
        for (int i = 0; i < psi.length; i++) {
			System.out.println(Arrays.toString(psi[i]) + ",");
		}
        
        for (int t = 0; t < T; t++) {
			System.out.println(Arrays.deepToString(outcome.getGa()[t]) + ",");
		}
    }    
}
