package edu.uiowa.javatm;

import java.util.Arrays;

public class ToTTestDrive {
    public static void main(String[] args) {

    				
    		String basePath, dataname;
   		//basePath = "/Users/zhiyzuo/OneDrive - University of Iowa/Data/journal-article-data/";
    		//dataname = "jasist";

		//basePath = "/Users/zhiyzuo/OneDrive - University of Iowa/Data/Stack-Exchange-Data-Dump/converted-data/";
    		//dataname = "chess";

		//basePath = "/Users/zhiyzuo/OneDrive - University of Iowa/Data/State-of-the-Union-Data/converted-data/";
    		//dataname = "sua";

    		basePath = "/Users/zhiyzuo/OneDrive - University of Iowa/Data/synthetic-data/";
    		dataname = "doc1";


    		String docFilePath = basePath + "WD-" + dataname + ".csv";
    		String vocFilePath = basePath + "vocabulary-" + dataname + ".csv";
    		String timeFilePath = basePath + "DK-" + dataname + ".csv";

    		int T = 2;	
    		int numIters = 1000;
    		double alpha = (double) T/50.;
    		double beta = .1;

        ToTGibbsSampler tot = new ToTGibbsSampler(T, numIters, 
												 alpha, beta, 
												 docFilePath, vocFilePath,
												 timeFilePath);

        tot.fit();
        ToTOutcome outcome = (ToTOutcome) tot.get_outcome();

        tot.showTopics(10);
        outcome.showTopicDistribution();
        System.out.println(Arrays.deepToString(outcome.getPsi()));
    }    

}
