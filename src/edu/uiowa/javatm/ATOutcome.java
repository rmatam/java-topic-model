package edu.uiowa.javatm;

import java.util.Arrays;
import java.util.HashMap;

public class ATOutcome extends LDAOutcome{

	protected String[] authArray;

	public ATOutcome(double[][] t, double[][] p, 
					HashMap<Integer, Double> tp, 
					String[][] tRep,
					String[] authArray) {
		super(t, p, tp, tRep);
		this.authArray = authArray;
	}
	
	public void showATDistribution(int top_t) {
		System.out.println("Name\tTopic");
		int j; // top topic indices
		int[] sortedIndices;
		double[] authorTheta;
		for(int a = 0; a < this.authArray.length; a++) {
			System.out.print(this.authArray[a]+": ");
			authorTheta = Arrays.copyOf(this.theta[a], this.theta[a].length);
			sortedIndices = Utils.sortIndex(authorTheta);
			for(int ii = 0; ii < top_t; ii++) {
				j = sortedIndices[ii];
				System.out.print("Topic " + j + ": (" + authorTheta[j] + "); ");
			}
			System.out.println("");
		}
	}
}

