package edu.uiowa.javatm;

import java.util.HashMap;

public class TIOTOutcome extends TMOutcome {

	private double[][] psi;
	private double[][][] ga;

    public TIOTOutcome(double[][] t, double[][] p, 
    					  double[][] k, double[][][] ga,
					  HashMap<Integer, Double> tp, 
					  String[][] tRep) {
        super(t, p, tp, tRep);
        this.psi = k;
        this.ga = ga;
    }

    public double[][] getPsi() {
        return psi;
    }

    public double[][][] getGa() {
        return ga;
    }
    
    public void printGa() {
    		for (int t = 0; t < this.ga.length; t++) {
    			System.out.println();
    			for (int k = 0; k < this.ga[k].length; k++) {
    				
			}
		}
    }

}
