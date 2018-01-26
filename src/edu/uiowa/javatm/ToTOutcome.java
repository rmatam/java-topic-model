package edu.uiowa.javatm;
import java.util.HashMap;

public class ToTOutcome extends TMOutcome{
	
	double[][] psi;

    public ToTOutcome(double[][] t, double[][] p, 
					 double[][] psi,
					 HashMap<Integer, Double> tp, 
					 String[][] tRep) {

        super(t, p, tp, tRep);

        this.psi = psi;
    }

    public double[][] getPsi() {
        return psi;
    }
}

