package edu.uiowa.javatm;

import java.util.HashMap;

public class ToTMultiOutcome extends ToTOutcome {

    public ToTMultiOutcome(double[][] t, double[][] p, double[][] k, 
						  HashMap<Integer, Double> tp, 
						  String[][] tRep) {
        super(t, p, k, tp, tRep);
    }

}
