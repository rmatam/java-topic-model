/*******************************************************************************
 * The MIT License (MIT)
 *
 * Copyright (c) Zhiya Zuo 2016
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 *******************************************************************************/
package edu.uiowa.javatm;

import java.io.IOException;

import org.apache.commons.math3.distribution.UniformRealDistribution;

class UtilsTestDrive {
    public static void main(String[] args) throws IOException {

    	
    		double tau = 0.01;
		UniformRealDistribution uniform = new UniformRealDistribution(-1*tau, tau);
		System.out.println(uniform.sample());

        /*
         * 
    		int a = +0;
    		int b = +0;
    		System.out.println(a/b);

    		String asFile = "/Users/zhiyzuo/Documents/java-topic-model/data/AS.csv";
    		ArrayList<int[]> AS = Utils.readAS(asFile);
    		System.out.println(Arrays.deepToString(AS.toArray()));

		// random number genrator
        final Random random = new Random();
        	System.out.println(random.ints(1, 0, 10).toArray()[0]);

    		String asFile = "/Users/zhiyzuo/Documents/java-topic-model/data/AS.csv";
    		ArrayList<int[]> AS = Utils.readAS(asFile);
    		System.out.println(Arrays.deepToString(AS.toArray()));
        /*
         * 
		double[][] mydoubleArray = new double [][] { { 1, 2.},
													{ 0.02, 0.12},
													{ 4.97, 3.56} };
		Utils.write2DArray(mydoubleArray, "/Users/zhiyzuo/Desktop/test.csv");

		String[][] myStringArray = new String [][] { { "X0", "Y0"},
													{ "X1", "Y1"},
													{ "X2", "Y2"},
													{ "X3", "Y3"},
													{ "X4", "Y4"} };
		Utils.write2DArray(myStringArray, "/Users/zhiyzuo/Desktop/test.csv");

        int[][] array1 = {{1,2,3},{4,5,6},{7,8,9}};
        System.out.println("Test getColumn()");
        System.out.println(Arrays.toString(Utils.getColumn(array1, 0)));

        int[] array2 = {1,1,2,3,2,3,2};
        System.out.println("Test find()");
        System.out.println(Utils.find(array2, 2).toString());

        double[] array3 = {0.25, 0.3, 0.15, 0.2, 0.1};
        System.out.println("Test sortIndex(double[])");
        int[] sortedIndices = Utils.sortIndex(array3);
        double[] array3_sorted = new double[array3.length];
        for (int i = 0; i < array3.length; i++) {
            array3_sorted[i] = array3[sortedIndices[i]];
        }
        System.out.println(Arrays.toString(array3_sorted));
        */

        /*
        System.out.println("Test _sum()");
        double sum_ = Utils._sum(array3);
        System.out.println("Sum of array3 is: " + sum_);

        double[] array4 = {0.5, 0.7, 0.8};
        System.out.println("Test _normalize()");
        double[] array4_normed = Utils._normalize(array4);
        System.out.println("Normalized array4 is: " + Arrays.toString(array4_normed));
        */

        /*
        double[] array4 = {0.4, 0.05, 0.1};
        System.out.println("Test sample()");
        System.out.println("Probablity List: " + Arrays.toString(array4));
        int[] choices = new int[10];
        for (int i = 0; i < 10; i++) {
            choices[i] = Utils.sample(array4);
        }
        System.out.println(Arrays.toString(choices));
        /*


        .*
        String vocFileName = "../data/test_mat.csv.clabel";
        System.out.println("Test readVocabulary()");
        String[] vocArray = Utils.readVocabulary(vocFileName);
        System.out.println(Arrays.toString(vocArray));
        */

    	/*
        String docFileName = "../data/WD.csv";
        System.out.println("Test readDocFile()");
        int[][] WD = Utils.readDocFile(docFileName);
        int[] DS = Utils.getColumn(WD, 0);
        int[] WS = Utils.getColumn(WD, 1);
        System.out.println(Arrays.deepToString(WD));
        System.out.println(Arrays.toString(DS));
        System.out.println(Arrays.toString(WS));
        System.out.println(DS[DS.length-1] + 1);
        */
    	
        /*
    	int[][] DA = Utils.readDA("data/DA.csv");
    	System.out.println(Arrays.deepToString(DA));
    	
    	for (int j = 0; j < DA.length; j++) {
    		System.out.println(j + " " + Arrays.toString(Utils.nonZeroIndex(DA[j])));
        }
        */

    	/*
        double[][] citationMatrix = Utils.readCitationFile("../data/doc1Citation.csv");
        System.out.println(Arrays.deepToString(citationMatrix));
        
        double[] doubleArray = new double[] {1, 2, 3, 4};
        System.out.println("Mean of " + Arrays.toString(doubleArray) + " is: " + Utils.mean(doubleArray));
        System.out.println("Variance of " + Arrays.toString(doubleArray) + " is: " + Utils.variance(doubleArray));
        */

    }
}
