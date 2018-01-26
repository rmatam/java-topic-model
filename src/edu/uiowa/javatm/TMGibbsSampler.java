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

import edu.uiowa.javatm.Utils;
import edu.uiowa.javatm.TMOutcome;

import java.util.Arrays;
import java.util.HashMap;

import java.text.DecimalFormat;

class TMGibbsSampler {

    /* private variables */
    // {{{
    protected int T = 100; // Number of topics
    protected int K; // Numer of timestamps
    protected int D; // Number of iterations for sampling
    protected int W; // Number of iterations for sampling
    protected int NNZ; // Total number of tokens
    protected int numIters = 1000; // Number of iterations for sampling

    // Dirichlet priors to control the shapes of document-topic and topic-word , respectively
    protected double alpha = 0.01;
    protected double beta = 0.01;

    // Vocaublary
    protected String[] vocArray;
    // Array of word indices and vocabulary indices
    protected int[] WS, DS;
    protected int[] ND;

    // Results of Gibbs Sampling
    protected double[][] theta;// theta: document/author over topic distributions - numDocs x  T
    protected double[][] phi;// phi: topic over word distributions - T x numWords

    // Samples
	int[] Z, n_t; // Z: topic assignment for each term; n_t: number of each topic
	int[][] n_t_w; // n_t_w: topic t with term w;

    // topics: a key-value pair with topic # as keys and topic proportion as values
    protected HashMap<Integer, Double> topics;
    // topicsRep: an array of strings for sorted keywords with 
    // corresponding probabilities in each topic
    protected String[][] topicsRep;
    // Result object using a self-defined class LDAOutcome
    protected TMOutcome outcome;
    
    
    // }}}

    public TMGibbsSampler() {}
    
    public String getModelType() {
		return "Generic TM";
    }

    public String getDistributionType() {
		return "N/A";
    }

    /* Getter methods for T, numIters, alpha, beta, and burnin */
    // {{{
    public int getNumTopics() {
        return T;
    }

    public int getNumIters() {
        return numIters;
    }

    public double getAlpha() {
        return alpha;
    }

    public double getBeta() {
        return beta;
    }

    public double[][] getTheta() {
        return this.theta;
    }

    public double[][] getPhi() {
        return this.phi;
    }
    // }}}

    /* Setter methods for T, numIters, alpha, beta, burnin, docFileName, and vocFileName */
    // {{{
    public void setNumTopics(int nt) {
        this.T = nt;
    }

    public void setNumIters(int ni) {
        this.numIters = ni;
    }

    public void setAlpha(double a) {
        this.alpha = a;
    }

    public void setBeta(double b) {
        this.beta = b;
    }


    // }}}

    /* Printing topic proportions and definitions (distribution over words)*/
    // Show top k words of (with regard to proportions) topics
    public void showTopics(int k) {
        for (int t = 0; t < T; t++) {
            System.out.println("Topic #" + t);  
            System.out.println(Arrays.toString(
                        Arrays.copyOfRange(topicsRep[t], 0, k)));
            System.out.println("---------------------------");  
        }
    }
    // Show all topics (in sorted order)
    public void showTopics() {
        showTopics(W);
    }

    /* Obtain overall topic distirbution  */
    protected void _getTopicProportion() {
    // {{{
        topics = new HashMap<Integer, Double>();

        for (int t = 0; t < T; t++) {
            topics.put(t, (double) this.n_t[t]/this.NNZ);
        }
     //}}}
    }

    /* Obtain topic representation: from most probable word to least */
    protected void _getTopicRepresentation() {
    // {{{
        DecimalFormat df = new DecimalFormat("#.####");
        this.topicsRep = new String[this.T][this.W];
        double[] topicOverWordArray;
        int[] sortedIndices;
        for (int t = 0; t < T; t++) {
            topicOverWordArray = Arrays.copyOf(this.phi[t], this.phi[t].length);
            sortedIndices = Utils.sortIndex(topicOverWordArray);
            for (int ii = 0; ii < this.W; ii++) {
                this.topicsRep[t][ii] = vocArray[sortedIndices[ii]] + "("
                    + df.format(this.phi[t][sortedIndices[ii]]) + ")";
            }
        }
    // }}}
    }
    
    /* Count the number of words for each document*/
    protected void _count_doc_size() {
    		this.ND = new int[this.D];
    		for (int di : DS) {
    			this.ND[di]++;
		}
	}
    
    /* Run Gibbs sampling on the given input for one iteration*/
    protected void _sample() {}

    /* Initialize before MCMC sampling*/
    protected void _init() {}

    /* Run Gibbs sampling */
	protected void _fit() {
		
		System.out.print(this.D + " documents; ");
		System.out.print(this.NNZ + " tokens; ");
		System.out.print(this.W + " unique tokens; ");
		if (this.getModelType().equals("ToTMulti") | this.getModelType().equals("TIOT")) {
			System.out.print(this.K + " unique timestamps; ");
		}
		System.out.print(this.T + " topics; ");
		if (this.getModelType().equals("TIOT")) {
			System.out.println("Using " + this.getDistributionType() + " distribution");
		}
		//System.out.println(this.numIters + " iterations; ");
		System.out.println("Start fitting " + this.getModelType() + " model for " + this.numIters + " iterations.");
		
		long endTime;
        long startTime = System.nanoTime();
		this._init();
        for(int iter = 1; iter <= this.numIters; iter++) {
			this._sample();
			if(iter == 1 | iter % 10 == 0) {
				endTime = System.nanoTime();
				System.out.println("Iteration " + iter + " (elapsed time: " + (endTime - startTime)/(1e9) + "s)");
			}
        }
        endTime = System.nanoTime();
        System.out.println("MCMC fit elapsed time: " + (endTime - startTime)/(1e9) + "s");
	}

    void _est_multinomial() {}

    public void fit() {
        this._fit();
        this._est_multinomial();
    }

    public TMOutcome get_outcome() {return null;}

}


