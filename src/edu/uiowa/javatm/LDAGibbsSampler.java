package edu.uiowa.javatm;

import edu.uiowa.javatm.Utils;
import edu.uiowa.javatm.LDAOutcome;

import java.util.Arrays;
import java.util.Random;
import java.util.stream.IntStream;

import org.apache.commons.math3.distribution.EnumeratedDistribution;
import org.apache.commons.math3.distribution.EnumeratedIntegerDistribution;
import org.apache.commons.math3.util.FastMath;

import com.sun.corba.se.impl.javax.rmi.CORBA.Util;

public class LDAGibbsSampler extends TMGibbsSampler {

	 protected int[][] n_d_t;
	
    /**
     * Minimal constructor for LDAGibbsSampler Class
     * @param doc Filename of the corpus of interest. The first column of doc should be the total number of tokens. 
     * Starting from the second row, the first column should be the current token's document index;
     * the second being its word index in the vocabulary array.
     * @param voc Filename of the vocabulary array. Each row contains a vocabulary term
     */
	
    public LDAGibbsSampler(String doc, String voc) {
        int[][] WD = Utils.readDocFile(doc);
        // 1st Column: Document Index DS
        this.DS = Utils.getColumn(WD, 0);
        // 2nd Column: Word Index (as in the vocabulary) WS
        this.WS = Utils.getColumn(WD, 1);
        // # of docs
		this.D = Arrays.stream(this.DS).max().getAsInt()+1;
        // Read vocabulary from text file
        this.vocArray = Utils.readVocabulary(voc);
        this.W = vocArray.length;
        // The total number of tokens
        this.NNZ = DS.length;
    }
    
    /**
     * Constructor for LDAGibbsSampler Class with detailed parameters
     * @param T The number of topics wanted.
     * @param numIters The number of iterations to run.
     * @param alpha Dirichlet prior for document over topic multinomial
     * @param beta Dirichlet prior for topic over word multinomial
     * @param doc Filename of the corpus of interest. The first column of doc should be the total number of tokens. 
     * Starting from the second row, the first column should be the current token's document index;
     * the second being its word index in the vocabulary array.
     * @param voc Filename of the vocabulary array. Each row contains a vocabulary term
     */
    public LDAGibbsSampler(int T, int numIters,
						  double alpha, double beta, 
						  String doc, String voc) {
        this(doc, voc);

        this.setNumTopics(T);
        this.setNumIters(numIters);
        this.setAlpha(alpha);
        this.setBeta(beta);
    }
    
    @Override
    public String getModelType() {
		return "LDA";
    }
    
    @Override
    public TMOutcome get_outcome() {
        LDAOutcome outcome = new LDAOutcome(this.theta, this.phi, this.topics, this.topicsRep);
        return outcome;
    }

    @Override
    void _est_multinomial() {
        this.theta = new double[this.D][this.T];
        this.phi = new double[this.T][this.W];

        for (int t = 0; t < this.T; t++) {
            for (int w = 0; w < this.W; w++) {
                this.phi[t][w] = (this.beta + this.n_t_w[t][w]) /
								(this.W*this.beta + Utils.sum(this.n_t_w[t]));
            }
            for (int d = 0; d < this.D; d++) {
                this.theta[d][t] = (this.alpha + this.n_d_t[d][t]) /
								  (this.T*this.alpha + Utils.sum(n_d_t[d]));
            }

			_getTopicProportion();
			_getTopicRepresentation();
        }
	}

    /* Run Gibbs Sampling on the given input */
    @Override
    protected void _sample() {
        double[] pi;
        double pi_cum, u;
        int di, wi, zi;
        int zi_sample = -1;
		// Sequentially sample for each token
		for (int i = 0; i < NNZ; i++) {
			di = this.DS[i]; // doc index
			wi = this.WS[i]; // voc index
			zi = this.Z[i]; // topic index
			// discount
			this.n_d_t[di][zi]--;
			this.n_t_w[zi][wi]--;
			this.n_t[zi]--;

			// pi is cumulative;
			pi = new double[T];
			pi_cum = 0;
			for (int t = 0; t < this.T; t++) {
				pi_cum += (this.n_d_t[di][t] + this.alpha) * 
					      (this.n_t_w[t][wi] + this.beta) / (n_t[t] + this.W * this.beta);
			    pi[t] = pi_cum;
			}
			
			zi_sample = Utils.sample1d(pi);
			
	        assert zi_sample != -1;
	        
			this.Z[i] = zi_sample;
			this.n_d_t[di][zi_sample]++;
			this.n_t_w[zi_sample][wi]++;
			this.n_t[zi_sample]++;
		}
	}
        

    @Override
	protected void _init() {
        this.n_d_t = new int[this.D][this.T];
        this.n_t_w = new int[this.T][this.W];
        this.n_t = new int[this.T];
        //this.randGen = new Random();

        // random init
        final Random random = new Random();
        this.Z = random.ints(0, T).limit(NNZ).toArray();
        
        int di, wi, zi;
		for (int i = 0; i < this.NNZ; i++) {
			di = this.DS[i]; // doc index
			wi = this.WS[i]; // voc index
			zi = this.Z[i]; // topic index

			this.n_d_t[di][zi]++;
			this.n_t_w[zi][wi]++;
			this.n_t[zi]++;
		}
	}

}


