package edu.uiowa.javatm;

import edu.uiowa.javatm.Utils;
import edu.uiowa.javatm.ATOutcome;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Random;

public class ATGibbsSampler extends TMGibbsSampler {

	protected int A;
	protected int[] Q; // latent author Q, as with latent topic array Z
	protected int[] n_a; // # of times an author got selected
	protected int[][] n_a_t;
	protected ArrayList<int[]> AS;
	protected String[] authArray;
	protected final Random random;
	
    /**
     * Minimal constructor for ATGibbsSampler Class
     * @param doc Filename of the corpus of interest. 
     * The first column of doc should be the total number of tokens. 
     * Starting from the second row, the first column should be the current token's document index;
     * the second being its word index in the vocabulary array.
     * @param voc Filename of the vocabulary array. Each row contains a vocabulary term
     * @param auth File name of author indices for each word
     * @param authArray File name of the author array. Each row contains a author
     */
	
    public ATGibbsSampler(String doc, String voc, 
						 String auth, String authArray) {
        int[][] WD = Utils.readDocFile(doc);
        // 1st Column: Document Index DS
        this.DS = Utils.getColumn(WD, 0);
        // 2nd Column: Word Index (as in the vocabulary) WS
        this.WS = Utils.getColumn(WD, 1);
        // Author Index (as in the author array) AS
        this.AS = Utils.readAS(auth);
        // # of docs
		this.D = Arrays.stream(this.DS).max().getAsInt()+1;
		//this.D = Arrays.stream(this.DS).max().getAsInt()+1;
        // Read vocabulary from text file
        this.vocArray = Utils.readVocabulary(voc);
        this.W = vocArray.length;
		// # of authors
        this.authArray = Utils.readAuthor(authArray);
		this.A = this.authArray.length;
        // The total number of tokens
        this.NNZ = DS.length;
        
        // init a random number generator for author sampler
        this.random = new Random();
    }
    
    /**
     * Constructor for ATGibbsSampler Class with detailed parameters
     * @param T The number of topics wanted.
     * @param numIters The number of iterations to run.
     * @param alpha Dirichlet prior for author over topic multinomial
     * @param beta Dirichlet prior for topic over word multinomial
     * @param doc Filename of the corpus of interest. The first column of doc should be the total number of tokens. 
     * Starting from the second row, the first column should be the current token's document index;
     * the second being its word index in the vocabulary array.
     * @param voc Filename of the vocabulary array. Each row contains a vocabulary term
     */
    public ATGibbsSampler(int T, int numIters,
						 double alpha, double beta, 
						 String doc, String voc,
						 String auth, String authArray) {
        this(doc, voc, auth, authArray);

        this.setNumTopics(T);
        this.setNumIters(numIters);
        this.setAlpha(alpha);
        this.setBeta(beta);
    }

    @Override
    public String getModelType() {
		return "AT";
    }
    
    @Override
    public TMOutcome get_outcome() {
        ATOutcome outcome = new ATOutcome(this.theta, this.phi, 
										 this.topics, this.topicsRep,
										 this.authArray);
        return outcome;
    }

    @Override
    void _est_multinomial() {
        this.theta = new double[this.A][this.T];
        this.phi = new double[this.T][this.W];

		double n_t_sum;

        for (int t = 0; t < T; t++) {
			n_t_sum = Utils.sum(this.n_t_w[t]);
            for (int w = 0; w < this.W; w++) {
                this.phi[t][w] = (this.beta + this.n_t_w[t][w]) /
								(this.W*this.beta + n_t_sum);
            }
            for (int a = 0; a < this.A; a++) {
                this.theta[a][t] = (this.alpha + this.n_a_t[a][t]) /
								  (this.T*this.alpha + Utils.sum(n_a_t[a]));
            }
        }
        
        _getTopicProportion();
        _getTopicRepresentation();
	}

    
    /* Run Gibbs Sampling on the given input */
    @Override
    protected void _sample() {
        double[] pi;
        int ai, wi, zi;
        int[] awi;
		// Sequentially sample for each token
		for (int i = 0; i < NNZ; i++) {
			wi = this.WS[i]; // voc index
			zi = this.Z[i]; // topic index
			ai = this.Q[i]; // author index
			// discount
			this.n_a_t[ai][zi]--;
			this.n_a[ai]--;
			this.n_t_w[zi][wi]--;
			this.n_t[zi]--;

			awi = this.AS.get(i); // list of authors to choose from
			// choose from the author list uniformally
			ai = awi[random.nextInt(awi.length)];
			this.Q[i] = ai;
			this.n_a[ai]++;

			// pi is cumulative;
			pi = new double[T];
			double pi_cum = 0;
			for (int t = 0; t < this.T; t++) {
				pi_cum += (this.n_a_t[ai][t] + this.alpha) / (n_a[ai] + this.T * this.alpha) * 
					      (this.n_t_w[t][wi] + this.beta) / (n_t[t] + this.W * this.beta);
			    pi[t] = pi_cum;
			}

			// Try to use cumulative method for sampling
	        double u = Math.random() * pi[this.T - 1];
	        for (zi = 0; zi < pi.length; zi++) {
	            if (u <= pi[zi])
	                break;
	        }
			this.Z[i] = zi;
			this.n_a_t[ai][zi]++;
			this.n_t_w[zi][wi]++;
			this.n_t[zi]++;
		}
	}
        

    @Override
	protected void _init() {
        this.n_a_t = new int[this.A][this.T];
        this.n_t_w = new int[this.T][this.W];
        this.n_t = new int[this.T];
        this.n_a = new int[this.A];
        //this.randGen = new Random();

        // random init
        this.Z = random.ints(0, T).limit(NNZ).toArray();
        this.Q = new int[this.NNZ];
        
        int ai, wi, zi;
        int[] awi;
		for (int i = 0; i < this.NNZ; i++) {
			wi = this.WS[i]; // voc index
			zi = this.Z[i]; // topic index
			awi = this.AS.get(i); // list of authors to choose from
			// choose from the author list uniformally
			ai = awi[random.nextInt(awi.length)];
			this.Q[i] = ai;
			this.n_a_t[ai][zi]++;
			this.n_t_w[zi][wi]++;
			this.n_t[zi]++;
			this.n_a[ai]++;
		}
	}

}



