package edu.uiowa.javatm;

import edu.uiowa.javatm.Utils;
import edu.uiowa.javatm.ToTOutcome;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Random;

import org.apache.commons.math3.util.FastMath;

//import org.apache.commons.math3.special.Gamma;
//import org.apache.commons.math3.distribution.UniformRealDistribution;


public class ToTGibbsSampler extends LDAGibbsSampler {

	private double kmin=1e-4, kmax=1-1e-4;

	protected double tau;
	protected double[][] psi; // beta distribution parameters
	//protected BetaDistribution[] betaDist; // array of beta distribution objects

	protected double[] DK; // integer timestamp array
	// array for saving timestamps for each topic. Should be of size T
	protected ArrayList<GlueList<Double>> topicTimeArr; 
	//protected UniformRealDistribution uniform;
	
	/**
	 * Minimal constructor for ToTGibbsSampler Class
	 * @param doc Filename of the corpus of interest. The first column of doc should be the total number of tokens. 
	 * Starting from the second row, the first column should be the current token's document index;
	 * the second being its word index in the vocabulary array.
	 * @param voc Filename of the vocabulary array. Each row contains a vocabulary term
	 * @param dkArray Filename of document time array;
	 */

	public ToTGibbsSampler(String doc, String voc, String dkArray) {
		super(doc, voc);
		
		// Doc-Time
		String[] dkStrArray = Utils.readVocabulary(dkArray);
		this.DK = new double[dkStrArray.length];
		for (int d = 0; d < dkStrArray.length; d++) {
			this.DK[d] = Double.parseDouble(dkStrArray[d]);
		}
		
		// normalize to (kmin, kmax) to avoid infinity and zero
		this.DK = Utils.minMaxScale(this.DK, kmin, kmax);

		// small amount of random noise to perturb topic timestamp noise.
		this.tau = 0;

		//System.out.println("this.tau set to " + this.tau );
		//this.uniform = new UniformRealDistribution(-1.*this.tau, this.tau);
		// Add random permutation before sampling
		/*
		for (int i = 0; i < this.DK.length; i++) {
			this.DK[i] = this.DK[i] + this.uniform.sample();
			if (this.DK[i] < kmin) {
				this.DK[i] = kmin;
			} else if (this.DK[i] > kmax) {
				this.DK[i] = kmax;
			}
				
		}*/
		//System.out.println(Arrays.toString(this.DK));
		//System.out.println(Utils.variance(this.DK));
		this._count_doc_size();
	}
	
	
	/**
	 * Constructor for ToTGibbsSampler Class with detailed parameters
	 * @param T The number of topics wanted.
	 * @param numIters The number of iterations to run.
	 * @param alpha Dirichlet prior for document over topic multinomial
	 * 
	 * @param beta Dirichlet prior for topic over word multinomial
	 * @param doc Filename of the corpus of interest. The first column of doc should be the total number of tokens. 
	 * Starting from the second row, the first column should be the current token's document index;
	 * the second being its word index in the vocabulary array.
	 * @param voc Filename of the vocabulary array. Each row contains a vocabulary term
	 * @param dkArray Filename of document time array;
	 */
	public ToTGibbsSampler(int T, int numIters,
						  double alpha, double beta,
						  String doc, String voc, String dkArray) {
		this(doc, voc, dkArray);

        this.setNumTopics(T);
        this.setNumIters(numIters);
        this.setAlpha(alpha);
        this.setBeta(beta);
	}
	
    @Override
    public String getModelType() {
		return "ToT";
    }
	
    @Override
    public TMOutcome get_outcome() {
        ToTOutcome outcome = new ToTOutcome(this.theta, this.phi, this.psi, this.topics, this.topicsRep);
        return outcome;
    }

	
	/* Run Gibbs Sampling on the given input */
	@Override
	protected void _sample() {
		double[] pi;
		int di, wi, zi;
		int zi_sample = -1;
		double pi_cum, normalized_density, u;
		// Sequentially sample for each token
		topicTimeArr = new ArrayList<GlueList<Double>>();
		for (int t = 0; t < this.T; t++) {
			topicTimeArr.add(new GlueList<Double>());
		}
		for (int i = 0; i < NNZ; i++) {
			di = this.DS[i]; // doc index
			wi = this.WS[i]; // voc index
			zi = this.Z[i]; // topic index
			// discount
			this.n_d_t[di][zi]--;
			this.n_t_w[zi][wi]--;
			this.n_t[zi]--;

			// pi is cumulative;
			pi = new double[this.T];
			pi_cum = 0;
			//System.out.println("*****************");
			for (int t = 0; t < this.T; t++) {
				normalized_density = FastMath.pow(Utils.betaPDF(this.psi[t][0], this.psi[t][1], this.DK[di]), 
												 1./this.ND[di]);
				pi_cum += (this.n_d_t[di][t] + this.alpha) * 
						  (this.n_t_w[t][wi] + this.beta) / (n_t[t] + this.W * this.beta) * 
						  normalized_density;
				//System.out.println(pi_cum + "; " + normalized_density + "; " + this.ND[di] + "; " + Utils.betaPDF(this.psi[t][0], this.psi[t][1], this.DK[di]));
				pi[t] = pi_cum;
			}
			//System.out.println("*****************");

			// Try to use cumulative method for sampling

			for (int t = 0; t < this.T; t++) {
				pi[t] = pi[t]/pi[this.T-1];
			}

			zi_sample = Utils.sample1d(pi);
			assert zi_sample >= 0;
			
			this.Z[i] = zi_sample;
			this.n_d_t[di][zi_sample]++;
			this.n_t_w[zi_sample][wi]++;
			this.n_t[zi_sample]++;

			// timestamp
			topicTimeArr.get(zi_sample).add(this.DK[di]);
		}

		
		// update psi: timestamp distribution given topics
		// method of moments
		//System.out.println(Arrays.deepToString(topicTimeArr.get(0).toArray()));
		this.psi = Utils.betaMoM(topicTimeArr, this.psi, this.tau);
		
	}
	
		

	@Override
	protected void _init() {
		this.n_d_t = new int[this.D][this.T];
		this.n_t_w = new int[this.T][this.W];
		this.n_t = new int[this.T];
		
		// random init
		final Random random = new Random();
		this.Z = random.ints(0, this.T).limit(this.NNZ).toArray();
		
		// init psi: T x 2 matrix (alpha and beta for Beta distribution)
		this.psi = new double[this.T][2];

		// array for recording topics' time distribution
		topicTimeArr = new ArrayList<GlueList<Double>>();
		for (int t = 0; t < this.T; t++) {
			topicTimeArr.add(new GlueList<Double>());
			psi[t][0] = psi[t][1] = 1;
		}
		
		int di, wi, zi;
		for (int i = 0; i < this.NNZ; i++) {
			di = this.DS[i]; // doc index
			wi = this.WS[i]; // voc index
			zi = this.Z[i]; // topic index

			this.n_d_t[di][zi]++;
			this.n_t_w[zi][wi]++;
			this.n_t[zi]++;
			topicTimeArr.get(zi).add(this.DK[di]);
			
		}
		//System.out.println(Arrays.toString(this.n_t) + " " + Utils.sum(this.n_t));
		//this.betaMoM(topicTimeArr);

	}

}
