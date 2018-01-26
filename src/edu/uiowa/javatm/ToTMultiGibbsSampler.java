package edu.uiowa.javatm;

import java.util.Random;

public class ToTMultiGibbsSampler extends LDAGibbsSampler {

    protected double pi; // topic over time dirichlet prior
	protected double[][] psi; // topic over time multinomial
	protected int[] dkArray; // document timestamp index array (integer)

	int[][] n_t_k; // n_t_k: topic t with time k;
	
	/**
	 * Constructor for ToTMultiGibbsSampler Class with detailed parameters
	 * @param T The number of topics wanted.
	 * @param K The number of total timestamps wanted.
	 * @param numIters The number of iterations to run.
	 * @param alpha Dirichlet prior for document over topic multinomial
	 * 
	 * @param beta Dirichlet prior for topic over word multinomial
	 * @param pi Dirichlet prior for topic over time multinomial
	 * @param doc Filename of the corpus of interest. The first column of doc should be the total number of tokens. 
	 * Starting from the second row, the first column should be the current token's document index;
	 * the second being its word index in the vocabulary array.
	 * @param voc Filename of the vocabulary array. Each row contains a vocabulary term
	 * @param time Filename of the document timestamp index array. Each row contains a timestamp index for each document.
	 */

	public ToTMultiGibbsSampler(int T, int K, int numIters,
							   double alpha, double beta, double pi,
							   String doc, String voc, String time) {
		super(doc, voc);
		this.dkArray = Utils.readTimestamp(time);

        this.setNumTopics(T);
        this.setNumTimestamp(K);
        this.setNumIters(numIters);
        this.setAlpha(alpha);
        this.setBeta(beta);
        this.setPi(pi);
	}
	
	public void setNumTimestamp(int K) {
		this.K = K;
	}

	public void setPi(double pi) {
		this.pi = pi;
	}
	
    @Override
    public String getModelType() {
		return "ToTMulti";
    }
	
    @Override
    public TMOutcome get_outcome() {
        ToTMultiOutcome outcome = new ToTMultiOutcome(this.theta, this.phi, this.psi, 
													 this.topics, this.topicsRep);
        return outcome;
    }

    @Override
    void _est_multinomial() {
        this.theta = new double[this.D][this.T];
        this.phi = new double[this.T][this.W];
        this.psi = new double[this.T][this.K];

        double k_sum, w_sum;
        for (int t = 0; t < this.T; t++) {
        		k_sum = Utils.sum(this.n_t_k[t]);
       		w_sum = Utils.sum(this.n_t_w[t]);
            for (int k = 0; k < this.K; k++) {
                this.psi[t][k] = (this.pi + this.n_t_k[t][k]) /
								(this.K*this.pi + k_sum);
            }
            for (int w = 0; w < this.W; w++) {
                this.phi[t][w] = (this.beta + this.n_t_w[t][w]) /
								(this.W*this.beta + w_sum);
                assert this.phi[t][w] >= 0;
            }
            for (int d = 0; d < this.D; d++) {
                this.theta[d][t] = (this.alpha + this.n_d_t[d][t]) /
								  (this.T*this.alpha + Utils.sum(n_d_t[d]));
                assert this.theta[d][t] >= 0;
            }

			_getTopicProportion();
			_getTopicRepresentation();
        }
	}

	
	/* Run Gibbs Sampling on the given input */
	@Override
	protected void _sample() {
		// instead of using a 1-d array.
		// I will 2-d for both topic and time dimensions.
		double[] pi; 
		int di, wi, zi, ki;
		int zi_sample;

		for (int i = 0; i < NNZ; i++) {
			di = this.DS[i]; // doc index
			wi = this.WS[i]; // voc index
			zi = this.Z[i]; // topic index
			ki = this.dkArray[di]; // time index
			// discount
			this.n_d_t[di][zi]--;
			this.n_t_w[zi][wi]--;
			this.n_t_k[zi][ki]--;
			this.n_t[zi]--;

			// pi is cumulative;
			pi = new double[this.T];
			double pi_cum = 0, den, num;
			for (int t = 0; t < this.T; t++) {
				den = (n_t[t] + this.W * this.beta) * (n_t[t] + this.K * this.pi);
				num = (this.n_d_t[di][t] + this.alpha) * (this.n_t_w[t][wi] + this.beta) * 
					  (this.n_t_k[t][ki] + this.pi);
				pi_cum += num / den; 
				pi[t] = pi_cum;
			}
			//System.out.println(Arrays.deepToString(pi_matrix));

			// Try to use cumulative method for sampling
			double u = Math.random() * pi[this.T-1];
			zi_sample = -1;
			for (int t = 0; t < this.T; t++) {
				if (u <= pi[t]) {
					zi_sample = t;
					break;
				}
			}

			this.Z[i] = zi_sample;
			this.n_d_t[di][zi_sample]++;
			this.n_t_w[zi_sample][wi]++;
			this.n_t_k[zi_sample][ki]++;
			this.n_t[zi_sample]++;
		}
	}


	@Override
	protected void _init() {
		this.n_d_t = new int[this.D][this.T];
		this.n_t_w = new int[this.T][this.W];
		this.n_t_k = new int[this.T][this.K];
		this.n_t = new int[this.T];
		
		// random init
		final Random random = new Random();
		this.Z = random.ints(0, this.T).limit(this.NNZ).toArray();
		
		int di, wi, zi, ki;
		for (int i = 0; i < this.NNZ; i++) {
			di = this.DS[i]; // doc index
			wi = this.WS[i]; // voc index
			zi = this.Z[i]; // topic index
			ki = this.dkArray[di]; // time index

			this.n_d_t[di][zi]++;
			this.n_t_w[zi][wi]++;
			this.n_t_k[zi][ki]++;
			this.n_t[zi]++;
		}

	}
	
}
