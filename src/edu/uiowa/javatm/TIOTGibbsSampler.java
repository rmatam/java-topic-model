package edu.uiowa.javatm;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Random;

import org.apache.commons.math3.distribution.NormalDistribution;
import org.apache.commons.math3.util.FastMath;

public class TIOTGibbsSampler extends ToTMultiGibbsSampler {
	
	protected int[] KZ; // latent variable for time
	protected double[][][] ga; // topic time beta params
	protected double[][] citationMat; // citation matrix: row as doc and col as time; entry: normalized citation

	protected ArrayList<ArrayList<GlueList<Double>>> citationSampleArr;  // citation matrix samples
	
	private String distributionType; // distribution to use when sampling for (normalized) citations
	private final NormalDistribution stdNormal = new NormalDistribution(0, 1);
			
	/*
	 * Note that, for each word (or their belonging doc),
	 * their timestamp can range from t_di to the largest timestamp in file
	 * Operationalized as `impact timestamp`.
	 * */
	
	/* Besides Beta, let me try Gaussian distribution as well (strictly unimodal)
	 * */
	
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
	 * @param citation Filename of the document citation matrix. Each row contains a doc's citation over each timestamp.
	 */

	public TIOTGibbsSampler(int T, int numIters,
						   double alpha, double beta, double pi,
						   String doc, String voc, String time, String citation) {
		super(T, 0, numIters, alpha, beta, pi, 
			  doc, voc, time);
		this.citationMat = Utils.readCitationFile(citation);
		//assert this.K == this.citationMat.length;
		this.setNumTimestamp(this.citationMat.length);

		this.ga = new double[this.T][this.K][2];
		//System.out.println(Arrays.deepToString(this.CitationMat));
		for (int t = 0; t < this.T; t++) {
			for (int k = 0; k < this.K; k++) {
				Arrays.fill(this.ga[t][k], 1);
			}
		}
		this.distributionType = "BETA";
		this._count_doc_size();
	}

	public TIOTGibbsSampler(int T, int numIters,
						   double alpha, double beta, double pi,
						   String doc, String voc, String time, String citation, 
						   String distribution) {
		this(T, numIters, alpha, beta, pi, doc, voc, time, citation);
		if (distribution.equals("BETA") | distribution.equals("GAUSSIAN")) {
			this.distributionType = distribution;
		} else {
			System.err.println(distribution + " is not a valid distribution type for TIOT.\n Choose from BETA or GAUSSIAN.");
			System.exit(1);
		}
	}
	
    @Override
    public String getModelType() {
		return "TIOT";
    }

    @Override
    public String getDistributionType() {
		return this.distributionType;
    }
	
    @Override
    public TMOutcome get_outcome() {
        TIOTOutcome outcome = new TIOTOutcome(this.theta, this.phi, this.psi, this.ga,
											 this.topics, this.topicsRep);
        return outcome;
    }

	
	/* Run Gibbs Sampling on the given input */
	@Override
	protected void _sample() {
		// instead of using a 1-d array.
		// I will 2-d for both topic and time dimensions.
		double[][] pi_matrix; 
		int[] zi_ki;
		double pi_cum = 0, den, num, u, unnormallized_density;
		int di, wi, zi, ki, kzi, kindex;
		int zi_sample, ki_sample;
		
		citationSampleArr = new ArrayList<ArrayList<GlueList<Double>>>();
		for (int t = 0; t < this.T; t++) {
			citationSampleArr.add(new ArrayList<GlueList<Double>>());
			for (int k = 0; k < this.K; k++) {
				citationSampleArr.get(t).add(new GlueList<Double>());
			}
		}

		for (int i = 0; i < NNZ; i++) {
			di = this.DS[i]; // doc index
			wi = this.WS[i]; // voc index
			zi = this.Z[i]; // topic index
			ki = this.dkArray[di]; // time index
			kzi = this.KZ[i]; // impact time index
			// discount
			this.n_d_t[di][zi]--;
			this.n_t_w[zi][wi]--;
			this.n_t_k[zi][kzi]--;
			this.n_t[zi]--;

			// pi is cumulative;
			pi_matrix = new double[this.T][this.K-ki];
			pi_cum = 0;
			for (int t = 0; t < this.T; t++) {
				// all timestep before ki will be assigned a zero value
				for (int k = ki; k < this.K; k++) {
					// denominator
					den = (n_t[t] + this.W * this.beta) * (n_t[t] + this.K * this.pi);
					// numerator
					num = (this.n_d_t[di][t] + this.alpha) * (this.n_t_w[t][wi] + this.beta) * 
						  (this.n_t_k[t][k] + this.pi);
					// times beta distribution on time
					if (this.distributionType.contains("BETA")) {
						unnormallized_density = Utils.betaPDF(this.ga[t][k][0], this.ga[t][k][1], this.citationMat[kzi][di]);
					} else {
						unnormallized_density = Utils.truncatedGaussianPDF(this.ga[t][k][0], this.ga[t][k][1], 
																		  0., 1.,
																		  this.citationMat[kzi][di],
																		  this.stdNormal);
					}

					num = num * FastMath.pow(unnormallized_density, 1./this.ND[di]);
					//num = num * FastMath.pow(1./(this.K - ki), unnormallized_density);
					//System.out.println("$$ " + FastMath.pow(1./this.ND[di], unnormallized_density));
					
					pi_cum += num / den; 

					kindex = k - ki;
					pi_matrix[t][kindex] = pi_cum;
				}
			}

			// Try to use cumulative method for sampling
			//u = Math.random() * pi_matrix[pi_matrix.length-1][pi_matrix[0].length-1];
			//u = Math.random() * pi[pi.length-1];
			
			zi_ki = Utils.sample2d(pi_matrix, ki, this.K);
			assert zi_ki != null;
			zi_sample = zi_ki[0];
			ki_sample = zi_ki[1];

			citationSampleArr.get(zi_sample).get(ki_sample).add(this.citationMat[ki_sample][di]);
			//citationSampleArr.get(zi_sample).get(ki).add(this.citationMat[ki][di]);
			
			this.Z[i] = zi_sample;
			this.KZ[i] = ki_sample;
			this.n_d_t[di][zi_sample]++;
			this.n_t_w[zi_sample][wi]++;
			this.n_t_k[zi_sample][ki_sample]++;
			this.n_t[zi_sample]++;
		}
			
		if(this.distributionType.equals("BETA")) {
			for (int t = 0; t < this.T; t++) {
				this.ga[t] = Utils.betaMoM(citationSampleArr.get(t), this.ga[t], 1e-4);
			}
		} 
		/*else {
			for (int t = 0; t < this.T; t++) {
				this.ga[t] = Utils.gaussianMoM(citationSampleArr.get(t), this.ga[t], 1e-4);
			}
		}*/
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
		this.KZ = new int[this.NNZ];
		
		int di, wi, zi, ki;
		for (int i = 0; i < this.NNZ; i++) {
			di = this.DS[i]; // doc index
			wi = this.WS[i]; // voc index
			zi = this.Z[i]; // topic index
			// min should be ki, and max be this.K
			ki = this.dkArray[di];
			ki = random.nextInt(this.K-ki)+ki; 
			assert ki < this.K;
			assert ki >= this.dkArray[di];
			this.KZ[i] = ki;

			this.n_d_t[di][zi]++;
			this.n_t_w[zi][wi]++;
			this.n_t_k[zi][ki]++;
			this.n_t[zi]++;
		}

	}
}
