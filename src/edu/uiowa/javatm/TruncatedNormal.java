package edu.uiowa.javatm;

import org.apache.commons.math3.distribution.NormalDistribution;

public class TruncatedNormal {
	
	private double mu, sigma, lo, hi, alpha, beta, Z;
	private final NormalDistribution stdNormal = new NormalDistribution(1, 0);
	
	public TruncatedNormal(double mu, double sigma, double lo, double hi) {
		this.hi = hi;
		this.lo = lo;
		this.mu = mu;
		this.sigma = sigma;
		
		this.alpha = (this.lo - this.mu)/(this.sigma);
		this.beta = (this.hi - this.mu)/(this.sigma);
		
		this.Z = stdNormal.cumulativeProbability(beta) - stdNormal.cumulativeProbability(alpha);
	}
	
	public double density(double x) {
		double den = this.sigma * this.Z;
		double num = stdNormal.density((x-this.mu)/this.sigma);
		return (num /den);
	}

}
