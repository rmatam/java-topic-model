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

import edu.uiowa.javatm.TMGibbsSampler;
import edu.uiowa.javatm.LDAGibbsSampler;
import edu.uiowa.javatm.ToTGibbsSampler;
import edu.uiowa.javatm.TMOutcome;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;

/**
 * Main class for reading and fitting data
 * @author Zhiya Zuo
 */
public class JavaTM {
	/**
	 * @param args First one indicates which topic model to use
	 */
	public static void main(String[] args) {
		TMGibbsSampler tmGibbsSampler = null;

		Option modelType = Option.builder("model").longOpt("model-type")
				.desc("Type of topic models to use")
				.hasArg().required().build();

		Option dataName = Option.builder("name").longOpt("data-name")
				.desc("Data name: used for saving outputs")
				.hasArg().required().build();

		Option alpha = Option.builder("a").longOpt("alpha")
				.desc("Dirichlet prior for document (author) over topic multinomial")
				.hasArg().required().build();

		Option beta = Option.builder("b").longOpt("beta")
				.desc("Dirichlet prior for topic over word multinomial")
				.hasArg().required().build();

		Option pi = Option.builder("p").longOpt("pi")
				    .desc("Dirichlet prior for topic over time multinomial")
					.hasArg().build();

		Option K = Option.builder("K").longOpt("K")
				.desc("The number of timestamp indices")
				.hasArg().build();
		
		/*Option tau = Option.builder("tau").longOpt("tau")
				.desc("Smoothing constant for topic time")
				.hasArg().build();*/

		Option doc = Option.builder("doc").longOpt("document-file")
				.desc("WD matrix to use")
				.hasArg().required().build();

		Option voc = Option.builder("voc").longOpt("vocabulary-file")
				.desc("Vocabulary file of the corpus of interest")
				.hasArg().required().build();

		Option auth = Option.builder("auth").longOpt("auth-file")
				.desc("Author indices for each token")
				.hasArg().build();

		Option authArray = Option.builder("authArray").longOpt("author-list-file")
				.desc("Author list")
				.hasArg().build();

		Option dkArray = Option.builder("dk").longOpt("document-time-file")
				.desc("Document timestamp file")
				.hasArg().build();
		
		Option citationMat = Option.builder("cm").longOpt("citation-matrix")
				.desc("Citation overtime for the corpus")
				.hasArg().build();

		Option numTopics = Option.builder("topic").longOpt("num-topics")
				.desc("The total number of topics")
				.hasArg().required().build();

		Option numIters = Option.builder("iter").longOpt("num-iters")
				.desc("The total number of iterations")
				.hasArg().required().build();

		Option outputDir = Option.builder("odir").longOpt("output-dir")
				.desc("Output directory")
				.hasArg().required().build();

		
		Options options = new Options();
		options.addOption(modelType).addOption(alpha).addOption(beta)
			   .addOption(numTopics).addOption(K).addOption(pi).addOption(citationMat)
			   .addOption(numIters).addOption(doc).addOption(voc).addOption(dkArray)
			   .addOption(outputDir).addOption(auth).addOption(authArray).addOption(dataName);

		CommandLineParser parser = new DefaultParser();
	    try {
	        // parse the command line arguments
	        CommandLine line = parser.parse( options, args );
	        String model = line.getOptionValue("model");
	        String name = line.getOptionValue("name");
	        String docFile = line.getOptionValue("doc");
	        String vocFile = line.getOptionValue("voc");
	        int topics = Integer.parseInt(line.getOptionValue("topic"));
	        int iters = Integer.parseInt(line.getOptionValue("iter"));
	        double a = Double.parseDouble(line.getOptionValue("a"));
	        double b = Double.parseDouble(line.getOptionValue("b"));

	        String modelLower = model.toLowerCase();
	        if (modelLower.equals("lda")) {
				tmGibbsSampler = new LDAGibbsSampler(topics, iters, a, b, docFile, vocFile);
			} else if (modelLower.equals("at")) {
				String authFile = line.getOptionValue("auth");
				String authArrayFile = line.getOptionValue("authArray");
				//double tau_val = Double.parseDouble(line.getOptionValue("tau"));
				tmGibbsSampler = new ATGibbsSampler(topics, iters, a, b, docFile, vocFile, 
												   authFile, authArrayFile);
			} else if (modelLower.equals("tot")) {
				String dkFile = line.getOptionValue("dk");
				//double tau_val = Double.parseDouble(line.getOptionValue("tau"));
				tmGibbsSampler = new ToTGibbsSampler(topics, iters, a, b, docFile, vocFile, dkFile);
			} else if (modelLower.equals("tiot")) {
				String timeFile = line.getOptionValue("dk");
				String citationFile = line.getOptionValue("cm");
				double p = Double.parseDouble(line.getOptionValue("p"));
				//int k = Integer.parseInt(line.getOptionValue("K"));
				tmGibbsSampler = new TIOTGibbsSampler(topics, iters, a, b, p, 
													 docFile, vocFile, timeFile, citationFile);			
			} else { 
				System.err.println("Invalid model type selection. Must be lda, at, tot or atot.");
				System.exit(ExitStatus.ILLEGAL_ARGUMENT);
			}

			long startTime = System.nanoTime();
			tmGibbsSampler.fit();
			TMOutcome outcome =  tmGibbsSampler.get_outcome();
			long endTime = System.nanoTime();
			long duration = (endTime - startTime); 
			System.out.println("Overall elapsed time: " + duration/1000000000.+" seconds");

	        tmGibbsSampler.showTopics(10);
			outcome.showTopicDistribution();

			String oDir = line.getOptionValue("odir");
			if (!oDir.endsWith("/")) {
				oDir = oDir + "/";
			}
			// append name to `oDir`
			oDir = oDir + name + "-";

			if (modelLower.contains("tot")) {
				// topic over time (tot and atot) has beta distribution parameters to write
				Utils.write2DArray(((ToTOutcome)outcome).getPsi(), oDir+"psi-"+modelLower+".csv");
			}

			if (modelLower.contains("tiot")) {
				// topic over time (tot and atot) has beta distribution parameters to write
				Utils.write2DArray(((TIOTOutcome) outcome).getPsi(), oDir+"psi-"+modelLower+".csv");
				double[][][] ga = ((TIOTOutcome) outcome).getGa();
				for (int t = 0; t < ga.length; t++) {
					Utils.write2DArray(ga[t], oDir+"gamma-"+t+"-"+modelLower+".csv");
				}
			}

			Utils.write2DArray(outcome.getPhi(), oDir+"phi-"+modelLower+".csv");
			Utils.write2DArray(outcome.getTheta(), oDir+"theta-"+modelLower+".csv");
			
			System.out.println("Output files saved to " + oDir);
	    }
	    catch( ParseException exp ) {
	        // oops, something went wrong
	        System.err.println( "Parsing failed. Reason: " + exp.getMessage() );
	    }
		
	}

}
