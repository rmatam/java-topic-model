package edu.uiowa.javatm;

import java.util.Arrays;
import java.util.ArrayList;
import java.util.Collections;
import java.util.function.DoubleToLongFunction;

import org.apache.commons.csv.CSVFormat;
import org.apache.commons.csv.CSVPrinter;
import org.apache.commons.math3.distribution.NormalDistribution;
import org.apache.commons.math3.stat.StatUtils;
import org.apache.commons.math3.util.FastMath;

import com.sun.media.sound.AlawCodec;

import java.io.*;

class Utils {
    /* * * * * * * * * * * * * * * * * *
     * This class contains static methods
     * to help preprocessing and facilitate
     * Gibbs Sampling process
     * * * * * * * * * * * * * * * * * */

    public static int[][] readDocFile(String docFileName) {
    // {{{ 
    // Convert doc2mat resulting file to a matrix NNZ x 2
    // 1st Column: Document Index DS
    // 2nd Column: Word Index (as in the vocabulary) WS
        int[][] WD = null;
        try {
            FileReader docReader = new FileReader(new File(docFileName));
            BufferedReader reader = new BufferedReader(docReader);
            String line = null;
            String[] wd = null;
            int index = -1;

            while ((line = reader.readLine()) != null) {

                if (index == -1) {
                    WD = new int[Integer.parseInt(line)][2];
                    index++;
                    continue;
                }

                wd = line.split(",");
                WD[index][0] = Integer.parseInt(wd[0]);
                WD[index][1] = Integer.parseInt(wd[1]);
                index++;
            }

            reader.close();
            docReader.close();
        } catch (Exception ex) {
            System.out.print("Error when reading WD file... ");
            ex.printStackTrace();
        } 
    //}}}
        return WD;
    }

    public static String[] readVocabulary(String vocFileName) {
    //{{{
        // Read vocabulary file
        // A word/token per line
        ArrayList<String> vocArrayList;
        String[] vocArray = null;
        try {
            FileReader vocabReader = new FileReader(new File(vocFileName));
            BufferedReader reader = new BufferedReader(vocabReader);
            String line = null;
            vocArrayList = new ArrayList<String>();

            while ((line = reader.readLine()) != null) {
                vocArrayList.add(line); 
            }

            reader.close();
            vocabReader.close();

            vocArray = vocArrayList.toArray(new String[vocArrayList.size()]);

        } catch (Exception ex) {
            System.out.print("Error when reading the vocabulary/author file... ");
            ex.printStackTrace();
        }
    //}}}
        return vocArray;
    }

    public static String[] readAuthor(String authorFileName) {
        return readVocabulary(authorFileName);
    }

    public static int[] readTimestamp(String timeFileName) {
        String[] dkStrArray = readVocabulary(timeFileName);
        	int[] dkArray = new int[dkStrArray.length];
        	for (int i = 0; i < dkArray.length; i++) {
        		dkArray[i] = Integer.valueOf(dkStrArray[i]);
		}
        	return dkArray;
    }

    public static int[][] readDA(String daMatrixFileName){
    // {{{
    	int[][] DA = null;
        try {
            FileReader daReader = new FileReader(new File(daMatrixFileName));
            BufferedReader reader = new BufferedReader(daReader);
            String line = null;
            String[] da = null;
            int row = -1;

            while ((line = reader.readLine()) != null) {

                if (row == -1) {
                    da = line.split(",");
                    DA = new int[Integer.parseInt(da[0])][Integer.parseInt(da[1])];
                    row++;
                    continue;
                }

                da = line.split(",");
                //System.out.println(Arrays.toString(da) + ": " + row);
                for (int col = 0; col < da.length; col++) {
                    DA[row][col] = Integer.parseInt(da[col]);
                }
                row++;
            }

            reader.close();
            daReader.close();
        } catch (IOException e) {
            System.err.print("Error redaing" + daMatrixFileName + "---- ");
            e.printStackTrace();
        }
    //}}}
        return DA;
    }
    
    public static ArrayList<int[]> readAS(String wordAuthorFileName) {
        ArrayList<int[]> AS = new ArrayList<int[]>();
        try {
            FileReader docReader = new FileReader(new File(wordAuthorFileName));
            BufferedReader reader = new BufferedReader(docReader);
            String line = null;
            String[] as = null;
            while ((line = reader.readLine()) != null) {
                as = line.split(",");
                AS.add(Arrays.stream(as).mapToInt(Integer::parseInt).toArray());
            }
            reader.close();
            docReader.close();
        } catch (Exception ex) {
            System.out.print("Error when reading AS file... ");
            ex.printStackTrace();
        } 
        return AS;
    }

    public static double[][] readCitationFile(String citationMatrixFile) {
    // {{{
        int[][] citationIntMatrix = readDA(citationMatrixFile);

        // Min-max for each timestamp index
        double[][] citationMatrix = new double[citationIntMatrix.length][citationIntMatrix[0].length];

		double[] tmp = new double[citationMatrix[0].length];
		for (int k = 0; k < citationMatrix.length; k++) {
			for (int d = 0; d < citationMatrix[k].length; d++) {
				tmp[d] = citationIntMatrix[k][d] * 1.0;
			}
			citationMatrix[k] = minMaxScale(tmp, 1e-2, 1-1e-2);
		}
    // }}}
        return citationMatrix;
    }
    
    public static Integer[] nonZeroIndex(int[] intArray) {
    	ArrayList<Integer> nnzIndices = new ArrayList<Integer>();
    	
    	for (int i = 0; i < intArray.length; i++) {
            if (intArray[i] > 0) {
                nnzIndices.add(i);
            }
        }
    	return nnzIndices.toArray(new Integer[nnzIndices.size()]);
    }

    // Get a certain column from a 2d int array
    public static int[] getColumn(int[][] int2DArray, int column) {
        int nrow = int2DArray.length;
        int[] columnArray = new int[nrow];
        for (int row = 0; row < nrow; row++) {
            columnArray[row] = int2DArray[row][column];
        }
        return columnArray;
    }

    public static double mean(double[] doubleArray) {
    // mean of a double array
        return sum(doubleArray) / doubleArray.length; 
    }

    public static double mean(ArrayList<Double> doubleArrayList) {
    // mean of a double array
        return mean(doubleArrayList.stream().mapToDouble(Double::doubleValue).toArray()); 
    }

    public static double mean(GlueList<Double> doubleGlueList) {
    		double avg = 0;
    		for (int i = 0; i < doubleGlueList.size; i++) {
    			avg += doubleGlueList.get(i);
		}
    		return avg/doubleGlueList.size;
    }
    
    public static double variance(GlueList<Double> doubleGluelist) {
        double sumSquare = 0;
        double avg = mean(doubleGluelist);
        for (int i = 0; i < doubleGluelist.size; i++) {
            sumSquare += (doubleGluelist.get(i) - avg) * (doubleGluelist.get(i) - avg);
        }
    //}}}
        return sumSquare / (doubleGluelist.size - 1);
    }

    public static double variance(ArrayList<Double> doubleArrayList) {
    // mean of a double array
        return variance(doubleArrayList.stream().mapToDouble(Double::doubleValue).toArray()); 
    }
    
    public static double variance(double[] doubleArray) {
    //{{{ sample variance of a double array
        double sumSquare = 0;
        double avg = mean(doubleArray);
        for (int i = 0; i < doubleArray.length; i++) {
            sumSquare += (doubleArray[i] - avg) * (doubleArray[i] - avg);
        }
    //}}}
        return sumSquare / (doubleArray.length - 1);
    }
    
    public static ArrayList<Integer> find(int[] intArray, int val) {
    //{{{ Find the indices of val in an integer array
        ArrayList<Integer> positions = new ArrayList<Integer> ();
        for (int i = 0; i < intArray.length; i++) {
            if(intArray[i] == val) {
                positions.add(i);
            }
        }
    //}}}
        return positions;
    }

    public static int[] sortIndex(double[] doubleArray) {
    // Get the index of a sorted double array
        return sortIndex(_arrayToArrayList(doubleArray));
    }

    public static int[] sortIndex(ArrayList<Double> doubleArrayList) {
    // Get the index of a sorted Double ArrayList
        ArrayList<Double> sortedArrayList = new ArrayList<Double>(doubleArrayList);
        Collections.sort(sortedArrayList);
        ArrayList<Double> copyArrayList = new ArrayList<Double>(doubleArrayList);
        int length = sortedArrayList.size();
        int[] sortedIndices = new int[length];
        Arrays.fill(sortedIndices, -1);

        for (int i = 0; i < length; i++) {
            // Traverse the value in descending value
            Double ithValue = sortedArrayList.get(length - 1 - i);
            while (_contains(sortedIndices, copyArrayList.indexOf(ithValue))) {
                copyArrayList.set(copyArrayList.indexOf(ithValue), -1.);
            }
            sortedIndices[i] = copyArrayList.indexOf(ithValue);
        }
        return sortedIndices;
    }

    static double[] normalize(double[] doubleArray) {
    // Normalize a double array
        double arraySum = sum(doubleArray);
        double[] arrayNormalized = new double[doubleArray.length];
        for (int i = 0; i < doubleArray.length; i++) {
            arrayNormalized[i] = doubleArray[i] / arraySum;
        }
        return arrayNormalized;
    }

    static double[] minMaxScale(double[] doubleArray) {
    		return minMaxScale(doubleArray, 0., 1.);
    }

    static double[] minMaxScale(double[] doubleArray, double low, double high) {
    // Normalize a double array to [0, 1] by min-max scaling
        double max_ = StatUtils.max(doubleArray);
        double min_ = StatUtils.min(doubleArray);
        double range_ = max_ - min_;
        double range_hl = high - low;
        double[] arrayNormalized = new double[doubleArray.length];
        for (int i = 0; i < doubleArray.length; i++) {
            arrayNormalized[i] = range_hl * (doubleArray[i]-min_) / range_ + low;
        }
        return arrayNormalized;
    }

    public static double sum(int[] intArray) {
		double[] doubleArray = Arrays.stream(intArray).asDoubleStream().toArray();
		return sum(doubleArray);
    }
    public static double sum(double[] doubleArray) {
    // Sum the elements in a double array
        double arraySum = 0;
        for (double d : doubleArray) {
			arraySum += d;
		}
        return arraySum;
    }

    // Sum the elements in a double matrix
    /*
    private static double _sum2d(double[][] doubleMatrix) {
        double matrixSum = 0;
        for (int i = 0; i < doubleMatrix.length; i++) {
            matrixSum += sum(doubleMatrix[i]);
        }
        return matrixSum;
    }
    private static double _sum2d(int[][] intMatrix) {
        return _sum2d(_intMatrix2doubleMatrix(intMatrix));
    }*/

    // Convert array into ArrayList
    private static ArrayList<Double> _arrayToArrayList(double[] doubleArray) {
        ArrayList<Double> doubleArrayList = new ArrayList<Double>();
        for (int i = 0; i < doubleArray.length; i++) {
            doubleArrayList.add(doubleArray[i]);
        }
        return doubleArrayList;
    }

    // Check if an int array contains an element
    private static boolean _contains(int[] intArray, int val) {
        for (int intVal : intArray) {
            if (intVal == val) {
                return true;
            }
        }
        return false;
    }

    /*
    private static double[][] _intMatrix2doubleMatrix(int[][] intMatrix) {

        double[][] doubleMatrix = new double[intMatrix.length][intMatrix[0].length];
        for (int row = 0; row < doubleMatrix.length; row++) {
            for (int column = 0; column < doubleMatrix[0].length; column++) {
                doubleMatrix[row][column] = (double) intMatrix[row][column];
            }
        }
        return doubleMatrix;
    }*/
    
    public static void write2DArray(String[][] arr, 
								   String csvFilePath) {
		FileWriter fWriter;
		String joinedString;
		try {
			fWriter = new FileWriter(csvFilePath);
			CSVPrinter csvPrinter = new CSVPrinter(fWriter, 
												  CSVFormat.DEFAULT.withQuote(null));
			for (String[] str : arr) {
				joinedString = String.join(",", str); 
				csvPrinter.printRecord(joinedString);
				//csvPrinter.printRecords(str);
			}
			csvPrinter.flush();
			csvPrinter.close();
			fWriter.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
    }
    
    public static void write2DArray(double[][] arr, 
								   String csvFilePath) {
    		String[][] stringArr = new String[arr.length][arr[0].length];
    		for (int i = 0; i < arr.length; i++) {
    			for(int j = 0; j < arr[0].length; j++) {
    				stringArr[i][j] = Double.toString(arr[i][j]);
    			}
		}
    		
    		write2DArray(stringArr, csvFilePath);
    }
    
    /*gamma functions: https://introcs.cs.princeton.edu/java/91float/Gamma.java.html*/
    private static double logGamma(double x) {
        double tmp = (x - 0.5) * Math.log(x + 4.5) - (x + 4.5);
        double ser = 1.0 + 76.18009173    / (x + 0)   - 86.50532033    / (x + 1)
                         + 24.01409822    / (x + 2)   -  1.231739516   / (x + 3)
                         +  0.00120858003 / (x + 4)   -  0.00000536382 / (x + 5);
        return tmp + Math.log(ser * Math.sqrt(2 * Math.PI));
     }
     public static double gamma(double x) { return Math.exp(logGamma(x)); }
     
	
	public static double betaPDF(double a, double b, double x) {
		// density of beta distribution with (a, b) given x
		double val = _betaFunc(a, b) / _betaNum(a, b, x);
		return _checkBounded(val);
	}
	
	private static double _betaNum(double a, double b, double x) {
		return Math.pow(x, a-1) * Math.pow(1-x, b-1);
	}
	
	private static double _betaFunc(double a, double b) {
		double val = (gamma(a)*gamma(b)/gamma(a+b));
		return _checkBounded(val);
	}
	
	private static double _checkBounded(double x) {
		// Check if a value `x` is bounded within finite or not
		if (Double.isInfinite(x)) {
			return 1e6;
		} else if (Double.isNaN(x)) {
			return 0;
		} else {
			return x;
		}
	}

	// method of moment to update beta distributions: topic over time
	public static double[][] betaMoM(ArrayList<GlueList<Double>> sampleArr, 
									double[][] psi, double tau) {
		double kz_mean, kz_var, common_val;
		GlueList<Double> arr;
		double[][] psiUpdated = new double[psi.length][psi[0].length];

		for (int t = 0; t < sampleArr.size(); t++) {
			arr = sampleArr.get(t);
			if(arr.size < 1) {
				// if nothing, do not update
				psiUpdated[t][0] = psi[t][0];
				psiUpdated[t][1] = psi[t][0];
				//psiUpdated[t][0] = 1.;
				//psiUpdated[t][1] = 1.;
				continue;
			} else {
				kz_mean = mean(arr)+tau;
				kz_var = _checkBounded(variance(arr))+tau;
			}

			common_val = (1-kz_mean)*kz_mean/kz_var - 1;

			//assert kz_mean * common_val > 0;
			//assert (1-kz_mean) * common_val > 0;
			
			if (common_val < 0) {
				psiUpdated[t][0] = psi[t][0];
				psiUpdated[t][1] = psi[t][0];
				continue;
				//System.err.println(Arrays.toString(arr.toArray()));
				//System.err.println(kz_mean + " " + kz_var + " " + common_val);
				//System.exit(1);
			} else {
				//System.out.println(Arrays.toString(arr.toArray()));
				//System.out.println(kz_mean + " " + kz_var + " " + common_val);
			}

			psiUpdated[t][0] = _checkBounded(kz_mean * common_val);
			psiUpdated[t][1] = _checkBounded((1-kz_mean) * common_val);
		}

		return psiUpdated;
	}

	public static double gaussianPDF(double mu, double sigma, double x) {
		// density of beta distribution with (a, b) given x
		double val_exp = -0.5 * ((x-mu)/sigma) * ((x-mu)/sigma);
		double val = 1/sigma * Math.exp(val_exp);
		return _checkBounded(val);
	}
	
	public static double truncatedGaussianPDF(double mu, double sigma,  
											 double lo, double hi, double x, 
											 NormalDistribution stdNormal) {
		double alpha = (lo - mu)/(sigma);
		double beta = (hi - mu)/(sigma);
		double Z = stdNormal.cumulativeProbability(beta) - stdNormal.cumulativeProbability(alpha);
		double den = sigma * Z;
		double num = stdNormal.density((x-mu)/sigma);
		return (num /den);
	}

	// method of moment to update gaussian distributions: topic over time

	public static double[][] truncatedGaussianMoM(ArrayList<GlueList<Double>> sampleArr, 
												 double lo, double hi,
												 double[][] psi, double tau) {
		/*
		GlueList<Double> arr;
		double kz_mean, kz_var, Z, ap, bt, phi_ap, phi_bt;
		double[][] psiUpdated = new double[psi.length][psi[0].length];

		for (int t = 0; t < sampleArr.size(); t++) {
			arr = sampleArr.get(t);
			if(arr.size < 1) {
				// if nothing, do not update
				psiUpdated[t][0] = psi[t][0];
				psiUpdated[t][1] = psi[t][0];
				//psiUpdated[t][0] = 1.;
				//psiUpdated[t][1] = 1.;
				continue;
			} else {
				kz_mean = mean(arr)+tau;
				kz_var = _checkBounded(variance(arr))+tau;
			}
			
			ap = 
			
			//psiUpdated[t][0] = kz_mean - 
		}

		return psiUpdated;
		*/
		return null;
	}
	
	
	// sample from unnormalized cumulative distribution
	public static int sample1d(double[] pArr) {
		for (int i = 0; i < pArr.length; i++) {
			pArr[i] = pArr[i]/pArr[pArr.length-1];
		}
		double u = FastMath.random();
		//System.out.println("+++++");
		//System.out.println(u);
		//System.out.println(Arrays.toString(pArr));
		for (int zi_tmp = 0; zi_tmp < pArr.length; zi_tmp++) {
			if (u <= pArr[zi_tmp]) {
				//System.out.println(u + "; " + zi_tmp);
				return zi_tmp;
			}
		}
		return -1;
	
	}

	public static int[] sample2d(double[][] pMat, int start, int end) {
		assert pMat[0].length == (end-start);
		//System.out.println("++++++++++++++++++++++++++");
		//System.out.println(Arrays.deepToString(pMat));
		for (int i = 0; i < pMat.length; i++) {
			for (int j = 0; j < pMat[0].length; j++) {
				pMat[i][j] = pMat[i][j]/pMat[pMat.length-1][pMat[0].length-1];
			}
		}
		double u = FastMath.random();
		for (int zi_tmp = 0; zi_tmp < pMat.length; zi_tmp++) {
			for (int j = start; j < end; j++) {
				if (u <= pMat[zi_tmp][j-start]) {
					//System.out.println(u + "; " + (j-start) + "; " + zi_tmp);
					//System.out.println(Arrays.deepToString(pMat));
					return new int[] {zi_tmp, j};
				}
			}
		}


		return null;
	
	}

}
