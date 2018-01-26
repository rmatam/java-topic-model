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

import java.util.HashMap;
import java.util.Iterator;
import java.util.Map;

class TMOutcome {

    /* * * * * * * * * * * * * * * * * *
     * This class is only for returning 
     * multiple values from LDAGibbsSampler.
     *
     * Only Getter methods are provided.
     * * * * * * * * * * * * * * * * * * */

    // Multinomial parameters
    protected double[][] theta, phi;
    // Topic proportions
    protected HashMap<Integer, Double> topics;
    // Keywords in topics in descending order
    protected String[][] topicsRep;

    public TMOutcome() {}

    public TMOutcome(double[][] t, double[][] p, HashMap<Integer, Double> tp, String[][] tRep) {
        this.theta = t;
        this.phi = p;
        this.topics = tp;
        this.topicsRep = tRep;
    }

    /* Getter methods to retrieve results */
    public double[][] getTheta() {
        return this.theta;
    }

    public double[][] getPhi() {
        return this.phi;
    }

    public HashMap<Integer, Double> getTopics() {
        return this.topics;
    }

    public void showTopicDistribution() {
        Map<Integer, Double>  map = this.topics;
        Iterator<Map.Entry<Integer, Double>> entries = map.entrySet().iterator();

        while(entries.hasNext()) {
            Map.Entry<Integer, Double> entry = entries.next();
            System.out.println("Topic #" + entry.getKey() + ": " + entry.getValue());
        }
    }

    public String[][] getTopicsRep() {
        return topicsRep;
    }

}

