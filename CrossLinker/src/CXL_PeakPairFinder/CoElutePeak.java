/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package CXL_PeakPairFinder;

import MSUmpire.PeakDataStructure.PeakCluster;

/**
 *
 * @author Chih-Chiang Tsou
 */
public class CoElutePeak{
        public PeakCluster peakpair;
        public float Correlation;
        public float PPM;
        
        public CoElutePeak(PeakCluster peakpair, float correlation, float ppm){
            this.peakpair=peakpair;
            this.Correlation=correlation;
            this.PPM=ppm;
        }
    }
