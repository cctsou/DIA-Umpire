/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package CXL_PeakPairFinder;

import MSUmpire.PeakDataStructure.PeakCluster;
import java.util.HashMap;

/**
 *
 * @author Chih-Chiang Tsou
 */
public class PairGroup {
    public PeakCluster lowMassPeak;
    public HashMap<Integer, CoElutePeak> highMassPeak;
    public HashMap<Integer,CoElutePeak> DeadEndpairs;
    private CoElutePeak BestPair=null; 
    
    public PairGroup(PeakCluster peakClusterA){
        this.lowMassPeak=peakClusterA;
    }
    public CoElutePeak GetBestPeakPair() {
        if (BestPair == null && highMassPeak!=null) {
            for (CoElutePeak coElutePeak : highMassPeak.values()) {
                if (BestPair==null || (coElutePeak.Correlation > BestPair.Correlation || (coElutePeak.Correlation == BestPair.Correlation && coElutePeak.PPM < BestPair.PPM))) {
                    BestPair = coElutePeak;
                }
            }
        }
        return BestPair;
    }
}
