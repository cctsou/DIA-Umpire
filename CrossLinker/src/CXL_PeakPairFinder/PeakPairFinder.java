/* 
 * Author: Chih-Chiang Tsou <chihchiang.tsou@gmail.com>
 *             Nesvizhskii Lab, Department of Computational Medicine and Bioinformatics, 
 *             University of Michigan, Ann Arbor
 *
 * Copyright 2014 University of Michigan, Ann Arbor, MI
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *      http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */
package CXL_PeakPairFinder;

import MSUmpire.BaseDataStructure.InstrumentParameter;
import MSUmpire.PeakDataStructure.PeakCluster;
import MSUmpire.PeptidePeakClusterDetection.PeakCurveCorrCalc;
import crosslinker.Linker;
import java.io.IOException;
import java.util.HashMap;
import net.sf.javaml.core.kdtree.KDTree;
import net.sf.javaml.core.kdtree.KeySizeException;
import org.apache.log4j.Logger;

/**
 *
 * @author Chih-Chiang Tsou <chihchiang.tsou@gmail.com>
 */
public class PeakPairFinder implements Runnable {

    private PeakCluster TargetPeak;
    InstrumentParameter parameter;
    Linker linker;
    private final KDTree PeakClusterKDTree;
    float lowrt = 0f;
    float highrt = 0f;
    public PairGroup pairgroup;

    public PeakPairFinder(PeakCluster peakClusterA, KDTree PeakClusterKDTree, InstrumentParameter parameter, Linker linker) {
        this.TargetPeak = peakClusterA;
        this.PeakClusterKDTree = PeakClusterKDTree;
        this.parameter = parameter;
        this.linker=linker;
        lowrt = peakClusterA.PeakHeightRT[0] - parameter.RTtol;
        if (lowrt < 0) {
            lowrt = 0f;
        }
        highrt = peakClusterA.PeakHeightRT[0] + parameter.RTtol;
    }

    @Override
    public void run() {
        //linked pair : x + Linker.Core
        pairgroup=new PairGroup(TargetPeak);
        
        //calculate MW of peak pair
        float pairmw= TargetPeak.NeutralMass()+linker.Core;
        pairgroup.highMassPeak=FindTargetMWPeak(pairmw);        
        float DeadEndpairMW = TargetPeak.NeutralMass()+ linker.Core + linker.H2O + linker.Arm;
        pairgroup.DeadEndpairs=FindTargetMWPeak(DeadEndpairMW);
    }

    private HashMap<Integer,CoElutePeak> FindTargetMWPeak(float targetmw) {

        HashMap<Integer,CoElutePeak> peaks=new HashMap<>();
        
        float lowmw = InstrumentParameter.GetMzByPPM(targetmw, 1, parameter.MS1PPM);
        float highmw = InstrumentParameter.GetMzByPPM(targetmw, 1, -parameter.MS1PPM);

        Object[] found = null;
        try {
            found = PeakClusterKDTree.range(new double[]{lowrt, lowmw}, new double[]{highrt, highmw});
        } catch (KeySizeException ex) {
        }
        if (found == null || found.length == 0) {
            return null;
        }
   
        for (Object foundpeak : found) {
            PeakCluster peakB = (PeakCluster) foundpeak;
            float ppm = InstrumentParameter.CalcPPM(targetmw, peakB.NeutralMass());
            if (ppm < parameter.MS1PPM) {
                float corr = 0f;
                try {
                    corr = PeakCurveCorrCalc.CalPeakCorr(TargetPeak.MonoIsotopePeak, peakB.MonoIsotopePeak, parameter.NoPeakPerMin);
                } catch (IOException ex) {
                    Logger.getRootLogger().error(ex.getMessage());
                }
                if (Float.isNaN(corr)) {
                    corr = 0f;
                }

                //if (corr > maxcorr) {
                if (corr>0.5f && (peaks.get(peakB.Charge)==null || (ppm < peaks.get(peakB.Charge).PPM))) {
                    if (!parameter.DetectSameChargePairOnly || peakB.Charge == TargetPeak.Charge) {
                        CoElutePeak peak = new CoElutePeak(peakB, corr, ppm);
                        peaks.put(peakB.Charge, peak);
                    }
                }
            }
        }        
        return peaks;
    }
}
