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
import MSUmpire.SpectralProcessingModule.IsotopePeakGroup;
import MSUmpire.SpectralProcessingModule.ScanPeakGroup;
import crosslinker.Linker;
import java.util.ArrayList;
import net.sf.javaml.core.kdtree.KDTree;
import net.sf.javaml.core.kdtree.KeyDuplicateException;
import net.sf.javaml.core.kdtree.KeySizeException;
import org.apache.log4j.Logger;

/**
 *
 * @author Chih-Chiang Tsou <chihchiang.tsou@gmail.com>
 */
public class MS2PeakPairFinder implements Runnable {

    public ScanPeakGroup ScanPeak;
    InstrumentParameter parameter;
    Linker linker;
    private KDTree PeakClusterKDTree;
    public ArrayList<PairGroupMS2> PeakPairGroupList;
    public ArrayList<PeakPair> peakPairs = new ArrayList<>();

    public class PeakPair {
        public PairGroupMS2 LowMassPeak;
        public PairGroupMS2 HighMassPeak;
    }

    private IsotopePeakGroup FindTargetMWPeak(IsotopePeakGroup peakCluster, float targetmw) {
        IsotopePeakGroup BestPeak = null;
        float BestPPM = Float.MAX_VALUE;
        float BestIntensity=0f;
        float lowmw = InstrumentParameter.GetMzByPPM(targetmw, 1, parameter.MS2PPM);
        float highmw = InstrumentParameter.GetMzByPPM(targetmw, 1, -parameter.MS2PPM);

        Object[] found = null;
        try {
            found = PeakClusterKDTree.range(new double[]{1, lowmw}, new double[]{1, highmw});
        } catch (KeySizeException ex) {
        }
        if (found == null || found.length == 0) {
            return null;
        }

        for (Object foundpeak : found) {
            IsotopePeakGroup peakB = (IsotopePeakGroup) foundpeak;
            float ppm = InstrumentParameter.CalcPPM(targetmw, peakB.NeutralMass());
            //if (ppm < parameter.MS2PPM && ppm < BestPPM) {
            if (ppm < parameter.MS2PPM && peakB.GetPeakXYPointByPeakidx(0).getY() > BestIntensity) {
                if (!parameter.DetectSameChargePairOnly || peakB.Charge == peakCluster.Charge) {
                    BestPeak = peakB;
                    BestIntensity=peakB.GetPeakXYPointByPeakidx(0).getY();
                    //BestPPM = ppm;
                }
            }
        }
        return BestPeak;
    }

    public MS2PeakPairFinder(ScanPeakGroup ScanPeak, InstrumentParameter parameter, Linker linker) {
        this.ScanPeak = ScanPeak;
        this.parameter = parameter;
        this.linker=linker;
    }

    
    
    @Override
    public void run() {
        PeakPairGroupList = new ArrayList<>();
        PeakClusterKDTree = new KDTree(2);
        
        //For each detected isotope peak group
        for (int i = 0; i < ScanPeak.peakGroupList.size(); i++) {
            IsotopePeakGroup peakGroup = ScanPeak.peakGroupList.get(i);
            try {
                PeakClusterKDTree.insert(new double[]{1, peakGroup.NeutralMass()}, peakGroup);
            } catch (KeyDuplicateException | KeySizeException ex) {
                Logger.getRootLogger().error(ex.getMessage());
            }
        }
        //Get TopN intensity threshold 
         float threshold = ScanPeak.Scan.GetTopNIntensity(parameter.MS2PairTopN);

        for (int i = 0; i < ScanPeak.peakGroupList.size(); i++) {
            //For each isotope peak group, initialize a PairGroupMS2 data structure
            IsotopePeakGroup PeakCluster = ScanPeak.peakGroupList.get(i);
            PairGroupMS2 Pairgroup = new PairGroupMS2(PeakCluster);
            //Calculate MW of peak pair
            float pairmw = PeakCluster.NeutralMass() + linker.Core;
            Pairgroup.HighMassPeak = FindTargetMWPeak(PeakCluster,pairmw);
            float DeadEndpairMW = PeakCluster.NeutralMass() + linker.Core + linker.H2O + linker.Arm;
            Pairgroup.DeadEndpairs = FindTargetMWPeak(PeakCluster,DeadEndpairMW);
            PeakPairGroupList.add(Pairgroup);
        }

        for (int i = 0; i < PeakPairGroupList.size(); i++) {
            if (PeakPairGroupList.get(i).HighMassPeak != null) {
                PairGroupMS2 lowmass = PeakPairGroupList.get(i);
                for (int j = 0; j < PeakPairGroupList.size(); j++) {
                    if (PeakPairGroupList.get(j).HighMassPeak != null) {
                        PairGroupMS2 highmass = PeakPairGroupList.get(j);
                        if (lowmass.LowMassPeak.NeutralMass() <= highmass.LowMassPeak.NeutralMass()) {
                            
                            float IntactMW = lowmass.LowMassPeak.NeutralMass() + highmass.LowMassPeak.NeutralMass() + linker.Core;
                            
                            //Check if the precursor mass value of the MS2 spectrum is the summation of MW values of low and high peak pairs
                            //If yes, then the two peak pairs are potentially crosslinking peptides
                            if (InstrumentParameter.CalcPPM(IntactMW, ScanPeak.Scan.PrecursorMass()) < parameter.MS2PPM) {
                                PeakPair pair = new PeakPair();
                                pair.HighMassPeak = highmass;
                                pair.LowMassPeak = lowmass;           
                                boolean LowPeakIntensity = false;
                                boolean HighPeakIntensity = false;
                           
                                //Check if any of peak of low mass peak pair higher than TopN intensity threshold
                                if (lowmass.LowMassPeak.GetPeakXYPointByPeakidx(0).getY() > threshold || lowmass.HighMassPeak.GetPeakXYPointByPeakidx(0).getY() > threshold) {
                                    LowPeakIntensity = true;
                                }
                                //Check if any of peak of high mass peak pair higher than TopN intensity threshold
                                if (highmass.LowMassPeak.GetPeakXYPointByPeakidx(0).getY() > threshold || highmass.HighMassPeak.GetPeakXYPointByPeakidx(0).getY() > threshold) {
                                    HighPeakIntensity = true;
                                }
                                
                                //criterion to determine whether the peak pair is valid
                                if (LowPeakIntensity && HighPeakIntensity) {
                                        peakPairs.add(pair);
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}
