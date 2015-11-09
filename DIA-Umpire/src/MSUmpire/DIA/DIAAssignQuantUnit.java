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
package MSUmpire.DIA;

import MSUmpire.BaseDataStructure.InstrumentParameter;
import MSUmpire.BaseDataStructure.SpectralDataType;
import MSUmpire.LCMSBaseStructure.LCMSPeakMS1;
import MSUmpire.PSMDataStructure.FragmentPeak;
import MSUmpire.PSMDataStructure.PepIonID;
import MSUmpire.PeakDataStructure.PeakCluster;
import MSUmpire.PeakDataStructure.PrecursorFragmentPairEdge;
import com.compomics.util.experiment.biology.Ion;
import com.compomics.util.experiment.biology.ions.ElementaryIon;
import com.compomics.util.experiment.biology.ions.PeptideFragmentIon;
import org.apache.commons.io.FilenameUtils;
import org.apache.log4j.Logger;

/**
 *
 * @author Chih-Chiang Tsou <chihchiang.tsou@gmail.com>
 */
public class DIAAssignQuantUnit implements Runnable {

    PepIonID pepIonID;
    LCMSPeakMS1 ms1lcms;
    double protonMass = ElementaryIon.proton.getTheoreticMass();
    InstrumentParameter parameter;

    public DIAAssignQuantUnit(PepIonID pepIonID, LCMSPeakMS1 ms1lcms, InstrumentParameter parameter) {
        this.pepIonID = pepIonID;
        this.ms1lcms = ms1lcms;
        this.parameter = parameter;
    }

    @Override
    public void run() {        
        PeakCluster targetCluster = null;
        for (PeakCluster peakCluster : pepIonID.MS1PeakClusters) {
            if (targetCluster == null || targetCluster.PeakHeight[0] < peakCluster.PeakHeight[0]) {
                targetCluster = peakCluster;
            }
        }
        pepIonID.CreateQuantInstance(ms1lcms.MaxNoPeakCluster);

        if (targetCluster != null) {
            pepIonID.PeakArea = targetCluster.PeakArea;
            pepIonID.PeakHeight = targetCluster.PeakHeight;
            pepIonID.PeakClusterScore = targetCluster.MS1Score;
            pepIonID.PeakRT = targetCluster.PeakHeightRT[0];
            pepIonID.ObservedMz = targetCluster.mz[0];
        } else {
            for (PeakCluster peakCluster : pepIonID.MS2UnfragPeakClusters) {
                if (targetCluster == null || targetCluster.PeakHeight[0] < peakCluster.PeakHeight[0]) {
                    targetCluster = peakCluster;
                }
            }
            if (targetCluster != null) {
                pepIonID.PeakRT = targetCluster.PeakHeightRT[0];
                pepIonID.ObservedMz = targetCluster.mz[0];
            }
        }
        //pepIonID.PeakClusterIndex = targetCluster.Index;

        if (targetCluster!=null && ms1lcms.datattype != SpectralDataType.DataType.pSMART) {
            MatchFragmentByTargetCluster(targetCluster);
            pepIonID.RemoveRedundantFrag();
            if (pepIonID.FragmentPeaks.isEmpty() && Math.max(pepIonID.MaxProbability, pepIonID.TargetedProbability()) > 0.8f) {
                Logger.getRootLogger().warn("Warning: " + pepIonID.ModSequence + "(MaxProb: " + pepIonID.MaxProbability + ") does not have matched fragment in "+FilenameUtils.getBaseName(ms1lcms.ParentmzXMLName));
                //MatchFragment();
            }
            pepIonID.ClearPepFragFactory();
        }
    }

    private void MatchFragmentByTargetCluster(PeakCluster peakCluster) {
        for (Ion frag : pepIonID.GetFragments()) {
            PrecursorFragmentPairEdge bestfragment = null;
            float targetmz = (float) frag.getTheoreticMz(1);
            for (PrecursorFragmentPairEdge fragmentClusterUnit : peakCluster.GroupedFragmentPeaks) {
                if (InstrumentParameter.CalcPPM(targetmz, fragmentClusterUnit.FragmentMz) <= parameter.MS2PPM) {
                    if (bestfragment == null || fragmentClusterUnit.Correlation > bestfragment.Correlation) {
                        bestfragment = fragmentClusterUnit;
                    }
                }
            }
           
            if (bestfragment != null) {
                FragmentPeak fragmentpeak = new FragmentPeak();
                fragmentpeak.ObservedMZ = bestfragment.FragmentMz;
                fragmentpeak.FragMZ = targetmz;
                fragmentpeak.corr = bestfragment.Correlation;
                fragmentpeak.intensity = bestfragment.Intensity;
                fragmentpeak.ApexDelta = bestfragment.ApexDelta;
                fragmentpeak.RTOverlapP = bestfragment.RTOverlapP;
                fragmentpeak.Charge = 1;
                fragmentpeak.ppm = InstrumentParameter.CalcSignedPPM(bestfragment.FragmentMz, targetmz);
                fragmentpeak.IonType = frag.getSubTypeAsString() + ((PeptideFragmentIon) frag).getNumber();
                pepIonID.FragmentPeaks.add(fragmentpeak);
            }

            targetmz = (float) frag.getTheoreticMz(2);
            bestfragment = null;
            for (PrecursorFragmentPairEdge fragmentClusterUnit : peakCluster.GroupedFragmentPeaks) {
                if (InstrumentParameter.CalcPPM(targetmz, fragmentClusterUnit.FragmentMz) <= parameter.MS2PPM) {
                    if (bestfragment == null || fragmentClusterUnit.Correlation > bestfragment.Correlation) {
                        bestfragment = fragmentClusterUnit;
                    }
                }
            }

            if (bestfragment != null) {
                FragmentPeak fragmentpeak = new FragmentPeak();
                fragmentpeak.ObservedMZ = bestfragment.FragmentMz;
                fragmentpeak.FragMZ = targetmz;
                fragmentpeak.corr = bestfragment.Correlation;
                fragmentpeak.intensity = bestfragment.Intensity;
                fragmentpeak.ApexDelta = bestfragment.ApexDelta;
                fragmentpeak.RTOverlapP = bestfragment.RTOverlapP;
                fragmentpeak.Charge = 2;
                fragmentpeak.ppm = InstrumentParameter.CalcSignedPPM(bestfragment.FragmentMz, targetmz);
                fragmentpeak.IonType = frag.getSubTypeAsString() + ((PeptideFragmentIon) frag).getNumber();
                pepIonID.FragmentPeaks.add(fragmentpeak);
            }
        }
    }
    
    private void MatchFragment() {
        for (Ion frag : pepIonID.GetFragments()) {
            PrecursorFragmentPairEdge bestfragment = null;
            float targetmz = (float) frag.getTheoreticMz(1);
            for (PeakCluster peakCluster : pepIonID.MS1PeakClusters) {
                for (PrecursorFragmentPairEdge fragmentClusterUnit : peakCluster.GroupedFragmentPeaks) {
                    if (InstrumentParameter.CalcPPM(targetmz, fragmentClusterUnit.FragmentMz) <= parameter.MS2PPM) {
                        if (bestfragment == null || fragmentClusterUnit.Correlation > bestfragment.Correlation) {
                            bestfragment = fragmentClusterUnit;
                        }
                    }
                }
            }
            for (PeakCluster peakCluster : pepIonID.MS2UnfragPeakClusters) {
                for (PrecursorFragmentPairEdge fragmentClusterUnit : peakCluster.GroupedFragmentPeaks) {
                    if (InstrumentParameter.CalcPPM(targetmz, fragmentClusterUnit.FragmentMz) <= parameter.MS2PPM) {
                        if (bestfragment == null || fragmentClusterUnit.Correlation > bestfragment.Correlation) {
                            bestfragment = fragmentClusterUnit;
                        }
                    }
                }
            }
            if (bestfragment != null) {
                FragmentPeak fragmentpeak = new FragmentPeak();
                fragmentpeak.ObservedMZ = bestfragment.FragmentMz;
                fragmentpeak.FragMZ = targetmz;
                fragmentpeak.corr = bestfragment.Correlation;
                fragmentpeak.intensity = bestfragment.Intensity;
                fragmentpeak.ApexDelta = bestfragment.ApexDelta;
                fragmentpeak.RTOverlapP = bestfragment.RTOverlapP;
                fragmentpeak.Charge = 1;
                fragmentpeak.ppm = InstrumentParameter.CalcSignedPPM(bestfragment.FragmentMz, targetmz);
                fragmentpeak.IonType = frag.getSubTypeAsString() + ((PeptideFragmentIon) frag).getNumber();
                pepIonID.FragmentPeaks.add(fragmentpeak);
            }

            targetmz = (float) frag.getTheoreticMz(2);
            bestfragment = null;
            for (PeakCluster peakCluster : pepIonID.MS1PeakClusters) {
                for (PrecursorFragmentPairEdge fragmentClusterUnit : peakCluster.GroupedFragmentPeaks) {
                    if (InstrumentParameter.CalcPPM(targetmz, fragmentClusterUnit.FragmentMz) <= parameter.MS2PPM) {
                        if (bestfragment == null || fragmentClusterUnit.Correlation > bestfragment.Correlation) {
                            bestfragment = fragmentClusterUnit;
                        }
                    }
                }
            }
            for (PeakCluster peakCluster : pepIonID.MS2UnfragPeakClusters) {
                for (PrecursorFragmentPairEdge fragmentClusterUnit : peakCluster.GroupedFragmentPeaks) {
                    if (InstrumentParameter.CalcPPM(targetmz, fragmentClusterUnit.FragmentMz) <= parameter.MS2PPM) {
                        if (bestfragment == null || fragmentClusterUnit.Correlation > bestfragment.Correlation) {
                            bestfragment = fragmentClusterUnit;
                        }
                    }
                }
            }
            if (bestfragment != null) {
                FragmentPeak fragmentpeak = new FragmentPeak();
                fragmentpeak.ObservedMZ = bestfragment.FragmentMz;
                fragmentpeak.FragMZ = targetmz;
                fragmentpeak.corr = bestfragment.Correlation;
                fragmentpeak.intensity = bestfragment.Intensity;
                fragmentpeak.ApexDelta = bestfragment.ApexDelta;
                fragmentpeak.RTOverlapP = bestfragment.RTOverlapP;
                fragmentpeak.Charge = 2;
                fragmentpeak.ppm = InstrumentParameter.CalcSignedPPM(bestfragment.FragmentMz, targetmz);
                fragmentpeak.IonType = frag.getSubTypeAsString() + ((PeptideFragmentIon) frag).getNumber();
                pepIonID.FragmentPeaks.add(fragmentpeak);
            }
        }
    }
}
