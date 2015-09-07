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

import MSUmpire.PeakDataStructure.PrecursorFragmentPairEdge;
import MSUmpire.BaseDataStructure.InstrumentParameter;
import MSUmpire.PeakDataStructure.PeakCurve;
import MSUmpire.PeakDataStructure.SortedCurveCollectionApexRT;
import Utility.UpdateProcess;
import com.compomics.util.experiment.biology.ions.ElementaryIon;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;

/**
 *
 * @author Chih-Chiang Tsou <chihchiang.tsou@gmail.com>
 */
public class FragDirectedGroupingUnit implements Runnable {

    public PeakCurve SeedFragment;
    private SortedCurveCollectionApexRT PeakCurveSortedListRT;
    InstrumentParameter parameter;
    UpdateProcess update;
    public ArrayList<PrecursorFragmentPairEdge> ResultList = new ArrayList<PrecursorFragmentPairEdge>();
    float RTtol = 0.5f;

    public FragDirectedGroupingUnit(PeakCurve SeedFragment, SortedCurveCollectionApexRT PeakCurveSortedListRT, InstrumentParameter parameter, UpdateProcess update) {
        this.SeedFragment = SeedFragment;
        this.PeakCurveSortedListRT = PeakCurveSortedListRT;
        this.parameter = parameter;
        this.update = update;
        //this.IsotopePatternMap = IsotopePatternMap;
    }

    @Override
    public void run() {

        int StartRTindex = PeakCurveSortedListRT.BinarySearchHigher(SeedFragment.ApexRT - RTtol);
        int EndRTindex = PeakCurveSortedListRT.BinarySearchLower(SeedFragment.ApexRT + RTtol);
        for (int i = StartRTindex; i <= EndRTindex; i++) {
//                        if (!ExcludeIndex.contains(SeedFragment.Index)) {
//                            PeakCurve groupedfragment = LCMSPeakBase.PeakCurveListRT.get(i);
//                            if (Math.abs(groupedfragment.ApexRT - fragment.ApexRT) < RTtol) {
//
//                            }
//                        }
//                    }
//        float ms1rtrange = targetMS1Curve.EndRT() - targetMS1Curve.StartRT();
//
//        for (int idx = startRTidx; idx <= endRTidx; idx++) {
//            PeakCurve peakCurve = PeakCurveSortedListRT.get(idx);
//            float peakcurvertrange = peakCurve.EndRT() - peakCurve.StartRT();            
//            float OverlapP = 0f;
//            if (targetMS1Curve.StartRT() >= peakCurve.StartRT() && targetMS1Curve.StartRT() <= peakCurve.EndRT() && targetMS1Curve.EndRT()>=peakCurve.EndRT()) {                
//                OverlapP = (peakCurve.EndRT() - targetMS1Curve.StartRT()) / ms1rtrange;
//            } else if (targetMS1Curve.EndRT() >= peakCurve.StartRT() && targetMS1Curve.EndRT() <= peakCurve.EndRT() && targetMS1Curve.StartRT() <= peakCurve.StartRT()) {                
//                OverlapP = (targetMS1Curve.EndRT() - peakCurve.StartRT()) / ms1rtrange;
//            } else if (targetMS1Curve.StartRT()<=peakCurve.StartRT() && targetMS1Curve.EndRT()>=peakCurve.EndRT()) {                
//                OverlapP = peakcurvertrange / ms1rtrange;
//            } else if (targetMS1Curve.StartRT()>=peakCurve.StartRT() && targetMS1Curve.EndRT()<=peakCurve.EndRT()) {                
//                OverlapP = 1;
//            }
//            if (OverlapP>0.3f && targetMS1Curve.ApexRT >= peakCurve.StartRT() && targetMS1Curve.ApexRT <= peakCurve.EndRT() && peakCurve.ApexRT >= targetMS1Curve.StartRT() && peakCurve.ApexRT <= targetMS1Curve.EndRT()) {
//                float corr = 0f;
//                float ApexDiff = Math.abs(targetMS1Curve.ApexRT - peakCurve.ApexRT);
//                
//                try {
//                    corr = PeakCurveCorrCalc.CalPeakCorr(targetMS1Curve, peakCurve, parameter.NoPeakPerMin);
//                } catch (IOException ex) {
//                    Logger.getLogger(FragDirectedGroupingUnit.class.getName()).log(Level.SEVERE, null, ex);
//                }
//                if (!Float.isNaN(corr) && corr > 0.2f) {
//                    PrecursorFragmentPairEdge PrecursorFragmentPair = new PrecursorFragmentPairEdge();
//                    PrecursorFragmentPair.Correlation = corr;
//                    PrecursorFragmentPair.PeakCurveIndexA = MS1PeakCluster.Index;
//                    float intensity = peakCurve.GetMaxIntensityByRegionRange(targetMS1Curve.StartRT(), targetMS1Curve.EndRT());
//                    PrecursorFragmentPair.PeakCurveIndexB = peakCurve.Index;
//                    PrecursorFragmentPair.FragmentMz = peakCurve.TargetMz;
//                    PrecursorFragmentPair.Intensity = peakCurve.ApexInt;
//                    PrecursorFragmentPair.RTOverlapP = OverlapP;
//                    PrecursorFragmentPair.ApexDelta = ApexDiff;
//                    //FragmentPeaks.put(peakCurve.Index, peakCurve);
//                    ResultList.add(PrecursorFragmentPair);
//                }
//            }
        }
        if (ResultList.size() > 1) {
            Collections.sort(ResultList, new Comparator<PrecursorFragmentPairEdge>() {
                @Override
                public int compare(PrecursorFragmentPairEdge o1, PrecursorFragmentPairEdge o2) {
                    return Float.compare(o1.FragmentMz, o2.FragmentMz);
                }
            });
            DeisotopingForPeakClusterFragment(ResultList);
            Collections.sort(ResultList, new Comparator<PrecursorFragmentPairEdge>() {
                @Override
                public int compare(PrecursorFragmentPairEdge o1, PrecursorFragmentPairEdge o2) {
                    return Float.compare(o1.FragmentMz, o2.FragmentMz);
                }
            });
        }
        if (update != null) {
            update.Update();
        }
    }

    public void DeisotopingForPeakClusterFragment(ArrayList<PrecursorFragmentPairEdge> fragments) {
        ArrayList<PrecursorFragmentPairEdge> newfragments = new ArrayList<>();
        boolean[] fragmentmarked = new boolean[fragments.size()];
        Arrays.fill(fragmentmarked, Boolean.TRUE);
        PrecursorFragmentPairEdge currentmaxfragment = fragments.get(0);
        int currentmaxindex = 0;
        for (int i = 1; i < fragments.size(); i++) {
            if (InstrumentParameter.CalcPPM(fragments.get(i).FragmentMz, currentmaxfragment.FragmentMz) > parameter.MS2PPM) {
                fragmentmarked[currentmaxindex] = false;
                currentmaxindex = i;
                currentmaxfragment = fragments.get(i);
            } else if (fragments.get(i).Intensity > currentmaxfragment.Intensity) {
                currentmaxindex = i;
                currentmaxfragment = fragments.get(i);
            }
        }
        fragmentmarked[currentmaxindex] = false;
        for (int i = 0; i < fragments.size(); i++) {
            if (!fragmentmarked[i]) {
                fragmentmarked[i] = true;
                PrecursorFragmentPairEdge startfrag = fragments.get(i);

                boolean groupped = false;
                for (int charge = 2; charge >= 1; charge--) {
                    float lastint = startfrag.Intensity;
                    boolean found = false;
                    for (int pkidx = 1; pkidx < 5; pkidx++) {
                        float targetmz = startfrag.FragmentMz + (float) pkidx / charge;
                        for (int j = i + 1; j < fragments.size(); j++) {
                            if (!fragmentmarked[j]) {
                                PrecursorFragmentPairEdge targetfrag = fragments.get(j);
                                if (InstrumentParameter.CalcPPM(targetfrag.FragmentMz, targetmz) < parameter.MS2PPM * (pkidx * 0.5 + 1)) {
                                    if (targetfrag.Intensity < lastint) {
                                        fragmentmarked[j] = true;
                                        lastint = targetfrag.Intensity;
                                        found = true;
                                        break;
                                    }
                                } else if (targetfrag.FragmentMz > targetmz) {
                                    break;
                                }
                            }
                        }
                        if (!found) {
                            break;
                        }
                    }
                    if (found) {
                        groupped = true;
                        //convert to charge 1 m/z
                        startfrag.FragmentMz = startfrag.FragmentMz * charge - (charge - 1) * (float) ElementaryIon.proton.getTheoreticMass();
                        newfragments.add(startfrag);
                    }
                }
                if (!groupped) {
                    newfragments.add(startfrag);
                }
            }
        }
        fragments.clear();
        for (PrecursorFragmentPairEdge fragment : newfragments) {
            fragments.add(fragment);
        }
    }
}
