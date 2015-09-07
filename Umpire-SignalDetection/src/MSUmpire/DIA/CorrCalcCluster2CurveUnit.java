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
import MSUmpire.PeakDataStructure.PeakCluster;
import MSUmpire.PeakDataStructure.PeakCurve;
import MSUmpire.PeakDataStructure.SortedCurveCollectionApexRT;
import MSUmpire.PeptidePeakClusterDetection.PeakCurveCorrCalc;
import Utility.UpdateProcess;
import java.io.IOException;
import java.util.ArrayList;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 *
 * @author Chih-Chiang Tsou <chihchiang.tsou@gmail.com>
 */
public class CorrCalcCluster2CurveUnit implements Runnable {

    public PeakCluster MS1PeakCluster;
    private SortedCurveCollectionApexRT PeakCurveSortedListApexRT;
    InstrumentParameter parameter;
    UpdateProcess update;
    //HashMap<Integer, PeakCurve> FragmentPeaks = new HashMap<>();
    public ArrayList<PrecursorFragmentPairEdge> GroupedFragmentList = new ArrayList<>();
    //TreeMap<Float, XYPoint>[] IsotopePatternMap;

    public CorrCalcCluster2CurveUnit(PeakCluster MS1PeakCluster, SortedCurveCollectionApexRT PeakCurveSortedListApexRT, InstrumentParameter parameter, UpdateProcess update) {
        this.MS1PeakCluster = MS1PeakCluster;
        this.PeakCurveSortedListApexRT = PeakCurveSortedListApexRT;
        this.parameter = parameter;
        this.update = update;
        //this.IsotopePatternMap = IsotopePatternMap;
    }

    @Override
    public void run() {

        int startRTidx = PeakCurveSortedListApexRT.BinarySearchLower(MS1PeakCluster.PeakHeightRT[0] - parameter.ApexDelta);
        int endRTidx = PeakCurveSortedListApexRT.BinarySearchHigher(MS1PeakCluster.PeakHeightRT[0] + parameter.ApexDelta);
        //XYPoint[] PatternRange = MS1PeakCluster.GetPatternRange(IsotopePatternMap);
        PeakCurve targetMS1Curve = MS1PeakCluster.MonoIsotopePeak;
//        if (PatternRange[0].Y > 1) {
//            targetMS1Curve = MS1PeakCluster.IsoPeaksCurves[1];
//        }
        float ms1rtrange = targetMS1Curve.EndRT() - targetMS1Curve.StartRT();
        int highCorrCnt = 0;
        
        for (int idx = startRTidx; idx <= endRTidx; idx++) {
            PeakCurve peakCurve = PeakCurveSortedListApexRT.get(idx);
            if(peakCurve.TargetMz>MS1PeakCluster.NeutralMass()){
                continue;
            }
            float peakcurvertrange = peakCurve.EndRT() - peakCurve.StartRT();
            float OverlapP = 0f;
            if (targetMS1Curve.StartRT() >= peakCurve.StartRT() && targetMS1Curve.StartRT() <= peakCurve.EndRT() && targetMS1Curve.EndRT() >= peakCurve.EndRT()) {
                OverlapP = (peakCurve.EndRT() - targetMS1Curve.StartRT()) / ms1rtrange;
            } else if (targetMS1Curve.EndRT() >= peakCurve.StartRT() && targetMS1Curve.EndRT() <= peakCurve.EndRT() && targetMS1Curve.StartRT() <= peakCurve.StartRT()) {
                OverlapP = (targetMS1Curve.EndRT() - peakCurve.StartRT()) / ms1rtrange;
            } else if (targetMS1Curve.StartRT() <= peakCurve.StartRT() && targetMS1Curve.EndRT() >= peakCurve.EndRT()) {
                OverlapP = peakcurvertrange / ms1rtrange;
            } else if (targetMS1Curve.StartRT() >= peakCurve.StartRT() && targetMS1Curve.EndRT() <= peakCurve.EndRT()) {
                OverlapP = 1;
            }
            if (OverlapP > parameter.RTOverlapThreshold && targetMS1Curve.ApexRT >= peakCurve.StartRT() && targetMS1Curve.ApexRT <= peakCurve.EndRT() && peakCurve.ApexRT >= targetMS1Curve.StartRT() && peakCurve.ApexRT <= targetMS1Curve.EndRT()) {
                float corr = 0f;
                float ApexDiff = Math.abs(targetMS1Curve.ApexRT - peakCurve.ApexRT);

                try {
                    corr = PeakCurveCorrCalc.CalPeakCorr(targetMS1Curve, peakCurve, parameter.NoPeakPerMin);
                } catch (IOException ex) {
                    Logger.getLogger(CorrCalcCluster2CurveUnit.class.getName()).log(Level.SEVERE, null, ex);                    
                }
                if (!Float.isNaN(corr) && corr > parameter.CorrThreshold) {
                    PrecursorFragmentPairEdge PrecursorFragmentPair = new PrecursorFragmentPairEdge();
                    PrecursorFragmentPair.Correlation = corr;
                    PrecursorFragmentPair.PeakCurveIndexA = MS1PeakCluster.Index;
                    //float intensity = peakCurve.GetMaxIntensityByRegionRange(targetMS1Curve.StartRT(), targetMS1Curve.EndRT());
                    PrecursorFragmentPair.PeakCurveIndexB = peakCurve.Index;
                    PrecursorFragmentPair.FragmentMz = peakCurve.TargetMz;
                    PrecursorFragmentPair.Intensity = peakCurve.ApexInt;
                    PrecursorFragmentPair.RTOverlapP = OverlapP;
                    PrecursorFragmentPair.ApexDelta = ApexDiff;
                    //FragmentPeaks.put(peakCurve.Index, peakCurve);
                    GroupedFragmentList.add(PrecursorFragmentPair);
                    if(PrecursorFragmentPair.Correlation>parameter.HighCorrThreshold){
                        highCorrCnt++;
                    }
                }
            }
        }       
        if(highCorrCnt<parameter.MinHighCorrCnt){
            GroupedFragmentList.clear();
        }
        if (update != null) {
            update.Update();
        }
    }
}
