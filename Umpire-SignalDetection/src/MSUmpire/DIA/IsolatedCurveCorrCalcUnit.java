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
import MSUmpire.PeptidePeakClusterDetection.PeakCurveCorrCalc;
import Utility.UpdateProcess;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 *
 * @author Chih-Chiang Tsou <chihchiang.tsou@gmail.com>
 */
public class IsolatedCurveCorrCalcUnit implements Runnable {

    public PeakCurve targetMS1Curve;
    private SortedCurveCollectionApexRT PeakCurveSortedListRT;
    InstrumentParameter parameter;
    UpdateProcess update;
    HashMap<Integer, PeakCurve> FragmentPeaks = new HashMap<>();
    public ArrayList<PrecursorFragmentPairEdge> ResultList = new ArrayList<>();
    //TreeMap<Float, XYPoint>[] IsotopePatternMap;

    public IsolatedCurveCorrCalcUnit(PeakCurve MS1Peakcurve, SortedCurveCollectionApexRT PeakCurveSortedListRT, InstrumentParameter parameter, UpdateProcess update) {
        this.targetMS1Curve = MS1Peakcurve;
        this.PeakCurveSortedListRT = PeakCurveSortedListRT;
        this.parameter = parameter;
        this.update = update;
        //this.IsotopePatternMap = IsotopePatternMap;
    }

    @Override
    public void run() {

        int startRTidx = PeakCurveSortedListRT.BinarySearchLower(targetMS1Curve.ApexRT - parameter.ApexDelta);
        int endRTidx = PeakCurveSortedListRT.BinarySearchHigher(targetMS1Curve.ApexRT + parameter.ApexDelta);
        //XYPoint[] PatternRange = MS1PeakCluster.GetPatternRange(IsotopePatternMap);

//        if (PatternRange[0].Y > 1) {
//            targetMS1Curve = MS1PeakCluster.IsoPeaksCurves[1];
//        }
        for (int idx = startRTidx; idx <= endRTidx; idx++) {
            PeakCurve peakCurve = PeakCurveSortedListRT.get(idx);

            boolean overlap = false;
            if (targetMS1Curve.StartRT() >= peakCurve.StartRT() && targetMS1Curve.StartRT() <= peakCurve.EndRT()) {
                overlap = true;
            } else if (targetMS1Curve.EndRT() >= peakCurve.StartRT() && targetMS1Curve.EndRT() <= peakCurve.EndRT()) {
                overlap = true;
            } else if (peakCurve.StartRT() >= targetMS1Curve.StartRT() && peakCurve.StartRT() <= targetMS1Curve.EndRT()) {
                overlap = true;
            } else if (peakCurve.EndRT() >= targetMS1Curve.StartRT() && peakCurve.EndRT() <= targetMS1Curve.EndRT()) {
                overlap = true;
            }
            if (overlap) {
                float corr = 0f;

                try {
                    corr = PeakCurveCorrCalc.CalPeakCorr(targetMS1Curve, peakCurve, parameter.NoPeakPerMin);
                } catch (IOException ex) {
                    Logger.getLogger(IsolatedCurveCorrCalcUnit.class.getName()).log(Level.SEVERE, null, ex);
                }
                if (!Float.isNaN(corr) && corr > parameter.CorrThreshold) {
                    PrecursorFragmentPairEdge overlapRegion = new PrecursorFragmentPairEdge();
                    overlapRegion.Correlation = corr;
                    overlapRegion.PeakCurveIndexA = targetMS1Curve.Index;
                    overlapRegion.PeakCurveIndexB = peakCurve.Index;
                    overlapRegion.FragmentMz = peakCurve.TargetMz;
                    overlapRegion.Intensity = peakCurve.ApexInt;
                    FragmentPeaks.put(peakCurve.Index, peakCurve);
                    ResultList.add(overlapRegion);
                }
            }
        }
        if (update != null) {
            update.Update();
        }
    }
}
