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

import MSUmpire.PeptidePeakClusterDetection.PeakCurveCorrCalc;
import MSUmpire.BaseDataStructure.InstrumentParameter;
import MSUmpire.PeakDataStructure.PeakCluster;
import MSUmpire.PeakDataStructure.PrecursorFragmentPairEdge;
import MSUmpire.PeakDataStructure.PeakCurve;
import MSUmpire.PeakDataStructure.SortedClusterCollectionClassApexRT;
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
public class CorrCalcCluster2ClusterUnit implements Runnable {

    public PeakCluster MS1PeakCluster;
    SortedClusterCollectionClassApexRT RTSortedClusters;
    InstrumentParameter parameter;
    UpdateProcess update;
    HashMap<Integer, PeakCurve> FragmentPeaks = new HashMap<>();
    public ArrayList<PrecursorFragmentPairEdge> ResultList = new ArrayList<>();

    public CorrCalcCluster2ClusterUnit(PeakCluster MS1PeakCluster, SortedClusterCollectionClassApexRT RTSortedClusters, InstrumentParameter parameter, UpdateProcess update) {
        this.MS1PeakCluster = MS1PeakCluster;
        this.RTSortedClusters = RTSortedClusters;
        this.parameter = parameter;
        this.update = update;
    }

    @Override
    public void run() {

        int startRTidx = RTSortedClusters.BinarySearchLower(MS1PeakCluster.PeakHeightRT[0] - parameter.ApexDelta);
        int endRTidx = RTSortedClusters.BinarySearchHigher(MS1PeakCluster.PeakHeightRT[0]+ parameter.ApexDelta);

        for (int idx = startRTidx; idx <= endRTidx; idx++) {
            PeakCluster peakCluster = RTSortedClusters.get(idx);
//            if(peakCluster.Corrs[0]<0.3f )
//                continue;            

            boolean overlap = false;
            if (MS1PeakCluster.startRT >= peakCluster.startRT && MS1PeakCluster.startRT <= peakCluster.endRT) {
                overlap = true;
            } else if (MS1PeakCluster.endRT >= peakCluster.startRT && MS1PeakCluster.endRT <= peakCluster.endRT) {
                overlap = true;
            } else if (peakCluster.startRT >= MS1PeakCluster.startRT && peakCluster.startRT <= MS1PeakCluster.endRT) {
                overlap = true;
            } else if (peakCluster.endRT >= MS1PeakCluster.startRT && peakCluster.endRT <= MS1PeakCluster.endRT) {
                overlap = true;
            }
            if (overlap) {

                float corr = 0f;

                try {
                    corr = PeakCurveCorrCalc.CalPeakCorr(MS1PeakCluster.MonoIsotopePeak, peakCluster.MonoIsotopePeak, parameter.NoPeakPerMin);
                } catch (IOException ex) {
                    Logger.getLogger(CorrCalcCluster2ClusterUnit.class.getName()).log(Level.SEVERE, null, ex);
                }
                if (!Float.isNaN(corr) && corr > parameter.CorrThreshold) {
                    PrecursorFragmentPairEdge overlapRegion = new PrecursorFragmentPairEdge();
                    overlapRegion.Correlation = corr;
                    overlapRegion.PeakCurveIndexA = MS1PeakCluster.Index;
                    overlapRegion.PeakCurveIndexB = peakCluster.Index;
                    overlapRegion.FragmentMz = peakCluster.TargetMz();
                    overlapRegion.Intensity = peakCluster.PeakHeight[0];
                    ResultList.add(overlapRegion);
                }
            }
        }

        if (update != null) {
            update.Update();
        }
    }
}
