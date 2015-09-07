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
package MSUmpire.PeptidePeakClusterDetection;

import MSUmpire.BaseDataStructure.ScanCollection;
import MSUmpire.BaseDataStructure.XYData;
import MSUmpire.LCMSBaseStructure.LCMSPeakBase;
import java.io.*;
import java.util.ArrayList;
import java.util.concurrent.ExecutionException;
import org.apache.log4j.Logger;

/**
 *
 * @author Chih-Chiang Tsou <chihchiang.tsou@gmail.com>
 */
public class PDHandlerMS1 extends PDHandlerBase {

    public PDHandlerMS1(LCMSPeakBase lcmspeak, int NoCPUs, float PPM) {
        this.NoCPUs = NoCPUs;
        this.PPM = PPM;
        this.LCMSPeakBase = lcmspeak;
        this.parameter = lcmspeak.parameter;
        //this.connectionManager = lcmspeak.connectionManager;
    }

    public void DetectPeakClusters(ArrayList<ScanCollection> scanCollections) throws InterruptedException, ExecutionException, IOException {        
        DetectSingleMZTraces(scanCollections);
        PeakCurveCorrClustering_V2(new XYData(Float.NEGATIVE_INFINITY, Float.POSITIVE_INFINITY));
    }

    public void DetectSingleMZTraces(ArrayList<ScanCollection> scanCollections) throws IOException {
        ReadPepIsoMS1PatternMap();
        LCMSPeakBase.UnSortedPeakCurves = new ArrayList<>();
        for (ScanCollection scanCollection : scanCollections) {
            FindAllPeakCurve(scanCollection);
        }
        Logger.getRootLogger().info("Inclusion mz values found: "+InclusionCheckInfo());        
        WaveletDetectMax();      
        //CreateSortedPeakCurveList();
        ClearRawPeaks();
    }
    
}
