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
import MSUmpire.LCMSBaseStructure.LCMSPeakMS1Meta;
import java.io.*;
import java.sql.SQLException;
import java.util.ArrayList;
import java.util.TreeMap;
import java.util.concurrent.ExecutionException;

/**
 *
 * @author Chih-Chiang Tsou <chihchiang.tsou@gmail.com>
 */
public class PDHandlerMS1Meta extends PDHandlerBase {

    public PDHandlerMS1Meta(LCMSPeakMS1Meta lCMSQuant, int NoCPUs, float PPM) throws SQLException {
        this.NoCPUs = NoCPUs;
        this.PPM = PPM;
        this.LCMSPeakBase = lCMSQuant;
        this.parameter = lCMSQuant.parameter;
        //this.connectionManager = lCMSQuant.connectionManager;
    }

    public void DetectPeakCurves(ScanCollection scanCollection) throws InterruptedException, ExecutionException, IOException, SQLException {
        ReadMetaIsoPatternMap();
        LCMSPeakBase.UnSortedPeakCurves = new ArrayList<>();
        FindAllPeakCurve(scanCollection);
        WaveletDetectMax();
        //CreateSortedPeakCurveList();
        PeakCurveCorrClustering_V2(new XYData(Float.NEGATIVE_INFINITY, Float.POSITIVE_INFINITY));
    }

    private void ReadMetaIsoPatternMap() throws FileNotFoundException, IOException {

        InputStream is = this.getClass().getClassLoader().getResourceAsStream("resource/MetaIsotopicPatternRange.csv");
        BufferedReader reader = new BufferedReader(new InputStreamReader(is));

        reader.readLine();
        IsotopePatternMap = new TreeMap[LCMSPeakBase.MinNoPeakCluster];
        IsotopePatternMap[0] = new TreeMap<>();
        IsotopePatternMap[1] = new TreeMap<>();
        String line = "";
        while ((line = reader.readLine()) != null) {
            float MW = Float.parseFloat(line.split(",")[0]);
            float MeanSecond = Float.parseFloat(line.split(",")[1]);
            float SDSecond = Float.parseFloat(line.split(",")[2]);
            float MeanThird = Float.parseFloat(line.split(",")[3]);
            float SDThird = Float.parseFloat(line.split(",")[4]);

            if (!Float.isNaN(MeanSecond)) {
                //IsoMapSecond.put(MW,new Normal(MeanSecond, SDSecond));
                //IsoMapThird.put(MW,new Normal(MeanThird, SDThird));                
                IsotopePatternMap[0].put(MW, new XYData(MeanSecond + 3.3f * SDSecond, MeanSecond - 3.3f * SDSecond));
                IsotopePatternMap[1].put(MW, new XYData(MeanThird + 3.3f * SDThird, MeanThird - 3.3f * SDThird));
            }
        }
        reader.close();
    }

}
