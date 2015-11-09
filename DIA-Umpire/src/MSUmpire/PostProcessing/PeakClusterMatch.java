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
package MSUmpire.PostProcessing;
import MSUmpire.BaseDataStructure.InstrumentParameter;
import MSUmpire.BaseDataStructure.XYPointCollection;
import MSUmpire.PeakDataStructure.PeakCluster;
import MSUmpire.SpectralProcessingModule.ScoreFunction;
import org.apache.log4j.Logger;

/**
 *
 * @author Chih-Chiang Tsou
 */
public class PeakClusterMatch implements Runnable{

    PeakCluster peakClusterA;
    PeakCluster peakClusterB;
    InstrumentParameter parameter;
    public float similiarity=0f;
    public PeakClusterMatch(PeakCluster peakClusterA, PeakCluster peakClusterB, InstrumentParameter parameter){
        this.peakClusterA=peakClusterA;
        this.peakClusterB=peakClusterB;
        this.parameter=parameter;
    }
    
    @Override
    public void run() {
        try {
            similiarity=0f;            
            XYPointCollection Scan1=peakClusterA.GetNormalizedFragmentScan();
            XYPointCollection Scan2=peakClusterB.GetNormalizedFragmentScan();
            if(Scan1.PointCount()>2 && Scan2.PointCount()>2){
                similiarity=ScoreFunction.CalcDotProductForScan(Scan1,Scan2);
            }            
        } catch (InterruptedException ex) {
            Logger.getRootLogger().error("error");
        }
    }

}
