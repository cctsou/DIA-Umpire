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
package MSUmpire.LCMSBaseStructure;


import MSUmpire.MSMSDBSearch.DBSearchParam;
import MSUmpire.PSMDataStructure.LCMSID;
import MSUmpire.PSMDataStructure.LCMSID;
import MSUmpire.PSMDataStructure.PepIonID;
import MSUmpire.PSMDataStructure.PepIonID;
import java.io.IOException;
import java.util.ArrayList;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;
import org.apache.log4j.Logger;

/**
 *
 * @author Chih-Chiang Tsou
 */
public class FDR_DataSetLevel_DDA {
    ExecutorService executorPool = null;
    public LCMSID combineID=null;
    
    public void GeneratePepIonList(ArrayList<LCMSPeakMS1> DDAFilenameList, DBSearchParam param, String combineIDPath) throws IOException{
        
        executorPool = Executors.newFixedThreadPool(param.NoCPUs);
        for (LCMSPeakMS1 DDA : DDAFilenameList) {
            DDA_ParsePepXMLThreadUnit thread = new DDA_ParsePepXMLThreadUnit(DDA, param);
            executorPool.execute(thread);
        }
        executorPool.shutdown();
        try {
            executorPool.awaitTermination(Long.MAX_VALUE, TimeUnit.NANOSECONDS);
        } catch (InterruptedException e) {
            Logger.getRootLogger().info("interrupted..");
        }
        
        //Estimate peptide level PepFDR in whole dataset
        combineID = new LCMSID(combineIDPath,param.DecoyPrefix,param.FastaPath);
        for (LCMSPeakMS1 Diafile : DDAFilenameList) {
            LCMSID lcms=Diafile.IDsummary;
            for (PepIonID pepIonID : lcms.GetPepIonList().values()) {
                if (!combineID.GetPepIonList().containsKey(pepIonID.GetKey())) {
                    PepIonID newpep = pepIonID.ClonePepIonID();
                    if (pepIonID.IsDecoy(param.DecoyPrefix)) {
                        newpep.IsDecoy = 1;
                    } else {
                        newpep.IsDecoy = 0;
                    }
                    combineID.AddPeptideID(newpep);
                }
                if (combineID.GetPepIonList().get(pepIonID.GetKey()).MaxProbability < pepIonID.MaxProbability) {
                    combineID.GetPepIonList().get(pepIonID.GetKey()).MaxProbability = pepIonID.MaxProbability;
                }
            }
        }
        //combineID.GeneratePepSeqList();
        combineID.DecoyTag = param.DecoyPrefix;
        combineID.FDR = param.PepFDR;
        combineID.ROCPepByMaxIniProb(param.DecoyPrefix);
        combineID.FindPepProbThresholdByFDR();
        //combineID.FindPepProbThresholdByFDRAtPepSeq();
        combineID.RemoveDecoyPep();
        combineID.RemoveLowProbPep();
        //combineID.GeneratePepSeqList();
        ////////////////////////////
        
    }
    
}
