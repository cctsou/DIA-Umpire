/*
 * Copyright 2014 Chih-Chiang Tsou.
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

package diaumpire;

import MSUmpire.DIA.DIAPack;
import MSUmpire.MSMSDBSearch.TandemParam;
import MSUmpire.PSMDataStructure.LCMSID;
import org.apache.commons.lang.exception.ExceptionUtils;
import org.apache.log4j.Logger;

/**
 *
 * @author Chih-Chiang Tsou
 */
public class QuantThread implements Runnable{
        DIAPack DiaFile;
        TandemParam tandemPara;
        LCMSID combineID;
        public QuantThread(DIAPack DIAfile, TandemParam searchpara, LCMSID combineID) throws CloneNotSupportedException{
            this.DiaFile=DIAfile;
            this.tandemPara=(TandemParam) searchpara.clone();
            this.combineID=combineID;
        }

        @Override
       public void run() {
           String tag="";
           if(combineID!=null){
               tag="DataSetFDR";
           }
            try {
                if (DiaFile.ReadSerializedLCMSID(tag)) {        
                    DiaFile.IDsummary.ClearMappedPep();
                    DiaFile.IDsummary.ReduceMemoryUsage();
                    DiaFile.IDsummary.FastaPath=tandemPara.FastaPath;
                    return;
                }
                                
                DiaFile.ParseSearchEngineResultiProphet(tandemPara, combineID);
                //DiaFile.ParseSearchEngineResult(tsearch);
                //DiaFile.ParsePepXML(tandemPara);                                   
                //DiaFile.IDsummary.CheckRT("step 1");
                DiaFile.BuildStructure();
                DiaFile.ms1lcms.ReadPeakCluster();                
                DiaFile.GenerateClusterScanNomapping();     
                //DiaFile.IDsummary.CheckRT("step 2");
                DiaFile.AssignQuant(false);
                //DiaFile.IDsummary.CheckRT("step 3");
                DiaFile.ExportID(tag);
                DiaFile.ms1lcms.ExportPeakClusterResultCSV();
                //DiaFile.IDsummary.CheckRT("step 4");
                DiaFile.ClearStructure();
            } catch (Exception ex) {
                Logger.getRootLogger().error(ExceptionUtils.getStackTrace(ex));
                System.exit(2);
            }
       }
    }
