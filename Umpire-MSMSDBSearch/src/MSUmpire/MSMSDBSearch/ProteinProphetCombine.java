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
package MSUmpire.MSMSDBSearch;

import Utility.PrintThread;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import org.apache.log4j.Logger;

/**
 *
 * @author Chih-Chiang Tsou <chihchiang.tsou@gmail.com>
 */
public class ProteinProphetCombine {
    
    public void ProteinProphetCombineSearch(DBSearchParam searchpara, ArrayList<String> PepXMLs, String OutputProtXML) throws IOException, InterruptedException {

        if(new File(OutputProtXML).exists()){
            return;
        }
        String[] ProteinProphetcmd = new String[PepXMLs.size() + 2];
        ProteinProphetcmd[0] = searchpara.xinteractpath.replace("xinteract", "ProteinProphet");

        for (int i = 0; i < PepXMLs.size(); i++) {
            ProteinProphetcmd[i + 1] = PepXMLs.get(i);
        }
        ProteinProphetcmd[ProteinProphetcmd.length - 1] = OutputProtXML;

        Process p = Runtime.getRuntime().exec(ProteinProphetcmd);
        Logger.getRootLogger().info("processing ProteinProphet...." + OutputProtXML);
        Logger.getRootLogger().debug("Command: " + Arrays.toString(ProteinProphetcmd));
        
        PrintThread printThread = new PrintThread(p);
        printThread.start();
        p.waitFor();
        if (p.exitValue() != 0) {
            Logger.getRootLogger().info("ProteinProphet...." + OutputProtXML + " failed");
        }
    }
    
    public void ProteinProphetCombineSearchiProphet(DBSearchParam searchpara, ArrayList<String> PepXMLs, String OutputProtXML) throws IOException, InterruptedException {

        if(new File(OutputProtXML).exists()){
            return;
        }
        String[] ProteinProphetcmd = new String[PepXMLs.size() + 4];
        ProteinProphetcmd[0] = searchpara.xinteractpath.replace("xinteract", "ProteinProphet");

        for (int i = 0; i < PepXMLs.size(); i++) {
            ProteinProphetcmd[i + 1] = PepXMLs.get(i);
        }
        ProteinProphetcmd[ProteinProphetcmd.length - 3] = OutputProtXML;
        ProteinProphetcmd[ProteinProphetcmd.length - 2] = "IPROPHET";
        ProteinProphetcmd[ProteinProphetcmd.length - 1] = "MINPROB0.5";

        Process p = Runtime.getRuntime().exec(ProteinProphetcmd);
        Logger.getRootLogger().info("processing ProteinProphet...." + OutputProtXML);
        Logger.getRootLogger().debug("Command: " + Arrays.toString(ProteinProphetcmd));
        
        PrintThread printThread = new PrintThread(p);
        printThread.start();
        p.waitFor();
        if (p.exitValue() != 0) {
            Logger.getRootLogger().info("ProteinProphet...." + OutputProtXML + " failed");
        }
    }
}
