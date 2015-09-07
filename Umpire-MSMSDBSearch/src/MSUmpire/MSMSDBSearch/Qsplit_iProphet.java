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
import org.apache.commons.io.FilenameUtils;
import org.apache.log4j.Logger;

/**
 *
 * @author Chih-Chiang Tsou
 */
public class Qsplit_iProphet implements Runnable{

    ArrayList<MSMSDBSearch> searches;
    String ipropepxml;
    String iprophetpath;
    String proteinprophetpath;
    
    public Qsplit_iProphet(String ipropepxml, ArrayList<MSMSDBSearch> searches) {
        this.searches=searches;
        this.ipropepxml=ipropepxml;
        this.iprophetpath=FilenameUtils.getFullPath(searches.get(0).GetParameter().xinteractpath)+"InterProphetParser";
        this.proteinprophetpath=FilenameUtils.getFullPath(searches.get(0).GetParameter().xinteractpath)+"ProteinProphet";
    }

    public void DoiProphet() throws InterruptedException, IOException{
        if (!(new File(ipropepxml)).exists()) {            
            Process p = null;
            ArrayList<String> cmdlist = new ArrayList<>();
            cmdlist.add(iprophetpath);
            cmdlist.add("THREADS="+searches.get(0).GetParameter().NoCPUs);
            for (MSMSDBSearch search : searches) {
                cmdlist.add(search.GetParameter().CombinedPepXML);
            }
            cmdlist.add(ipropepxml);
            String[] iprohcmd = new String[cmdlist.size()];
            iprohcmd = cmdlist.toArray(iprohcmd);
            p = Runtime.getRuntime().exec(iprohcmd);

            Logger.getRootLogger().info("InterProphetParser iProphet analysis...." + Arrays.toString(iprohcmd));
            Logger.getRootLogger().debug("Command: " + Arrays.toString(iprohcmd));

            PrintThread printThread=new PrintThread(p);
            printThread.start();
            p.waitFor();
            if (p.exitValue() != 0) {
                Logger.getRootLogger().info("InterProphetParser iProphet analysis. : " + Arrays.toString(iprohcmd) + " failed");
                //PrintOutput(p);
                return;
            }
            
            String[] protcmd ={ proteinprophetpath, "IPROPHET" , ipropepxml, ipropepxml.replace("pep.xml", "prot.xml")};
            p = Runtime.getRuntime().exec(protcmd);
            
            Logger.getRootLogger().info("ProteinProphet for iProphet pepXML analysis...." + Arrays.toString(protcmd));
            Logger.getRootLogger().debug("Command: " + Arrays.toString(protcmd));

            printThread=new PrintThread(p);
            printThread.start();
            p.waitFor();
            if (p.exitValue() != 0) {
                Logger.getRootLogger().info("ProteinProphet for iProphet pepXML analysis. : " + Arrays.toString(protcmd) + " failed");
                //PrintOutput(p);
                return;
            }
        }
    }
    
    @Override
    public void run() {
        try {
            DoiProphet();
        } catch (InterruptedException ex) {
            Logger.getRootLogger().debug(ex.getMessage());
        } catch (IOException ex) {
            Logger.getRootLogger().debug(ex.getMessage());
        }
    }
    
}
