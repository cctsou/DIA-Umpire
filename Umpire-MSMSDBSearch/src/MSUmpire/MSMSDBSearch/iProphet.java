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
public class iProphet implements Runnable {

    ArrayList<MSMSDBSearch> searches;
    String ipropepxml;
    String iprophetpath;
    String proteinprophetpath;

    public iProphet(String ipropepxml, ArrayList<MSMSDBSearch> searches) {
        this.searches = searches;
        this.ipropepxml = ipropepxml;
        this.iprophetpath = FilenameUtils.getFullPath(searches.get(0).GetParameter().xinteractpath) + "InterProphetParser";
        this.proteinprophetpath = FilenameUtils.getFullPath(searches.get(0).GetParameter().xinteractpath) + "ProteinProphet";
    }

    public void DoiProphetCombineProteinXML(DBSearchParam searchpara, ArrayList<String> PepXMLs, String OutputProtXML) throws IOException, InterruptedException {

        if (new File(OutputProtXML).exists()) {
            return;
        }
        String[] ProteinProphetcmd = new String[PepXMLs.size() + 3];
        ProteinProphetcmd[0] = proteinprophetpath;
        ProteinProphetcmd[1] = "IPROPHET";

        for (int i = 0; i < PepXMLs.size(); i++) {
            ProteinProphetcmd[i + 2] = PepXMLs.get(i);
        }
        ProteinProphetcmd[ProteinProphetcmd.length - 1] = OutputProtXML;

        Logger.getRootLogger().info("processing ProteinProphet for iProphet pepXML analysis...." + OutputProtXML);
        Logger.getRootLogger().debug("Command: " + Arrays.toString(ProteinProphetcmd));
        Process p = Runtime.getRuntime().exec(ProteinProphetcmd);
        

        PrintThread printThread = new PrintThread(p);
        printThread.start();
        p.waitFor();
        if (p.exitValue() != 0) {
            Logger.getRootLogger().info("ProteinProphet for iProphet pepXML analysis...." + OutputProtXML + " failed");
        }
    }

    public void DoiProphetProteinXML() throws IOException, InterruptedException {
        if (new File(ipropepxml.replace("pep.xml", "prot.xml")).exists()) {
            return;
        }
        Process p = null;
        String[] protcmd = {proteinprophetpath, "IPROPHET", ipropepxml, ipropepxml.replace("pep.xml", "prot.xml")};
        Logger.getRootLogger().info("ProteinProphet for iProphet pepXML analysis...." + Arrays.toString(protcmd));
        Logger.getRootLogger().debug("Command: " + Arrays.toString(protcmd));
        
        p = Runtime.getRuntime().exec(protcmd);

        PrintThread printThread = new PrintThread(p);
        printThread.start();
        p.waitFor();
        if (p.exitValue() != 0) {
            Logger.getRootLogger().info("ProteinProphet for iProphet pepXML analysis. : " + Arrays.toString(protcmd) + " failed");
            //PrintOutput(p);
        }
    }

    public void DoiProphetPepXML() throws InterruptedException, IOException {
        if (!(new File(ipropepxml)).exists()) {
            Process p = null;
            ArrayList<String> cmdlist = new ArrayList<>();
            cmdlist.add(iprophetpath);
            cmdlist.add("THREADS=" + searches.get(0).GetParameter().NoCPUs);
            cmdlist.add("NONSP");
            for (MSMSDBSearch search : searches) {
                cmdlist.add(search.GetParameter().InteractPepXMLPath);
            }
            cmdlist.add(ipropepxml);
            String[] iprohcmd = new String[cmdlist.size()];
            iprohcmd = cmdlist.toArray(iprohcmd);
            Logger.getRootLogger().info("InterProphetParser iProphet analysis...." + Arrays.toString(iprohcmd));
            Logger.getRootLogger().debug("Command: " + Arrays.toString(iprohcmd));
            
            p = Runtime.getRuntime().exec(iprohcmd);
            PrintThread printThread = new PrintThread(p);
            printThread.start();
            p.waitFor();
            if (p.exitValue() != 0) {
                Logger.getRootLogger().info("InterProphetParser iProphet analysis. : " + Arrays.toString(iprohcmd) + " failed");
                //PrintOutput(p);
                return;
            }
        }
    }

    public void DoiProphetByCombinePepXML() throws InterruptedException, IOException {
        if ((new File(ipropepxml)).exists()) {
            Logger.getRootLogger().info("iProphet pepXML file exists : " + ipropepxml);
            return;
        }
        Process p = null;
        ArrayList<String> cmdlist = new ArrayList<>();
        cmdlist.add(iprophetpath);
        cmdlist.add("THREADS=" + searches.get(0).GetParameter().NoCPUs);
        cmdlist.add("NONSP");
        for (MSMSDBSearch search : searches) {
            cmdlist.add(search.GetParameter().CombinedPepXML);
        }
        cmdlist.add(ipropepxml);
        String[] iprohcmd = new String[cmdlist.size()];
        iprohcmd = cmdlist.toArray(iprohcmd);
        Logger.getRootLogger().info("InterProphetParser iProphet analysis...." + Arrays.toString(iprohcmd));
        Logger.getRootLogger().debug("Command: " + Arrays.toString(iprohcmd));
        p = Runtime.getRuntime().exec(iprohcmd);

        PrintThread printThread = new PrintThread(p);
        printThread.start();
        p.waitFor();
        if (p.exitValue() != 0) {
            Logger.getRootLogger().info("InterProphetParser iProphet analysis. : " + Arrays.toString(iprohcmd) + " failed");
            //PrintOutput(p);
            return;
        }

        String[] protcmd = {proteinprophetpath, "IPROPHET", ipropepxml, ipropepxml.replace("pep.xml", "prot.xml")};

        Logger.getRootLogger().info("ProteinProphet for iProphet pepXML analysis...." + Arrays.toString(protcmd));
        Logger.getRootLogger().debug("Command: " + Arrays.toString(protcmd));
        
        p = Runtime.getRuntime().exec(protcmd);

        printThread = new PrintThread(p);
        printThread.start();
        p.waitFor();
        if (p.exitValue() != 0) {
            Logger.getRootLogger().info("ProteinProphet for iProphet pepXML analysis. : " + Arrays.toString(protcmd) + " failed");
            //PrintOutput(p);
            return;
        }

    }

    @Override
    public void run() {
        try {
            DoiProphetByCombinePepXML();
        } catch (InterruptedException ex) {
            Logger.getRootLogger().debug(ex.getMessage());
        } catch (IOException ex) {
            Logger.getRootLogger().debug(ex.getMessage());
        }
    }
}
