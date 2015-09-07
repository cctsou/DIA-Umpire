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
import org.apache.commons.lang.exception.ExceptionUtils;
import org.apache.log4j.Logger;

/**
 *
 * @author Chih-Chiang Tsou <chihchiang.tsou@gmail.com>
 */
public abstract class MSMSDBSearch implements Runnable {

    public boolean Resume = true;

    public abstract void DBSearch();
    public abstract DBSearchParam GetParameter();
    public abstract MSMSDBSearch Clone();
    public abstract void GeneratePBS();
    
    protected void PBSHeader(StringBuilder sb) {
        sb.append("#PBS -N " + FilenameUtils.getBaseName(GetParameter().PepXMLPath) + "\n");
        sb.append("#PBS -A nesvi_fluxod\n");
        sb.append("#PBS -q fluxod\n");
        sb.append("#PBS -l qos=flux\n");
        sb.append("#PBS -l nodes=1:ppn="+GetParameter().NoCPUs+",walltime=8:00:00,pmem=3gb\n");
        sb.append("#PBS -j oe\n");
        sb.append("#PBS -M tsouc@umich.edu\n");
        sb.append("#PBS -m abe\n");
        sb.append("#PBS -V\n");
    }
    
    public void GenerateParamFile() {
        GetParameter().GenerateParamFile();
    }
    public void SetResultFilePath(String mzXMLfile) {
        GetParameter().SetResultFilePath(mzXMLfile);
    }
    
    public void SetCombineFileName(String filename){
        SetCombineFileName(filename,"");
    }
    
    public void SetCombineFileName(String filename,String tag){
        GetParameter().SetCombineFileName(filename,tag);
    }
    
    public void CombineProteinProphet(ArrayList<String> PepXMLs) throws IOException, InterruptedException{             
        ProteinProphetCombine proteinProphetCombine=new ProteinProphetCombine();
        proteinProphetCombine.ProteinProphetCombineSearch(GetParameter(), PepXMLs, GetParameter().CombinedProt);
    }
    
    public void CombinePepXML(ArrayList<String> PepXMLs) {
        if (new File(GetParameter().CombinedPepXML).exists()) {
            return;
        }
        try {
            Process p = null;
            ArrayList<String> cmdlist = new ArrayList<>();
            cmdlist.add(GetParameter().xinteractpath);
            for (String op : GetParameter().xinteractpara.split(" ")) {
                cmdlist.add(op);
            }
            for (String pepxml : PepXMLs) {
                cmdlist.add(pepxml);
            }
            cmdlist.add("-N" + GetParameter().CombinedPepXML);
            String[] xinteractcmd = new String[cmdlist.size()];
            xinteractcmd = cmdlist.toArray(xinteractcmd);
            p = Runtime.getRuntime().exec(xinteractcmd);
            
            Logger.getRootLogger().info("xinteract combined analysis...." + GetParameter().CombinedPepXML + " Option:" + GetParameter().xinteractpara );

            PrintThread printThread = new PrintThread(p);
            printThread.start();
            p.waitFor();            
        } catch (Exception ex) {
            Logger.getRootLogger().error(ExceptionUtils.getStackTrace(ex));
        }
    }
    
    public void ConvertFile() throws InterruptedException, IOException{
         if (GetParameter().SpectrumPath.toLowerCase().endsWith(".mgf")) {
            if (!new File(GetParameter().SpectrumPath.replace(".mgf", ".mzXML")).exists()) {
                String[] msconvertcmd = {GetParameter().msconvertpath, "--mzXML", "--32", "-z", GetParameter().SpectrumPath, "-o", FilenameUtils.getFullPath(GetParameter().SpectrumPath)};
                Process p = Runtime.getRuntime().exec(msconvertcmd);
                Logger.getRootLogger().info("MGF file coversion by msconvert.exe...." + GetParameter().SpectrumPath);
                Logger.getRootLogger().debug("Command: " + Arrays.toString(msconvertcmd));
                PrintThread printThread = new PrintThread(p);
                printThread.start();
                p.waitFor();
                if (p.exitValue() != 0) {
                    Logger.getRootLogger().info("msconvert : " + GetParameter().SpectrumPath + " failed");
                    //PrintOutput(p);
                    return;
                }
                //PrintOutput(p);
            }
            GetParameter().SpectrumPath = GetParameter().SpectrumPath.replace(".mgf", ".mzXML");
        }
    }

    @Override
    public void run() {
        DBSearch();
    }
}
