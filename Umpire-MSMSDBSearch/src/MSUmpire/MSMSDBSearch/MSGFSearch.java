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
import java.io.FileWriter;
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
public class MSGFSearch extends MSMSDBSearch {

    MSGFParam parameter;

    public MSGFSearch(MSGFParam parameter) {        
        this.parameter = parameter;
    }

    public void BuildDatabaseIndex() throws IOException, InterruptedException {
        //Usage: java -Xmx3500M -cp MSGFPlus.jar edu.ucsd.msjava.msdbsearch.BuildSA -d DatabaseFile (*.fasta or *.fa) [-tda 0/1/2] (0: target only, 1: target-decoy database only, 2: both)
        if (!new File(parameter.FastaPath.replace("fasta", "").replace("fa", "") + ".cnlcp").exists()) {
            String[] buildcmd = {parameter.JavaPath, "-Xmx5G", "-cp", parameter.msgfpath, "edu.ucsd.msjava.msdbsearch.BuildSA", "-d", parameter.FastaPath, "-tda", "0"};
            Process p = Runtime.getRuntime().exec(buildcmd);           
            Logger.getRootLogger().info("Building index of fasta: " + parameter.FastaPath + " for MSGF search");
            Logger.getRootLogger().debug("Command: " + Arrays.toString(buildcmd));
            PrintThread printThread=new PrintThread(p);
            printThread.start();
            p.waitFor();             
            if (p.exitValue() != 0) {
                Logger.getRootLogger().info("Building index of fasta: " + parameter.FastaPath + " for MSGF search failed");
                return;
            }   
        }
    }
    
    public void GeneratePBS() {
        FileWriter writer = null;
        try {
            StringBuilder sb = new StringBuilder();
            PBSHeader(sb);
            sb.append("if [ ! -f " + parameter.SpectrumPath.replace(".mgf", ".mzXML") + " ]; then\n");
            sb.append(parameter.msconvertpath + " --mzXML --32 -z " + parameter.SpectrumPath + " -o " + FilenameUtils.getFullPath(parameter.SpectrumPath) + "\n\n");
            sb.append("fi\n");
            parameter.SpectrumPath = parameter.SpectrumPath.replace(".mgf", ".mzXML");
            parameter.GenerateParamFile();
            String isotopeerror = "0,0";
            int EnzymeID = 1;
            if (parameter.IsotopeError) {
                isotopeerror = "0,1";
            }
            if (parameter.NonSpecificCleavage) {
                EnzymeID = 0;
            }
            sb.append(parameter.JavaPath + " -jar -Xmx" + parameter.NoCPUs * 3 + "G " + parameter.msgfpath + " -s " + parameter.SpectrumPath + " -d " + parameter.FastaPath + " -o " + parameter.RawSearchResult + " -t " + String.valueOf(Math.round(parameter.PrecursorPPM)) + "ppm -thread " + String.valueOf(parameter.NoCPUs) + " -tda 0 -m 0 -ti " + isotopeerror + " -inst " + String.valueOf(parameter.MSGFInstrumentID) + " -e " + EnzymeID + " -mod " + parameter.parameterPath + "\n\n");
            sb.append(parameter.idconvert + " " + parameter.RawSearchResult + " --pepXML -o " + FilenameUtils.getFullPath(parameter.PepXMLPath) + "\n\n");
            //sb.append(parameter.xinteractpath + " " + parameter.xinteractpara + " " + parameter.PepXMLPath + " -N" + parameter.InteractPepXMLPath + "\n\n");
            sb.append("exit\n");
            writer = new FileWriter(parameter.parameterPath.replace("param", "pbs"));
            writer.write(sb.toString());
            writer.close();
        } catch (Exception ex) {
            Logger.getRootLogger().error(ExceptionUtils.getStackTrace(ex));
        }
    }
    public void RunMSGF() throws IOException, InterruptedException {        
        ConvertFile();
        parameter.GenerateParamFile();
        
        if (parameter.Overwrite || !new File(parameter.RawSearchResult).exists() || new File(parameter.RawSearchResult).length() / 1024 < 2) {
            String isotopeerror="0,0";
            if(parameter.IsotopeError){
                isotopeerror="0,1";
            }            
            if(parameter.NonSpecificCleavage){
                 parameter.MSGFEnzymeID=0;
            }
            //-jar -Xmx8G MSGFPlus.jar -s 18484_REP3_1ug_Ecoli_NewStock2_SWATH_1_Q2.mzXML -d Ecoli_PlusRev.fa -t 30ppm -thread 10 -tda 0 -m 1 -inst 2 -mod Mods.txt
            //String[] msgfcmd = {"-s", parameter.SpectrumPath, "-d", parameter.FastaPath,"-o",parameter.RawSearchResult, "-t", String.valueOf(Math.round(parameter.PrecursorPPM)) + "ppm", "-thread", String.valueOf(parameter.NoCPUs), "-tda", "0", "-m", String.valueOf(parameter.MissCleavage), "-inst", "2", "-mod", parameter.parameterPath};
            String[] msgfcmd = {parameter.JavaPath, "-jar","-Xmx"+parameter.NoCPUs*3+"G", parameter.msgfpath, "-s", parameter.SpectrumPath, "-d", parameter.FastaPath,"-o",parameter.RawSearchResult, "-t", String.valueOf(Math.round(parameter.PrecursorPPM)) + "ppm", "-thread", String.valueOf(parameter.NoCPUs), "-tda", "0", "-m",  String.valueOf(parameter.MSGFFragmentMethodID),"-ti",isotopeerror, "-inst", String.valueOf(parameter.MSGFInstrumentID), "-e", String.valueOf(parameter.MSGFEnzymeID), "-mod", parameter.parameterPath};
            Logger.getRootLogger().info("processing MSGF+ search " + parameter.SpectrumPath);
            Logger.getRootLogger().debug("Command: " + Arrays.toString(msgfcmd));
            
            Process p = Runtime.getRuntime().exec(msgfcmd);                                    
            PrintThread printThread = new PrintThread(p);
            printThread.start();
            p.waitFor();
            if (p.exitValue() != 0) {
                Logger.getRootLogger().info("msgf : " + parameter.RawSearchResult + " failed");                  
                return;
            }
        }

        if (parameter.Overwrite || !new File(parameter.PepXMLPath).exists() || new File(parameter.PepXMLPath).length() / 1024 < 2) {
            //idconvert 18484_REP3_1ug_Ecoli_NewStock2_SWATH_1_Q2.mzid --pepXML
            String[] idconvertcmd = {parameter.idconvert, parameter.RawSearchResult, "--pepXML","-o",FilenameUtils.getFullPath(parameter.PepXMLPath)};
            Logger.getRootLogger().info("Converting .mzid to .pep.xml by idconvert...." + parameter.PepXMLPath);
            Logger.getRootLogger().debug("Command: " + Arrays.toString(idconvertcmd));
            Process p = Runtime.getRuntime().exec(idconvertcmd);                                    
            PrintThread printThread=new PrintThread(p);
            printThread.start();
            p.waitFor();
            if (p.exitValue() != 0) {
                System.out.print("idconvert : " + parameter.RawSearchResult + " failed");                
                return;
            }
        }

        if (parameter.Overwrite || !(new File(parameter.InteractPepXMLPath)).exists()) {
            new File(parameter.InteractPepXMLPath).delete();
            new File(parameter.ProtXMLPath).delete();

            Process p = null;
            ArrayList<String> cmdlist = new ArrayList<>();
            cmdlist.add(parameter.xinteractpath);
            for (String op : parameter.xinteractpara.split(" ")) {
                cmdlist.add(op);
            }
            cmdlist.add(parameter.PepXMLPath);
            cmdlist.add("-N" + parameter.InteractPepXMLPath);
            String[] xinteractcmd = new String[cmdlist.size()];
            xinteractcmd = cmdlist.toArray(xinteractcmd);
            p = Runtime.getRuntime().exec(xinteractcmd);

            Logger.getRootLogger().info("xinteract PeptideProphet/ProteinProphet analysis...." + parameter.InteractPepXMLPath + " Option:" + parameter.xinteractpara );
            Logger.getRootLogger().debug("Command: " + Arrays.toString(xinteractcmd));

            PrintThread printThread=new PrintThread(p);
            printThread.start();
            p.waitFor();
            if (p.exitValue() != 0) {
                Logger.getRootLogger().info("xinteract : " + parameter.PepXMLPath + " failed");
                return;
            }        
        }
    }

    @Override
    public void DBSearch() {
        try {
            RunMSGF();
        } catch (Exception ex) {
            Logger.getRootLogger().error(ExceptionUtils.getStackTrace(ex));
        } 
    }

//    @Override
//    public void GenerateParamFile() {
//        this.parameter.GenerateParamFile();
//    }
//    
//     @Override
//    public void SetResultFilePath(String mzXMLfile) {
//        parameter.SetResultFilePath(mzXMLfile);
//    }
    
     @Override
    public DBSearchParam GetParameter() {
       return parameter;
    }
    
    public MSMSDBSearch Clone() {
        try {
            MSGFParam newparameter = (MSGFParam) parameter.clone();
            MSGFSearch newsearch = new MSGFSearch(newparameter);
            return newsearch;
        } catch (Exception ex) {
            Logger.getRootLogger().error(ExceptionUtils.getStackTrace(ex));       
        }
        return null;
    }
}
