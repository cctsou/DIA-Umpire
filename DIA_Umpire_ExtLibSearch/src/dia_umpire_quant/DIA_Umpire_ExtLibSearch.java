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
package DIA_Umpire_Quant;

/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
import MSUmpire.BaseDataStructure.UmpireInfo;
import MSUmpire.DIA.DIAPack;
import MSUmpire.DIA.RTMappingExtLib;
import MSUmpire.FragmentLib.FragmentLibManager;
import MSUmpire.BaseDataStructure.DBSearchParam;
import MSUmpire.BaseDataStructure.TandemParam;
import MSUmpire.PSMDataStructure.PTMManager;
import MSUmpire.Utility.ConsoleLogger;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import org.apache.commons.io.FilenameUtils;
import org.apache.log4j.Level;
import org.apache.log4j.Logger;

/**
 *
 * @author Chih-Chiang Tsou
 */
public class DIA_Umpire_ExtLibSearch {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) throws FileNotFoundException, IOException, Exception {
        System.out.println("=================================================================================================");
        System.out.println("DIA-Umpire targeted re-extraction analysis using external library (version: " + UmpireInfo.GetInstance().Version + ")");
        if (args.length != 1) {
            System.out.println("command format error, the correct format should be: java -jar -Xmx10G DIA_Umpire_ExtLibSearch.jar diaumpire_module.params");
            return;
        }
        try {
            ConsoleLogger.SetConsoleLogger(Level.INFO);
            ConsoleLogger.SetFileLogger(Level.DEBUG, FilenameUtils.getFullPath(args[0]) + "diaumpire_extlibsearch.log");
        } catch (Exception e) {
        }

        Logger.getRootLogger().info("Version: "+UmpireInfo.GetInstance().Version);
        Logger.getRootLogger().info("Parameter file:" + args[0]);

        BufferedReader reader = new BufferedReader(new FileReader(args[0]));
        String line = "";
        String WorkFolder = "";
        int NoCPUs = 2;

        String ExternalLibPath = "";
        String ExternalLibDecoyTag = "DECOY";
        
        float ExtProbThreshold =0.99f;
        float RTWindow_Ext=-1f;
                
        TandemParam tandemPara = new TandemParam(DBSearchParam.SearchInstrumentType.TOF5600);
        HashMap<String, File> AssignFiles = new HashMap<>();

        //<editor-fold defaultstate="collapsed" desc="Reading parameter file">
        while ((line = reader.readLine()) != null) {
            line=line.trim();
            Logger.getRootLogger().info(line);
            if (!"".equals(line) && !line.startsWith("#")) {
                //System.out.println(line);
                if (line.equals("==File list begin")) {
                    do {
                        line = reader.readLine();
                        line = line.trim();
                        if (line.equals("==File list end")) {
                            continue;
                        } else if (!"".equals(line)) {
                            File newfile = new File(line);
                            if (newfile.exists()) {
                                AssignFiles.put(newfile.getAbsolutePath(), newfile);
                            } else {
                                Logger.getRootLogger().info("File: " + newfile + " does not exist.");
                            }
                        }
                    } while (!line.equals("==File list end"));
                }
                if (line.split("=").length < 2) {
                    continue;
                }
                String type = line.split("=")[0].trim();
                String value = line.split("=")[1].trim();
                switch (type) {                    

                    case "Path": {
                        WorkFolder = value;
                        break;
                    }
                    case "path": {
                        WorkFolder = value;
                        break;
                    }
                    case "Thread": {
                        NoCPUs = Integer.parseInt(value);
                        break;
                    }
                    case "Fasta": {
                        tandemPara.FastaPath = value;
                        break;
                    }
                    case "DecoyPrefix": {
                        if (!"".equals(value)) {
                            tandemPara.DecoyPrefix = value;
                        }
                        break;
                    }
                    case "ExternalLibPath": {
                        ExternalLibPath = value;
                        break;
                    }
                    case "ExtProbThreshold":{
                        ExtProbThreshold = Float.parseFloat(value);
                        break;
                    }
                      case "RTWindow_Ext": {
                        RTWindow_Ext = Float.parseFloat(value);
                        break;
                    }
                    case "ExternalLibDecoyTag": {
                        ExternalLibDecoyTag = value;
                        if(ExternalLibDecoyTag.endsWith("_")){
                           ExternalLibDecoyTag=ExternalLibDecoyTag.substring(0, ExternalLibDecoyTag.length()-1);
                        }
                        break;
                    }     
                }
            }
        }
//</editor-fold>

        //Initialize PTM manager using compomics library
        PTMManager.GetInstance();
        
        
        //Check if the fasta file can be found
        if (!new File(tandemPara.FastaPath).exists()) {
            Logger.getRootLogger().info("Fasta file :"+tandemPara.FastaPath + " cannot be found, the process will be terminated, please check.");
            System.exit(1);
        }               
                
        //Generate DIA file list
        ArrayList<DIAPack> FileList = new ArrayList<>();
       
            File folder = new File(WorkFolder);
            for (final File fileEntry : folder.listFiles()) {
                if (fileEntry.isFile() && (fileEntry.getAbsolutePath().toLowerCase().endsWith(".mzxml") | fileEntry.getAbsolutePath().toLowerCase().endsWith(".mzml"))
                        && !fileEntry.getAbsolutePath().toLowerCase().endsWith("q1.mzxml")
                        && !fileEntry.getAbsolutePath().toLowerCase().endsWith("q2.mzxml")
                        && !fileEntry.getAbsolutePath().toLowerCase().endsWith("q3.mzxml")) {
                    AssignFiles.put(fileEntry.getAbsolutePath(), fileEntry);
                }
                if (fileEntry.isDirectory()) {
                    for (final File fileEntry2 : fileEntry.listFiles()) {
                        if (fileEntry2.isFile() && (fileEntry2.getAbsolutePath().toLowerCase().endsWith(".mzxml") | fileEntry2.getAbsolutePath().toLowerCase().endsWith(".mzml"))
                                && !fileEntry2.getAbsolutePath().toLowerCase().endsWith("q1.mzxml")
                                && !fileEntry2.getAbsolutePath().toLowerCase().endsWith("q2.mzxml")
                                && !fileEntry2.getAbsolutePath().toLowerCase().endsWith("q3.mzxml")) {
                            AssignFiles.put(fileEntry2.getAbsolutePath(), fileEntry2);
                        }
                    }
                }
            }

            Logger.getRootLogger().info("No. of files assigned :" + AssignFiles.size());
            for (File fileEntry : AssignFiles.values()) {
                Logger.getRootLogger().info(fileEntry.getAbsolutePath());
            }
            
             for (File fileEntry : AssignFiles.values()) {
                String mzXMLFile = fileEntry.getAbsolutePath();
                if (mzXMLFile.toLowerCase().endsWith(".mzxml") | mzXMLFile.toLowerCase().endsWith(".mzml")) {
                    DIAPack DiaFile = new DIAPack(mzXMLFile, NoCPUs);
                    Logger.getRootLogger().info("=================================================================================================");
                    Logger.getRootLogger().info("Processing " + mzXMLFile);
                    if (!DiaFile.LoadDIASetting()) {
                        Logger.getRootLogger().info("Loading DIA setting failed, job is incomplete");
                        System.exit(1);
                    }
                    if (!DiaFile.LoadParams()) {
                        Logger.getRootLogger().info("Loading parameters failed, job is incomplete");
                        System.exit(1);
                    }
                    Logger.getRootLogger().info("Loading identification results " + mzXMLFile + "....");

                    //If the serialization file for ID file existed
                    if (DiaFile.ReadSerializedLCMSID()) {
                        DiaFile.IDsummary.ReduceMemoryUsage();
                        DiaFile.IDsummary.FastaPath=tandemPara.FastaPath;
                        FileList.add(DiaFile);
                    }
                }
            }
            
            //<editor-fold defaultstate="collapsed" desc="Targeted re-extraction using external library">
            
            //External library search
            
                Logger.getRootLogger().info("Targeted extraction using external library");
                
                //Read exteranl library
                FragmentLibManager ExlibManager = FragmentLibManager.ReadFragmentLibSerialization(WorkFolder, FilenameUtils.getBaseName(ExternalLibPath));
                if (ExlibManager == null) {
                    ExlibManager = new FragmentLibManager(FilenameUtils.getBaseName(ExternalLibPath));
                    
                    //Import traML file
                    ExlibManager.ImportFragLibByTraML(ExternalLibPath, ExternalLibDecoyTag);
                    //Check if there are decoy spectra
                    ExlibManager.CheckDecoys();
                    //ExlibManager.ImportFragLibBySPTXT(ExternalLibPath);
                    ExlibManager.WriteFragmentLibSerialization(WorkFolder);
                }
                Logger.getRootLogger().info("No. of peptide ions in external lib:" + ExlibManager.PeptideFragmentLib.size());
                for (DIAPack diafile : FileList) {
                    if (diafile.IDsummary == null) {
                        diafile.ReadSerializedLCMSID();
                    }
                    //Generate RT mapping
                    RTMappingExtLib RTmap = new RTMappingExtLib(diafile.IDsummary, ExlibManager, diafile.GetParameter());
                    RTmap.GenerateModel();
                    RTmap.GenerateMappedPepIon();
                    
                    diafile.BuildStructure();
                    diafile.MS1FeatureMap.ReadPeakCluster();
                    diafile.GenerateMassCalibrationRTMap();
                    //Perform targeted re-extraction
                    diafile.TargetedExtractionQuant(false, ExlibManager,ExtProbThreshold,RTWindow_Ext);
                    diafile.MS1FeatureMap.ClearAllPeaks();
                    diafile.IDsummary.ReduceMemoryUsage();
                    //Remove target IDs below the defined probability threshold
                    diafile.IDsummary.RemoveLowProbMappedIon(ExtProbThreshold);                    
                    diafile.ExportID();
                    diafile.ClearStructure();
                    Logger.getRootLogger().info("Peptide ions: " + diafile.IDsummary.GetPepIonList().size() + " Mapped ions: " + diafile.IDsummary.GetMappedPepIonList().size());
                }
            
            //</editor-fold>
            
            Logger.getRootLogger().info("Job done");
            Logger.getRootLogger().info("=================================================================================================");
        
    }
}
