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

import MSUmpire.BaseDataStructure.UmpireInfo;
import MSUmpire.DIA.DIAPack;
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
import java.util.concurrent.TimeUnit;
import org.apache.commons.io.FilenameUtils;
import org.apache.log4j.Level;
import org.apache.log4j.Logger;

/**
 *
 * @author Chih-Chiang Tsou
 */
public class DIA_Umpire_LCMSIDGen {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) throws FileNotFoundException, IOException, Exception {
        System.out.println("=================================================================================================");
        System.out.println("DIA-Umpire LCMSID geneartor (version: " + UmpireInfo.GetInstance().Version + ")");
        if (args.length != 1) {
            System.out.println("command format error, the correct format should be: java -jar -Xmx10G DIA_Umpire_LCMSIDGen.jar diaumpire_module.params");
            return;
        }
        try {
            ConsoleLogger.SetConsoleLogger(Level.INFO);
            ConsoleLogger.SetFileLogger(Level.DEBUG, FilenameUtils.getFullPath(args[0]) + "diaumpire_lcmsidgen.log");
        } catch (Exception e) {
        }

        Logger.getRootLogger().info("Version: " + UmpireInfo.GetInstance().Version);
        Logger.getRootLogger().info("Parameter file:" + args[0]);

        BufferedReader reader = new BufferedReader(new FileReader(args[0]));
        String line = "";
        String WorkFolder = "";
        int NoCPUs = 2;

        TandemParam tandemPara = new TandemParam(DBSearchParam.SearchInstrumentType.TOF5600);
        HashMap<String, File> AssignFiles = new HashMap<>();

        //<editor-fold defaultstate="collapsed" desc="Reading parameter file">
        while ((line = reader.readLine()) != null) {
            line = line.trim();
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
                    case "DecoyPrefix": {
                        if (!"".equals(value)) {
                            tandemPara.DecoyPrefix = value;
                        }
                        break;
                    }
                    case "PeptideFDR": {
                        tandemPara.PepFDR = Float.parseFloat(value);
                        break;
                    }                                        
                }
            }
        }
//</editor-fold>

        //Initialize PTM manager using compomics library
        PTMManager.GetInstance();
        
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

        //process each DIA file to genearate untargeted identifications
        for (File fileEntry : AssignFiles.values()) {
            String mzXMLFile = fileEntry.getAbsolutePath();
            if (mzXMLFile.toLowerCase().endsWith(".mzxml") | mzXMLFile.toLowerCase().endsWith(".mzml")) {
                long time = System.currentTimeMillis();

                DIAPack DiaFile = new DIAPack(mzXMLFile, NoCPUs);
                FileList.add(DiaFile);
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

                DiaFile.ParsePepXML(tandemPara);
                DiaFile.BuildStructure();
                if (!DiaFile.MS1FeatureMap.ReadPeakCluster()) {
                    Logger.getRootLogger().info("Loading peak and structure failed, job is incomplete");
                    System.exit(1);
                }
                DiaFile.MS1FeatureMap.ClearMonoisotopicPeakOfCluster();
                //Generate mapping between index of precursor feature and pseudo MS/MS scan index 
                DiaFile.GenerateClusterScanNomapping();
                //Doing quantification
                DiaFile.AssignQuant();
                DiaFile.ClearStructure();

                DiaFile.IDsummary.ReduceMemoryUsage();
                time = System.currentTimeMillis() - time;
                Logger.getRootLogger().info(mzXMLFile + " processed time:" + String.format("%d hour, %d min, %d sec", TimeUnit.MILLISECONDS.toHours(time), TimeUnit.MILLISECONDS.toMinutes(time) - TimeUnit.HOURS.toMinutes(TimeUnit.MILLISECONDS.toHours(time)), TimeUnit.MILLISECONDS.toSeconds(time) - TimeUnit.MINUTES.toSeconds(TimeUnit.MILLISECONDS.toMinutes(time))));
            }
            Logger.getRootLogger().info("Job done");
            Logger.getRootLogger().info("=================================================================================================");
        }
    }
}
