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
import MSUmpire.DIA.RTAlignedPepIonMapping;
import MSUmpire.DIA.RTMappingExtLib;
import MSUmpire.FragmentLib.FragmentLibManager;
import MSUmpire.MSMSDBSearch.DBSearchParam;
import MSUmpire.MSMSDBSearch.TandemParam;
import MSUmpire.PSMDataStructure.FragmentPeak;
import MSUmpire.PSMDataStructure.LCMSID;
import MSUmpire.PSMDataStructure.PTMManager;
import MSUmpire.PSMDataStructure.ProtID;
import MSUmpire.PSMDataStructure.FragmentSelection;
import MSUmpire.Utility.ExportTable;
import MSUmpire.SearchResultParser.ProtXMLParser;
import MSUmpire.Utility.ConsoleLogger;
import MSUmpire.Utility.DateTimeTag;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.TimeUnit;
import java.util.zip.DataFormatException;
import javax.xml.parsers.ParserConfigurationException;
import org.apache.commons.io.FilenameUtils;
import org.apache.commons.lang.exception.ExceptionUtils;
import org.apache.log4j.Level;
import org.apache.log4j.Logger;
import org.xml.sax.SAXException;

/**
 *
 * @author Chih-Chiang Tsou
 */
public class DIA_Umpire_Quant {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) throws FileNotFoundException, IOException, Exception {
        System.out.println("=================================================================================================");
        System.out.println("DIA-Umpire quantitation with targeted re-extraction analysis (version: " + UmpireInfo.GetInstance().Version + ")");
        if (args.length != 1) {
            System.out.println("command format error, it should be like: java -jar -Xmx10G DIA_Umpire_Quant.jar diaumpire.quant_params");
            return;
        }
        try {
            ConsoleLogger.SetConsoleLogger(Level.DEBUG);
            ConsoleLogger.SetFileLogger(Level.DEBUG, FilenameUtils.getFullPath(args[0]) + "diaumpire_quant.log");
        } catch (Exception e) {
        }

        Logger.getRootLogger().info("Version: "+UmpireInfo.GetInstance().Version);
        Logger.getRootLogger().info("Parameter file:" + args[0]);

        BufferedReader reader = new BufferedReader(new FileReader(args[0]));
        String line = "";
        String WorkFolder = "";
        int NoCPUs = 2;

        String UserMod = "";
        String Combined_Prot = "";
        String InternalLibID = "";
        String ExternalLibPath = "";
        String ExternalLibDecoyTag = "DECOY";
        
        float ProbThreshold = 0.9f;
        float ReSearchProb =0.8f;
        float Freq = 0f;
        int TopNPep = 6;
        int TopNFrag = 6;
        String FilterWeight = "GW";
        float MinWeight = 0.9f;
        float RTWindow_Int=-1f;
        float RTWindow_Ext=-1f;
        
        TandemParam tandemPara = new TandemParam(DBSearchParam.SearchInstrumentType.TOF5600);
        HashMap<String, File> AssignFiles = new HashMap<>();
        boolean InternalLibSearch = false;
        boolean ExternalLibSearch = false;

        boolean ExportSaint = false;
        boolean SAINT_MS1 = false;
        boolean SAINT_MS2 = true;
//        String SaintBaitFile="";
//        String SaintPreyFile="";
//        String SaintInteractionFile="";
        HashMap<String, String[]> BaitList = new HashMap<>();
        HashMap<String, String> BaitName = new HashMap<>();
        HashMap<String, String[]> ControlList = new HashMap<>();
        HashMap<String, String> ControlName = new HashMap<>();

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
                    case "TargetedExtraction": {
                        InternalLibSearch = Boolean.parseBoolean(value);
                        break;
                    }
                    case "InternalLibSearch": {
                        InternalLibSearch = Boolean.parseBoolean(value);
                        break;
                    }
                    case "ExternalLibSearch": {
                        ExternalLibSearch = Boolean.parseBoolean(value);
                        break;
                    }

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
                    case "Combined_Prot": {
                        Combined_Prot = value;
                        break;
                    }
                    case "DecoyPrefix": {
                        if (!"".equals(value)) {
                            tandemPara.DecoyPrefix = value;
                        }
                        break;
                    }
                    case "UserMod": {
                        UserMod = value;
                        break;
                    }
                    case "ProteinFDR": {
                        tandemPara.ProtFDR = Float.parseFloat(value);
                        break;
                    }
                    case "PeptideFDR": {
                        tandemPara.PepFDR = Float.parseFloat(value);
                        break;
                    }
                    case "InternalLibID": {
                        InternalLibID = value;
                        break;
                    }
                    case "ExternalLibPath": {
                        ExternalLibPath = value;
                        break;
                    }
                     case "RTWindow_Int": {
                        RTWindow_Int = Float.parseFloat(value);
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
                    case "ProbThreshold": {
                        ProbThreshold = Float.parseFloat(value);
                        break;
                    }
                    case "ReSearchProb": {
                        ReSearchProb = Float.parseFloat(value);
                        break;
                    }                    
                    case "FilterWeight": {
                        FilterWeight = value;
                        break;
                    }
                    case "MinWeight": {
                        MinWeight = Float.parseFloat(value);
                        break;
                    }
                    case "TopNFrag": {
                        TopNFrag = Integer.parseInt(value);
                        break;
                    }
                    case "TopNPep": {
                        TopNPep = Integer.parseInt(value);
                        break;
                    }
                    case "Freq": {
                        Freq = Float.parseFloat(value);
                        break;
                    }

                    //<editor-fold defaultstate="collapsed" desc="SaintOutput">
                    case "ExportSaintInput": {
                        ExportSaint = Boolean.parseBoolean(value);
                        break;
                    }
                    case "QuantitationType": {
                        switch (value) {
                            case "MS1": {
                                SAINT_MS1 = true;
                                SAINT_MS2 = false;
                                break;
                            }
                            case "MS2": {
                                SAINT_MS1 = false;
                                SAINT_MS2 = true;
                                break;
                            }
                            case "BOTH": {
                                SAINT_MS1 = true;
                                SAINT_MS2 = true;
                                break;
                            }
                        }
                        break;
                    }
//                    case "BaitInputFile": {
//                        SaintBaitFile = value;
//                        break;
//                    }
//                    case "PreyInputFile": {
//                        SaintPreyFile = value;
//                        break;
//                    }
//                    case "InterationInputFile": {
//                        SaintInteractionFile = value;
//                        break;
//                    }
                    default: {
                        if (type.startsWith("BaitName_")) {
                            BaitName.put(type.substring(9), value);
                        }
                        if (type.startsWith("BaitFile_")) {
                            BaitList.put(type.substring(9), value.split("\t"));
                        }
                        if (type.startsWith("ControlName_")) {
                            ControlName.put(type.substring(12), value);
                        }
                        if (type.startsWith("ControlFile_")) {
                            ControlList.put(type.substring(12), value.split("\t"));
                        }
                        break;
                    }
//</editor-fold>                    
                }
            }
        }
//</editor-fold>

        PTMManager.GetInstance();
        if (!UserMod.equals("")) {
            PTMManager.GetInstance().ImportUserMod(UserMod);
        }

        if (!new File(tandemPara.FastaPath).exists()) {
            Logger.getRootLogger().info("Fasta file :"+tandemPara.FastaPath + " cannot be found, the process will be terminated, please check.");
            System.exit(1);
        }               
        
        if (!new File(Combined_Prot).exists()) {
            Logger.getRootLogger().info("ProtXML file: " +Combined_Prot + " cannot be found, the export protein summary table will be empty.");
        }
        LCMSID protID = null;

        if (Combined_Prot != null && !Combined_Prot.equals("")) {
            protID = LCMSID.ReadLCMSIDSerialization(Combined_Prot);
            if (!"".equals(Combined_Prot) && protID == null) {
                protID = new LCMSID(Combined_Prot, tandemPara.DecoyPrefix, tandemPara.FastaPath);
                ProtXMLParser protxmlparser = new ProtXMLParser(protID, Combined_Prot, 0f);
                protID.RemoveLowLocalPWProtein(0.8f);
                protID.RemoveLowMaxIniProbProtein(0.9f);
                protID.FilterByProteinDecoyFDRUsingMaxIniProb(tandemPara.DecoyPrefix, tandemPara.ProtFDR);
                protID.LoadSequence();
                protID.WriteLCMSIDSerialization(Combined_Prot);
            }
            Logger.getRootLogger().info("Protein No.:" + protID.ProteinList.size());
        }
        HashMap<String, HashMap<String, FragmentPeak>> IDSummaryFragments = new HashMap<>();

        ArrayList<DIAPack> FileList = new ArrayList<>();
        try {
            File folder = new File(WorkFolder);
            for (final File fileEntry : folder.listFiles()) {
                if (fileEntry.isFile() && fileEntry.getAbsolutePath().toLowerCase().endsWith(".mzxml")
                        && !fileEntry.getAbsolutePath().toLowerCase().endsWith("q1.mzxml")
                        && !fileEntry.getAbsolutePath().toLowerCase().endsWith("q2.mzxml")
                        && !fileEntry.getAbsolutePath().toLowerCase().endsWith("q3.mzxml")) {
                    AssignFiles.put(fileEntry.getAbsolutePath(), fileEntry);
                }
                if (fileEntry.isDirectory()) {
                    for (final File fileEntry2 : fileEntry.listFiles()) {
                        if (fileEntry2.isFile() && fileEntry2.getAbsolutePath().toLowerCase().endsWith(".mzxml")
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
                ProcessDIA(fileEntry, NoCPUs, tandemPara, FileList, IDSummaryFragments, protID);
            }

            Logger.getRootLogger().info("=================================================================================================");
            if (InternalLibSearch && FileList.size() > 1) {
                Logger.getRootLogger().info("Module C: Targeted extraction (internal library search)");

                FragmentLibManager libManager = FragmentLibManager.ReadFragmentLibSerialization(WorkFolder, InternalLibID);
                if (libManager == null) {
                    Logger.getRootLogger().info("Building internal spectral library");
                    libManager = new FragmentLibManager(InternalLibID);
                    ArrayList<LCMSID> LCMSIDList = new ArrayList<>();
                    for (DIAPack dia : FileList) {
                        LCMSIDList.add(dia.IDsummary);
                    }
                    //libManager.ImportFragLib(LCMSIDList);
                    libManager.ImportFragLibTopFrag(LCMSIDList, Freq, TopNFrag);
                    libManager.WriteFragmentLibSerialization(WorkFolder);
                }
                libManager.ReduceMemoryUsage();

                Logger.getRootLogger().info("Building retention time prediction model and generate candiate peptide list");
                for (int i = 0; i < FileList.size(); i++) {
                    FileList.get(i).IDsummary.ClearMappedPep();
                }
                for (int i = 0; i < FileList.size(); i++) {
                    for (int j = i + 1; j < FileList.size(); j++) {
                        RTAlignedPepIonMapping alignment = new RTAlignedPepIonMapping(WorkFolder, FileList.get(i).GetParameter(), FileList.get(i).IDsummary, FileList.get(j).IDsummary);
                        alignment.GenerateModel();
                        alignment.GenerateMappedPepIon();
                        //alignment.ExportMappedPepIon(connectionManager);
                    }
                    FileList.get(i).ExportID();
                    FileList.get(i).IDsummary = null;
                }

                Logger.getRootLogger().info("Targeted matching........");
                for (DIAPack diafile : FileList) {
                    if (diafile.IDsummary == null) {
                        diafile.ReadSerializedLCMSID();
                    }
                    if (!diafile.IDsummary.GetMappedPepIonList().isEmpty()) {
                        diafile.UseMappedIon = true;
                        diafile.FilterMappedIonByProb = false;
                        diafile.BuildStructure();
                        diafile.ms1lcms.ReadPeakCluster();
                        diafile.ms1lcms.ClearMonoisotopicPeakOfCluster();
                        diafile.GenerateMassCalibrationRTMap();
                        diafile.AssignMappedPepQuant(false, libManager,1.1f,RTWindow_Int);
                        diafile.ms1lcms.ClearAllPeaks();
                        diafile.IDsummary.ReduceMemoryUsage();
                        diafile.IDsummary.RemoveLowProbMappedIon(ProbThreshold);
                        if (protID != null) {
                            diafile.IDsummary.GenerateProteinByRefIDByPepSeq(protID, true);
                            diafile.IDsummary.ReMapProPep();
                        }
                        diafile.ExportID();
                        Logger.getRootLogger().info("Peptide ions: " + diafile.IDsummary.GetPepIonList().size() + " Mapped ions: " + diafile.IDsummary.GetMappedPepIonList().size());
                        diafile.ClearStructure();
                    }
                    diafile.IDsummary = null;
                    System.gc();
                }
                Logger.getRootLogger().info("=================================================================================================");
            }

            if (ExternalLibSearch) {
                Logger.getRootLogger().info("Module C: Targeted extraction (external library search)");
                FragmentLibManager ExlibManager = FragmentLibManager.ReadFragmentLibSerialization(WorkFolder, FilenameUtils.getBaseName(ExternalLibPath));
                if (ExlibManager == null) {
                    ExlibManager = new FragmentLibManager(FilenameUtils.getBaseName(ExternalLibPath));
                    ExlibManager.ImportFragLibByTraML(ExternalLibPath, ExternalLibDecoyTag);
                    ExlibManager.CheckDecoys();
                    //ExlibManager.ImportFragLibBySPTXT(ExternalLibPath);
                    ExlibManager.WriteFragmentLibSerialization(WorkFolder);
                }
                Logger.getRootLogger().info("No. of peptide ions in external lib:" + ExlibManager.PeptideFragmentLib.size());
                for (DIAPack diafile : FileList) {
                    if (diafile.IDsummary == null) {
                        diafile.ReadSerializedLCMSID();
                    }
                    RTMappingExtLib RTmap = new RTMappingExtLib(diafile.IDsummary, ExlibManager, diafile.GetParameter());
                    RTmap.GenerateModel();
                    RTmap.GenerateMappedPepIon();

                    diafile.BuildStructure();
                    diafile.ms1lcms.ReadPeakCluster();
                    diafile.GenerateMassCalibrationRTMap();
                    diafile.AssignMappedPepQuant(false, ExlibManager,ReSearchProb,RTWindow_Ext);
                    diafile.ms1lcms.ClearAllPeaks();
                    diafile.IDsummary.ReduceMemoryUsage();
                    diafile.IDsummary.RemoveLowProbMappedIon(ProbThreshold);
                    if (protID != null) {
                        diafile.IDsummary.GenerateProteinByRefIDByPepSeq(protID, true);
                        diafile.IDsummary.ReMapProPep();
                    }
                    diafile.ExportID();
                    diafile.ClearStructure();
                    Logger.getRootLogger().info("Peptide ions: " + diafile.IDsummary.GetPepIonList().size() + " Mapped ions: " + diafile.IDsummary.GetMappedPepIonList().size());
                }
            }
            Logger.getRootLogger().info("Peptide and fragment selection across the whole dataset");
            ArrayList<LCMSID> SummaryList = new ArrayList<>();
            for (DIAPack diafile : FileList) {
                if (diafile.IDsummary == null) {
                    diafile.ReadSerializedLCMSID();
                    ////////////////////////////
                    diafile.IDsummary.GenerateProteinByRefIDByPepSeq(protID, true);
                    diafile.IDsummary.ReMapProPep();
                    /////////////////////////////
                    diafile.IDsummary.ClearAssignPeakCluster();
                    //diafile.IDsummary.ClearPSMs();
                }
                if ("GW".equals(FilterWeight)) {
                    diafile.IDsummary.SetFilterByGroupWeight();
                } else if ("PepW".equals(FilterWeight)) {
                    diafile.IDsummary.SetFilterByWeight();
                }
                SummaryList.add(diafile.IDsummary);
            }
            FragmentSelection fragselection = new FragmentSelection(SummaryList);
            fragselection.freqPercent = Freq;
            fragselection.GeneratePepFragScoreMap();
            fragselection.GenerateTopFragMap(TopNFrag);
            fragselection.GenerateProtPepScoreMap(MinWeight);
            fragselection.GenerateTopPepMap(TopNPep);

            //<editor-fold defaultstate="collapsed" desc="Writing general reports">                 
            ExportTable export = new ExportTable(WorkFolder, SummaryList, IDSummaryFragments, protID, fragselection);
            export.Export(TopNPep, TopNFrag, Freq);
            //</editor-fold>

            //<editor-fold defaultstate="collapsed" desc="//<editor-fold defaultstate="collapsed" desc="Generate SAINT input files">
            if (ExportSaint && protID != null) {
                HashMap<String, DIAPack> Filemap = new HashMap<>();
                for (DIAPack DIAfile : FileList) {
                    Filemap.put(DIAfile.GetBaseName(), DIAfile);
                }

                FileWriter baitfile = new FileWriter(WorkFolder + "SAINT_Bait_" + DateTimeTag.GetTag() + ".txt");
                FileWriter preyfile = new FileWriter(WorkFolder + "SAINT_Prey_" + DateTimeTag.GetTag() + ".txt");
                FileWriter interactionfileMS1 = null;
                FileWriter interactionfileMS2 = null;
                if (SAINT_MS1) {
                    interactionfileMS1 = new FileWriter(WorkFolder + "SAINT_Interaction_MS1_" + DateTimeTag.GetTag() + ".txt");
                }
                if (SAINT_MS2) {
                    interactionfileMS2 = new FileWriter(WorkFolder + "SAINT_Interaction_MS2_" + DateTimeTag.GetTag() + ".txt");
                }
                HashMap<String, String> PreyID = new HashMap<>();

                for (String samplekey : ControlName.keySet()) {
                    String name = ControlName.get(samplekey);
                    for (String file : ControlList.get(samplekey)) {
                        baitfile.write(FilenameUtils.getBaseName(file) + "\t" + name + "\t" + "C\n");
                        LCMSID IDsummary = Filemap.get(FilenameUtils.getBaseName(file)).IDsummary;
                        if (SAINT_MS1) {
                            SaintOutput(protID, IDsummary, fragselection, interactionfileMS1, file, name, PreyID, 1);
                        }
                        if (SAINT_MS2) {
                            SaintOutput(protID, IDsummary, fragselection, interactionfileMS2, file, name, PreyID, 2);
                        }
                    }
                }
                for (String samplekey : BaitName.keySet()) {
                    String name = BaitName.get(samplekey);
                    for (String file : BaitList.get(samplekey)) {
                        baitfile.write(FilenameUtils.getBaseName(file) + "\t" + name + "\t" + "T\n");
                        LCMSID IDsummary = Filemap.get(FilenameUtils.getBaseName(file)).IDsummary;
                        if (SAINT_MS1) {
                            SaintOutput(protID, IDsummary, fragselection, interactionfileMS1, file, name, PreyID, 1);
                        }
                        if (SAINT_MS2) {
                            SaintOutput(protID, IDsummary, fragselection, interactionfileMS2, file, name, PreyID, 2);
                        }
                    }
                }
                baitfile.close();
                if (SAINT_MS1) {
                    interactionfileMS1.close();
                }
                if (SAINT_MS2) {
                    interactionfileMS2.close();
                }
                for (String AccNo : PreyID.keySet()) {
                    preyfile.write(AccNo + "\t" + PreyID.get(AccNo) + "\n");
                }
                preyfile.close();
            }

//</editor-fold>
            
            Logger.getRootLogger().info("Job done");
            Logger.getRootLogger().info("=================================================================================================");

        } catch (Exception e) {
            Logger.getRootLogger().error(ExceptionUtils.getStackTrace(e));
            throw e;
        }
    }

    private static void SaintOutput(LCMSID protID, LCMSID IDsummary, FragmentSelection fragselection, FileWriter interactionfile, String filename, String samplename, HashMap<String, String> PreyID, int quanttype) throws IOException {
        for (String key : protID.ProteinList.keySet()) {
            if (IDsummary.ProteinList.containsKey(key)) {
                ProtID protein = IDsummary.ProteinList.get(key);
                float abundance = 0f;

                if (quanttype == 1) {
                    abundance = protein.GetAbundanceByMS1_IBAQ();
                } else if (quanttype == 2) {
                    abundance = protein.GetAbundanceByTopCorrFragAcrossSample(fragselection.TopPeps.get(protein.getAccNo()), fragselection.TopFrags);
                }
                if (abundance > 0) {
                    interactionfile.write(FilenameUtils.getBaseName(filename) + "\t" + samplename + "\t" + protein.getAccNo() + "\t" + abundance + "\n");
                    if (!PreyID.containsKey(protein.getAccNo())) {
                        PreyID.put(protein.getAccNo(), /*protein.Sequence.length()+"\t"+*/ protein.GetGeneName());
                    }
                }
            }
        }
    }

    private static void ProcessDIA(final File fileEntry, int NoCPUs, TandemParam tandemPara, ArrayList<DIAPack> FileList, HashMap<String, HashMap<String, FragmentPeak>> IDSummaryFragments, LCMSID protID) throws IOException, FileNotFoundException, DataFormatException, InterruptedException, ExecutionException, ParserConfigurationException, SAXException, Exception {
        String mzXMLFile = fileEntry.getAbsolutePath();
        if (mzXMLFile.toLowerCase().endsWith(".mzxml")) {
            long time = System.currentTimeMillis();

            DIAPack DiaFile = new DIAPack(mzXMLFile, NoCPUs);
            FileList.add(DiaFile);
            HashMap<String, FragmentPeak> FragMap = new HashMap<>();
            IDSummaryFragments.put(FilenameUtils.getBaseName(mzXMLFile), FragMap);
            if (!new File(FilenameUtils.getFullPath(DiaFile.Filename) + DiaFile.GetQ1Name() + ".mzXML").exists()
                    | !new File(FilenameUtils.getFullPath(DiaFile.Filename) + DiaFile.GetQ2Name() + ".mzXML").exists()
                    | !new File(FilenameUtils.getFullPath(DiaFile.Filename) + DiaFile.GetQ3Name() + ".mzXML").exists()) {
                return;
            }
            Logger.getRootLogger().info("=================================================================================================");
            Logger.getRootLogger().info("Processing " + mzXMLFile);
            if (!DiaFile.LoadDIASetting()) {
//                DiaFile.dIA_Setting=new DIA_Setting();
//                DiaFile.LoadParams();
//                DiaFile.SetDataType(SpectralDataType.DataType.DIA_F_Window);
//                DiaFile.SetWindowSize(25);
//                DiaFile.GetSpectrumParser();
                Logger.getRootLogger().info("Loading DIA setting failed, job is incomplete");
                System.exit(1);
            }
            if (!DiaFile.LoadParams()) {
                Logger.getRootLogger().info("Loading parameters failed, job is incomplete");
                System.exit(1);
            }
            Logger.getRootLogger().info("Loading identification results " + mzXMLFile + "....");

            if (!DiaFile.ReadSerializedLCMSID()) {
                DiaFile.ParsePepXML(tandemPara);
                DiaFile.BuildStructure();
                if (!DiaFile.ms1lcms.ReadPeakCluster()) {
                    Logger.getRootLogger().info("Loading peak and structure failed, job is incomplete");
                    System.exit(1);
                }
                DiaFile.IDsummary.GenerateProteinByRefIDByPepSeq(protID, true);
                DiaFile.IDsummary.ReMapProPep();
                DiaFile.ms1lcms.ClearMonoisotopicPeakOfCluster();
                DiaFile.GenerateClusterScanNomapping();
                DiaFile.AssignQuant();
                DiaFile.ClearStructure();
            }            
            DiaFile.IDsummary.ReduceMemoryUsage();   
            time = System.currentTimeMillis() - time;
            Logger.getRootLogger().info(mzXMLFile + " processed time:" + String.format("%d hour, %d min, %d sec", TimeUnit.MILLISECONDS.toHours(time), TimeUnit.MILLISECONDS.toMinutes(time) - TimeUnit.HOURS.toMinutes(TimeUnit.MILLISECONDS.toHours(time)), TimeUnit.MILLISECONDS.toSeconds(time) - TimeUnit.MINUTES.toSeconds(TimeUnit.MILLISECONDS.toMinutes(time))));
        }
    }

}
