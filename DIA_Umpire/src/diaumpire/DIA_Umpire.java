/* 
 * Copyright 2014 Chih-Chiang Tsou <chihchiang.tsou@gmail.com>.
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

import MSUmpire.BaseDataStructure.InstrumentParameter;
import MSUmpire.BaseDataStructure.SpectralDataType;
import MSUmpire.BaseDataStructure.UmpireInfo;
import MSUmpire.BaseDataStructure.XYData;
import MSUmpire.DIA.DIAPack;
import MSUmpire.DIA.FDR_DataSetLevel;
import MSUmpire.FragmentLib.FragmentLibManager;
import MSUmpire.MSMSDBSearch.CometParam;
import MSUmpire.MSMSDBSearch.CometSearch;
import MSUmpire.MSMSDBSearch.DBSearchParam;
import MSUmpire.MSMSDBSearch.MSGFSearch;
import MSUmpire.MSMSDBSearch.MSGFParam;
import MSUmpire.MSMSDBSearch.MSMSDBSearch;
import MSUmpire.MSMSDBSearch.ProteinProphetCombine;
import MSUmpire.MSMSDBSearch.TandemParam;
import MSUmpire.MSMSDBSearch.TandemSearch;
import MSUmpire.PSMDataStructure.EnzymeManager;
import MSUmpire.PSMDataStructure.FragmentPeak;
import MSUmpire.PSMDataStructure.FragmentSelection;
import MSUmpire.PSMDataStructure.LCMSID;
import MSUmpire.PSMDataStructure.PTMManager;
import MSUmpire.PostProcessing.ExportTable;
import MSUmpire.QuantModule.RTAlignedPepIonMapping;
import MSUmpire.QuantModule.RTMappingExtLib;
import MSUmpire.SearchResultParser.ProtXMLParser;
import Utility.ConsoleLogger;
import java.io.*;;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;
import java.util.zip.DataFormatException;
import javax.xml.parsers.ParserConfigurationException;
import org.apache.commons.io.FilenameUtils;
import org.apache.commons.lang.exception.ExceptionUtils;
import org.apache.log4j.Level;
import org.apache.log4j.Logger;
import org.xml.sax.SAXException;
import org.xmlpull.v1.XmlPullParserException;

/**
 *
 * @author Chih-Chiang Tsou <chihchiang.tsou@gmail.com>
 */
public class DIA_Umpire {

    /**
     * @param args the command line arguments DIA_Umpire parameterfile
     */
    public static void main(String[] args) throws InterruptedException, FileNotFoundException, ExecutionException, IOException, ParserConfigurationException, DataFormatException, SAXException, XmlPullParserException, Exception {
        System.out.println("=================================================================================================");
        System.out.println("DIA-Umpire analysis (version: "+UmpireInfo.GetInstance().Version+")");
        if (args.length != 1) {
            System.out.println("command: java -jar -Xmx10G DIA_Umpire.jar diaumpire.paras");
            return;
        }
        ConsoleLogger consoleLogger = new ConsoleLogger();
        consoleLogger.SetConsoleLogger(Level.DEBUG);
        consoleLogger.SetFileLogger(Level.DEBUG, FilenameUtils.getFullPath(args[0]) + "diaumpire_debug.log");
        Logger logger = Logger.getRootLogger();
        Logger.getRootLogger().info("Version: "+UmpireInfo.GetInstance().Version);
        logger.info("Parameter file:" + args[0]);
        
        BufferedReader reader = new BufferedReader(new FileReader(args[0]));
        String line = "";
        String WorkFolder = "";

        InstrumentParameter para = new InstrumentParameter(InstrumentParameter.InstrumentType.TOF5600);
        para.DetermineBGByID = false;
        para.EstimateBG = true;
        int NoCPUs = 2;
        int NoOfConcurrentFile = 2;
        int NoThreadPerFile = 1;

        boolean SignalExtraction = false;
        boolean DBSerach = false;
        boolean GeneratePBS=false;
        boolean Quantitation = false;
        boolean TargetedExtraction = false;
        boolean ExternalLibSearch=false;
        boolean BuildPeptideCandiate = false;
        boolean ExportTable = false;

        SpectralDataType.DataType dataType = SpectralDataType.DataType.DIA_F_Window;
        String SearchEngine = "";
        String MGFtag = "";
        String UserMod = "";
        String InternalLibID = "";
        String ExternalLibPath = "";
        String ExternalLibDecoyTag = "DECOY";
        float ProbThreshold = 0.9f;
        float ReSearchProb =0.8f;
        String CombineProtXML = "";
        String WindowType = "";
        String tag="";
        float WindowSize = 25;
        //ConnectionManager connectionManager = null;
        TandemParam tandemparam = new TandemParam(DBSearchParam.SearchInstrumentType.TOF5600);
        TandemSearch tandemsearch = new TandemSearch(tandemparam);
        MSGFParam msgfparam = new MSGFParam(DBSearchParam.SearchInstrumentType.TOF5600);
        MSGFSearch msgfsearch = new MSGFSearch(msgfparam);
        CometParam cometparam = new CometParam(DBSearchParam.SearchInstrumentType.TOF5600);
        CometSearch cometsearch = new CometSearch(cometparam);
        ArrayList<MSMSDBSearch> SearchEngines = new ArrayList<>();
        ArrayList<XYData> WindowList = new ArrayList<>();
        ArrayList<XYData> MS1WindowList = new ArrayList<>();

        boolean ExportPrecursorPeak = false;
        boolean ExportFragmentPeak = false;
        boolean DataSetLevelPepFDR=false;
        boolean Mysql = false;
        String Address = "";
        String Port = "";
        String DBName = "";
        String user = "";
        String pwd = "";

        //<editor-fold defaultstate="collapsed" desc="Read parameter file">
        while ((line = reader.readLine()) != null) {
            Logger.getRootLogger().info(line);
            if (!"".equals(line) && !line.startsWith("#")) {
                //System.out.println(line);
                if (line.equals("==window setting begin")) {
                    while (!(line = reader.readLine()).equals("==window setting end")) {
                        if (!"".equals(line)) {
                            WindowList.add(new XYData(Float.parseFloat(line.split("\t")[0]), Float.parseFloat(line.split("\t")[1])));
                        }
                    }
                    continue;
                }
                if (line.equals("==MS1 window setting begin")) {
                    while (!(line = reader.readLine()).equals("==MS1 window setting end")) {
                        if (!"".equals(line)) {
                            MS1WindowList.add(new XYData(Float.parseFloat(line.split("\t")[0]), Float.parseFloat(line.split("\t")[1])));
                        }
                    }
                    continue;
                }
                if (line.split("=").length < 2) {
                    continue;
                }
                String type = line.split("=")[0].trim();
                String value = line.split("=")[1].trim();
                
                if (type.startsWith("para.")) {
                    type = type.replace("para.", "SE.");
                }
                switch (type) {
                    case "A-SignalExtraction": {
                        SignalExtraction = Boolean.parseBoolean(value);
                        break;
                    }
                    case "B-DBSerach": {
                        DBSerach = Boolean.parseBoolean(value);
                        break;
                    }
                    case "B2-GeneratePBS": {
                        GeneratePBS = Boolean.parseBoolean(value);
                        break;
                    }                    
                    case "C-TargetedExtraction": {
                        TargetedExtraction = Boolean.parseBoolean(value);
                        break;
                    }
                    case "C2-ExternalLibSearch": {
                        ExternalLibSearch = Boolean.parseBoolean(value);
                        break;
                    }                    
                    case "D-Quantitation": {
                        Quantitation = Boolean.parseBoolean(value);
                        break;
                    }
                    case "E-Export": {
                        ExportTable = Boolean.parseBoolean(value);
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
                    
                    //<editor-fold defaultstate="collapsed" desc="MySQL">
                    case "Mysql": {
                        Mysql = Boolean.valueOf(value);
                        break;
                    }
                    case "Mysql.Address": {
                        Address = value;
                        break;
                    }
                    case "Mysql.Port": {
                        Port = value;
                        break;
                    }
                    case "Mysql.DBName": {
                        DBName = value;
                        break;
                    }
                    case "Mysql.User": {
                        user = value;
                        break;
                    }
                    case "Mysql.Pwd": {
                        pwd = value;
                        break;
                    }
//</editor-fold>
                    
                    case "ExportPrecursorPeak": {
                        ExportPrecursorPeak = Boolean.parseBoolean(value);
                        break;
                    }
                    case "ExportFragmentPeak": {
                        ExportFragmentPeak = Boolean.parseBoolean(value);
                        break;
                    }

                    //<editor-fold defaultstate="collapsed" desc="instrument parameters">
                    case "PrecursorRank": {
                        para.PrecursorRank = Integer.parseInt(value);
                        break;
                    }
                    case "FragmentRank": {
                        para.FragmentRank = Integer.parseInt(value);
                        break;
                    }
                    case "CorrThreshold": {
                        para.CorrThreshold = Float.parseFloat(value);
                        break;
                    }
                    case "DeltaApex": {
                        para.ApexDelta = Float.parseFloat(value);
                        break;
                    }
                    case "SE.MS1PPM": {
                        para.MS1PPM = Float.parseFloat(value);
                        break;
                    }
                    case "SE.MS2PPM": {
                        para.MS2PPM = Float.parseFloat(value);
                        break;
                    }
                    case "SE.SN": {
                        para.SNThreshold = Float.parseFloat(value);
                        break;
                    }
                    case "SE.MS2SN": {
                        para.MS2SNThreshold = Float.parseFloat(value);
                        break;
                    }
                    case "SE.MinMSIntensity": {
                        para.MinMSIntensity = Float.parseFloat(value);
                        break;
                    }
                    case "SE.MinMSMSIntensity": {
                        para.MinMSMSIntensity = Float.parseFloat(value);
                        break;
                    }
                    case "SE.MinRTRange": {
                        para.MinRTRange = Float.parseFloat(value);
                        break;
                    }
                    case "SE.MaxNoPeakCluster": {
                        para.MaxNoPeakCluster = Integer.parseInt(value);
                        break;
                    }
                    case "SE.MinNoPeakCluster": {
                        para.MinNoPeakCluster = Integer.parseInt(value);
                        break;
                    }
                    case "SE.MinMS2NoPeakCluster": {
                        para.MinMS2NoPeakCluster = Integer.parseInt(value);
                        break;
                    }
                    case "SE.MaxCurveRTRange": {
                        para.MaxCurveRTRange = Float.parseFloat(value);
                        break;
                    }
                    case "SE.Resolution": {
                        para.Resolution = Integer.parseInt(value);
                        break;
                    }
                    case "SE.RTtol": {
                        para.RTtol = Float.parseFloat(value);
                        break;
                    }
                    case "SE.NoPeakPerMin": {
                        para.NoPeakPerMin = Integer.parseInt(value);
                        break;
                    }
                    case "SE.StartCharge": {
                        para.StartCharge = Integer.parseInt(value);
                        break;
                    }
                    case "SE.EndCharge": {
                        para.EndCharge = Integer.parseInt(value);
                        break;
                    }
                    case "SE.MS2StartCharge": {
                        para.MS2StartCharge = Integer.parseInt(value);
                        break;
                    }
                    case "SE.MS2EndCharge": {
                        para.MS2EndCharge = Integer.parseInt(value);
                        break;
                    }
                    case "SE.NoMissedScan": {
                        para.NoMissedScan = Integer.parseInt(value);
                        break;
                    }
                    case "SE.Denoise": {
                        para.Denoise = Boolean.valueOf(value);
                        break;
                    }
                    case "SE.EstimateBG": {
                        para.EstimateBG = Boolean.valueOf(value);
                        break;
                    }
                    case "SE.RemoveGroupedPeaks": {
                        para.RemoveGroupedPeaks = Boolean.valueOf(value);
                        break;
                    }
                    case "SE.RemoveGroupedPeaksRTOverlap": {
                        para.RemoveGroupedPeaksRTOverlap = Float.valueOf(value);
                        break;
                    }
                    case "SE.RemoveGroupedPeaksCorr": {
                        para.RemoveGroupedPeaksCorr = Float.valueOf(value);
                        break;
                    }
                    case "SE.MinMZ": {
                        para.MinMZ = Float.valueOf(value);
                        break;
                    }
                    case "SE.IsoCorrThreshold": {
                        para.IsoCorrThreshold = Float.valueOf(value);
                        break;
                    }
                    case "SE.IsoPattern": {
                        para.IsoPattern = Float.valueOf(value);
                        break;
                    }                    
                    case "SE.StartRT":{
                        para.startRT = Float.valueOf(value);
                    }
                    case "SE.EndRT":{
                        para.endRT = Float.valueOf(value);
                    }
                    case "SE.MassDefectFilter":{
                        para.MassDefectFilter = Boolean.parseBoolean(value);
                    }

//</editor-fold>
                    
                    //<editor-fold defaultstate="collapsed" desc="General search parameters">
                    case "SearchEngine": {
                        SearchEngine = value;
                        for (int i = 0; i < SearchEngine.split(",").length; i++) {
                            switch (SearchEngine.split(",")[i]) {
                                case "Tandem": {
                                    SearchEngines.add(tandemsearch);
                                    break;
                                }
                                case "Comet": {
                                    SearchEngines.add(cometsearch);
                                }
                                break;
                                case "MSGF": {
                                    SearchEngines.add(msgfsearch);
                                }
                                break;
                            }
                        }
                        break;
                    }
                    case "PrecursorPPM": {
                        tandemparam.PrecursorPPM = Float.parseFloat(value);
                        msgfparam.PrecursorPPM = Float.parseFloat(value);
                        cometparam.PrecursorPPM = Float.parseFloat(value);
                        break;
                    }
                    case "MissCleavage": {
                        tandemparam.MissCleavage = Integer.parseInt(value);
                        msgfparam.MissCleavage = Integer.parseInt(value);
                        cometparam.MissCleavage = Integer.parseInt(value);
                        break;
                    }
                    case "IsotopeError": {
                        tandemparam.IsotopeError = Boolean.parseBoolean(value);
                        msgfparam.IsotopeError = Boolean.parseBoolean(value);
                        cometparam.IsotopeError = Boolean.parseBoolean(value);
                        break;
                    }
                    case "FragPPM": {
                        tandemparam.FragPPM = Float.parseFloat(value);
                        msgfparam.FragPPM = Float.parseFloat(value);
                        cometparam.FragPPM = Float.parseFloat(value);
                        break;
                    }
                    case "NonSpecificCleavage": {
                        tandemparam.NonSpecificCleavage = Boolean.parseBoolean(value);
                        msgfparam.NonSpecificCleavage = Boolean.parseBoolean(value);
                        cometparam.NonSpecificCleavage = Boolean.parseBoolean(value);
                        break;
                    }
                    case "MGFtag": {
                        MGFtag = value;
                        break;
                    }
                    case "NoOfConcurrentFile": {
                        NoOfConcurrentFile = Integer.parseInt(value);
                        break;
                    }
                    case "SearchFileThread": {
                        NoOfConcurrentFile = Integer.parseInt(value);
                        break;
                    }
                    case "NoThreadPerFile": {
                        tandemparam.NoCPUs = Integer.parseInt(value);
                        cometparam.NoCPUs = Integer.parseInt(value);
                        msgfparam.NoCPUs = Integer.parseInt(value);
                        NoThreadPerFile = Integer.parseInt(value);
                        break;
                    }
                     case "SearchThreadPerFile": {
                        tandemparam.NoCPUs = Integer.parseInt(value);
                        cometparam.NoCPUs = Integer.parseInt(value);
                        msgfparam.NoCPUs = Integer.parseInt(value);
                        NoThreadPerFile = Integer.parseInt(value);
                        break;
                    }
                    case "PepFDR": {
                        tandemparam.PepFDR = Float.parseFloat(value);
                        cometparam.PepFDR = Float.parseFloat(value);
                        msgfparam.PepFDR = Float.parseFloat(value);
                        break;
                    }
                    case "ProtFDR": {
                        tandemparam.ProtFDR = Float.parseFloat(value);
                        cometparam.ProtFDR = Float.parseFloat(value);
                        msgfparam.ProtFDR = Float.parseFloat(value);
                        break;
                    }
                    case "DataSetLevelPepFDR": {
                        DataSetLevelPepFDR = Boolean.parseBoolean(value);
                    }
                    case "DecoyPrefix": {
                        if (!"".equals(value)) {
                            tandemparam.DecoyPrefix = value;
                            cometparam.DecoyPrefix = value;
                            msgfparam.DecoyPrefix = value;
                        }
                        break;
                    }
                    case "UserMod": {
                        UserMod = value;
                        break;
                    }
                    case "TPP":{
                        tandemparam.msconvertpath = value;
                        cometparam.msconvertpath = value;
                        msgfparam.msconvertpath = value;
                        tandemparam.tandempath = value;
                        tandemparam.tandem2XML = value;
                        tandemparam.xinteractpath = value;
                        cometparam.xinteractpath = value;
                        msgfparam.xinteractpath = value;
                        msgfparam.idconvert = value;
                    }
                    case "msconvertpath": {
                        tandemparam.msconvertpath = value;
                        cometparam.msconvertpath = value;
                        msgfparam.msconvertpath = value;
                        break;
                    }
                    case "tandempath": {
                        tandemparam.tandempath = value;
                        break;
                    }
                    case "tandem2XML": {
                        tandemparam.tandem2XML = value;
                        break;
                    }
                    case "xinteractpath": {
                        tandemparam.xinteractpath = value;
                        cometparam.xinteractpath = value;
                        msgfparam.xinteractpath = value;
                        break;
                    }
                    case "cometpath": {
                        cometparam.cometpath = value;
                        break;
                    }
                    case "msgfpath": {
                        msgfparam.msgfpath = value;
                        break;
                    }
                    case "idconvertpath": {
                        msgfparam.idconvert = value;
                    }
                    case "Fasta": {
                        tandemparam.FastaPath = value;
                        cometparam.FastaPath = value;
                        msgfparam.FastaPath = value;
                        break;
                    }
//</editor-fold>

                    //<editor-fold defaultstate="collapsed" desc="tandem parameters">
                    case "tandem.parampath": {
                        tandemparam.parameterPath = value;
                        break;
                    }
                    case "tandem.MinNoPeaksScoring": {
                        tandemparam.MinNoPeaksScoring = Integer.parseInt(value);
                        break;
                    }
                    case "tandem.MinNoPeaks": {
                        tandemparam.MinNoPeaks = Integer.parseInt(value);
                        break;
                    }
                    case "tandem.TotalPeaks": {
                        tandemparam.TotalPeaks = Integer.parseInt(value);
                        break;
                    }
                    case "tandem.SemiCleavage": {
                        tandemparam.SemiCleavage = Boolean.parseBoolean(value);
                        break;
                    }
                    case "tandem.VarMod": {
                        tandemparam.PotentialModification = value;
                        break;
                    }
                    case "tandem.FixMod": {
                        tandemparam.FixModification = value;
                        break;
                    }
                    case "tandem.ModMotif": {
                        tandemparam.PotentialModMotif = value;
                        break;
                    }
                    case "tandem.xinteractpara": {
                        tandemparam.xinteractpara = value;
                        break;
                    }
                    case "tandem.SpectrumConditioning": {
                        tandemparam.SpectrumConditioning = Boolean.parseBoolean(value);
                        break;
                    }
                    case "tandem.Scoring": {
                        tandemparam.Scoring = value;
                        break;
                    }
//</editor-fold>

                    //<editor-fold defaultstate="collapsed" desc="comet paramters">
                    case "comet.parampath": {
                        cometparam.parameterPath = value;
                        break;
                    }
                    case "comet.xinteractpara": {
                        cometparam.xinteractpara = value;
                        break;
                    }
                    
//</editor-fold>
                    
                    //<editor-fold defaultstate="collapsed" desc="msgf parameters">
                    case "msgf.parampath":{
                        msgfparam.parameterPath=value;
                    }
                    case "msgf.xinteractpara": {
                        msgfparam.xinteractpara = value;
                        break;
                    }
                    case "msgf.insID":{
                        msgfparam.MSGFInstrumentID=Integer.parseInt(value);
                    }
                    case "msgf.fragID":{
                        msgfparam.MSGFFragmentMethodID=Integer.parseInt(value);
                    }
                    case "msgf.enzID":{
                        msgfparam.MSGFEnzymeID=Integer.parseInt(value);
                    }
//</editor-fold>

                    case "BuildPeptideCandiate": {
                        BuildPeptideCandiate = Boolean.parseBoolean(value);
                        break;
                    }
                    case "CombineProtXML": {
                        CombineProtXML = value;
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
                    case "ExternalLibDecoyTag": {
                        ExternalLibDecoyTag = value;
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
                    case "WindowType": {
                        WindowType = value;
                        switch (WindowType) {
                            case "SWATH": {
                                dataType = SpectralDataType.DataType.DIA_F_Window;
                                break;
                            }
                            case "V_SWATH": {
                                dataType = SpectralDataType.DataType.DIA_V_Window;
                                break;
                            }
                            case "MSX": {
                                dataType = SpectralDataType.DataType.MSX;
                                break;
                            }
                            case "MSE": {
                                dataType = SpectralDataType.DataType.MSe;
                                break;
                            }
                            case "pSMART": {
                                dataType = SpectralDataType.DataType.pSMART;
                                break;
                            }
                            case "WiSIM": {
                                dataType = SpectralDataType.DataType.WiSIM;
                                break;
                            }
                        }
                        break;
                    }
                    case "WindowSize": {
                        WindowSize = Float.parseFloat(value);
                        break;
                    }
                }
            }
        }
//</editor-fold>

//        if (Mysql) {
//            try {
//                connectionManager = new ConnectionManager(Address, Integer.parseInt(Port), DBName, user, pwd);
//                Connection connection = connectionManager.GetConnection();
//            } catch (NumberFormatException | SQLException e) {
//                logger.info("MySQL Connection failed, please check");
//                logger.trace(e.getMessage());
//                return;
//            }
//        }

        PTMManager.GetInstance();
        EnzymeManager.GetInstance();
        if (!UserMod.equals("")) {
            try {
                PTMManager.GetInstance().ImportUserMod(UserMod);
            } catch (Exception e) {
                logger.info("Importing " + UserMod + " failed, please check");
                logger.error(ExceptionUtils.getStackTrace(e));
            }
        }

        msgfsearch.BuildDatabaseIndex();
        ArrayList<DIAPack> FileList = new ArrayList<>();
        ArrayList<LCMSID> LCMSIDList = new ArrayList<>();
        ArrayList<String> PepXMLs = new ArrayList<>();
        HashMap<String, HashMap<String, FragmentPeak>> IDSummaryFragments = new HashMap<>();

        ExecutorService executorPool = null;
        //executorPool = Executors.newFixedThreadPool(NoOfConcurrentFile);
        try {
            File folder = new File(WorkFolder);
            for (final File fileEntry : folder.listFiles()) {
                if (fileEntry.isDirectory()) {
                    String SWATH_mzXML = fileEntry.getPath() + "/" + fileEntry.getName().replace(".mzXML", "") + ".mzXML";
                    if (new File(SWATH_mzXML).exists()) {
                        //long time = System.currentTimeMillis();
                        logger.info("=================================================================================================");
                        logger.info("Processing " + SWATH_mzXML + "....");
                        DIAPack DiaFile = new DIAPack(SWATH_mzXML, NoCPUs);                        
                        FileList.add(DiaFile);
                        if (SignalExtraction) {
                            DiaFile.SetDataType(dataType);
                            DiaFile.SetParameter(para);
                            //DiaFile.SetMySQLConnection(connectionManager);
                            if (dataType == SpectralDataType.DataType.DIA_F_Window) {
                                DiaFile.SetWindowSize(WindowSize);
                            } else if (dataType == SpectralDataType.DataType.DIA_V_Window) {
                                for (XYData window : WindowList) {
                                    DiaFile.AddVaribleWindow(window);
                                }
                            }
                            if (dataType == SpectralDataType.DataType.WiSIM) {
                                for (XYData window : MS1WindowList) {
                                    DiaFile.AddMS1Window(window);
                                }
                            }
                            DiaFile.SaveDIASetting();
                            DiaFile.SaveParams();
                            DiaFile.ExportPrecursorPeak = ExportPrecursorPeak;
                            DiaFile.ExportFragmentPeak = ExportFragmentPeak;

                            //SEThread thread=new SEThread(DiaFile);
                            //executorPool.execute(thread);
                            logger.info("Module A: Signal extraction");
                            DiaFile.process();
                            DiaFile.ms1lcms.ClearAllPeaks();
                            DiaFile.ClearStructure();
                        }
                        else {
                            DiaFile.LoadDIASetting();
                            DiaFile.LoadParams();
                        }
                        //time = System.currentTimeMillis() - time;
                        //logger.info(SWATH_mzXML + " processed time:" + String.format("%d hour, %d min, %d sec", TimeUnit.MILLISECONDS.toHours(time), TimeUnit.MILLISECONDS.toMinutes(time) - TimeUnit.HOURS.toMinutes(TimeUnit.MILLISECONDS.toHours(time)), TimeUnit.MILLISECONDS.toSeconds(time) - TimeUnit.MINUTES.toSeconds(TimeUnit.MILLISECONDS.toMinutes(time))));
                    }
                }
            }

//            executorPool.shutdown();
//            try {
//                executorPool.awaitTermination(Long.MAX_VALUE, TimeUnit.NANOSECONDS);
//            } catch (InterruptedException e) {
//                logger.info("interrupted..");
//            }
            if(GeneratePBS){
                for (DIAPack DiaFile : FileList) {
                    DiaFile.GenerateSEJob(args[0], "");
                    DiaFile.GetPSBSaerchJobs(SearchEngines, MGFtag);
                }
            }            

            if (!DBSerach && !TargetedExtraction && !ExternalLibSearch  && !Quantitation) {
                return;
            }

            if (DBSerach) {
                executorPool = Executors.newFixedThreadPool(NoOfConcurrentFile);
                logger.info("Module B: Untargeted spectrum-centric search, No of files:" + FileList.size() + ", No of search engines:" + SearchEngines.size());
                logger.info("MS/MS database searching");
                for (DIAPack DiaFile : FileList) {
                    ArrayList<MSMSDBSearch> jobs = DiaFile.GetAllDBSearchJobs(SearchEngines, MGFtag);
                    for (MSMSDBSearch db : jobs) {
                        executorPool.execute(db);
                    }
                }
                executorPool.shutdown();
//                while (!executorPool.isTerminated()) {
//                }
                try {
                    executorPool.awaitTermination(Long.MAX_VALUE, TimeUnit.NANOSECONDS);
                } catch (InterruptedException e) {
                    logger.info("interrupted..");
                }
                
                executorPool = Executors.newFixedThreadPool(NoOfConcurrentFile);
                logger.info("Performing ProteinProphet");
                for (DIAPack DiaFile : FileList) {
                    //System.out.println("Add "+DiaFile.Filename);
                    TPPThread thread = new TPPThread(DiaFile, SearchEngines, MGFtag);
                    executorPool.execute(thread);                    
                }
                executorPool.shutdown();
//                while (!executorPool.isTerminated()) {
//                }
                try {
                    executorPool.awaitTermination(Long.MAX_VALUE, TimeUnit.NANOSECONDS);
                } catch (InterruptedException e) {
                    logger.info("interrupted..");
                }
            }
            if(DataSetLevelPepFDR){
                tag = "DataSetFDR";
            }
            if (Quantitation) {
                executorPool = Executors.newFixedThreadPool(NoOfConcurrentFile);
                logger.info("Module D: Quantitation for identified peptide ions");
                
                LCMSID combinePepID=null;
                if(DataSetLevelPepFDR){
                    combinePepID=LCMSID.ReadLCMSIDSerialization(WorkFolder+"combinePepID.SerFS");
                    if (combinePepID == null) {
                        FDR_DataSetLevel fdr = new FDR_DataSetLevel();
                        fdr.GeneratePepIonList(FileList, tandemparam, WorkFolder + "combinePepID.SerFS");
                        //fdr.combineID.CheckRT("combine");
                        combinePepID = fdr.combineID;
                        combinePepID.WriteLCMSIDSerialization(WorkFolder + "combinePepID.SerFS");
                    }
                }
                
                for (DIAPack DiaFile : FileList) {
                    DiaFile.SetNoCPUs(NoThreadPerFile);
                    QuantThread td = new QuantThread(DiaFile, tandemparam,combinePepID);
                    executorPool.execute(td);
                    //time = System.currentTimeMillis() - time;
                    //System.out.println(DiaFile.Filename + " processed time:" + String.format("%d hour, %d min, %d sec\n", TimeUnit.MILLISECONDS.toHours(time), TimeUnit.MILLISECONDS.toMinutes(time) - TimeUnit.HOURS.toMinutes(TimeUnit.MILLISECONDS.toHours(time)), TimeUnit.MILLISECONDS.toSeconds(time) - TimeUnit.MINUTES.toSeconds(TimeUnit.MILLISECONDS.toMinutes(time))));
                    HashMap<String, FragmentPeak> FragMap = new HashMap<>();
                    IDSummaryFragments.put(DiaFile.Filename, FragMap);                    
                }
                executorPool.shutdown();
//                while (!executorPool.isTerminated()) {
//                }
                  try {
                    executorPool.awaitTermination(Long.MAX_VALUE, TimeUnit.NANOSECONDS);
                } catch (InterruptedException e) {
                    logger.info("interrupted..");
                }                
            }
            
            if (TargetedExtraction) {
                logger.info("=================================================================================================");
                logger.info("Module C: Targeted extraction");

                FragmentLibManager libManager = new FragmentLibManager(InternalLibID);
                libManager = FragmentLibManager.ReadFragmentLibSerialization(WorkFolder, InternalLibID);
                
                if (libManager == null) {
                    libManager = new FragmentLibManager(InternalLibID);
                    logger.info("Building internal spectral library");
                    LCMSIDList = new ArrayList<>();
                    for (DIAPack dia : FileList) {
                        LCMSIDList.add(dia.IDsummary);
                    }
                    //libManager.ImportFragLib(LCMSIDList);
                    libManager.ImportFragLibTopFrag(LCMSIDList,0.5f,6);
                    libManager.WriteFragmentLibSerialization(WorkFolder);
                }
                libManager.ReduceMemoryUsage();
                
                if (BuildPeptideCandiate) {
                    logger.info("Building retention time prediction model and generate candidate peptide list");
                    for (int i = 0; i < FileList.size(); i++) {
                        for (int j = i + 1; j < FileList.size(); j++) {
                            RTAlignedPepIonMapping alignment = new RTAlignedPepIonMapping(WorkFolder, FileList.get(i).GetParameter(), FileList.get(i).IDsummary, FileList.get(j).IDsummary);
                            alignment.GenerateModel();
                            alignment.GenerateMappedPepIon();
                            //alignment.ExportMappedPepIon(connectionManager);
                        }
                        FileList.get(i).ExportID(tag);
                        FileList.get(i).IDsummary = null;
                    }
                }

                for (int i = 0; i < FileList.size(); i++) {
                    DIAPack diafile = FileList.get(i);
                    diafile.IDsummary = null;
                    diafile.ReadSerializedLCMSID(tag);
                    diafile.IDsummary.FastaPath = tandemparam.FastaPath;
                    diafile.UseMappedIon = true;
                    diafile.FilterMappedIonByProb = false;
                    diafile.BuildStructure();
                    diafile.ms1lcms.ReadPeakCluster();
                    diafile.GenerateMassCalibrationRTMap();
                    diafile.AssignMappedPepQuant(false, libManager);
                    PepXMLs.add(diafile.GetiProphExtPepxml(libManager.LibID));
                    diafile.ms1lcms.ClearAllPeaks();
                    diafile.IDsummary.RemoveLowProbMappedIon(0.1f);
                    diafile.ExportID("DIA_lib"+tag);
                    diafile.ClearStructure();
                    diafile.IDsummary = null;
                    System.gc();
                }
            }
            
            if (ExternalLibSearch) {
                Logger.getRootLogger().info("Module C: Targeted extraction (external library search)");
                FragmentLibManager ExlibManager = FragmentLibManager.ReadFragmentLibSerialization(FilenameUtils.getFullPath(ExternalLibPath), FilenameUtils.getBaseName(ExternalLibPath));
                if (ExlibManager == null) {
                    ExlibManager = new FragmentLibManager(FilenameUtils.getBaseName(ExternalLibPath), null);
                    ExlibManager.ImportFragLibByTraML(ExternalLibPath, ExternalLibDecoyTag);
                    //ExlibManager.ImportFragLibBySPTXT(ExternalLibPath);
                    ExlibManager.WriteFragmentLibSerialization(WorkFolder);
                }
                Logger.getRootLogger().info("No. of peptide ions in external lib:" + ExlibManager.PeptideFragmentLib.size());
                for (DIAPack diafile : FileList) {
                    if (diafile.IDsummary == null) {
                        diafile.ReadSerializedLCMSID("DIA_lib"+tag);
                    }
                    RTMappingExtLib RTmap = new RTMappingExtLib(diafile.IDsummary, ExlibManager, diafile.GetParameter());
                    RTmap.GenerateModel();
                    RTmap.GenerateMappedPepIon();

                    diafile.BuildStructure();
                    diafile.ms1lcms.ReadPeakCluster();
                    diafile.GenerateMassCalibrationRTMap();
                    diafile.AssignMappedPepQuant(false, ExlibManager,ReSearchProb);
                     PepXMLs.add(diafile.GetiProphExtPepxml(ExlibManager.LibID));
                    diafile.ms1lcms.ClearAllPeaks();
                    diafile.IDsummary.ReduceMemoryUsage();
                    diafile.IDsummary.RemoveLowProbMappedIon(ProbThreshold);                    
                    diafile.ExportID(FilenameUtils.getBaseName(ExternalLibPath));
                    diafile.ClearStructure();
                    Logger.getRootLogger().info("Peptide ions: " + diafile.IDsummary.GetPepIonList().size() + " Mapped ions: " + diafile.IDsummary.GetMappedPepIonList().size());
                }
            }
            LCMSID protID = null;
            if (!"".equals(CombineProtXML)) {
                if (!new File(CombineProtXML).exists()) {
                    for (DIAPack DiaFile : FileList) {
                        DiaFile.SetIPROPHETPepXML();
                        for (String pepxml : DiaFile.iProphPepXMLs) {
                            PepXMLs.add(pepxml);
                        }
                    }
                    try {
                        ProteinProphetCombine proteinprophet = new ProteinProphetCombine();
                        proteinprophet.ProteinProphetCombineSearchiProphet(tandemparam, PepXMLs, CombineProtXML);
                    } catch (Exception e) {
                        logger.error("ProteinProphet failed");
                    }
                }

                protID = LCMSID.ReadLCMSIDSerialization(CombineProtXML);
                if (!"".equals(CombineProtXML) && protID == null) {
                    protID = new LCMSID(CombineProtXML, tandemparam.DecoyPrefix, tandemparam.FastaPath);
                    ProtXMLParser protxmlparser = new ProtXMLParser(protID, CombineProtXML, 0f);
                    protID.RemoveLowLocalPWProtein(0.2f);
                    protID.FilterByProteinDecoyFDRUsingMaxIniProb(tandemparam.DecoyPrefix, tandemparam.ProtFDR);
                    protID.RemoveLowWeightPep(0.9f);
                    protID.GenerateIndisProtMap();
                    protID.LoadSequence();
                    protID.WriteLCMSIDSerialization(CombineProtXML);
                    System.out.println("Protein No in combined file:" + protID.ProteinList.size());
                }
            }
            if (ExportTable) {
                for (int i = 0; i < FileList.size(); i++) {
                    DIAPack diafile = FileList.get(i);
                    if (diafile.IDsummary == null) {
                        if(!diafile.ReadSerializedLCMSID(FilenameUtils.getBaseName(ExternalLibPath))){
                            if(!diafile.ReadSerializedLCMSID("DIA_lib")){
                                diafile.ReadSerializedLCMSID();
                            }
                        }
                    }
                    diafile.UseMappedIon = true;
                    diafile.IDsummary.RemoveLowProbMappedIon(0.99f);
                    diafile.IDsummary.GenerateProteinByRefIDByPepSeq(protID, true);
                    diafile.IDsummary.ReMapProPep();
                    //diafile.IDsummary.ClearPSMs();
                    diafile.IDsummary.ClearAssignPeakCluster();
                }
                FragmentSelection fragselection = new FragmentSelection(LCMSIDList);
                fragselection.freqPercent = 0.5f;
                fragselection.GeneratePepFragScoreMap();
                fragselection.GenerateTopFragMap(6);
                fragselection.GenerateProtPepScoreMap(0.9f);
                fragselection.GenerateTopPepMap(6);
                ExportTable export = new ExportTable(WorkFolder, LCMSIDList, IDSummaryFragments, protID, fragselection);
                export.Export(6, 6, 0.5f);
            }

        } catch (Exception e) {
            Logger.getRootLogger().error(ExceptionUtils.getStackTrace(e));
            throw e;
        }
    }
}
