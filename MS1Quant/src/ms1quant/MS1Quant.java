/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package ms1quant;

import MSUmpire.BaseDataStructure.InstrumentParameter;
import MSUmpire.BaseDataStructure.UmpireInfo;
import MSUmpire.LCMSBaseStructure.LCMSPeakMS1;
import MSUmpire.MSMSDBSearch.DBSearchParam;
import MSUmpire.MSMSDBSearch.TandemParam;
import MSUmpire.PSMDataStructure.LCMSID;
import MSUmpire.PSMDataStructure.PTMManager;
import MSUmpire.SearchResultParser.PepXMLParser;
import MSUmpire.Utility.ConsoleLogger;
import java.io.*;
import java.util.Arrays;
import java.util.HashMap;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;
import org.apache.commons.io.FilenameUtils;
import org.apache.commons.lang.exception.ExceptionUtils;
import org.apache.log4j.Level;
import org.apache.log4j.Logger;

/**
 *
 * @author Chih-Chiang Tsou <chihchiang.tsou@gmail.com>
 */
public class MS1Quant {

    /**
     * @param args the command line arguments MS1Quant parameterfile
     */
    public static void main(String[] args) throws Exception {

        BufferedReader reader = null;
        try {
            System.out.println("=================================================================================================");
            System.out.println("Umpire MS1 quantification and feature detection analysis (version: " + UmpireInfo.GetInstance().Version + ")");
            if (args.length < 3 || !args[1].startsWith("-mode")) {
                System.out.println("command : java -jar -Xmx10G MS1Quant.jar ms1quant.params -mode[1 or 2] [Option]");
                System.out.println("\n-mode");
                System.out.println("\t1:Single file mode--> mzXML_file PepXML_file");
                System.out.println("\t\tEx: -mode1 file1.mzXML file1.pep.xml");
                System.out.println("\t2:Folder mode--> mzXML_Folder PepXML_Folder, all generated csv tables will be merged into a single csv file");
                System.out.println("\t\tEx: -mode2 /data/mzxml/ /data/pepxml/");
                System.out.println("\nOptions");
                System.out.println("\t-C\tNo of concurrent files to be processed (only for folder mode), Ex. -C5, default:1");
                System.out.println("\t-p\tMinimum probability, Ex. -p0.9, default:0.9");
                System.out.println("\t-ID\tDetect identified feature only");
                System.out.println("\t-O\toutput folder, Ex. -O/data/");
                return;
            }
            ConsoleLogger consoleLogger = new ConsoleLogger();
            consoleLogger.SetConsoleLogger(Level.DEBUG);
            consoleLogger.SetFileLogger(Level.DEBUG, FilenameUtils.getFullPath(args[0]) + "ms1quant_debug.log");
            Logger logger = Logger.getRootLogger();
            logger.debug("Command: " + Arrays.toString(args));
            logger.info("MS1Quant version: " + UmpireInfo.GetInstance().Version);

            String parameterfile = args[0];
            logger.info("Parameter file: " + parameterfile);
            File paramfile = new File(parameterfile);
            if (!paramfile.exists()) {
                logger.error("Parameter file " + paramfile.getAbsolutePath() + " cannot be found. The program will exit.");
            }

            reader = new BufferedReader(new FileReader(paramfile.getAbsolutePath()));
            String line = "";
            InstrumentParameter param = new InstrumentParameter(InstrumentParameter.InstrumentType.TOF5600);
            int NoCPUs = 2;
            int NoFile = 1;
            param.DetermineBGByID = false;
            param.EstimateBG = true;

            //<editor-fold defaultstate="collapsed" desc="Read parameter file">
            while ((line = reader.readLine()) != null) {
                if (!"".equals(line) && !line.startsWith("#")) {
                    logger.info(line);
                    //System.out.println(line);
                    if (line.split("=").length < 2) {
                        continue;
                    }
                    if (line.split("=").length < 2) {
                        continue;
                    }
                    String type = line.split("=")[0].trim();
                    if (type.startsWith("para.")) {
                        type = type.replace("para.", "SE.");
                    }
                    String value = line.split("=")[1].trim();
                    switch (type) {
                        case "Thread": {
                            NoCPUs = Integer.parseInt(value);
                            break;
                        }
                        //<editor-fold defaultstate="collapsed" desc="instrument parameters">

                        case "SE.MS1PPM": {
                            param.MS1PPM = Float.parseFloat(value);
                            break;
                        }
                        case "SE.MS2PPM": {
                            param.MS2PPM = Float.parseFloat(value);
                            break;
                        }
                        case "SE.SN": {
                            param.SNThreshold = Float.parseFloat(value);
                            break;
                        }
                        case "SE.MS2SN": {
                            param.MS2SNThreshold = Float.parseFloat(value);
                            break;
                        }
                        case "SE.MinMSIntensity": {
                            param.MinMSIntensity = Float.parseFloat(value);
                            break;
                        }
                        case "SE.MinMSMSIntensity": {
                            param.MinMSMSIntensity = Float.parseFloat(value);
                            break;
                        }
                        case "SE.MinRTRange": {
                            param.MinRTRange = Float.parseFloat(value);
                            break;
                        }
                        case "SE.MaxNoPeakCluster": {
                            param.MaxNoPeakCluster = Integer.parseInt(value);
                            param.MaxMS2NoPeakCluster = Integer.parseInt(value);
                            break;
                        }
                        case "SE.MinNoPeakCluster": {
                            param.MinNoPeakCluster = Integer.parseInt(value);
                            param.MinMS2NoPeakCluster = Integer.parseInt(value);
                            break;
                        }
                        case "SE.MinMS2NoPeakCluster": {
                            param.MinMS2NoPeakCluster = Integer.parseInt(value);
                            break;
                        }
                        case "SE.MaxCurveRTRange": {
                            param.MaxCurveRTRange = Float.parseFloat(value);
                            break;
                        }
                        case "SE.Resolution": {
                            param.Resolution = Integer.parseInt(value);
                            break;
                        }
                        case "SE.RTtol": {
                            param.RTtol = Float.parseFloat(value);
                            break;
                        }
                        case "SE.NoPeakPerMin": {
                            param.NoPeakPerMin = Integer.parseInt(value);
                            break;
                        }
                        case "SE.StartCharge": {
                            param.StartCharge = Integer.parseInt(value);
                            break;
                        }
                        case "SE.EndCharge": {
                            param.EndCharge = Integer.parseInt(value);
                            break;
                        }
                        case "SE.MS2StartCharge": {
                            param.MS2StartCharge = Integer.parseInt(value);
                            break;
                        }
                        case "SE.MS2EndCharge": {
                            param.MS2EndCharge = Integer.parseInt(value);
                            break;
                        }
                        case "SE.NoMissedScan": {
                            param.NoMissedScan = Integer.parseInt(value);
                            break;
                        }
                        case "SE.Denoise": {
                            param.Denoise = Boolean.valueOf(value);
                            break;
                        }
                        case "SE.EstimateBG": {
                            param.EstimateBG = Boolean.valueOf(value);
                            break;
                        }
                        case "SE.RemoveGroupedPeaks": {
                            param.RemoveGroupedPeaks = Boolean.valueOf(value);
                            break;
                        }
                        case "SE.MinFrag": {
                            param.MinFrag = Integer.parseInt(value);
                            break;
                        }
                        case "SE.IsoPattern": {
                            param.IsoPattern = Float.valueOf(value);
                            break;
                        }
                        case "SE.StartRT": {
                            param.startRT = Float.valueOf(value);
                        }
                        case "SE.EndRT": {
                            param.endRT = Float.valueOf(value);
                        }

//</editor-fold>
                    }
                }
            }
//</editor-fold>

            int mode = 1;
            if (args[1].equals("-mode2")) {
                mode = 2;
            } else if (args[1].equals("-mode1")) {
                mode = 1;
            } else {
                logger.error("-mode number not recongized. The program will exit.");
            }

            String mzXML = "";
            String pepXML = "";
            String mzXMLPath = "";
            String pepXMLPath = "";
            File mzXMLfile = null;
            File pepXMLfile = null;
            File mzXMLfolder = null;
            File pepXMLfolder = null;
            int idx = 0;
            if (mode == 1) {
                mzXML = args[2];
                logger.info("Mode1 mzXML file: " + mzXML);
                mzXMLfile = new File(mzXML);
                if (!mzXMLfile.exists()) {
                    logger.error("Mode1 mzXML file " + mzXMLfile.getAbsolutePath() + " cannot be found. The program will exit.");
                    return;
                }
                pepXML = args[3];
                logger.info("Mode1 pepXML file: " + pepXML);
                pepXMLfile = new File(pepXML);
                if (!pepXMLfile.exists()) {
                    logger.error("Mode1 pepXML file " + pepXMLfile.getAbsolutePath() + " cannot be found. The program will exit.");
                    return;
                }
                idx = 4;
            } else if (mode == 2) {
                mzXMLPath = args[2];
                logger.info("Mode2 mzXML folder: " + mzXMLPath);
                mzXMLfolder = new File(mzXMLPath);
                if (!mzXMLfolder.exists()) {
                    logger.error("Mode2 mzXML folder " + mzXMLfolder.getAbsolutePath() + " does not exist. The program will exit.");
                    return;
                }
                pepXMLPath = args[3];
                logger.info("Mode2 pepXML folder: " + pepXMLPath);
                pepXMLfolder = new File(pepXMLPath);
                if (!pepXMLfolder.exists()) {
                    logger.error("Mode2 pepXML folder " + pepXMLfolder.getAbsolutePath() + " does not exist. The program will exit.");
                    return;
                }
                idx = 4;
            }

            String outputfolder = "";
            float MinProb = 0f;
            for (int i = idx; i < args.length; i++) {
                if (args[i].startsWith("-")) {
                    if (args[i].equals("-ID")) {
                        param.TargetIDOnly = true;
                        logger.info("Detect ID feature only: true");
                    }
                    if (args[i].startsWith("-O")) {
                        outputfolder = args[i].substring(2);
                        logger.info("Output folder: " + outputfolder);

                        File outputfile = new File(outputfolder);
                        if (!outputfolder.endsWith("\\") | outputfolder.endsWith("/")) {
                            outputfolder += "/";
                        }
                        if (!outputfile.exists()) {
                            outputfile.mkdir();
                        }
                    }
                    if (args[i].startsWith("-C")) {
                        try {
                            NoFile = Integer.parseInt(args[i].substring(2));
                            logger.info("No of concurrent files: " + NoFile);
                        } catch (Exception ex) {
                            logger.error(args[i] + " is not a correct integer format, will process only one file at a time.");
                        }
                    }
                    if (args[i].startsWith("-p")) {
                        try {
                            MinProb = Float.parseFloat(args[i].substring(2));
                            logger.info("probability threshold: " + MinProb);
                        } catch (Exception ex) {
                            logger.error(args[i] + " is not a correct format, will use 0 as threshold instead.");
                        }
                    }
                }
            }

            reader.close();
            TandemParam tandemparam = new TandemParam(DBSearchParam.SearchInstrumentType.TOF5600);
            PTMManager.GetInstance();

            if (param.TargetIDOnly) {
                param.EstimateBG = false;
                param.ApexDelta = 1.5f;
                param.NoMissedScan = 10;
                param.MiniOverlapP = 0.2f;
                param.RemoveGroupedPeaks = false;
                param.CheckMonoIsotopicApex = false;
                param.DetectByCWT = false;
                param.FillGapByBK = false;
                param.IsoCorrThreshold = -1f;
                param.SmoothFactor = 3;
            }

            if (mode == 1) {
                logger.info("Processing " + mzXMLfile.getAbsolutePath() + "....");
                long time = System.currentTimeMillis();
                LCMSPeakMS1 LCMS1 = new LCMSPeakMS1(mzXMLfile.getAbsolutePath(), NoCPUs);
                LCMS1.SetParameter(param);

                LCMS1.Resume = false;
                if (!param.TargetIDOnly) {
                    LCMS1.CreatePeakFolder();
                }
                LCMS1.ExportPeakClusterTable = true;

                if (pepXMLfile.exists()) {
                    tandemparam.InteractPepXMLPath = pepXMLfile.getAbsolutePath();
                    LCMS1.ParsePepXML(tandemparam, MinProb);
                    logger.info("No. of PSMs included: " + LCMS1.IDsummary.PSMList.size());
                    logger.info("No. of Peptide ions included: " + LCMS1.IDsummary.GetPepIonList().size());
                }

                if (param.TargetIDOnly) {
                    LCMS1.SaveSerializationFile = false;
                }

                if (param.TargetIDOnly || !LCMS1.ReadPeakCluster()) {
                    LCMS1.PeakClusterDetection();
                }

                if (pepXMLfile.exists()) {
                    LCMS1.AssignQuant(false);
                    LCMS1.IDsummary.ExportPepID(outputfolder);
                }
                time = System.currentTimeMillis() - time;
                logger.info(LCMS1.ParentmzXMLName + " processed time:" + String.format("%d hour, %d min, %d sec", TimeUnit.MILLISECONDS.toHours(time), TimeUnit.MILLISECONDS.toMinutes(time) - TimeUnit.HOURS.toMinutes(TimeUnit.MILLISECONDS.toHours(time)), TimeUnit.MILLISECONDS.toSeconds(time) - TimeUnit.MINUTES.toSeconds(TimeUnit.MILLISECONDS.toMinutes(time))));
                LCMS1.BaseClearAllPeaks();
                LCMS1.SetmzXML(null);
                LCMS1.IDsummary = null;
                LCMS1 = null;
                System.gc();
            } else if (mode == 2) {

                LCMSID IDsummary = new LCMSID("", "", "");
                logger.info("Parsing all pepXML files in " + pepXMLPath + "....");
                for (File file : pepXMLfolder.listFiles()) {
                    if (file.getName().toLowerCase().endsWith("pep.xml") || file.getName().toLowerCase().endsWith("pepxml")) {
                        PepXMLParser pepXMLParser = new PepXMLParser(IDsummary, file.getAbsolutePath(), MinProb);
                    }
                }
                HashMap<String, LCMSID> LCMSIDMap = IDsummary.GetLCMSIDFileMap();

                ExecutorService executorPool = null;
                executorPool = Executors.newFixedThreadPool(NoFile);

                logger.info("Processing all mzXML files in " + mzXMLPath + "....");
                for (File file : mzXMLfolder.listFiles()) {
                    if (file.getName().toLowerCase().endsWith("mzxml")) {
                        LCMSID id = LCMSIDMap.get(FilenameUtils.getBaseName(file.getName()));
                        if (!id.PSMList.isEmpty()) {
                            MS1TargetQuantThread thread = new MS1TargetQuantThread(file, id, NoCPUs, outputfolder, param);
                            executorPool.execute(thread);
                        }
                    }
                }
                LCMSIDMap.clear();
                LCMSIDMap = null;
                IDsummary = null;
                executorPool.shutdown();
                try {
                    executorPool.awaitTermination(Long.MAX_VALUE, TimeUnit.NANOSECONDS);
                } catch (InterruptedException e) {
                    logger.info("interrupted..");
                }

                if (outputfolder == null | outputfolder.equals("")) {
                    outputfolder = mzXMLPath;
                }

                logger.info("Merging PSM files..");
                File output = new File(outputfolder);
                FileWriter writer = new FileWriter(output.getAbsolutePath() + "/PSM_merge.csv");
                boolean header = false;
                for (File csvfile : output.listFiles()) {
                    if (csvfile.getName().toLowerCase().endsWith("_psms.csv")) {
                        BufferedReader outreader = new BufferedReader(new FileReader(csvfile));
                        String outline = outreader.readLine();
                        if (!header) {
                            writer.write(outline + "\n");
                            header = true;
                        }
                        while ((outline = outreader.readLine()) != null) {
                            writer.write(outline + "\n");
                        }
                        outreader.close();
                        csvfile.delete();
                    }
                }
                writer.close();
            }
            logger.info("MS1 quant module is complete.");
        } catch (Exception e) {
            Logger.getRootLogger().error(ExceptionUtils.getStackTrace(e));
            throw e;
        }
    }
}
