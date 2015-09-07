/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package MSUmpire;

import MSUmpire.BaseDataStructure.InstrumentParameter;
import MSUmpire.BaseDataStructure.UmpireInfo;
import MSUmpire.LCMSBaseStructure.LCMSPeakMS1;
import MSUmpire.MSMSDBSearch.DBSearchParam;
import MSUmpire.MSMSDBSearch.TandemParam;
import MSUmpire.PSMDataStructure.PTMManager;
import Utility.ConsoleLogger;
import java.io.*;
import java.util.concurrent.TimeUnit;
import org.apache.commons.io.FilenameUtils;
import org.apache.commons.lang.exception.ExceptionUtils;
import org.apache.log4j.Level;
import org.apache.log4j.Logger;
import org.xmlpull.v1.XmlPullParserException;

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
            System.out.println("Umpire MS1 quantitation and feature detection analysis (version: " + UmpireInfo.GetInstance().Version + ")");
            if (args.length < 2 || args.length > 5) {
                System.out.println("command : java -jar -Xmx8G MS1Quant.jar mzMXL_file ms1quant.params PepXML_file MinProb [Option]");

                System.out.println("");
                System.out.println("\nOptions");
                System.out.println("\t-ID\tDetect identified feature only");
                return;
            }
            ConsoleLogger consoleLogger = new ConsoleLogger();
            consoleLogger.SetConsoleLogger(Level.DEBUG);
            consoleLogger.SetFileLogger(Level.DEBUG, FilenameUtils.getFullPath(args[0]) + "ms1quant_debug.log");
            Logger logger = Logger.getRootLogger();
            
            String parameterfile = args[1];
            logger.info("MS1Quant version: "+UmpireInfo.GetInstance().Version);            
            logger.info("parameter file: "+parameterfile);            
            File paramfile = new File(parameterfile);
            if (!paramfile.exists()) {
                logger.error("parameter file "+paramfile.getAbsolutePath() + " does not exist. The program will exit.");
            }
            
            String mzXML = args[0];
            logger.info("mzXML file: "+mzXML);
            File mzXMLfile = new File(mzXML);
            if (!mzXMLfile.exists()) {
                logger.error("mzXML file "+mzXMLfile.getAbsolutePath() + " does not exist. The program will exit.");                
            }
                        
            String pepXML="";
            if (args.length > 2) {
                pepXML = args[2];
                logger.info("pepXML file: " + pepXML);
            }
            File pepXMLfile = new File(pepXML);
            if (!pepXMLfile.exists()) {
                logger.error("pepXML file "+pepXMLfile.getAbsolutePath() + " does not exist. Will only do feature detection.");                
            }
            float MinProb=0f;
            if (args.length > 3) {
                try {
                    MinProb = Float.parseFloat(args[3]);
                    logger.info("probability threshold: "+MinProb);
                } catch (Exception ex) {
                    logger.error(args[3] + " is not a correct format, will use 0 as threshold instead.");
                }
            }
          
            reader = new BufferedReader(new FileReader(paramfile.getAbsolutePath()));
            String line = "";
            InstrumentParameter param = new InstrumentParameter(InstrumentParameter.InstrumentType.TOF5600);
            int NoCPUs = 2;            
            String UserMod = "";
            param.DetermineBGByID = false;
            param.EstimateBG = true;
            if (args.length > 4) {
                if (args[4].equals("-ID")) {
                    param.TargetIDOnly = true;
                    logger.info("Detect ID feature only: true");
                }
            }
            logger.info("Parameter file:" + paramfile.getAbsolutePath());
            
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
                        param.MaxMS2NoPeakCluster= Integer.parseInt(value);
                        break;
                    }
                    case "SE.MinNoPeakCluster": {
                        param.MinNoPeakCluster = Integer.parseInt(value);
                        param.MinMS2NoPeakCluster= Integer.parseInt(value);
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
            
            reader.close();
            TandemParam tandemparam = new TandemParam(DBSearchParam.SearchInstrumentType.TOF5600);
            PTMManager.GetInstance();
            if (!UserMod.equals("")) {
                try {
                    PTMManager.GetInstance().ImportUserMod(UserMod);
                } catch (IOException | XmlPullParserException e) {
                    logger.info("Importing " + UserMod + " failed, please check");
                    logger.error(ExceptionUtils.getStackTrace(e));
                }
            }
            logger.info("Processing " + mzXMLfile.getAbsolutePath() + "....");
            long time = System.currentTimeMillis();
            LCMSPeakMS1 LCMS1 = new LCMSPeakMS1(mzXMLfile.getAbsolutePath(), NoCPUs);
            LCMS1.SetParameter(param);
                        
            if (param.TargetIDOnly) {
                param.EstimateBG = false;
                param.ApexDelta = 3f;
                param.NoMissedScan = 10;
                param.MiniOverlapP = 0.1f;
                param.RemoveGroupedPeaks = false;
                param.CheckMonoIsotopicApex = false;
                param.DetectByCWT = false;
                param.FillGapByBK = false;
                param.IsoCorrThreshold = -1f;
            }
            
            LCMS1.Resume = false;
            LCMS1.CreatePeakFolder();
            LCMS1.ExportPeakClusterTable = true;
        
            if (pepXMLfile.exists()) {
                tandemparam.InteractPepXMLPath = pepXMLfile.getAbsolutePath();
                LCMS1.ParsePepXML(tandemparam, MinProb);
            }

            if (!LCMS1.ReadPeakCluster()) {
                LCMS1.PeakClusterDetection();
            }
            if (pepXMLfile.exists()) {
                LCMS1.AssignQuant(false);
                LCMS1.IDsummary.ExportPepID();
            }
            time = System.currentTimeMillis() - time;
            logger.info(LCMS1.ParentmzXMLName + " processed time:" + String.format("%d hour, %d min, %d sec", TimeUnit.MILLISECONDS.toHours(time), TimeUnit.MILLISECONDS.toMinutes(time) - TimeUnit.HOURS.toMinutes(TimeUnit.MILLISECONDS.toHours(time)), TimeUnit.MILLISECONDS.toSeconds(time) - TimeUnit.MINUTES.toSeconds(TimeUnit.MILLISECONDS.toMinutes(time))));
            LCMS1.BaseClearAllPeaks();
            LCMS1.SetmzXML(null);
            LCMS1.IDsummary = null;
            LCMS1 = null;
            System.gc();
        } catch (Exception e) {
            Logger.getRootLogger().error(ExceptionUtils.getStackTrace(e));
            throw e;      
        }
    }
}
