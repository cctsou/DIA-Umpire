/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package Exec;

import CXL_PeakPairFinder.CoElutePeak;
import CXL_PeakPairFinder.CrosslinkerPepFinder;
import CXL_PeakPairFinder.MS2PeakPairFinder;
import CXL_PeakPairFinder.MS2PeakPairFinder.PeakPair;
import CXL_PeakPairFinder.PeakPairFinder;
import CXL_PeakPairFinder.PrecursorCrossPepFinder;
import MSUmpire.BaseDataStructure.InstrumentParameter;
import MSUmpire.BaseDataStructure.ScanCollection;
import MSUmpire.BaseDataStructure.ScanData;
import MSUmpire.LCMSPeakStructure.LCMSPeakMS1;
import MSUmpire.SpectralProcessingModule.ScanIsotopePeakDetectionThread;
import MSUmpire.SpectralProcessingModule.ScanPeakGroup;
import MSUmpire.Utility.ConsoleLogger;
import crosslinker.Linker;
import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
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
public class CrossLinkerDetection {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) throws IOException, InterruptedException, ExecutionException {
        System.out.println("=================================================================================================");
        System.out.println("CXL peak pair finder");
        if (args.length < 2 || args.length > 3) {
            System.out.println("command format error, the correct format is: java -jar -Xmx8G CrossLinkerDetection.jar mzMXL_file cxl.params");
            return;
        }

        try {
            ConsoleLogger.SetConsoleLogger(Level.INFO);
            ConsoleLogger.SetFileLogger(Level.DEBUG, FilenameUtils.getFullPath(args[0]) + "cxl.log");
        } catch (Exception e) {
        }

        Logger logger = Logger.getRootLogger();

        InstrumentParameter parameter = null;
        Linker linker=new Linker();
        parameter = new InstrumentParameter(InstrumentParameter.InstrumentType.Orbitrap);
        parameter.MS1PPM = 20;
        parameter.MS2PPM = 600f;
        parameter.StartCharge = 1;
        parameter.EndCharge = 2;
        parameter.MS2StartCharge = 1;
        parameter.MS2EndCharge = 2;
        parameter.MinMS2NoPeakCluster = 1;
        parameter.MaxMS2NoPeakCluster = 4;
        parameter.RemoveGroupedPeaks = true;
        parameter.RemoveGroupedPeaksCorr = 0.6f;
        parameter.RemoveGroupedPeaksRTOverlap = 0.6f;
        parameter.MS2PairTopN = 5;
        parameter.DetectSameChargePairOnly = true;
        
        parameter.MinPeakPerPeakCurve = 3;
        parameter.MaxCurveRTRange = 5f;
        parameter.NoMissedScan = 3;
        parameter.MinMZ = 150f;
        parameter.DetectByCWT = true;
        parameter.FillGapByBK = false;
        parameter.IsoCorrThreshold = -1f;
        parameter.SmoothFactor = 3;
        parameter.IsoPattern = 0.6f;
        parameter.startRT = 15f;
        parameter.endRT = 150f;
        parameter.RTtol = 0.8f;
        parameter.MassDefectFilter = true;
        int NoCPUs = 10;
        String parameterfile = args[1];
        String FilePath = args[0];
        Logger.getRootLogger().info("Parameter file:" + parameterfile);
        Logger.getRootLogger().info("mzXML file:" + FilePath);
        BufferedReader reader = new BufferedReader(new FileReader(parameterfile));

        String line = "";

        //<editor-fold defaultstate="collapsed" desc="Read parameter file">
        while ((line = reader.readLine()) != null) {
            Logger.getRootLogger().info(line);
            if (!"".equals(line) && !line.startsWith("#")) {
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
                        parameter.MS1PPM = Float.parseFloat(value);
                        break;
                    }
                    case "SE.MS2PPM": {
                        parameter.MS2PPM = Float.parseFloat(value);
                        break;
                    }
                    case "SE.SN": {
                        parameter.SNThreshold = Float.parseFloat(value);
                        break;
                    }
                    case "SE.MS2SN": {
                        parameter.MS2SNThreshold = Float.parseFloat(value);
                        break;
                    }
                    case "SE.MinMSIntensity": {
                        parameter.MinMSIntensity = Float.parseFloat(value);
                        break;
                    }
                    case "SE.MinMSMSIntensity": {
                        parameter.MinMSMSIntensity = Float.parseFloat(value);
                        break;
                    }
                    case "SE.MinRTRange": {
                        parameter.MinRTRange = Float.parseFloat(value);
                        break;
                    }
                    case "SE.MaxNoPeakCluster": {
                        parameter.MaxNoPeakCluster = Integer.parseInt(value);
                        parameter.MaxMS2NoPeakCluster = Integer.parseInt(value);
                        break;
                    }
                    case "SE.MinNoPeakCluster": {
                        parameter.MinNoPeakCluster = Integer.parseInt(value);
                        parameter.MinMS2NoPeakCluster = Integer.parseInt(value);
                        break;
                    }
                    case "SE.MinMS2NoPeakCluster": {
                        parameter.MinMS2NoPeakCluster = Integer.parseInt(value);
                        break;
                    }
                    case "SE.MaxCurveRTRange": {
                        parameter.MaxCurveRTRange = Float.parseFloat(value);
                        break;
                    }
                    case "SE.Resolution": {
                        parameter.Resolution = Integer.parseInt(value);
                        break;
                    }
                    case "SE.RTtol": {
                        parameter.RTtol = Float.parseFloat(value);
                        break;
                    }
                    case "SE.NoPeakPerMin": {
                        parameter.NoPeakPerMin = Integer.parseInt(value);
                        break;
                    }
                    case "SE.StartCharge": {
                        parameter.StartCharge = Integer.parseInt(value);
                        break;
                    }
                    case "SE.EndCharge": {
                        parameter.EndCharge = Integer.parseInt(value);
                        break;
                    }
                    case "SE.MS2StartCharge": {
                        parameter.MS2StartCharge = Integer.parseInt(value);
                        break;
                    }
                    case "SE.MS2EndCharge": {
                        parameter.MS2EndCharge = Integer.parseInt(value);
                        break;
                    }
                    case "SE.NoMissedScan": {
                        parameter.NoMissedScan = Integer.parseInt(value);
                        break;
                    }
                    case "SE.Denoise": {
                        parameter.Denoise = Boolean.valueOf(value);
                        break;
                    }
                    case "SE.EstimateBG": {
                        parameter.EstimateBG = Boolean.valueOf(value);
                        break;
                    }
                    case "SE.RemoveGroupedPeaks": {
                        parameter.RemoveGroupedPeaks = Boolean.valueOf(value);
                        break;
                    }
                    case "SE.MinFrag": {
                        parameter.MinFrag = Integer.parseInt(value);
                        break;
                    }
                    case "SE.IsoPattern": {
                        parameter.IsoPattern = Float.valueOf(value);
                        break;
                    }
                    case "SE.StartRT": {
                        parameter.startRT = Float.valueOf(value);
                        break;
                    }
                    case "SE.EndRT": {
                        parameter.endRT = Float.valueOf(value);
                        break;
                    }
                    case "SE.RemoveGroupedPeaksRTOverlap": {
                        parameter.RemoveGroupedPeaksRTOverlap = Float.valueOf(value);
                        break;
                    }
                    case "SE.RemoveGroupedPeaksCorr": {
                        parameter.RemoveGroupedPeaksCorr = Float.valueOf(value);
                        break;
                    }
                    case "SE.MinMZ": {
                        parameter.MinMZ = Float.valueOf(value);
                        break;
                    }
                    case "SE.IsoCorrThreshold": {
                        parameter.IsoCorrThreshold = Float.valueOf(value);
                        break;
                    }
                    case "SE.MassDefectFilter": {
                        parameter.MassDefectFilter = Boolean.parseBoolean(value);
                        break;
                    }
                    case "SE.DetectSameChargePairOnly": {
                        parameter.DetectSameChargePairOnly = Boolean.parseBoolean(value);
                        break;
                    }
                    case "MS2PairTopN": {
                        parameter.MS2PairTopN = Integer.parseInt(value);
                        break;
                    }
                    case "SE.MS2Pairing": {
                        parameter.MS2Pairing = Boolean.parseBoolean(value);
                        break;
                    }
                    case "CXL.CoreMass":{
                        linker.Core=Float.parseFloat(value);
                        break;
                    }
                    case "CXL.ArmMass":{
                        linker.Arm=Float.parseFloat(value);
                        break;
                    }                            

//</editor-fold>//</editor-fold>                   
                }
            }
        }
//</editor-fold>

        String basename = FilenameUtils.getBaseName(FilePath);

        //Initialize an LCMSPeakMS1 data structure to do MS1 isotope peak cluster detection
        LCMSPeakMS1 lCMSPeakBase = null;
        ArrayList<ScanPeakGroup> MS2ScanPeakGroup = null;

        lCMSPeakBase = new LCMSPeakMS1(FilePath, parameter, NoCPUs);
        lCMSPeakBase.Resume=false;
        try {
            //Detect MS1 peak clusters
            lCMSPeakBase.PeakClusterDetection();
        } catch (FileNotFoundException | ParserConfigurationException | SAXException | DataFormatException | XmlPullParserException ex) {
            logger.error("Peak detection failed:");
            logger.error(ExceptionUtils.getStackTrace(ex));
        }

        //Remove detected peak clusters which appear in more than 90% of LC.
        lCMSPeakBase.RemoveContaminantPeaks(0.9f);

        //Sort peak clusters based on peak apex intensity of monoisotope peak
        float[] numArray = new float[lCMSPeakBase.PeakClusters.size()];
        for (int i = 0; i < lCMSPeakBase.PeakClusters.size(); i++) {
            numArray[i] = lCMSPeakBase.PeakClusters.get(i).PeakHeight[0];
        }
        Arrays.sort(numArray);

        //Get median of peak intensity of all peak clusters
        float medianIntensity;
        if (numArray.length % 2 == 0) {
            medianIntensity = ((float) numArray[numArray.length / 2] + (float) numArray[numArray.length / 2 - 1]) / 2;
        } else {
            medianIntensity = (float) numArray[numArray.length / 2];
        }

        ///////////Detect MS2 peak clusters for all MS2 spectra
        if (parameter.MS2Pairing) {
            ScanCollection MS2Scans = lCMSPeakBase.GetSpectrumParser().GetAllScanCollectionByMSLabel(false, true, false, true);
            MS2ScanPeakGroup = new ArrayList<>();
            ArrayList<ScanIsotopePeakDetectionThread> ResultList = new ArrayList<>();

            ExecutorService executorPool = null;
            executorPool = Executors.newFixedThreadPool(NoCPUs);
            Logger.getRootLogger().info("MS2 peak detection");
            for (ScanData scan : MS2Scans.ScanHashMap.values()) {
                //Detect isotope peak cluster for each MS2 scan
                ScanIsotopePeakDetectionThread unit = new ScanIsotopePeakDetectionThread(scan, parameter);
                ResultList.add(unit);
                executorPool.execute(unit);
            }

            executorPool.shutdown();
            try {
                executorPool.awaitTermination(Long.MAX_VALUE, TimeUnit.NANOSECONDS);
            } catch (InterruptedException e) {
                Logger.getRootLogger().info("interrupted..");
            }

            //Extract detected MS2 isotope peak clusters
            for (ScanIsotopePeakDetectionThread scan : ResultList) {
                MS2ScanPeakGroup.add(scan.scanpeak);
            }
        //detection.DetectSingleMZTraces(scans);
            //lCMSPeakBase.ExportPeakCurveResult();
        }
        
        //Initialize crosslinker peak pair finder
        CrosslinkerPepFinder finder = new CrosslinkerPepFinder(lCMSPeakBase, medianIntensity, linker);

        //Find all peak pairs
        finder.FindAllPairPeaks(NoCPUs);

        //Find peak pairs for MS2 scans
        finder.FindPairPeakMS2(MS2ScanPeakGroup, NoCPUs);

        Logger.getRootLogger().info("Exporting tables");
        FileWriter writer = new FileWriter(FilenameUtils.getFullPath(FilePath) + basename + "_112_pair.xls");
        writer.write("A_MW\tA_mz\tA_Charge\tA_startRT\tA_endRT\tA_intensity\tB_MW\tB_mz\tB_Charge\tB_startRT\tB_endRT\tB_intensity\tCorr\tppm\n");
        for (PeakPairFinder pair : finder.PairList) {
            if (pair.pairgroup != null && pair.pairgroup.GetBestPeakPair() != null) {
                CoElutePeak peakB = pair.pairgroup.GetBestPeakPair();
                writer.write(pair.pairgroup.lowMassPeak.NeutralMass() + "\t" + pair.pairgroup.lowMassPeak.TargetMz() + "\t" + pair.pairgroup.lowMassPeak.Charge + "\t" + pair.pairgroup.lowMassPeak.MonoIsotopePeak.StartRT() + "\t" + pair.pairgroup.lowMassPeak.MonoIsotopePeak.EndRT() + "\t" + pair.pairgroup.lowMassPeak.PeakHeight[0] + "\t" + peakB.peakpair.NeutralMass() + "\t" + peakB.peakpair.TargetMz() + "\t" + peakB.peakpair.Charge + "\t" + peakB.peakpair.MonoIsotopePeak.StartRT() + "\t" + peakB.peakpair.MonoIsotopePeak.EndRT() + "\t" + peakB.peakpair.PeakHeight[0] + "\t" + peakB.Correlation + "\t" + peakB.PPM + "\n");
            }
        }
        writer.close();

        writer = new FileWriter(FilenameUtils.getFullPath(FilePath) + basename + "_PrecNPair.xls");
        writer.write("LowA_MW\tLowA_mz\tLowA_Charge\tLowA_startRT\tLowA_endRT\tLowA_Intensity\t"
                + "LowB_MW\tLowB_mz\tLowB_Charge\tLowB_startRT\tLowB_endRT\tLowB_Intensity\t"
                + "HighA_MW\tHighA_mz\tHighA_Charge\tHighA_startRT\tHighA_endRT\tHighA_Intensity\t"
                + "HighB_MW\tHighB_mz\tHighB_Charge\tHighB_startRT\tHighB_endRT\tHighB_Intensity\t"
                + "LowHighCorr\tPrecMW\tPrecMz\tPrecCharge\tPrecStartRT\tPrecEndRT\tPrecIntensity\tTotalCorr\tPPM\n");
        for (PrecursorCrossPepFinder pair : finder.IntactPepList) {
            if (pair.LowMassPeakGroup.GetBestPeakPair() != null && pair.HighMassPeakGroup.GetBestPeakPair() != null && pair.PrecursorCrossPepPeaks != null && pair.LowHighPeakCorr > 0.2f && pair.MaxHighMassPeakCorr > 0.2f && pair.MaxLowMassPeakCorr > 0.2f) {
                CoElutePeak peakB = pair.GetBestPrecursor();
                writer.write(pair.LowMassPeakGroup.lowMassPeak.NeutralMass() + "\t" + pair.LowMassPeakGroup.lowMassPeak.TargetMz() + "\t" + pair.LowMassPeakGroup.lowMassPeak.Charge + "\t" + pair.LowMassPeakGroup.lowMassPeak.MonoIsotopePeak.StartRT() + "\t" + pair.LowMassPeakGroup.lowMassPeak.MonoIsotopePeak.EndRT() + "\t" + pair.LowMassPeakGroup.lowMassPeak.PeakHeight[0]
                        + "\t" + pair.LowMassPeakGroup.GetBestPeakPair().peakpair.NeutralMass() + "\t" + pair.LowMassPeakGroup.GetBestPeakPair().peakpair.TargetMz() + "\t" + pair.LowMassPeakGroup.GetBestPeakPair().peakpair.Charge + "\t" + pair.LowMassPeakGroup.GetBestPeakPair().peakpair.MonoIsotopePeak.StartRT() + "\t" + pair.LowMassPeakGroup.GetBestPeakPair().peakpair.MonoIsotopePeak.EndRT() + "\t" + pair.LowMassPeakGroup.GetBestPeakPair().peakpair.PeakHeight[0]
                        + "\t" + pair.HighMassPeakGroup.lowMassPeak.NeutralMass() + "\t" + pair.HighMassPeakGroup.lowMassPeak.TargetMz() + "\t" + pair.HighMassPeakGroup.lowMassPeak.Charge + "\t" + pair.HighMassPeakGroup.lowMassPeak.MonoIsotopePeak.StartRT() + "\t" + pair.HighMassPeakGroup.lowMassPeak.MonoIsotopePeak.EndRT() + "\t" + pair.HighMassPeakGroup.lowMassPeak.PeakHeight[0]
                        + "\t" + pair.HighMassPeakGroup.GetBestPeakPair().peakpair.NeutralMass() + "\t" + pair.HighMassPeakGroup.GetBestPeakPair().peakpair.TargetMz() + "\t" + pair.HighMassPeakGroup.GetBestPeakPair().peakpair.Charge + "\t" + pair.HighMassPeakGroup.GetBestPeakPair().peakpair.MonoIsotopePeak.StartRT() + "\t" + pair.HighMassPeakGroup.GetBestPeakPair().peakpair.MonoIsotopePeak.EndRT() + "\t" + pair.HighMassPeakGroup.GetBestPeakPair().peakpair.PeakHeight[0]
                        + "\t" + pair.LowHighPeakCorr + "\t" + peakB.peakpair.NeutralMass() + "\t" + peakB.peakpair.TargetMz() + "\t" + peakB.peakpair.Charge + "\t" + peakB.peakpair.MonoIsotopePeak.StartRT() + "\t" + peakB.peakpair.MonoIsotopePeak.EndRT() + "\t" + peakB.peakpair.PeakHeight[0] + "\t" + peakB.Correlation + "\t" + peakB.PPM + "\n");
            }
        }
        writer.close();

        if (parameter.MS2Pairing) {
            writer = new FileWriter(FilenameUtils.getFullPath(FilePath) + basename + "_MS2PrecNPair.xls");
            writer.write("ScanNum\tRT\tPrecursorMW\tPrecursorMZ\tCharge\t"
                    + "LowA_MW\tLowA_mz\tLowA_Charge\tLowA_Intensity\t"
                    + "LowB_MW\tLowB_mz\tLowB_Charge\tLowB_Intensity\t"
                    + "HighA_MW\tHighA_mz\tHighA_Charge\tHighA_Intensity\t"
                    + "HighB_MW\tHighB_mz\tHighB_Charge\tHighB_Intensity\n");
            for (MS2PeakPairFinder linkerfinder : finder.IntactPepListMS2) {
                if (!linkerfinder.peakPairs.isEmpty()) {
                    for (PeakPair pair : linkerfinder.peakPairs) {
                        writer.write(linkerfinder.ScanPeak.Scan.ScanNum + "\t" + linkerfinder.ScanPeak.Scan.RetentionTime + "\t" + linkerfinder.ScanPeak.Scan.PrecursorMass() + "\t" + linkerfinder.ScanPeak.Scan.PrecursorMz + "\t" + linkerfinder.ScanPeak.Scan.PrecursorCharge
                                + "\t" + pair.LowMassPeak.LowMassPeak.NeutralMass() + "\t" + pair.LowMassPeak.LowMassPeak.PrecursorMz() + "\t" + pair.LowMassPeak.LowMassPeak.Charge + "\t" + pair.LowMassPeak.LowMassPeak.GetPeakXYPointByPeakidx(0).getY()
                                + "\t" + pair.LowMassPeak.HighMassPeak.NeutralMass() + "\t" + pair.LowMassPeak.HighMassPeak.PrecursorMz() + "\t" + pair.LowMassPeak.HighMassPeak.Charge + "\t" + pair.LowMassPeak.HighMassPeak.GetPeakXYPointByPeakidx(0).getY()
                                + "\t" + pair.HighMassPeak.LowMassPeak.NeutralMass() + "\t" + pair.HighMassPeak.LowMassPeak.PrecursorMz() + "\t" + pair.HighMassPeak.LowMassPeak.Charge + "\t" + pair.HighMassPeak.LowMassPeak.GetPeakXYPointByPeakidx(0).getY()
                                + "\t" + pair.HighMassPeak.HighMassPeak.NeutralMass() + "\t" + pair.HighMassPeak.HighMassPeak.PrecursorMz() + "\t" + pair.HighMassPeak.HighMassPeak.Charge + "\t" + pair.HighMassPeak.HighMassPeak.GetPeakXYPointByPeakidx(0).getY()
                                + "\n");
                    }
                }
            }
            writer.close();
        }
         Logger.getRootLogger().info("Done");
    }
}
