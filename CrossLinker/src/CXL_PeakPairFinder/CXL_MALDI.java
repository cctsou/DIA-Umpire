/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package CXL_PeakPairFinder;

import MSUmpire.BaseDataStructure.InstrumentParameter;
import MSUmpire.BaseDataStructure.ScanCollection;
import MSUmpire.LCMSPeakStructure.LCMSPeakBase;
import MSUmpire.PeptidePeakClusterDetection.PDHandlerMS1;
import MSUmpire.SpectrumParser.MALDIDataParser;
import MSUmpire.Utility.ConsoleLogger;
import java.io.FileWriter;
import java.io.IOException;;
import java.util.ArrayList;
import java.util.concurrent.ExecutionException;
import org.apache.log4j.Level;
import org.apache.log4j.Logger;

/**
 *
 * @author Chih-Chiang Tsou
 */
public class CXL_MALDI {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) throws IOException, InterruptedException, ExecutionException {
        ConsoleLogger cl = new ConsoleLogger();
        cl.SetConsoleLogger(Level.DEBUG);
        Logger logger = Logger.getRootLogger();
        //String FilePath = "F:\\Data\\CXL\\March 2015 samples\\MALDI MS Peaklists\\MALDI_35_SN\\";
        String FilePath = "F:\\Data\\CXL\\March 2015 samples\\MALDI MS Peaklists\\HSP-CHIP\\";
        int NoCPUs=10;
        InstrumentParameter parameter=null;
        LCMSPeakBase lCMSPeakBase = null;
        parameter = new InstrumentParameter(InstrumentParameter.InstrumentType.MALDI);
        MALDIDataParser parser = new MALDIDataParser(FilePath);
        parser.cycletime = 120f / 576f;
        parser.Parse();
        ArrayList<ScanCollection> scans = new ArrayList<>();
        parameter.IsoCorrThreshold=0.2f;
        parameter.MS1PPM=600;
        parameter.RemoveGroupedPeaks=true;
        parameter.RemoveGroupedPeaksCorr=0.7f;
        parameter.RemoveGroupedPeaksRTOverlap=0.7f;
        parameter.MinPeakPerPeakCurve=3;
        scans.add(parser.scanCollection);
        lCMSPeakBase = new LCMSPeakBase();
        lCMSPeakBase.parameter = parameter;
        lCMSPeakBase.ScanCollectionName = FilePath + "/HSP_CHIP.mzXML";
        lCMSPeakBase.ParentmzXMLName = FilePath + "/HSP_CHIP.mzXML";
        lCMSPeakBase.MinNoPeakCluster = 1;
        lCMSPeakBase.MaxNoPeakCluster = 3;
        lCMSPeakBase.StartCharge = 1;
        lCMSPeakBase.EndCharge = 1;
        lCMSPeakBase.CreatePeakFolder();
        lCMSPeakBase.ExportPeakClusterTable = true;
        lCMSPeakBase.ExportPeakCurveTable = true;
        PDHandlerMS1 detection = new PDHandlerMS1(lCMSPeakBase, NoCPUs, parameter.MS1PPM);
        detection.DetectPeakClusters(scans);
        //lCMSPeakBase.GenerateMZSortedClusterList(false);
        
        //detection.DetectSingleMZTraces(scans);
        //lCMSPeakBase.ExportPeakCurveResult();
        
        CrosslinkerPepFinder finder = new CrosslinkerPepFinder(lCMSPeakBase,0f);
        lCMSPeakBase.parameter.MS1PPM=200;
        lCMSPeakBase.parameter.RTtol=1f;
        finder.FindAllPairPeaks();
        
        FileWriter writer = new FileWriter(FilePath + "112_pair.xls");
        writer.write("A_MW\tA_startScan\tA_endScan\tA_intensity\tB_MW\tB_startScan\tB_endScan\tB_intensity\tCorr\n");
        for (PeakPairFinder pair : finder.PairList) {
            if (pair.pairgroup!=null && pair.pairgroup.highMassPeak != null) {
                for (CoElutePeak peakB : pair.pairgroup.highMassPeak.values()) {
                    writer.write(pair.pairgroup.lowMassPeak.NeutralMass() + "\t" + parser.scanCollection.GetScan(pair.pairgroup.lowMassPeak.MonoIsotopePeak.StartScan).MGFTitle + "\t" + parser.scanCollection.GetScan(pair.pairgroup.lowMassPeak.MonoIsotopePeak.EndScan).MGFTitle + "\t" + pair.pairgroup.lowMassPeak.PeakHeight[0] + "\t" + peakB.peakpair.NeutralMass() + "\t" + parser.scanCollection.GetScan(peakB.peakpair.MonoIsotopePeak.StartScan).MGFTitle + "\t" + parser.scanCollection.GetScan(peakB.peakpair.MonoIsotopePeak.EndScan).MGFTitle + "\t" + peakB.peakpair.PeakHeight[0] + "\t" + peakB.Correlation + "\n");
                }
            }
        }
        writer.close();
        writer = new FileWriter(FilePath + "PrecNPair.xls");
        writer.write("Low_MW\tLow_startScan\tLow_endScan\tLow_Intensity\t"
                + "High_MW\tHigh_startScan\tHigh_endScan\tHigh_Intensity\t"
                + "LowHighCorr\tPrecMW\tPrecStartScan\tPrecEndScan\tPrecIntensity\tTotalCorr\tPPM\n");
        for (PrecursorCrossPepFinder pair : finder.IntactPepList) {
            if (pair.PrecursorCrossPepPeaks!=null && pair.LowHighPeakCorr>0.5f && pair.MaxHighMassPeakCorr>0.5f && pair.MaxLowMassPeakCorr>0.5f) {
                for (CoElutePeak peakB : pair.PrecursorCrossPepPeaks.values()) {
                    writer.write(pair.LowMassPeakGroup.lowMassPeak.NeutralMass() + "\t" + parser.scanCollection.GetScan(pair.LowMassPeakGroup.lowMassPeak.MonoIsotopePeak.StartScan).MGFTitle + "\t" + parser.scanCollection.GetScan(pair.LowMassPeakGroup.lowMassPeak.MonoIsotopePeak.EndScan).MGFTitle + "\t" + pair.LowMassPeakGroup.lowMassPeak.PeakHeight[0] + "\t" + pair.HighMassPeakGroup.lowMassPeak.NeutralMass() + "\t" + parser.scanCollection.GetScan(pair.HighMassPeakGroup.lowMassPeak.MonoIsotopePeak.StartScan).MGFTitle + "\t" + parser.scanCollection.GetScan(pair.HighMassPeakGroup.lowMassPeak.MonoIsotopePeak.EndScan).MGFTitle + "\t" + pair.HighMassPeakGroup.lowMassPeak.PeakHeight[0] + "\t" + pair.LowHighPeakCorr + "\t" + peakB.peakpair.NeutralMass() + "\t" + parser.scanCollection.GetScan(peakB.peakpair.MonoIsotopePeak.StartScan).MGFTitle + "\t" + parser.scanCollection.GetScan(peakB.peakpair.MonoIsotopePeak.EndScan).MGFTitle + "\t" + peakB.peakpair.PeakHeight[0] + "\t" + peakB.Correlation + "\t" + peakB.PPM + "\n");
                }
            }
        }
        writer.close();
    }    
}
