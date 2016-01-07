/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package crosslinker;

import MSUmpire.BaseDataStructure.ScanData;
import MSUmpire.BaseDataStructure.XYData;
import MSUmpire.SpectrumParser.MGFParser;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import org.apache.commons.io.FilenameUtils;
import MSUmpire.SpectrumParser.PKLScanParser;

/**
 *
 * @author Chih-Chiang Tsou
 */
public class CrossLinkerHandler {

    public HashMap<String, MALDI_Spot> ResultMap = new HashMap<>();
    public ArrayList<CrossLinkerPep> CrossLinkingPeps = new ArrayList<>();

    private void ReadFile(File folder) throws IOException {
        File[] listOfFiles = folder.listFiles();
        if (listOfFiles == null) {
            return;
        }
        ArrayList<CrossLinkerScanResult> Results = new ArrayList<>();
        for (int i = 0; i < listOfFiles.length; i++) {
            if (listOfFiles[i].isFile()) {
                String files = listOfFiles[i].getName();
                if (files.toLowerCase().endsWith(".mgf")) {
                    MGFParser mgf = new MGFParser(folder + "/" + files);
                    for (int scanNo : mgf.scanCollection.GetMS2DescendingArray()) {
                        ScanData Scan = mgf.scanCollection.GetScan(scanNo);
                        Scan.MsLevel = 2;
                        Scan.MGFTitle = Scan.MGFTitle.split("_")[Scan.MGFTitle.split("_").length - 2];
                        Scan.Normalization();
                        CrossLinkerScanResult scanresult = new CrossLinkerScanResult(Scan);
                        Results.add(scanresult);
                    }
                } else if (files.toLowerCase().endsWith(".pkl")) {
                    PKLScanParser pkl = new PKLScanParser(folder + "/" + files);
                    pkl.scan.MsLevel = 1;
                    pkl.scan.MGFTitle = pkl.scan.MGFTitle.split("_")[pkl.scan.MGFTitle.split("_").length - 2];
                    pkl.scan.Normalization();
                    CrossLinkerScanResult scanresult = new CrossLinkerScanResult(pkl.scan);
                    Results.add(scanresult);
                }
            } else {
                ReadFile(listOfFiles[i]);
            }
        }

        //Process scans
        for (CrossLinkerScanResult scanresult : Results) {
            scanresult.FindAllPairPeaks();
            scanresult.MatchPairMz();
            for (FragmentPair pairresult : scanresult.CrossLinkingEvidence) {
                for (CrossLinkerScanResult result : Results) {
                    for (int i = 0; i < 2; i++) {
                        for (int j = 0; j < 2; j++) {
                            if (Math.abs(pairresult.FragmentPair[i][j].getX() - result.Scan.PrecursorMz) < Parameter.Tolerance) {
                                scanresult.PossibleFragmentScans.add(result);
                            }
                        }
                    }
                }
            }
            for (FragmentPair pairresult : scanresult.CrossLinkingEvidence) {
                boolean found = false;
                for (CrossLinkerPep crosspep : CrossLinkingPeps) {
                    if (crosspep.AddPair(pairresult.FragmentPair[0][0].getX(), pairresult.FragmentPair[1][0].getX())) {
                        crosspep.AddEvidence(scanresult);
                        float intactint = 0f;
                        if (scanresult.MSlevel == 1) {
                            crosspep.AddIntensityMS1(pairresult.FragmentPair[0][0].getY(), pairresult.FragmentPair[1][0].getY(), pairresult.IntactCrossPeptide.getY());
                        }
                        if (scanresult.MSlevel == 2) {
                            crosspep.AddIntensityMS2(pairresult.FragmentPair[0][0].getY(), pairresult.FragmentPair[1][0].getY());
                        }
                        if (scanresult.Fragment199 != null) {
                            crosspep.Frag199int += scanresult.Fragment199.getY();
                        }
                        found = true;
                    }
                }
                if (!found) {
                    CrossLinkerPep newpep = new CrossLinkerPep();
                    newpep.AddPair(pairresult.FragmentPair[0][0].getX(), pairresult.FragmentPair[1][0].getX());
                    CrossLinkingPeps.add(newpep);
                    newpep.AddEvidence(scanresult);
                    if (scanresult.MSlevel == 1) {
                        newpep.AddIntensityMS1(pairresult.FragmentPair[0][0].getY(), pairresult.FragmentPair[1][0].getY(), pairresult.IntactCrossPeptide.getY());
                    }
                    if (scanresult.MSlevel == 2) {
                        newpep.AddIntensityMS2(pairresult.FragmentPair[0][0].getY(), pairresult.FragmentPair[1][0].getY());
                    }
                    if (scanresult.Fragment199 != null) {
                        newpep.Frag199int += scanresult.Fragment199.getY();
                    }
                }
            }
        }

        for (CrossLinkerScanResult scanresult : Results) {
            if (ResultMap.containsKey(scanresult.Scan.MGFTitle)) {
                MALDI_Spot spot = ResultMap.get(scanresult.Scan.MGFTitle);
                if (spot.MS1result != null && scanresult.MSlevel == 1) {
                    spot.MS1result = scanresult;
                }
                if (scanresult.MSlevel == 2) {
                    spot.MS2result.add(scanresult);
                }
            } else {
                MALDI_Spot newspot = new MALDI_Spot();
                newspot.SpotTag = scanresult.Scan.MGFTitle;
                if (scanresult.MSlevel == 1) {
                    newspot.MS1result = scanresult;
                } else {
                    newspot.MS2result.add(scanresult);
                }
                ResultMap.put(newspot.SpotTag, newspot);
            }
        }
    }

    public void Process() throws IOException {

        String FilePath = "F:\\Data\\CXL\\March 2015 samples\\MALDI MS Peaklists\\HSP-CHIP-Tau\\";
        File folder = new File(FilePath);
        //Read scans
        ReadFile(folder);

        //Export result
        FileWriter writer = new FileWriter(FilenameUtils.getFullPath(FilePath) + "/CrossLinkerResult.csv");

        writer.write("ScanTitle,MS1Or2,PrecursorMZ,Charge,RT,Fragment199,CrossLinkingEvidence,DeadEndEvidence,PossibleFragMS2Scans\n");
        for (MALDI_Spot spot : ResultMap.values()) {
            WriterResult(writer, spot.MS1result);
            for (CrossLinkerScanResult result : spot.MS2result) {
                WriterResult(writer, result);
            }
        }
        writer.close();

        FileWriter writer2 = new FileWriter(FilenameUtils.getFullPath(FilePath) + "/CrossLinkerPeptides2.csv");
        writer2.write("MZ_A,MZ_B,ConScanMS1,ConScanMS2,IntAMS1,IntBMS1,IntCMS1,IntAMS2,IntBMS2,Frag199Cnt,MS1_spots, MS2_spots, MSMS_A, MSMS_B\n");
        for (CrossLinkerPep linkerpep : CrossLinkingPeps) {
            linkerpep.CountConsScans();
            writer2.write(linkerpep.GetAmz() + "," + linkerpep.GetBmz() + "," + linkerpep.MS1ConsScanCnt + "," + linkerpep.MS2ConsScanCnt + "," + linkerpep.maxAintMS1 + "," + linkerpep.maxBintMS1 + "," + linkerpep.maxCintMS1 + "," + linkerpep.maxAintMS2 + "," + linkerpep.maxBintMS2 + "," + linkerpep.Frag199int + "," + linkerpep.ExportOrderedSpotTag() + "\n");
        }
        writer2.close();
//        for (MALDI_Spot spot : ResultMap.values()) {
//            for (CrossLinkerScanResult result : spot.MS2result) {
//                for (XYPoint[][] pairresult : result.CrossLinkingEvidence) {
//                    writer.write("(A1:" + pairresult[0][0].X + "_ " + pairresult[0][0].Y + "-A2:" + pairresult[0][1].X + "_ " + pairresult[0][1].Y + ")(B1:" + pairresult[1][0].X + "_ " + pairresult[1][0].Y + "-B2:" + pairresult[1][1].X + "_ " + pairresult[1][1].Y + ")");
//                }
//            }
//        }

    }

    private void WriterResult(FileWriter writer, CrossLinkerScanResult result) throws IOException {
        writer.write(result.Scan.MGFTitle + "," + result.Scan.MsLevel + "," + result.Scan.PrecursorMz + "," + result.Scan.PrecursorCharge + "," + result.Scan.RetentionTime + "," + String.valueOf(result.Fragment199 != null) + ",");
        for (FragmentPair pairresult : result.CrossLinkingEvidence) {
            writer.write("(A1:" + pairresult.FragmentPair[0][0].getX() + "_ " + pairresult.FragmentPair[0][0].getY() + "-A2:" + pairresult.FragmentPair[0][1].getX() + "_ " + pairresult.FragmentPair[0][1].getY() + ")(B1:" + pairresult.FragmentPair[1][0].getX() + "_ " + pairresult.FragmentPair[1][0].getY() + "-B2:" + pairresult.FragmentPair[1][1].getX() + "_ " + pairresult.FragmentPair[1][1].getY() + ")");
        }
        writer.write(",");
        for (XYData[] pairresult : result.DeadEndEvidence) {
            writer.write("(A1:" + pairresult[0].getX() + "_ " + pairresult[0].getY() + "-A2:" + pairresult[1].getX() + "_ " + pairresult[1].getY() + ")");
        }
        String ScanString = ",";

        for (CrossLinkerScanResult fragresult : result.PossibleFragmentScans) {
            ScanString += fragresult.Scan.MGFTitle + ";";
        }
        writer.write(ScanString + "\n");
    }
}
