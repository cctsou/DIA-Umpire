/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package crosslinker;

import MSUmpire.BaseDataStructure.ScanData;
import MSUmpire.BaseDataStructure.XYData;
import MSUmpire.BaseDataStructure.XYPointCollection;
import java.util.ArrayList;

/**
 *
 * @author Chih-Chiang Tsou
 */
public class CrossLinkerScanResult {

    public ScanData Scan;
    public ArrayList<XYData[]> PairList = new ArrayList<>();
    public XYData Fragment199;
    public ArrayList<FragmentPair> CrossLinkingEvidence = new ArrayList<>();
    public ArrayList<XYData[]> DeadEndEvidence = new ArrayList<>();
    public int DeadEndLinker = 0;
    public int MSlevel;
    public ArrayList<CrossLinkerScanResult> PossibleFragmentScans = new ArrayList<>();

    public CrossLinkerScanResult(ScanData scan) {
        Scan = scan;
        MSlevel = scan.MsLevel;
    }

    public void FindAllPairPeaks() {

        for (int i = 0; i < Scan.PointCount(); i++) {
            XYData peak = Scan.Data.get(i);
            if (peak.getX() < DC4.Arm + 147) {
                continue;
            }
            if (Math.abs(Scan.GetPoinByXCloset(DC4.DABCO + DC4.Arm + DC4.H2O + DC4.H).getX() - (DC4.DABCO + DC4.Arm + DC4.H2O + DC4.H)) < Parameter.Tolerance) {
                Fragment199 = Scan.GetPoinByXCloset(DC4.DABCO + DC4.Arm + DC4.H2O + DC4.H);
            }
            XYPointCollection peakrange = Scan.GetSubSetByXRange(peak.getX() + DC4.DABCO - Parameter.Tolerance, peak.getX() + DC4.DABCO + Parameter.Tolerance);
            for (int pk = 0; pk < peakrange.PointCount(); pk++) {
                if (Math.abs(peakrange.Data.get(pk).getX() - (peak.getX() + DC4.DABCO)) < Parameter.Tolerance) {
                    XYData[] pair = new XYData[2];
                    pair[0] = peak;
                    pair[1] = peakrange.Data.get(pk);
                    PairList.add(pair);
                }
            }
        }
    }

    public void MatchPairMz() {

        for (int pairindex = 0; pairindex < PairList.size(); pairindex++) {
            XYData[] pair = PairList.get(pairindex);
            for (int pairindex2 = pairindex + 1; pairindex2 < PairList.size(); pairindex2++) {
                XYData[] pair2 = PairList.get(pairindex2);
                float matchmz = -1f;
                XYData IntactPeps = null;
                if (MSlevel == 2) {
                    matchmz = Scan.PrecursorMz;
                } else if (MSlevel == 1) {
                    matchmz = Scan.GetPoinByXCloset(pair[0].getX() + pair2[0].getX() + DC4.DABCO - DC4.H).getX();
                    IntactPeps = Scan.GetPoinByXCloset(pair[0].getX() + pair2[0].getX() + DC4.DABCO - DC4.H);
                }
                if (Math.abs(matchmz - (pair[0].getX() + pair2[0].getX() + DC4.DABCO - DC4.H)) < Parameter.Tolerance) {
                    FragmentPair fragpair = new FragmentPair();
                    fragpair.FragmentPair = new XYData[2][2];
                    fragpair.FragmentPair[0] = pair;
                    fragpair.FragmentPair[1] = pair2;
                    fragpair.IntactCrossPeptide = IntactPeps;
                    CrossLinkingEvidence.add(fragpair);
                }
            }
            float matchmz = -1f;
            if (MSlevel == 2) {
                matchmz = Scan.PrecursorMz;
            } else if (MSlevel == 1) {
                matchmz = Scan.GetPoinByXCloset(pair[0].getX() + DC4.DABCO + DC4.H2O + DC4.Arm).getX();
            }
            if (Math.abs(matchmz - (pair[0].getX() + DC4.DABCO + DC4.H2O + DC4.Arm)) < Parameter.Tolerance) {
                DeadEndEvidence.add(pair);
            }
        }
    }
}
