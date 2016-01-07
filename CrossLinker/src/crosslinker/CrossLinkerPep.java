/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package crosslinker;

import MSUmpire.BaseDataStructure.ScanData;
import java.util.ArrayList;
import java.util.HashMap;

/**
 *
 * @author Chih-Chiang Tsou
 */
public class CrossLinkerPep {

    private float Asum;
    private int Acount;
    private float Bsum;
    private int Bcount;
    public float maxAintMS1;
    public float maxBintMS1;
    public float maxAintMS2;
    public float maxBintMS2;
    public float maxCintMS1;
    public float Frag199int;
    public HashMap<String, ArrayList<CrossLinkerScanResult>> MS1EvidenceList = new HashMap<>();
    public HashMap<String, ArrayList<CrossLinkerScanResult>> MS2EvidenceList = new HashMap<>();
    private ArrayList<ScanData> MS2ScanA = new ArrayList<>();
    private ArrayList<ScanData> MS2ScanB = new ArrayList<>();

    public boolean AddPair(float A, float B) {
        if ((Acount == 0 && Bcount == 0) || CheckPairMZ(A, B)) {
            Asum += A;
            Acount++;
            Bsum += B;
            Bcount++;
            return true;
        }
        return false;
    }

    public boolean CheckPairMZ(float A, float B) {
        if (Math.abs(B - GetBmz()) < Parameter.Tolerance && Math.abs(A - GetAmz()) < Parameter.Tolerance) {
            return true;
        }
        return false;
    }

    public float GetAmz() {
        return Asum / Acount;
    }

    public float GetBmz() {
        return Bsum / Bcount;
    }

    public void AddEvidence(CrossLinkerScanResult scanresult) {

        for (CrossLinkerScanResult ms2 : scanresult.PossibleFragmentScans) {
            if (Math.abs(ms2.Scan.PrecursorMz - GetAmz()) < Parameter.Tolerance) {
                MS2ScanA.add(ms2.Scan);
            }
            if (Math.abs(ms2.Scan.PrecursorMz - GetBmz()) < Parameter.Tolerance) {
                MS2ScanB.add(ms2.Scan);
            }
        }

        HashMap<String, ArrayList<CrossLinkerScanResult>> EvidenceList = null;
        if (scanresult.Scan.MsLevel == 1) {
            EvidenceList = MS1EvidenceList;
        } else if (scanresult.Scan.MsLevel == 2) {
            EvidenceList = MS2EvidenceList;
        }
        if (EvidenceList.containsKey(scanresult.Scan.MGFTitle)) {
            EvidenceList.get(scanresult.Scan.MGFTitle).add(scanresult);
        } else {
            ArrayList<CrossLinkerScanResult> scanResults = new ArrayList<>();
            scanResults.add(scanresult);
            EvidenceList.put(scanresult.Scan.MGFTitle, scanResults);
        }
    }
    public int MS2ConsScanCnt = 0;
    public int MS1ConsScanCnt = 0;

    public void CountConsScans() {
        boolean ms1flag = false;
        boolean ms2flag = false;
        for (char letter = 'A'; letter <= 'H'; letter++) {
            for (int i = 1; i <= 24; i++) {
                String tag = String.valueOf(letter) + String.format("%02d", i);
                if (MS1EvidenceList.containsKey(tag)) {
                    if (MS1ConsScanCnt == 0) {
                        MS1ConsScanCnt++;
                    }
                    if (ms1flag) {
                        MS1ConsScanCnt++;
                    }
                    ms1flag = true;
                } else {
                    ms1flag = false;
                }
                if (MS2EvidenceList.containsKey(tag)) {
                    if (MS2ConsScanCnt == 0) {
                        MS2ConsScanCnt++;
                    }
                    if (ms2flag) {
                        MS2ConsScanCnt++;
                    }
                    ms2flag = true;
                }
            }
        }
    }

    public void AddIntensityMS1(float Aint, float Bint, float Cint) {
        if (Aint > maxAintMS1) {
            maxAintMS1 = Aint;
        }
        if (Bint > maxBintMS1) {
            maxBintMS1 = Bint;
        }
        if (Cint > maxCintMS1) {
            maxCintMS1 = Cint;
        }
    }

    public void AddIntensityMS2(float Aint, float Bint) {
        if (Aint > maxAintMS2) {
            maxAintMS2 = Aint;
        }
        if (Bint > maxBintMS1) {
            maxBintMS2 = Bint;
        }
    }

    public String ExportOrderedSpotTag() {
        String MS1string = "";
        String MS2string = "";
        String MSMSSpotsA = "";
        String MSMSSpotsB = "";
        for (char letter = 'A'; letter <= 'H'; letter++) {
            for (int i = 1; i <= 24; i++) {
                String tag = String.valueOf(letter) + String.format("%02d", i);
                if (MS1EvidenceList.containsKey(tag)) {
                    MS1string += tag + ";";
                }
                if (MS2EvidenceList.containsKey(tag)) {
                    MS2string += tag + ";";
                }

                for (ScanData msms : MS2ScanA) {
                    if (msms.MGFTitle.equals(tag)) {
                        MSMSSpotsA += tag + ";";
                        break;
                    }
                }
                for (ScanData msms : MS2ScanB) {
                    if (msms.MGFTitle.equals(tag)) {
                        MSMSSpotsB += tag + ";";
                        break;
                    }
                }
            }
        }
        return MS1string + "," + MS2string + "," + MSMSSpotsA + "," + MSMSSpotsB;
    }
}
