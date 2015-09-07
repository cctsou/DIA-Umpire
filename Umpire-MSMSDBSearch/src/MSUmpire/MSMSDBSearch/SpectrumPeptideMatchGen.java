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
package MSUmpire.MSMSDBSearch;

import MSUmpire.BaseDataStructure.InstrumentParameter;
import MSUmpire.BaseDataStructure.ScanData;
import MSUmpire.BaseDataStructure.XYData;
import MSUmpire.BaseDataStructure.XYPointCollection;
import MSUmpire.PSMDataStructure.PepIonID;
import com.compomics.util.experiment.biology.Ion;
import com.compomics.util.experiment.biology.ions.ElementaryIon;
import com.compomics.util.experiment.biology.ions.PeptideFragmentIon;
import java.io.IOException;
import java.util.ArrayList;
import org.xmlpull.v1.XmlPullParserException;

/**
 *
 * @author Chih-Chiang Tsou <chihchiang.tsou@gmail.com>
 */
public class SpectrumPeptideMatchGen {

    public PepIonID pepIonID;
    public ScanData Scan;
    public ArrayList<MatchUnit> MatchList;
    public int HighBMatch = -1;
    public int HighYMatch = -1;
    public int BMatch = -1;
    public int YMatch = -1;
    public float MatchScore = -1f;
    DBSearchParam parameter;
    public ArrayList<XYData> MatchFragMzInt;//X: mz, Y: intensity

    public int TotalMatch() {
        return HighBMatch + HighYMatch;
    }

    public SpectrumPeptideMatchGen(PepIonID pepIonID, ScanData Scan, DBSearchParam parameter) throws XmlPullParserException, IOException {
        this.pepIonID = pepIonID;
        this.Scan = Scan;
        if (this.Scan == null) {
            this.Scan = new ScanData();
            this.Scan.RetentionTime = -1f;
        }
        this.parameter = parameter;
        if (pepIonID != null && Scan != null) {
            Matching();
        }
    }
    
    private void Matching() {
        HighBMatch = 0;
        HighYMatch = 0;
        BMatch = 0;
        YMatch = 0;

        MatchList = new ArrayList<>();
        MatchFragMzInt = new ArrayList<>();
        if (Scan == null || Scan.Data.isEmpty()) {
            return;
        }

        Scan.GenerateTopPeakScanData(parameter.TotalPeaks);

        //ModificationMatch mod=new ModificationMatch(pt, true, BMatch)
        double protonMass = ElementaryIon.proton.getTheoreticMass();
        for (Ion frag : pepIonID.GetFragments()) {
            if ((frag.getSubType() == PeptideFragmentIon.B_ION || frag.getSubType() == PeptideFragmentIon.Y_ION) && "".equals(frag.getNeutralLossesAsString())) {
                float targetmz = (float) frag.getTheoreticMz(1);
                float targetmz2 = (float) frag.getTheoreticMz(2);
                XYData closetPeak = Scan.GetHighestPeakInMzWindow(targetmz, parameter.FragPPM);
                if (closetPeak == null) {
                    closetPeak = Scan.GetHighestPeakInMzWindow(targetmz2, parameter.FragPPM);
                }

                if (closetPeak != null) {
                    if (closetPeak.getY() > (Scan.MaxY / 20)) {
                        MatchList.add(new MatchUnit(frag, closetPeak));
                        MatchFragMzInt.add(new XYData(closetPeak.getX(), closetPeak.getY()));
                        MatchScore += closetPeak.getY();
                        if (frag.getSubType() == PeptideFragmentIon.B_ION) {
                            BMatch++;
                        }
                        if (frag.getSubType() == PeptideFragmentIon.Y_ION) {
                            YMatch++;
                        }
                    }
                }

                //Match uisng toppeaks                
                XYData closetPeaktop = Scan.TopPeakScan.GetHighestPeakInMzWindow(targetmz, parameter.FragPPM);
                if (closetPeaktop == null) {
                    closetPeaktop = Scan.TopPeakScan.GetHighestPeakInMzWindow(targetmz2, parameter.FragPPM);
                }
                if (closetPeaktop != null) {
                    if (frag.getSubType() == PeptideFragmentIon.B_ION) {
                        HighBMatch++;
                    }
                    if (frag.getSubType() == PeptideFragmentIon.Y_ION) {
                        HighYMatch++;
                    }
                }
            }
        }
        //MatchScore/=Scan.MaxY;
    }

    public class MatchUnit {

        public Ion frag;
        public XYData peak;

        public MatchUnit(Ion frag, XYData peak) {
            this.frag = frag;
            this.peak = peak;
        }
    }
}
