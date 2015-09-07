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
package MSUmpire.DIA;

import MSUmpire.BaseDataStructure.InstrumentParameter;
import MSUmpire.PeakDataStructure.PrecursorFragmentPairEdge;
import MSUmpire.UmpireSearchDataStructure.FragmentIntensity;
import MSUmpire.UmpireSearchDataStructure.PepIonCandidate;
import MSUmpire.UmpireSearchDataStructure.SearchUnit;
import com.compomics.util.experiment.biology.ions.ElementaryIon;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import MSUmpire.MSMSDBSearch.DBSearchParam;

/**
 *
 * @author Chih-Chiang Tsou <chihchiang.tsou@gmail.com>
 */
public class MS1PeakGroupSearch implements SearchUnit {

    public ArrayList<PrecursorFragmentPairEdge> GroupedFragments;
    DBSearchParam parameter;
    public HashMap<Integer, Double> FactorialTable;
    public int BMatch = 0;
    public int YMatch = 0;
    float Score;
    public float Bintensitysum = 0f;
    public float Bcorrsum = 0f;
    public float Yintensitysum = 0f;
    public float Ycorrsum = 0f;
    public int Bconverage;
    public int Yconverage;
    public float maxint = 0f;
    public PepIonCandidate candidate;

    public MS1PeakGroupSearch(HashMap<Integer, Double> FactorialTable) {
        this.FactorialTable = FactorialTable;
    }

    @Override
    public float CalScore() {
        //XTandem hyperscore
        Score = (float) ((Bintensitysum + Yintensitysum) * Math.pow(2, FactorialTable.get(BMatch) + FactorialTable.get(YMatch)));
        return Score;
    }

    @Override
    public void Matching(PepIonCandidate pepion, DBSearchParam parameter) {
        BMatch = 0;
        YMatch = 0;
        Bintensitysum = 0f;
        Bcorrsum = 0f;
        Yintensitysum = 0f;
        Ycorrsum = 0f;
        Bconverage = 0;
        Yconverage = 0;
        maxint = 0f;

        candidate = pepion;
        this.parameter = parameter;
        ArrayList<Float> FragmentKey = new ArrayList<>();
        for (PrecursorFragmentPairEdge frag : GroupedFragments) {
            FragmentKey.add(frag.FragmentMz);
            if (frag.Intensity > maxint) {
                maxint = frag.Intensity;
            }
        }

        //ModificationMatch mod=new ModificationMatch(pt, true, BMatch)
        double protonMass = ElementaryIon.proton.getTheoreticMass();
        FragmentIntensity[] Bions = pepion.GetBFragments();
        FragmentIntensity[] Yions = pepion.GetYFragments();
        Bconverage = 0;
        boolean gap = true;
        for (FragmentIntensity Bion : Bions) {
            float targetmz = (float) Bion.fragmentIon.getTheoreticMz(1);
            int closetidx = MatchBymz(targetmz, FragmentKey);
            if (closetidx == -1) {
                targetmz = (float) Bion.fragmentIon.getTheoreticMz(2);
                closetidx = MatchBymz(targetmz, FragmentKey);
            }
            if (closetidx != -1) {
                PrecursorFragmentPairEdge peak = GroupedFragments.get(closetidx);
                BMatch++;
                Bintensitysum += peak.Intensity;
                Bcorrsum += peak.Correlation;
                peak.MatchedFragMz = targetmz;
                if (!gap) {
                    Bconverage++;
                }
                gap = false;
            } else {
                gap = true;
            }
        }
        gap = true;
        for (FragmentIntensity Yion : Yions) {
            float targetmz = (float) (Yion.fragmentIon.getTheoreticMass() + protonMass);
            int closetidx = MatchBymz(targetmz, FragmentKey);
            if (closetidx == -1) {
                targetmz = (float) (Yion.fragmentIon.getTheoreticMass() + protonMass * 2) / 2;
                closetidx = MatchBymz(targetmz, FragmentKey);
            }
            if (closetidx != -1) {
                PrecursorFragmentPairEdge peak = GroupedFragments.get(closetidx);
                YMatch++;
                Yintensitysum += peak.Intensity;
                Ycorrsum += peak.Correlation;
                peak.MatchedFragMz = targetmz;
                if (!gap) {
                    Yconverage++;
                }
                gap = false;
            } else {
                gap = true;
            }
        }
    }

    private int MatchBymz(float targetmz, ArrayList<Float> FragmentKey) {
        float lowmz = InstrumentParameter.GetMzByPPM(targetmz, 1, parameter.FragPPM);
        int startidx = Collections.binarySearch(FragmentKey, lowmz);

        if (startidx < 0) {
            startidx = -(startidx + 1);
        }
        if (startidx > 0) {
            startidx--;
        }
        Float closetPeak = null;
        int returnidx = -1;
        for (int idx = startidx; idx < FragmentKey.size(); idx++) {
            Float peakmz = FragmentKey.get(idx);
            if (InstrumentParameter.CalcPPM(targetmz, peakmz) <= parameter.FragPPM) {
                if (closetPeak == null || peakmz > closetPeak) {
                    closetPeak = peakmz;
                    returnidx = idx;
                }
            } else if (peakmz > targetmz) {
                break;
            }
        }
        return returnidx;
    }
}
