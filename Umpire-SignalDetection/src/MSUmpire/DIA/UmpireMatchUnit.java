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

import MSUmpire.PeakDataStructure.PeakCluster;
import MSUmpire.PeakDataStructure.PrecursorFragmentPairEdge;
import MSUmpire.UmpireSearchDataStructure.PepIonCandidate;
import MSUmpire.UmpireSearchDataStructure.PepIonLib;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.logging.Level;
import java.util.logging.Logger;
import MSUmpire.MSMSDBSearch.DBSearchParam;
import org.xmlpull.v1.XmlPullParserException;

/**
 *
 * @author Chih-Chiang Tsou <chihchiang.tsou@gmail.com>
 */
public class UmpireMatchUnit implements Runnable {

    public PeakCluster ms1cluster;
    ArrayList<PrecursorFragmentPairEdge> fragments;
    PepIonLib IonLib;
    DBSearchParam searchpara;
    public ArrayList<MS1PeakGroupSearch> Ranking = new ArrayList<>();
    HashMap<Integer, Double> FactorialTable;
    public ArrayList<PepIonCandidate> candidates;

    public UmpireMatchUnit(PeakCluster ms1cluster, ArrayList<PrecursorFragmentPairEdge> fragments, PepIonLib IonLib, DBSearchParam searchpara, HashMap<Integer, Double> FactorialTable) {
        this.ms1cluster = ms1cluster;
        this.fragments = fragments;
        this.IonLib = IonLib;
        this.searchpara = searchpara;
        this.FactorialTable = FactorialTable;
    }

    @Override
    public void run() {
        try {
            candidates = IonLib.GetCandidates(ms1cluster.TargetMz(), ms1cluster.Charge, searchpara.PrecursorPPM);
            for (PepIonCandidate candidate : candidates) {
                MS1PeakGroupSearch search = new MS1PeakGroupSearch(FactorialTable);
                search.GroupedFragments = fragments;

                search.Matching(candidate, searchpara);
                if (search.BMatch + search.YMatch > 0) {
                    Ranking.add(search);
                }
//                if(ms1cluster.Index==210 && candidate.peptide.getSequence().equals(ms1cluster.AssignedPepIon)){
//                    search.Matching(candidate, searchpara);
//                }
            }
            Collections.sort(Ranking, new Comparator<MS1PeakGroupSearch>() {
                @Override
                public int compare(MS1PeakGroupSearch o1, MS1PeakGroupSearch o2) {
                    return -Float.compare(o1.CalScore(), o2.CalScore());
                }
            });
            
            ///remove low ranking candidates
            for (int i = Ranking.size() - 1; i > 4; i--) {
                Ranking.remove(i);
            }
        } catch (IOException ex) {
            Logger.getLogger(UmpireMatchUnit.class.getName()).log(Level.SEVERE, null, ex);
        } catch (XmlPullParserException ex) {
            Logger.getLogger(UmpireMatchUnit.class.getName()).log(Level.SEVERE, null, ex);
        }

    }

}
