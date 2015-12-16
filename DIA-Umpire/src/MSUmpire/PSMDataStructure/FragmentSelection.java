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
package MSUmpire.PSMDataStructure;

import java.util.ArrayList;
import java.util.HashMap;
import org.apache.commons.io.FilenameUtils;

/**
 * Fragment and peptide selection algorithms (as described in DIA-Umpire paper) 
 * @author Chih-Chiang Tsou <chihchiang.tsou@gmail.com>
 */
public class FragmentSelection {

    ArrayList<LCMSID> FileList;
    public HashMap<String, HashMap<String, Float>> PepFragScore = new HashMap<>();
    public HashMap<String, HashMap<String, Float>> ProtPepScore = new HashMap<>();
    public HashMap<String, ArrayList<String>> TopFrags = new HashMap<>();
    public HashMap<String, ArrayList<String>> TopPeps = new HashMap<>();
    public HashMap<String, HashMap<String, Float>> FragInt = new HashMap<>();
    public float freqPercent = 0f;
    public float CorrelationThreshold=0.6f;
    public int NoConsecutiveRun = 3;
    public Scoring scoring=Scoring.IntensityCorr;

    public enum Scoring{
        IntensityCorr,
        Intensity,
    }
    public FragmentSelection(ArrayList<LCMSID> FileList) {
        this.FileList = FileList;
    }

    public void GeneratePepFragScoreMap() {
        for (LCMSID IDSummary : FileList) {
            for (String key : IDSummary.GetPepIonList().keySet()) {
                if (!PepFragScore.containsKey(key)) {
                    PepFragScore.put(key, new HashMap<String, Float>());
                }
            }
            for (String key : IDSummary.GetMappedPepIonList().keySet()) {
                if (!PepFragScore.containsKey(key)) {
                    PepFragScore.put(key, new HashMap<String, Float>());
                }
            }
        }

        for (String PepKey : PepFragScore.keySet()) {
            int IDNo = 0;
            HashMap<String, Float> fragmentscore = new HashMap<>();
            HashMap<String, Integer> fragmentfreq = new HashMap<>();

            for (LCMSID IDSummary : FileList) {
                PepIonID pep = null;
                if (IDSummary.GetPepIonList().containsKey(PepKey)) {
                    pep = IDSummary.GetPepIonList().get(PepKey);
                } else if (IDSummary.GetMappedPepIonList().containsKey(PepKey)) {
                    pep = IDSummary.GetMappedPepIonList().get(PepKey);
                }
                if (pep != null) {
                    IDNo++;
                    for (FragmentPeak frag : pep.FragmentPeaks) {
                        if (frag.corr > CorrelationThreshold) {
                            if (!fragmentscore.containsKey(frag.GetFragKey())) {
                                fragmentscore.put(frag.GetFragKey(), 0f);
                                fragmentfreq.put(frag.GetFragKey(), 0);
                            }
                            if (scoring == Scoring.Intensity) {
                                fragmentscore.put(frag.GetFragKey(), fragmentscore.get(frag.GetFragKey()) + frag.intensity);
                            } else if (scoring == Scoring.IntensityCorr) {
                                fragmentscore.put(frag.GetFragKey(), fragmentscore.get(frag.GetFragKey()) + frag.corr * frag.intensity);
                            }
                            fragmentfreq.put(frag.GetFragKey(), fragmentfreq.get(frag.GetFragKey()) + 1);

                            if (!FragInt.containsKey(pep.GetKey() + "_" + frag.GetFragKey())) {
                                FragInt.put(pep.GetKey() + "_" + frag.GetFragKey(), new HashMap<String, Float>());
                            }
                            FragInt.get(pep.GetKey() + "_" + frag.GetFragKey()).put(FilenameUtils.getBaseName(IDSummary.Filename), frag.intensity);
                        }
                    }
                }
            }
            for (String fragkey : fragmentfreq.keySet()) {
                if (fragmentfreq.get(fragkey) > IDNo * freqPercent) {
                    PepFragScore.get(PepKey).put(fragkey, fragmentscore.get(fragkey));
                }
            }
        }
    }

    public void GenerateProtPepScoreMap(float pepweight) {
        for (LCMSID IDSummary : FileList) {
            for (String key : IDSummary.ProteinList.keySet()) {
                if (!ProtPepScore.containsKey(key)) {
                    ProtPepScore.put(key, new HashMap<String, Float>());
                }
            }
        }

        for (String ProteinKey : ProtPepScore.keySet()) {
            HashMap<String, Float> pepscore = new HashMap<>();
            HashMap<String, Integer> pepfreq = new HashMap<>();

            int IDNo = 0;
            for (LCMSID IDSummary : FileList) {
                if (IDSummary.ProteinList.containsKey(ProteinKey)) {
                    ProtID protein = IDSummary.ProteinList.get(ProteinKey);
                    if (protein.IDByDBSearch) {
                        IDNo++;
                    }
                    for (PepIonID pep : protein.PeptideID.values()) {
                        if (pep.FilteringWeight > pepweight) {
                            if (!pepscore.containsKey(pep.GetKey())) {
                                pepscore.put(pep.GetKey(), 0f);
                                pepfreq.put(pep.GetKey(), 0);
                            }
                            pepscore.put(pep.GetKey(), pepscore.get(pep.GetKey()) + pep.GetPepAbundanceByTopCorrFragAcrossSample(TopFrags.get(pep.GetKey())));
                            pepfreq.put(pep.GetKey(), pepfreq.get(pep.GetKey()) + 1);
                        }
                    }
                }
            }
            for (String pepkey : pepfreq.keySet()) {
                if (pepfreq.get(pepkey) > IDNo * freqPercent) {
                    ProtPepScore.get(ProteinKey).put(pepkey, pepscore.get(pepkey));
                }
            }
        }
    }

    public void FillMissingFragScoreMap() {
        for (String PepKey : PepFragScore.keySet()) {
            if (PepFragScore.get(PepKey).isEmpty()) {
                HashMap<String, Float> fragmentscore = new HashMap<>();
                for (LCMSID IDSummary : FileList) {
                    PepIonID pep = null;
                    if (IDSummary.GetPepIonList().containsKey(PepKey)) {
                        pep = IDSummary.GetPepIonList().get(PepKey);
                    } else if (IDSummary.GetMappedPepIonList().containsKey(PepKey)) {
                        pep = IDSummary.GetMappedPepIonList().get(PepKey);
                    }
                    if (pep != null) {
                        for (FragmentPeak frag : pep.FragmentPeaks) {
                            if (!fragmentscore.containsKey(frag.GetFragKey())) {
                                fragmentscore.put(frag.GetFragKey(), 0f);
                            }
                            fragmentscore.put(frag.GetFragKey(), fragmentscore.get(frag.GetFragKey()) + frag.corr * frag.intensity);

                            if (!FragInt.containsKey(pep.GetKey() + "_" + frag.GetFragKey())) {
                                FragInt.put(pep.GetKey() + "_" + frag.GetFragKey(), new HashMap<String, Float>());
                            }
                            FragInt.get(pep.GetKey() + "_" + frag.GetFragKey()).put(FilenameUtils.getBaseName(IDSummary.Filename), frag.intensity);
                        }
                    }
                }
                for (String fragkey : fragmentscore.keySet()) {
                    PepFragScore.get(PepKey).put(fragkey, fragmentscore.get(fragkey));
                }
            }
        }
    }

    public void GenerateTopFragMap(int topNFrag) {
        for (String PepKey : PepFragScore.keySet()) {
            HashMap<String, Float> Frags = PepFragScore.get(PepKey);
            ArrayList<String> frags = new ArrayList<>();
            for (int i = 0; i < topNFrag; i++) {
                float bestscore = 0f;
                String bestfrag = "";
                for (String frag : Frags.keySet()) {
                    if (!frags.contains(frag) && Frags.get(frag) > bestscore) {
                        bestscore = Frags.get(frag);
                        bestfrag = frag;
                    }
                }
                if (!"".equals(bestfrag)) {
                    frags.add(bestfrag);
                }
            }
            TopFrags.put(PepKey, frags);
        }
    }

    public void GenerateTopPepMap(int topNPep) {
        for (String ProteinKey : ProtPepScore.keySet()) {
            HashMap<String, Float> Peptides = ProtPepScore.get(ProteinKey);
            ArrayList<String> peps = new ArrayList<>();
            for (int i = 0; i < topNPep; i++) {
                float bestscore = 0f;
                String bestpep = "";
                for (String pep : Peptides.keySet()) {
                    if (!peps.contains(pep) && Peptides.get(pep) > bestscore) {
                        bestscore = Peptides.get(pep);
                        bestpep = pep;
                    }
                }
                if (!"".equals(bestpep)) {
                    peps.add(bestpep);
                }
            }
            TopPeps.put(ProteinKey, peps);
        }
    }
  
}
