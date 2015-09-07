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
package MSUmpire.QuantModule;

import MSUmpire.PSMDataStructure.PepIonID;
import MSUmpire.PSMDataStructure.LCMSID;
import MSUmpire.PSMDataStructure.ModStringConvert;
import MSUmpire.PSMDataStructure.ModificationInfo;
import com.compomics.util.experiment.identification.matches.ModificationMatch;
import java.sql.SQLException;
import java.util.ArrayList;

/**
 *
 * @author Chih-Chiang Tsou <chihchiang.tsou@gmail.com>
 */
public class LabelingPepIonMapping {

    private LCMSID LCMSID;
    private ArrayList<LabelingPair> LabelingPairs = new ArrayList<>();

    public class LabelingPair {

        private ModificationInfo light;
        private ModificationInfo heavy;

        public LabelingPair(ModificationInfo light, ModificationInfo heavy) {
            this.light = light;
            this.heavy = heavy;
        }
    }

    public LabelingPepIonMapping(LCMSID LCMSID) {
        this.LCMSID = LCMSID;
    }

    public void AddLabelingMod(ModificationInfo light, ModificationInfo heavy) {
        LabelingPairs.add(new LabelingPair(light, heavy));
    }

    ////////////////Need to check when one of the labeling is unmodified status......
    public void GenerateMappedPepIon() throws SQLException {
        for (PepIonID pepion : LCMSID.GetPepIonList().values()) {            
            for (LabelingPair modpair : LabelingPairs) {
                MapLabelingPair(pepion, modpair.light,modpair.heavy);
                MapLabelingPair(pepion, modpair.heavy,modpair.light);
            }
        }
    }

    private void MapLabelingPair(PepIonID pepion, ModificationInfo checkmod, ModificationInfo mapmod) {        
        boolean needToMap=true;
        if (mapmod.modification != null) {
            for (ModificationMatch mod : pepion.Modifications) {
                if (mod.getTheoreticPtm().equals(mapmod.modification.getName())) {
                    needToMap = false;
                    break;
                }
            }
        }
        ////////Whether it is a light peptide
        ////////light version is no-mod
        if (needToMap) {
            if (checkmod.massdiff == 0f) {
                needToMap = false;
                for (int i = 0; i < pepion.Sequence.length(); i++) {
                    String AA = String.valueOf(pepion.Sequence.charAt(i));
                    if (AA.equals(checkmod.site)) {
                        needToMap = true;
                        break;
                    }
                }
            } else {
                needToMap = false;
                for (ModificationMatch mod : pepion.Modifications) {
                    if (mod.getTheoreticPtm().equals(checkmod.modification.getName())) {
                        needToMap = true;
                        break;
                    }
                }
            }
            if (needToMap) {
                PepIonID predictedPepIon = pepion.ClonePepIonID();
                if (checkmod.modification != null) {
                    ArrayList<ModificationMatch> removelist = new ArrayList<>();
                    for (ModificationMatch mod : predictedPepIon.Modifications) {
                        if (mod.getTheoreticPtm().equals(checkmod.modification.getName())) {
                            removelist.add(mod);
                        }
                    }
                    for (ModificationMatch mod : removelist) {
                        predictedPepIon.Modifications.remove(mod);
                    }
                }
                predictedPepIon.ModSequence=predictedPepIon.ModSequence.replace(checkmod.GetKey(), "");                
                if (mapmod.modification != null) {
                    for (int i = 0; i < predictedPepIon.Sequence.length(); i++) {
                        String AA = String.valueOf(predictedPepIon.Sequence.charAt(i));
                        if (AA.equals(mapmod.site)) {
                            predictedPepIon.Modifications.add(new ModificationMatch(mapmod.modification.getName(), true, (i + 1)));
                            predictedPepIon.ModSequence = ModStringConvert.AddModIntoSeqBeforeSite(predictedPepIon.ModSequence, mapmod.GetKey(), i);
                        }
                    }
                }
                if (!LCMSID.GetPepIonList().containsKey(predictedPepIon.GetKey()) && !LCMSID.GetMappedPepIonList().containsKey(predictedPepIon.GetKey())) {
                    predictedPepIon.ReFreshPepFactory();
                    predictedPepIon.SetRT(pepion.GetIDRT());
                    LCMSID.GetMappedPepIonList().put(predictedPepIon.GetKey(), predictedPepIon);
                }
            }
        }
    }
}
