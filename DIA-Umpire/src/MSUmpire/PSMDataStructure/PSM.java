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

import com.compomics.util.experiment.identification.matches.ModificationMatch;
import com.compomics.util.general.IsotopicDistribution;
import com.compomics.util.protein.AASequenceImpl;
import com.compomics.util.protein.MolecularFormula;
import java.io.Serializable;
import java.util.ArrayList;

/**
 *
 * @author Chih-Chiang Tsou <chihchiang.tsou@gmail.com>
 */
public class PSM implements Serializable{
    private static final long serialVersionUID = 1098763624628L;

    public float Probability;
    public String Sequence;
    public int Charge;
    public String ModSeq = "";
    public String TPPModSeq = "";
    public int Rank;
    public int ScanNo;
    public ArrayList<String> ParentProtIDs;
    public String RawDataName;
    public String SpecNumber;
    public float RetentionTime=-1;
    public float NeighborMaxRetentionTime;
    public ArrayList<ModificationMatch> Modifications;
    public String PreAA;
    public String NextAA;
    public int MissedCleavage;
    public float MassError;
    public float hyperscore;
    public float nextscore;
    public float bscore;
    public float yscore;
    public float cscore;
    public float zscore;
    public float ascore;
    public float xscore;
    public float expect;
    public PepIonID pepIonID;
    public float LuciphorLFLR;
    public float LuciphorFLR;
    public float LuciphorScore;

    public String GetRawNameString() {
        return SpecNumber.substring(0, SpecNumber.indexOf("."));
    }

    public PSM() {
        Modifications = new ArrayList<>();
        ParentProtIDs = new ArrayList<>();
    }

    public String GetModificationString() {
        String ModificationString = "";
        for (ModificationMatch mod : Modifications) {
            ModificationString += mod.getTheoreticPtm() + "(" + mod.getModificationSite() + ");";
        }
        return ModificationString;
    }
    
    public String GetPepKey() {
        return ModSeq + "_" + Charge;
    }

    public void AddParentProtein(String protein) {
        if (!ParentProtIDs.contains(protein)) {
            ParentProtIDs.add(protein);
        }
    }

    public boolean IsDecoy(String decoytag) {
        boolean decoy = true;
        for (String pro : ParentProtIDs) {
            if (!pro.startsWith(decoytag)) {
                decoy = false;
            }
        }
        return decoy;
    }

    public float NeutralPrecursorMz() {
        return (NeutralPepMass + Charge * 1.00727f) / Charge;
    }

    public float GetObsrIsotopicMz(int pkdix) {
        //pkdix starts from 0
        return ObserPrecursorMz() + (float) (pkdix) / Charge;
    }

    public float NeutralPepMass;
    public float ObserPrecursorMass;

    public float ObserPrecursorMz() {
        return (ObserPrecursorMass + Charge * 1.00727f) / Charge;
    }

    private MolecularFormula GetMolecularFormula() {
        MolecularFormula formula = new MolecularFormula(GetAASequenceImpl());
        //formula.addMolecularFormula(GetModMolecularFormula());
        return formula;
    }
    
    AASequenceImpl AAimple;
    private AASequenceImpl GetAASequenceImpl() {
        if (AAimple == null) {
            AAimple = new AASequenceImpl(Sequence);
        }
        return AAimple;
    }
    
    

    float[] TheoIso;
    public float[] IsotopicDistrubtionRatio(int NoOfIsoPeaks) {
        if (TheoIso == null) {
            IsotopicDistribution calc = GetIsotopicDistribution();
            Double[] isopeak = calc.getPercMax();
            float firstPattern = (float) (double) isopeak[0];

            TheoIso = new float[NoOfIsoPeaks];
            for (int i = 0; i < NoOfIsoPeaks; i++) {
                TheoIso[i] = (float) (double) (isopeak[i] / firstPattern);
            }
        }
        return TheoIso;
    }

    IsotopicDistribution calc;
    private IsotopicDistribution GetIsotopicDistribution() {
        if (calc == null) {
            calc = new IsotopicDistribution(GetMolecularFormula());
            calc.calculate();
        }
        //AAimple.addModification(new ModificationImplementation(Modifications , Modifications, new HashMap(), 0));
        return calc;
    }
}
