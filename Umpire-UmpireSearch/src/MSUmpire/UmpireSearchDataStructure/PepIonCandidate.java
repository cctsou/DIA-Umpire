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
package MSUmpire.UmpireSearchDataStructure;

import com.compomics.util.experiment.biology.Ion;
import com.compomics.util.experiment.biology.IonFactory;
import com.compomics.util.experiment.biology.Peptide;
import com.compomics.util.experiment.biology.ions.ElementaryIon;
import com.compomics.util.experiment.biology.ions.PeptideFragmentIon;
import com.compomics.util.experiment.identification.matches.ModificationMatch;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.concurrent.locks.Lock;
import java.util.concurrent.locks.ReentrantLock;

/**
 *
 * @author Chih-Chiang Tsou <chihchiang.tsou@gmail.com>
 */
public class PepIonCandidate {

    //public String Sequence;
    public int Charge;
    private float mz = -1f;
    public Peptide peptide = null;
    public float RT;
    public String Sequence;

    private FragmentIntensity[] BFragments;
    private FragmentIntensity[] YFragments;
    public String ModString;

    public void SetPeptide() {
        peptide = new Peptide(Sequence, new ArrayList<ModificationMatch>());
    }

    public void SetMod() {
        if (ModString != null && !ModString.isEmpty()) {
            String[] modifications = ModString.split(";");
            for (int i = 0; i < modifications.length; i++) {
                String name = modifications[i].split("@")[0];
                int site = Integer.parseInt(modifications[i].split("@")[1]);
                peptide.addModificationMatch(new ModificationMatch(name, true, site));
            }
        }
    }

    public String GetModString() {
        ModString = "";
        for (ModificationMatch mod : peptide.getModificationMatches()) {
            ModString += mod.getTheoreticPtm() + "@" + mod.getModificationSite() + ";";
        }
        return ModString;
    }

    public void SetMz(float mz) {
        this.mz = mz;
    }

    public boolean DecoyPep() {
        for (String prot : peptide.getParentProteinsNoRemapping()) {
            if (!prot.startsWith("rev_") && !prot.endsWith("_REV") && !prot.endsWith("REVERSED")) {
                return false;
            }
        }
        return true;
    }

    public float NeutralPrecursorMz() {
        if (mz == -1f) {
            mz = (CalcNeutralPepMass() + Charge * (float) ElementaryIon.proton.getTheoreticMass()) / Charge;
        }
        return mz;
    }

    public float UpdateNeutralPrecursorMz() {
        mz = (CalcNeutralPepMass() + Charge * (float) ElementaryIon.proton.getTheoreticMass()) / Charge;
        return mz;
    }

    public float CalcNeutralPepMass() {
        return peptide.getMass().floatValue();
    }

    Lock fraglock = new ReentrantLock(true);
    boolean SetFragmentFinished = false;

    public FragmentIntensity[] GetBFragments() {
        if (!SetFragmentFinished) {
            fraglock.lock();
            try {
                if (!SetFragmentFinished) {
                    SetFragments();
                    SetFragmentFinished = true;
                }
            } finally {
                fraglock.unlock();
            }
        }
        return BFragments;
    }

    public FragmentIntensity[] GetYFragments() {
        if (!SetFragmentFinished) {
            fraglock.lock();
            try {
                if (!SetFragmentFinished) {
                    SetFragments();
                    SetFragmentFinished = true;
                }
            } finally {
                fraglock.unlock();
            }
        }
        return YFragments;
    }

    public void ReleaseFragments() {
        fraglock.lock();
        try {
            SetFragmentFinished = false;
            BFragments = null;
            YFragments = null;
        } finally {
            fraglock.unlock();
        }
    }

    private void SetFragments() {                        
        HashMap<Integer, HashMap<Integer, ArrayList<Ion>>> allfragment = IonFactory.getInstance().getFragmentIons(peptide);                
        ArrayList<Ion> fragments=new ArrayList<>();
        fragments.addAll(allfragment.get(Ion.IonType.PEPTIDE_FRAGMENT_ION.index).get(PeptideFragmentIon.B_ION));
        fragments.addAll(allfragment.get(Ion.IonType.PEPTIDE_FRAGMENT_ION.index).get(PeptideFragmentIon.Y_ION));
        
        BFragments = new FragmentIntensity[peptide.getSequence().length() - 1];
        YFragments = new FragmentIntensity[peptide.getSequence().length() - 1];
        for (Ion frag : fragments) {            
            if ("".equals(frag.getNeutralLossesAsString())) {
                if (frag.getSubType() == PeptideFragmentIon.B_ION) {
                    int index = ((PeptideFragmentIon) frag).getNumber() - 1;
                    BFragments[index] = new FragmentIntensity();
                    BFragments[index].fragmentIon = frag;
                } else if (frag.getSubType() == PeptideFragmentIon.Y_ION) {
                    int index = ((PeptideFragmentIon) frag).getNumber() - 1;
                    YFragments[index] = new FragmentIntensity();
                    YFragments[index].fragmentIon = frag;
                }
            }
        }
    }
}
