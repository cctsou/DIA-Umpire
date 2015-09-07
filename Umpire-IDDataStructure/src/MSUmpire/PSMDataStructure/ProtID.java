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

import java.io.IOException;
import java.io.Serializable;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.PriorityQueue;
import org.apache.commons.lang.exception.ExceptionUtils;
import org.apache.log4j.Logger;
import org.xmlpull.v1.XmlPullParserException;

/**
 *
 * @author Chih-Chiang Tsou <chihchiang.tsou@gmail.com>
 */
public class ProtID implements Serializable{
    private static final long serialVersionUID = 863236278237L;

    public transient boolean IDByDBSearch=false;
    public ArrayList<String> IndisProteins;
    public transient ArrayList<String> IndisProtDes;
    private String AccNo;
    public String UniProtID;
    public String ProteinGroup;
    public String Description;
    public String GeneName;
    public int ProteinLength;
    public float Mass;
    public float Probability;
    public float GroupProb;//added 0828, needs to make transent for older serialization file
    public float MaxLocalPW; // added 0828, needs to make transent for older serialization file
    public float MaxIniProb = 0f;
    public String Sequence;
    public HashMap<String, PepIonID> PeptideID;
    public HashMap<String, PepIonID> ProtPeptideID;
    public ArrayList<String> TheoPeptides = null;
    public ArrayList<String> ProtPepSeq = null;

    public void AddDisProteins(String proID) {
        if (!IndisProteins.contains(proID)) {
            IndisProteins.add(proID);
        }
    }

    public void SetDescription(String Des) {
        Description = Des;
        if (Description.contains("ID=") && Description.contains("GeneName")) {
            UniProtID = Description.substring(Description.indexOf("ID=") + 4, Description.indexOf("GeneName") - 2);
        }
        if (Description.contains("GeneName=") && Description.contains("OtherGeneNames")) {
            GeneName = Description.substring(Description.indexOf("GeneName=") + 9, Description.indexOf("OtherGeneNames") - 2).trim();
        }
    }
    
    public ProtID CloneProtein() {
        ProtID newprotein = new ProtID();
        try {
            if (Sequence != null) {
                newprotein.SetSequence(Sequence);
            } else {
                Logger.getRootLogger().error("Sequence of protein:"+getAccNo()+" is null");
            }
        } catch (Exception ex) {
            Logger.getRootLogger().error(ExceptionUtils.getStackTrace(ex));
        }
        newprotein.Description = Description;
        newprotein.Mass = Mass;
        newprotein.setAccNo(getAccNo());
        newprotein.UniProtID = UniProtID;
        newprotein.IndisProteins.add(getAccNo());
        return newprotein;
    }

    public float GetAbundanceByFragment_IBAQ(int topNfragment) {
        float totalabundance = 0f;
        for (PepIonID peptide : PeptideID.values()) {
            totalabundance += peptide.GetPepAbundanceByTopFragments(topNfragment);
        }
        return totalabundance / TheoPeptides.size();
    }

    public float GetAbundanceByTopPepFrag(int toppep, int topfrag) {
        return GetAbundanceByTopPepFrag(toppep,topfrag,-1f);
    }
    public float GetAbundanceByTopPepFrag(int toppep, int topfrag, float pepweight) {
         if(PeptideID.isEmpty()){
            return 0;
        }
        PriorityQueue<Float> TopQueue = new PriorityQueue<>(PeptideID.size(), Collections.reverseOrder());
        for (PepIonID peptide : PeptideID.values()) {
            if (peptide.FilteringWeight > pepweight) {
                TopQueue.add(peptide.GetPepAbundanceByTopFragments(topfrag));
            }
        }
        float totalabundance = 0f;
        int num = Math.min(toppep, TopQueue.size());

        for (int i = 0; i < num; i++) {
            totalabundance += TopQueue.poll();
        }
        return totalabundance / num;
    }

    public float GetAbundanceByAllFragment_IBAQ() {
        return GetAbundanceByAllFragment_IBAQ(-1);
    }
    
    public float GetAbundanceByAllFragment_IBAQ(float pepweight) {
        float totalabundance = 0f;
        for (PepIonID peptide : PeptideID.values()) {
            if (peptide.FilteringWeight>pepweight) {
                totalabundance += peptide.GetPepAbundanceByAllFragments();
            }
        }
        return totalabundance / TheoPeptides.size();
    }

     public int GetSILAC_H_Count(String Mod) {
        int spc = 0;
        for (PepIonID pep : PeptideID.values()) {
            if (pep.ModSequence.contains(Mod)) {
                spc += pep.GetSpectralCount();
            }
        }
        return spc;
    }
     
    public int GetSILAC_L_Count(String AA, String Mod) {
        int spc = 0;
        for (PepIonID pep : PeptideID.values()) {
            if (pep.ModSequence.contains(AA) && !pep.ModSequence.contains(Mod)) {
                spc += pep.GetSpectralCount();
            }
        }
        return spc;
    }
    
    
    public int GetSpectralCount() {
        int spc = 0;
        for (PepIonID pep : PeptideID.values()) {
            spc += pep.GetSpectralCount();            
        }
        return spc;
    }
    
      public int GetFragCount() {
        int frag = 0;
        for (PepIonID pep : PeptideID.values()) {
            frag += pep.GetFragCount();
        }
        return frag;
    }

    public float GetAbundanceByMS1_IBAQ_SILAC_H(String Mod) {
        float totalabundance = 0f;
        for (PepIonID peptide : PeptideID.values()) {
            if (peptide.PeakHeight != null && peptide.ModSequence.contains(Mod)) {
                totalabundance += peptide.PeakHeight[0];
            }
        }
        return totalabundance / TheoPeptides.size();
    }

    public float GetAbundanceByMS1_IBAQ_SILAC_L(String AA, String Mod) {
        float totalabundance = 0f;
        for (PepIonID peptide : PeptideID.values()) {
            if (peptide.PeakHeight != null && peptide.ModSequence.contains(AA) && !peptide.ModSequence.contains(Mod)) {
                totalabundance += peptide.PeakHeight[0];
            }
        }
        return totalabundance / TheoPeptides.size();
    }
      
    public float GetAbundanceByMS1_IBAQ() {
        return GetAbundanceByMS1_IBAQ(-1f);
    }
    
    public float GetAbundanceByMS1_IBAQ(float pepweight) {
        float totalabundance = 0f;
        for (PepIonID peptide : PeptideID.values()) {
            if (peptide.PeakHeight != null && peptide.FilteringWeight>pepweight) {
                totalabundance += peptide.PeakHeight[0];
            }
        }
        return totalabundance / TheoPeptides.size();
    }

    public float GetAbundanceByMS1_TopN(int topN) {
        return GetAbundanceByMS1_TopN(topN,-1f);
    }
    public float GetAbundanceByMS1_TopN(int topN, float pepweight) {
        if(PeptideID.isEmpty()){
            return 0;
        }
        PriorityQueue<Float> TopQueue =  new PriorityQueue<>(PeptideID.size(), Collections.reverseOrder());
        for (PepIonID peptide : PeptideID.values()) {
            if (peptide.PeakHeight != null && peptide.FilteringWeight>pepweight) {
                TopQueue.add(peptide.PeakHeight[0]);
            }
        }

        float totalabundance = 0f;
        int num = Math.min(topN, TopQueue.size());

        for (int i = 0; i < num; i++) {
            totalabundance += TopQueue.poll();
        }
        return totalabundance / num;
    }

     public float GetAbundanceByTopAcrossSample(ArrayList<String> ProPep) {
        float totalabundance = 0f;
        if (ProPep != null) {
            for (PepIonID pepIonID : PeptideID.values()) {
                if (ProPep.contains(pepIonID.GetKey())) {
                    totalabundance += pepIonID.PeakHeight[0];
                }
            }
        }
        return totalabundance;
    }
     
     public float GetAbundanceByTopAcrossSampleTopN(ArrayList<String> ProPep, int TopN) {
        float totalabundance = 0f;
        if (ProPep != null) {
            for (PepIonID pepIonID : PeptideID.values()) {
                if (ProPep.contains(pepIonID.GetKey())) {
                    totalabundance += pepIonID.PeakHeight[0];
                }
            }
        }
        return totalabundance;
    }
          
    
    public float GetAbundanceByTopCorrFragAcrossSample(ArrayList<String> ProPep, HashMap<String, ArrayList<String>> PepFrag) {
        float totalabundance = 0f;
        int count=0;
        if (ProPep != null) {
            for (PepIonID pepIonID : PeptideID.values()) {
                if (ProPep.contains(pepIonID.GetKey())) {
                    totalabundance += pepIonID.GetPepAbundanceByTopCorrFragAcrossSample(PepFrag.get(pepIonID.GetKey()));                    
                    count++;
                }
            }
        }
//        if(count<2){
//            totalabundance=0;
//        }
        return totalabundance;
    }
        
    public ProtID() {
        PeptideID = new HashMap<>();
        IndisProteins = new ArrayList<>();
        IndisProtDes=new ArrayList<>();
        ProtPeptideID = new HashMap<>();
        ProtPepSeq = new ArrayList<>();
    }

    public void AddPeptideID(PepIonID pepID) {
        if (!PeptideID.containsKey(pepID.GetKey())) {
            PeptideID.put(pepID.GetKey(), pepID);
        }
    }

    public void SetSequence(String Seq) throws IOException, XmlPullParserException {
        Sequence = Seq;
        InsilicosDigestion(1, 6, 30);
    }

    public void InsilicosDigestion(int missedcleave, int minlength, int maxlength) throws XmlPullParserException, IOException {
        TheoPeptides = EnzymeManager.GetInstance().GetTrypsin().digest(Sequence, missedcleave, minlength, maxlength);        
        if (String.valueOf(Sequence.charAt(0)).equals("M")) {
            int mc=0;
            for (int i = 1; i < Sequence.length(); i++) {
                if (String.valueOf(Sequence.charAt(i)).equals("K") || String.valueOf(Sequence.charAt(i)).equals("R")) {
                    mc++;
                    if(mc>missedcleave){
                        return;
                    }
                    String pep=Sequence.substring(1,i+1);
                    if (pep.length()>=minlength && pep.length()<=maxlength && !TheoPeptides.contains(pep)) {                        
                        TheoPeptides.add(pep);
                    }
                }
            }
        }
    }

    public String GetGeneName() {
        if (Description.contains("GN=")) {
            return Description.split("GN=")[1].split(" ")[0];
        }
        return getAccNo();
    }

    public boolean IsDecoy(String decoytag) {
        if (getAccNo().startsWith(decoytag)) {
            return true;
        }
        return false;
    }

    /**
     * @return the AccNo
     */
    public String getAccNo() {
        //return AccNo.replace("|", "_");
        return AccNo;
    }

    /**
     * @param AccNo the AccNo to set
     */
    public void setAccNo(String AccNo) {
        this.AccNo = AccNo;
    }
    
    
    public void UpdateMaxIniProb() {
        MaxIniProb = 0f;
        for (PepIonID pepion : PeptideID.values()) {
            if (pepion.MaxProbability > MaxIniProb) {
                MaxIniProb = pepion.MaxProbability;
            }
        }
    }
    
    public void RemoveLowWeightPepID(float threshold) {
        ArrayList<PepIonID> removelist=new ArrayList<>();
        MaxIniProb=0f;
        for(PepIonID pepion: PeptideID.values()){
            if(pepion.FilteringWeight<threshold){
                removelist.add(pepion);
            }
        }
        for(PepIonID pepIonID : removelist){
            PeptideID.remove(pepIonID.GetKey());
        }
         for(PepIonID pepion: ProtPeptideID.values()){
            if(pepion.FilteringWeight<threshold){
                removelist.add(pepion);
            }
        }
        for(PepIonID pepIonID : removelist){
            ProtPeptideID.remove(pepIonID.GetKey());
        }
    }

}
