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

import java.io.Serializable;
import java.util.ArrayList;

/**
 * Precursor-fragment group class
 * @author Chih-Chiang Tsou <chihchiang.tsou@gmail.com>
 */
public class FragmentPeakGroup extends PeptideFragment implements Serializable{
    private static final long serialVersionUID = 736478436L;

    public ArrayList<Float> IntensityGroup = new ArrayList<>();
    public ArrayList<Float> CorrGroup = new ArrayList<>();
    public ArrayList<Float> PPMGroup = new ArrayList<>();
    public ArrayList<Float> ApexDeltaGroup = new ArrayList<>();
    public ArrayList<Float> RTOverlapPGroup = new ArrayList<>();
    
    private float AvgInt=-1f;
    
    public void ClearGroups(){
        GetAvgInt();
        IntensityGroup=null;
        CorrGroup=null;
        PPMGroup=null;
        ApexDeltaGroup=null;
        RTOverlapPGroup=null;                
    }
    
    public float GetAvgInt(){
        if(AvgInt==-1){
            AvgInt=0;
            for(float intensity : IntensityGroup){
                AvgInt+=intensity;
            }
            AvgInt/=IntensityGroup.size();
            AvgInt=(float) Math.sqrt(AvgInt);
        }        
        return AvgInt;
    }
    
    public String GetCorrString() {
        String output = "";
        for (float corr : CorrGroup) {
            output += corr + ";";
        }
        return output;
    }

    public String GetPPMString() {
        String output = "";
        for (float ppm : PPMGroup) {
            output += ppm + ";";
        }
        return output;
    }

    public String GetIntString() {
        String output = "";
        for (float intensity : IntensityGroup) {
            output += intensity + ";";
        }
        return output;
    }
    
    public String GetApexDeltaString() {
        String output = "";
        for (float delta : ApexDeltaGroup) {
            output += delta + ";";
        }
        return output;
    }
    
    public String GetRTOverlapString() {
        String output = "";
        for (float overlap : RTOverlapPGroup) {
            output += overlap + ";";
        }
        return output;
    }
   
    
    @Override
    public FragmentPeakGroup clone(){
        FragmentPeakGroup fragmentPeakGroup=new FragmentPeakGroup();
        fragmentPeakGroup.AvgInt=AvgInt;
        fragmentPeakGroup.Charge=Charge;
        fragmentPeakGroup.CorrGroup=CorrGroup;
        fragmentPeakGroup.IntensityGroup=IntensityGroup;
        fragmentPeakGroup.ApexDeltaGroup=ApexDeltaGroup;
        fragmentPeakGroup.RTOverlapPGroup=RTOverlapPGroup;
        fragmentPeakGroup.PPMGroup=PPMGroup;
        fragmentPeakGroup.ObservedMZ=ObservedMZ;
        fragmentPeakGroup.IonType=IonType;
        return fragmentPeakGroup;
    }
}
