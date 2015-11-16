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
import java.io.Serializable;
import java.util.ArrayList;
import java.util.HashMap;

/**
 *
 * @author Chih-Chiang Tsou <chihchiang.tsou@gmail.com>
 */
public class PepFragmentLib implements Serializable{
    private static final long serialVersionUID = 8372548274L;

    public String ModSequence;
    public String Sequence;
    public String ModificationString;
    public int Charge;
    public float MaxProbability = 0f;
    public float PrecursorMz;
    public HashMap<String, FragmentPeakGroup> FragmentGroups = new HashMap<>();
    public ArrayList<ModificationMatch> Modifications = new ArrayList<>();
    public ArrayList<Float> RetentionTime=new ArrayList<>();
    public float MS1Score;

    public String GetKey() {
        return ModSequence + "_" + Charge;
    }
    
    public HashMap<String, FragmentPeakGroup> CloneFragmentGroup(){
        HashMap<String, FragmentPeakGroup> NewFragmentGroups=new HashMap<>();
        for(FragmentPeakGroup frag : FragmentGroups.values()){
            FragmentPeakGroup newfrag=frag.clone();
            NewFragmentGroups.put(newfrag.GetFragKey(), newfrag);
        }
        return NewFragmentGroups;
    }

    public void AddFragments(ArrayList<FragmentPeak> FragmentPeaks) {
        float topintensity = 0f;
        for (FragmentPeak frag : FragmentPeaks) {
            if (frag.intensity > topintensity) {
                topintensity = frag.intensity;
            }
        }
        for (FragmentPeak fragment : FragmentPeaks) {
            if (!FragmentGroups.containsKey(fragment.GetFragKey())) {
                FragmentPeakGroup frag = new FragmentPeakGroup();
                frag.IonType = fragment.IonType;
                frag.FragMZ = fragment.FragMZ;  
                frag.Charge = fragment.Charge;
                FragmentGroups.put(fragment.GetFragKey(), frag);
            }
            FragmentGroups.get(fragment.GetFragKey()).CorrGroup.add(fragment.corr);
            FragmentGroups.get(fragment.GetFragKey()).IntensityGroup.add(fragment.intensity / topintensity);
            FragmentGroups.get(fragment.GetFragKey()).PPMGroup.add(fragment.ppm);
            FragmentGroups.get(fragment.GetFragKey()).ApexDeltaGroup.add(fragment.ApexDelta);
            FragmentGroups.get(fragment.GetFragKey()).RTOverlapPGroup.add(fragment.RTOverlapP);
        }
    }
}
