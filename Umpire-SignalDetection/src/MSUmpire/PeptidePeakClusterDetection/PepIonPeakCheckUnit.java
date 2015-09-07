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
package MSUmpire.PeptidePeakClusterDetection;

import MSUmpire.PSMDataStructure.PepIonID;
import java.util.ArrayList;

/**
 *
 * @author Chih-Chiang Tsou <chihchiang.tsou@gmail.com>
 */
public class PepIonPeakCheckUnit {

    public PepIonPeakCheckUnit(PepIonID pepIonID) {
        this.pepIonID = pepIonID;
    }
    public PepIonID pepIonID;
    public CheckPeakUnit MS1CheckUnit;
    public ArrayList<CheckPeakUnit> MS2CheckUnits = new ArrayList<>();

    public int FoundMS2Peaks() {
        int count = 0;
        if (MS2CheckUnits != null) {
            for (CheckPeakUnit ms2 : MS2CheckUnits) {
                if (ms2.Stage[2]) {
                    count++;
                }
            }
        }
        return count;
    }
}
