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
package MSUmpire.PeakDataStructure;

import MSUmpire.SortedListLib.SortedList;
import java.util.Comparator;

/**
 *
 * @author Chih-Chiang Tsou <chihchiang.tsou@gmail.com>
 */
public class SortedCurveCollectionIntensity extends SortedList<PeakCurve> {
    private static final long serialVersionUID = 964997167414074601L;

    public SortedCurveCollectionIntensity() {
        super(new Comparator<PeakCurve>() {
            @Override
            public int compare(PeakCurve x, PeakCurve y) {
                if (x.ApexInt == y.ApexInt) {
                    return Float.compare(x.TargetMz, y.TargetMz);
                }
                return -Float.compare(x.ApexInt, y.ApexInt);
            }
        });
    }
}
