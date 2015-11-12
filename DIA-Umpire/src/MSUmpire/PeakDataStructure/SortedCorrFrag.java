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

import ExtPackages.SortedListLib.SortedList;
import java.util.Comparator;

/**
 *
 * @author Chih-Chiang Tsou <chihchiang.tsou@gmail.com>
 */
public class SortedCorrFrag extends SortedList<PrecursorFragmentPairEdge> {
    private static final long serialVersionUID = -2283085827334564318L;
    
    public SortedCorrFrag() {
        super(new Comparator<PrecursorFragmentPairEdge>() {
            @Override
            public int compare(PrecursorFragmentPairEdge x, PrecursorFragmentPairEdge y) {
                if (x.Correlation == y.Correlation) {
                    return Float.compare(x.FragmentMz, y.FragmentMz);
                }
                return -Float.compare(x.Correlation, y.Correlation);
            }
        });
    }
}
