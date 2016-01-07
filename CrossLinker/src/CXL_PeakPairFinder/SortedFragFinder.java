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
package CXL_PeakPairFinder;

import ExternalPackages.SortedListLib.SortedList;
import java.io.Serializable;
import java.util.Comparator;

/**
 *
 * @author Chih-Chiang Tsou <chihchiang.tsou@gmail.com>
 */
public class SortedFragFinder extends SortedList<PeakPairFinder> implements Serializable{

    public SortedFragFinder() {
        super(new Comparator<PeakPairFinder>() {
            @Override
            public int compare(PeakPairFinder x, PeakPairFinder y) {                
                return Float.compare(x.pairgroup.lowMassPeak.NeutralMass(), y.pairgroup.lowMassPeak.NeutralMass());
            }
        });
    }

    public int BinarySearchLower(float value) {
        if (isEmpty()) {
            return 0;
        }
        int lower = 0;
        int upper = size() - 1;

        if (value - get(upper).pairgroup.lowMassPeak.NeutralMass() >= 0) {
            return upper;
        }
        if (value - get(0).pairgroup.lowMassPeak.NeutralMass() <= 0) {
            return 0;
        }

        while (lower <= upper) {
            int middle = (lower + upper) / 2;
            float comparisonResult = value - get(middle).pairgroup.lowMassPeak.NeutralMass();
            if (comparisonResult == 0) {
                return middle;
            } else if (comparisonResult < 0) {
                upper = middle - 1;
            } else {
                lower = middle + 1;
            }
        }
        if (upper < 0) {
            return 0;
        }
        while (upper > 0 && get(upper).pairgroup.lowMassPeak.NeutralMass() >= value) {
            upper--;
        }
        return upper;
    }

    public int BinarySearchHigher(float value) {
        if (isEmpty()) {
            return 0;
        }
        int lower = 0;
        int upper = size() - 1;

        if (value - get(upper).pairgroup.lowMassPeak.NeutralMass() >= 0) {
            return upper;
        }
        if (value - get(0).pairgroup.lowMassPeak.NeutralMass() <= 0) {
            return 0;
        }

        while (lower <= upper) {
            int middle = (lower + upper) / 2;
            float comparisonResult = value - get(middle).pairgroup.lowMassPeak.NeutralMass();
            if (comparisonResult == 0) {
                return middle;
            } else if (comparisonResult < 0) {
                upper = middle - 1;
            } else {
                lower = middle + 1;
            }
        }
        if (lower > size() - 1) {
            return size() - 1;
        }
        while (upper < size() && get(upper).pairgroup.lowMassPeak.NeutralMass() <= value) {
            upper++;
        }
        return upper;
    }

    public int BinarySearchClosest(float value) {
        if (isEmpty()) {
            return 0;
        }
        int lower = 0;
        int upper = size() - 1;

        if (value - get(upper).pairgroup.lowMassPeak.NeutralMass() >= 0) {
            return upper;
        }
        if (value - get(0).pairgroup.lowMassPeak.NeutralMass() <= 0) {
            return 0;
        }

        while (lower <= upper) {
            int middle = (lower + upper) / 2;
            float comparisonResult = value - get(middle).pairgroup.lowMassPeak.NeutralMass();
            if (comparisonResult == 0) {
                return middle;
            } else if (comparisonResult < 0) {
                upper = middle - 1;
            } else {
                lower = middle + 1;
            }
        }

        if (Math.abs(value - get(lower).pairgroup.lowMassPeak.NeutralMass()) > Math.abs(value - get(upper).pairgroup.lowMassPeak.NeutralMass())) {
            return upper;
        } else {
            return lower;
        }
    }
}
