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
import java.io.Serializable;
import java.util.Comparator;

/**
 *
 * @author Chih-Chiang Tsou <chihchiang.tsou@gmail.com>
 */
public class SortedCurveCollectionMZ extends SortedList<PeakCurve> implements Serializable {
    private static final long serialVersionUID = 1506046406460464L;

    public SortedCurveCollectionMZ() {
        super(new Comparator<PeakCurve>() {

            @Override
            public int compare(PeakCurve x, PeakCurve y) {
                if (x.TargetMz == y.TargetMz) {
                    return Float.compare(x.StartRT(), y.StartRT());
                }
                return Float.compare(x.TargetMz, y.TargetMz);
            }
        });
    }

    public int BinarySearchLower(float value) {
        if (isEmpty()) {
            return 0;
        }
        int lower = 0;
        int upper = size() - 1;

        if (value - get(upper).TargetMz >= 0) {
            return upper;
        }
        if (value - get(0).TargetMz <= 0) {
            return 0;
        }

        while (lower <= upper) {
            int middle = (lower + upper) / 2;
            float comparisonResult = value - get(middle).TargetMz;
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
        while (upper > 0 && get(upper).TargetMz >= value) {
            upper--;
        }
        return upper;
    }

    public int BinarySearchClosest(float value) {
        if (isEmpty()) {
            return 0;
        }
        int lower = 0;
        int upper = size() - 1;

        if (value - get(upper).TargetMz >= 0) {
            return upper;
        }
        if (value - get(0).TargetMz <= 0) {
            return 0;
        }

        while (lower <= upper) {
            int middle = (lower + upper) / 2;
            float comparisonResult = value - get(middle).TargetMz;
            if (comparisonResult == 0) {
                return middle;
            } else if (comparisonResult < 0) {
                upper = middle - 1;
            } else {
                lower = middle + 1;
            }
        }

        if (Math.abs(value - get(lower).TargetMz) > Math.abs(value - get(upper).TargetMz)) {
            return upper;
        } else {
            return lower;
        }
    }

}
