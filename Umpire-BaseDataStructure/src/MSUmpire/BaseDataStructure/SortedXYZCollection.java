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
package MSUmpire.BaseDataStructure;

import MSUmpire.SortedListLib.SortedList;
import java.io.Serializable;
import java.util.Comparator;

/**
 *
 * @author Chih-Chiang Tsou <chihchiang.tsou@gmail.com>
 */
public class SortedXYZCollection extends SortedList<XYZData> implements Serializable {
    private static final long serialVersionUID = 836465863L;

    public SortedXYZCollection() {
        super(new Comparator<XYZData>() {
            @Override
            public int compare(XYZData x, XYZData y) {
                if (x.getX() == y.getX()) {
                    return 1;
                }
                return -Float.compare(y.getX(), x.getX());
            }
        });
    }

    public XYZData GetCloset(float value) {
        return get(BinarySearchClosest(value));
    }
    public int BinarySearchClosest(float value) {

        if (isEmpty()) {
            return 0;
        }
        int lower = 0;
        int upper = size() - 1;

        if (value - get(upper).getX() >= 0) {
            return upper;
        }
        if (value - get(0).getX() <= 0) {
            return 0;
        }

        while (lower <= upper) {
            int middle = (lower + upper) / 2;
            float comparisonResult = value - get(middle).getX();
            if (comparisonResult == 0) {
                return middle;
            } else if (comparisonResult < 0) {
                upper = middle - 1;
            } else {
                lower = middle + 1;
            }
        }

        if (Math.abs(value - get(lower).getX()) > Math.abs(value - get(upper).getX())) {
            return upper;
        } else {
            return lower;
        }
    }

}
