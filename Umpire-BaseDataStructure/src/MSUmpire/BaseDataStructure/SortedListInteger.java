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
/*
 * 
 */

import MSUmpire.SortedListLib.SortedList;
import java.io.Serializable;
import java.util.Comparator;

/**
 *
 * @author Chih-Chiang Tsou <chihchiang.tsou@gmail.com>
 */
public class SortedListInteger extends SortedList<Integer> implements Serializable{

    public SortedListInteger() {
        super(new Comparator<Integer>() {
            @Override
            public int compare(Integer x, Integer y) {
                return -x.compareTo(y);
            }
        });
    }

    public int BinarySearchLower(Integer value) {
        return BinarySearchLower(value);
    }

    public int BinarySearchHigher(float value) {
        if (isEmpty()) {
            return 0;
        }
        int lower = 0;
        int upper = size() - 1;

        if (value - get(upper) >= 0) {
            return upper;
        }
        if (value - get(0) <= 0) {
            return 0;
        }

        while (lower <= upper) {
            int middle = (lower + upper) / 2;
            float comparisonResult = value - get(middle);
            if (comparisonResult == 0) {
                while (middle - 1 >= 0 && get(middle - 1) == value) {
                    middle--;
                }
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
        while (lower < size() - 1 && get(lower) <= value) {
            lower++;
        }
        return lower;
    }

    public int BinarySearchLower(float value) {
        if (isEmpty()) {
            return 0;
        }
        int lower = 0;
        int upper = size() - 1;

        if (value - get(upper) >= 0) {
            return upper;
        }
        if (value - get(0) <= 0) {
            return 0;
        }

        while (lower <= upper) {
            int middle = (lower + upper) / 2;
            float comparisonResult = value - get(middle);
            if (comparisonResult == 0) {
                while (middle - 1 >= 0 && get(middle - 1) == value) {
                    middle--;
                }
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
        while (upper > 0 && get(upper) >= value) {
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

        if (value - get(upper) >= 0) {
            return upper;
        }
        if (value - get(0) <= 0) {
            return 0;
        }

        while (lower <= upper) {
            int middle = (lower + upper) / 2;
            float comparisonResult = value - get(middle);
            if (comparisonResult == 0) {
                return middle;
            } else if (comparisonResult < 0) {
                upper = middle - 1;
            } else {
                lower = middle + 1;
            }
        }

        if (Math.abs(value - get(lower)) > Math.abs(value - get(upper))) {
            return upper;
        } else {
            return lower;
        }
    }

}
