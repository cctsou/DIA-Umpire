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

import java.io.Serializable;
import java.util.Comparator;

/**
 *
 * @author Chih-Chiang Tsou <chihchiang.tsou@gmail.com>
 */
public class XYComparator implements Serializable{

    public static <XYPoint extends Comparable<? super XYPoint>> Comparator<? super XYPoint> ascending() {
        return new Comparator<XYPoint>() {
            public int compare(XYPoint a, XYPoint b) {
                return a.compareTo(b);
            }
        };
    }

    public static <T extends Comparable<? super T>> Comparator<? super T> descending() {
        return new Comparator<T>() {
            public int compare(T a, T b) {
                return b.compareTo(a);
            }
        };
    }

}
