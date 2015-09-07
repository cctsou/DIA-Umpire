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
package MSUmpire.SortedListLib;

import MSUmpire.SortedListLib.SortedList;
import java.util.Comparator;

/**
 * Provides a {@code SortedList} which sorts the elements by their natural
 * order.
 *
 * @author Mark Rhodes
 * @version 1.1
 * @see SortedList
 * @param <T> any {@code Comparable}
 */
public class NaturalSortedList<T extends Comparable<? super T>>
        extends SortedList<T> {

    private static final long serialVersionUID = -8834713008973648930L;

    /**
     * Constructs a new {@code NaturalSortedList} which sorts elements according
     * to their <i>natural order</i>.
     */
    public NaturalSortedList() {
        super(new Comparator<T>() {
            public int compare(T one, T two) {
                return one.compareTo(two);
            }
        });
    }
}
