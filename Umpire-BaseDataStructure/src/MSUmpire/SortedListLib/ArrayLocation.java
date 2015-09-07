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

/**
 * Simple class to define a location within the a patch work array.
 */
class ArrayLocation {

    final int backingListIndex;
    final int subListIndex; //should be -1 in the case that there is no sublist..

    ArrayLocation(int backingListIndex, int subListIndex) {
        this.backingListIndex = backingListIndex;
        this.subListIndex = subListIndex;
    }

    //location where there is no sub list index..
    ArrayLocation(int backingListIndex) {
        this.backingListIndex = backingListIndex;
        this.subListIndex = -1;
    }

    /**
     * Returns whether or not this location has a sub-list index.
     */
    public boolean hasSubListIndex() {
        return subListIndex != -1;
    }

    @Override
    public String toString() {
        return "ArrayLocation:[" + backingListIndex + ", " + subListIndex + "]";
    }

    @Override
    public int hashCode() {
        final int prime = 31;
        int result = 1;
        result = prime * result + backingListIndex;
        result = prime * result + subListIndex;
        return result;
    }

    @Override
    public boolean equals(Object obj) {
        if (this == obj) {
            return true;
        }
        if (obj == null) {
            return false;
        }
        if (getClass() != obj.getClass()) {
            return false;
        }
        ArrayLocation other = (ArrayLocation) obj;
        if (backingListIndex != other.backingListIndex) {
            return false;
        }
        if (subListIndex != other.subListIndex) {
            return false;
        }
        return true;
    }

}
