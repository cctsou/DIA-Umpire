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
package MSUmpire.SpectralProcessingModule;

import MSUmpire.BaseDataStructure.XYPointCollection;
import java.util.ArrayList;

/**
 *
 * @author Chih-Chiang Tsou
 */
public class Binning {
    
    public XYPointCollection Binning(XYPointCollection scan, float threshold, ArrayList<Integer> removeList) {

        float Binsize = 0.1f;
        XYPointCollection returnCollection = new XYPointCollection();
        int ArraySize = (int) Math.ceil((scan.Data.get(scan.PointCount() - 1).getX() - scan.Data.get(0).getX()) / Binsize) + 2;

        float[] mzarray = new float[ArraySize];
        float[] valueindex = new float[ArraySize];

        for (int i = 0; i < ArraySize; i++) {
            valueindex[i] = scan.Data.get(0).getX() + Binsize * i;
        }
        int arrayidx = 1;
        for (int i = 0; i < scan.PointCount(); i++) {
            if (scan.Data.get(i).getY() > threshold) {
                
                while (scan.Data.get(i).getX() > valueindex[arrayidx]) {
                    arrayidx++;
                }
                float intensity = scan.Data.get(i).getY();
                float mz = scan.Data.get(i).getX();
                float intenlow = intensity * (Binsize - (mz - valueindex[arrayidx - 1])) / Binsize;
                float intenup = intensity * (Binsize - (valueindex[arrayidx] - mz)) / Binsize;

                if (intenlow > mzarray[arrayidx - 1]) {
                    mzarray[arrayidx - 1] = intenlow;
                }
                if (intenup > mzarray[arrayidx]) {
                    mzarray[arrayidx] = intenup;
                }
            }
        }
        for (int i = 0; i < mzarray.length; i++) {
            if (removeList == null || removeList.isEmpty() || !removeList.contains(i)) {
                if (mzarray[i] > threshold) {
                    returnCollection.AddPoint(i, mzarray[i]);
                }
            }
        }
        return returnCollection;
    }
}
