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
package MSUmpire.MathPackage;

import java.util.concurrent.locks.ReadWriteLock;
import java.util.concurrent.locks.ReentrantReadWriteLock;
import jsc.distributions.ChiSquared;

/**
 *
 * @author Chih-Chiang Tsou <chihchiang.tsou@gmail.com>
 */
public class ChiSquareGOF {

    private static ChiSquareGOF models = null;
    public static ChiSquared[] chimodels;
    public static ReadWriteLock lock = new ReentrantReadWriteLock();    
    
    private ChiSquareGOF(int maxpeak) {        
        chimodels = new ChiSquared[maxpeak-1];
        for (int i = 1; i <= maxpeak; i++) {
            chimodels[i - 1] = new ChiSquared(i);
        }
    }

    public static ChiSquareGOF GetInstance(int maxpeak) {
        if (models == null || (maxpeak>1 && maxpeak >= chimodels.length)) {
            lock.writeLock().lock();
            try {
                if (models == null) {
                    models = new ChiSquareGOF(maxpeak);
                }
            } finally {
                lock.writeLock().unlock();
            }
        }
        return models;
    }

    public float GetGoodNessOfFitProb(float[] expected, float[] observed) {
        float gof = 0f;
        int nopeaks = 0;
        for (int i = 0; i < Math.min(observed.length,expected.length); i++) {
            if (observed[i] > 0) {
                float error = expected[i] - observed[i];
                gof += (error * error) / (expected[i] * expected[i]);
                nopeaks++;
            }
        }
        if (Float.isNaN(gof) || nopeaks<2){
            return 0f;
        }
        float prob = 1 - (float) chimodels[nopeaks-2].cdf(gof);
        return prob;
    }
}
