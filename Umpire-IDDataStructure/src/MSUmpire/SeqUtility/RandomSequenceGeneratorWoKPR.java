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
package MSUmpire.SeqUtility;

import java.util.Random;
import java.util.concurrent.locks.ReadWriteLock;
import java.util.concurrent.locks.ReentrantReadWriteLock;

public class RandomSequenceGeneratorWoKPR {

    private String RandomSeq;
    int RandomIdx = 0;
    /**
     * All possible characters
     */
    private final char[] CHARS = {
        'A', 'N', 'D', 'C', 'Q', 'E',
        'G', 'H', 'I', 'L', 'M', 'F',
        'S', 'T', 'W', 'Y', 'V', 'B', 'Z', 'X'};

    public static RandomSequenceGeneratorWoKPR GetInstance(){
        if(randomSequenceGeneratorWoKPR==null){
            randomSequenceGeneratorWoKPR=new RandomSequenceGeneratorWoKPR();
        }
        return randomSequenceGeneratorWoKPR;
    }
    
    private static RandomSequenceGeneratorWoKPR randomSequenceGeneratorWoKPR=null;
    
    private RandomSequenceGeneratorWoKPR(){
        RandomSeq = generate(50);
    }
    
    transient ReadWriteLock lock = new ReentrantReadWriteLock();
   
    public char GetNext(){
        lock.writeLock().lock();
        try {
            if (RandomIdx == 50) {
                RandomIdx %= 50;
            }
            return RandomSeq.charAt(RandomIdx++);
        } finally {
            lock.writeLock().unlock();
        }
    }
    
    /**
     * Number of possible characters
     */
    private final int NUMBER_OF_CHARS = CHARS.length;
   
    
    /**
     * Random generator
     */
    private Random random = new Random();
    
            
    /**
     * Returns random sequence
     *
     * @param length Size of the sequence
     * @return Random sequence
     */
    private String generate(int length) {
        StringBuffer buffer = new StringBuffer();
        char randomChar;
        int randomInt;
        for (int i = 0; i < length; i++) {
            randomInt = random.nextInt(NUMBER_OF_CHARS);
            randomChar = CHARS[randomInt];
            buffer.append(randomChar);
        }
        return buffer.toString();
    }
}
