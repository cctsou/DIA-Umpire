/* 
 * Author: Chih-Chiang Tsou <chihchiang.tsou@gmail.com>
 *             Nesvizhskii Lab, Department of computational medicine and bioinformatics, 
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
package MSUmpire.FastaUtility;

import com.compomics.util.experiment.identification.SequenceFactory;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;

/**
 *
 * @author Chih-Chiang Tsou <chihchiang.tsou@gmail.com>
 */
public class FastaTest {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) throws IOException, FileNotFoundException, StringIndexOutOfBoundsException, ClassNotFoundException, IllegalArgumentException, InterruptedException {
        String filename = "F:\\Data\\SWATH_Glyco\\swissprot_Hs_plusREV.2013Jan09.fa";
        //String filename="F:/fasta/UPS_1_2_Ecoli_PlusRevTag2.fa";
        
        System.out.println("Testing " + filename);
        File fastaFile = new File(filename);
        SequenceFactory sequenceFactory = SequenceFactory.getInstance();
        sequenceFactory.loadFastaFile(fastaFile);
        System.out.print(sequenceFactory.getAccessions().size());
                
    }
}
