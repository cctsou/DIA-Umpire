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

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;

/**
 *
 * @author Chih-Chiang Tsou <chihchiang.tsou@gmail.com>
 */
public class FastaCombine {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) throws FileNotFoundException, IOException {
        String filename = "F:/Data/OpenSWATH_GoldStandard/SGS.fasta";
        String filename2 = "F:/Data/OpenSWATH_GoldStandard/uniprot_yeast.fasta";
        String resultname = "F:/Data/OpenSWATH_GoldStandard/uniprot_yeast_SGS.fasta";
        BufferedReader reader = new BufferedReader(new FileReader(filename));
        BufferedReader reader2 = new BufferedReader(new FileReader(filename2));
        FileWriter writer = new FileWriter(resultname);

        String line = "";
        int count = 1;
        while ((line = reader.readLine()) != null) {
            writer.write(line + "\n");
        }
        while ((line = reader2.readLine()) != null) {
            writer.write(line + "\n");
        }
        writer.close();
        reader.close();
        reader2.close();
    }

}
