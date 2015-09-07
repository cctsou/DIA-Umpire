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
public class FastaManipulator {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) throws FileNotFoundException, IOException {
        String filename = "F:/fasta/Homo_sapiens.GRCh37.74.pep.all_plusREV.fa";
        //String filename="F:/fasta/UPS_1_2_Ecoli_PlusRevTag2.fa";
        String resultname = "F:/fasta/Homo_sapiens.GRCh37.74.pep.all_plusREV_mod.fa";
        BufferedReader reader = new BufferedReader(new FileReader(filename));

        FileWriter writer = new FileWriter(resultname);
        String line = "";
        while ((line = reader.readLine()) != null) {
            if (line.startsWith(">")) {
                line = line.substring(0, line.indexOf(" "));
            }
            writer.write(line + "\n");
        }

        writer.close();
        reader.close();
    }

}
