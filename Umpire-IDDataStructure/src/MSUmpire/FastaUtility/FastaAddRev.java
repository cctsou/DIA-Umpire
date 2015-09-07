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

import MSUmpire.FastaParser.FastaParser_V2;
import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;

/**
 *
 * @author Chih-Chiang Tsou <chihchiang.tsou@gmail.com>
 */
public class FastaAddRev {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) throws IOException, FileNotFoundException, StringIndexOutOfBoundsException, ClassNotFoundException, IllegalArgumentException, InterruptedException {
        String filename = "F:\\Data\\Gomez_QE\\uniprot_sprot_2014-12-11_HUMAN.fasta";
        //String filename="F:/fasta/UPS_1_2_Ecoli_PlusRevTag2.fa";
        String resultfile = "F:\\Data\\Gomez_QE\\uniprot_sprot_2014-12-11_HUMAN_rev.fasta";
        BufferedReader reader = new BufferedReader(new FileReader(filename));

        FileWriter resultfastaFile = new FileWriter(resultfile);

        FastaParser_V2 fastaparser=new FastaParser_V2(filename);
                
        System.out.println("No. of proteins: "+fastaparser.ProteinList.size());

        StringBuilder outputBuilder=new StringBuilder();
        
        for (String ProteinACC : fastaparser.ProteinList.keySet()) {
            String Sequence = fastaparser.ProteinList.get(ProteinACC).Seq;
            String Des = fastaparser.ProteinList.get(ProteinACC).Des;
            outputBuilder.append(">" + ProteinACC+ " "+Des+ "\n" + Sequence +"\n");
            outputBuilder.append(">rev_" + ProteinACC + "\n" + new StringBuilder(Sequence).reverse().toString() + "\n");
        }
        resultfastaFile.write(outputBuilder.toString());
        resultfastaFile.close();
    }
}
