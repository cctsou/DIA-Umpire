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

import com.compomics.util.experiment.biology.Protein;
import com.compomics.util.experiment.identification.SequenceFactory;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;

/**
 *
 * @author Chih-Chiang Tsou <chihchiang.tsou@gmail.com>
 */
public class Fasta_AccMod {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) throws IOException, FileNotFoundException, StringIndexOutOfBoundsException, ClassNotFoundException, IllegalArgumentException, InterruptedException {
        String filename = "F:/fasta/sprot_Drosophila_UPS_plusREV.2014Jan24.fa";
        //String filename="F:/fasta/UPS_1_2_Ecoli_PlusRevTag2.fa";
        String resultname = "F:/fasta/sprot_Drosophila_UPS_plusREV.2014Jan24.fa";
        BufferedReader reader = new BufferedReader(new FileReader(filename));

        FileWriter writer = new FileWriter(resultname);
        String line = "";
        while ((line = reader.readLine()) != null) {
//            if(line.startsWith(">rev_")){
//                line=line.replace("rev_sp|", "rev_").replace("rev_tr|", "rev_");
//            }
            if (line.startsWith(">")) {
                line = line.replace("sp|", "sp_").replace("tr|", "tr_");

//                if (!line.startsWith(">sp")) {
//                    if (line.indexOf(" ##") != -1) {
//                        line = line.replace(" ## ", " ID=|" + line.substring(1, line.indexOf(" ##")) + "|");
//                    } else {
//                        line = line.replaceFirst(" ", " ID=|" + line.substring(1, line.indexOf(" ")) + "|");
//                    }
//                }
//                if (!line.startsWith(">sp")) {
//                    if (line.indexOf(" ##") != -1) {
//                        line = line.replace(" ## ", "|"+ line.substring(1, line.indexOf(" ##")) + "_BOVIN ");
//                    } else {
//                        line = line.replaceFirst(" ", "|" + line.substring(1, line.indexOf(" ")) + "_BOVIN ");
//                    }
//                }
//            if(line.startsWith(">rev_")){
//                line=line.replace(">rev_", ">rev_sp_");
//            }
//            else{
//                line=line.replace(">", ">sp_");
//            }            
//                if (line.contains("| COMMON CONTAMINANT!")) {
//                    String pro = (String) line.subSequence(1, line.indexOf("|"));
//                    line = line.replace("| COMMON CONTAMINANT!", "_CC");
//                }
            }
            writer.write(line + "\n");
        }

        writer.close();
        reader.close();

        File fastaFile2 = new File(resultname);
        SequenceFactory sequenceFactory2 = SequenceFactory.getInstance();
        sequenceFactory2.loadFastaFile(fastaFile2);

        Protein protein = sequenceFactory2.getProtein("P09111");
        String seq = protein.getSequence();

        System.out.print(sequenceFactory2.getAccessions().size());

//       for(String acc : sequenceFactory.getAccessions()){
//           if(acc.contains("UBE2I_HUMAN")){
//               System.out.print(acc);
//           }
//       }
        //System.out.print(sequenceFactory.getProtein("sp_TRYP_PIG").getAccession());
    }

}
