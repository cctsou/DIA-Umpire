/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
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
public class AddReverse {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) throws FileNotFoundException, IOException {
        String filename = "F:/Data/OpenSWATH_GoldStandard/uniprot_human_SGS.fasta";
        String resultfile = "F:/Data/OpenSWATH_GoldStandard/uniprot_human_SGS_Rev.fasta";
        BufferedReader reader = new BufferedReader(new FileReader(filename));

        FileWriter writer = new FileWriter(resultfile);
        String ID = "";
        String RevID = "";
        String Seq = "";
        String RevSeq = "";
        String line = "";
        while ((line = reader.readLine()) != null) {
            if (line.startsWith(">")) {
                if (!ID.isEmpty()) {
                    writer.write(ID + "\n");
                    writer.write(Seq + "\n");
                    writer.write(RevID + "\n");
                    writer.write(new StringBuilder(RevSeq).reverse().toString() + "\n");
                    ID = "";
                    RevID = "";
                    Seq = "";
                    RevSeq = "";
                }
                RevID = line.replace(">", ">rev_");
                ID = line;
            } else {
                Seq += line;
                RevSeq += line;
            }
        }
        if (!ID.isEmpty()) {
            writer.write(ID + "\n");
            writer.write(Seq + "\n");
            writer.write(RevID + "\n");
            writer.write(new StringBuilder(RevSeq).reverse().toString() + "\n");
        }
        writer.close();
        reader.close();
    }

}
