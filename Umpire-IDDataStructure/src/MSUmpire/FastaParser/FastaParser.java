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
package MSUmpire.FastaParser;

import MSUmpire.PSMDataStructure.EnzymeManager;
import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.io.Serializable;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.logging.Level;
import java.util.logging.Logger;
import org.xmlpull.v1.XmlPullParserException;

/**
 *
 * @author Chih-Chiang Tsou <chihchiang.tsou@gmail.com>
 */
public class FastaParser implements Serializable{
    private static final long serialVersionUID = 19398249L;

    public HashMap<String, String[]> ProteinList;
    public transient HashMap<String, ArrayList<String>> PeptideList;
    
    
    public FastaParser(String filename){
        ProteinList=new HashMap<>();
        try {
            Parse(filename);
        } catch (IOException ex) {
            Logger.getLogger(FastaParser.class.getName()).log(Level.SEVERE, null, ex);
        }
    }
    
    private void Parse(String filename) throws FileNotFoundException, IOException {

        BufferedReader reader = new BufferedReader(new FileReader(filename));
        String line = "";
        String ACC = "";
        String des="";
        StringBuilder Seq = new StringBuilder();
        while ((line = reader.readLine()) != null) {
            if (line.startsWith(">")) {
                if (!"".equals(ACC) && !"".equals(Seq.toString())) {
                    ProteinList.put(ACC, new String[]{Seq.toString(),des});
                    Seq=new StringBuilder();
                }
                ACC = line.trim().split(" ")[0].substring(1);
                des=line.replace(">"+ACC+" ", "");
            } else {
                if (!"".equals(line.trim())) {
                    Seq.append(line);
                }
            }
        }
        if (!"".equals(ACC) && !"".equals(Seq.toString())) {
            ProteinList.put(ACC, new String[]{Seq.toString(),des});
        }
        reader.close();
    }
        
    public void digestion(int missedcleave, int minlength, int maxlength) throws XmlPullParserException, IOException {
        PeptideList=new HashMap<>();
        for (String acc : ProteinList.keySet()) {
            String Sequence = ProteinList.get(acc)[0];            
            ArrayList<String> TheoPeptides = EnzymeManager.GetInstance().GetTrypsinNoP().digest(Sequence, missedcleave, minlength, maxlength);
            AddFirstMetDroppedPep(Sequence, missedcleave, minlength, maxlength, TheoPeptides);
            for (String pep : TheoPeptides) {
                if (!PeptideList.containsKey(pep)) {
                    PeptideList.put(pep, new ArrayList<String>());
                }
                PeptideList.get(pep).add(acc);
            }
        }
    }

    public void AddFirstMetDroppedPep(String Sequence, int missedcleave, int minlength, int maxlength, ArrayList<String> TheoPeptides) {
        if (String.valueOf(Sequence.charAt(0)).equals("M")) {
            int mc = 0;
            for (int i = 1; i < Sequence.length(); i++) {
                if (String.valueOf(Sequence.charAt(i)).equals("K") || String.valueOf(Sequence.charAt(i)).equals("R")) {
                    mc++;
                    if (mc > missedcleave) {
                        return;
                    }
                    String pep = Sequence.substring(1, i + 1);
                    if (pep.length() >= minlength && pep.length() <= maxlength && !TheoPeptides.contains(pep)) {
                        TheoPeptides.add(pep);
                    }
                }
            }
        }
    }
    
    
}
