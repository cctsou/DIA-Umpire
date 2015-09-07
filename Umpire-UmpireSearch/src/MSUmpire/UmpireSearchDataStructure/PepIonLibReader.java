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
package MSUmpire.UmpireSearchDataStructure;

import MSUmpire.BaseDataStructure.XYData;
import com.compomics.util.experiment.biology.Peptide;
import com.compomics.util.experiment.identification.matches.ModificationMatch;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import org.apache.commons.io.FilenameUtils;
import org.xmlpull.v1.XmlPullParserException;

/**
 *
 * @author Chih-Chiang Tsou <chihchiang.tsou@gmail.com>
 */
public class PepIonLibReader {

    public HashMap<String, ArrayList<String>> PepLib = new HashMap<>();
    public SortedPepCandidate IonLib = new SortedPepCandidate();

    String fasta;
    public XYData LibMzRange;

    public void ReadPepLib(String fasta) throws IOException, XmlPullParserException {
        this.fasta = fasta;
        BufferedReader reader2 = new BufferedReader(new FileReader(FilenameUtils.getFullPath(fasta) + FilenameUtils.getBaseName(fasta) + ".pepprot"));
        String line = "";
        reader2.readLine();

        while ((line = reader2.readLine()) != null) {
            String peptideseq = line.split("\t")[0];
            String prots = line.split("\t")[1];

            ArrayList<String> proteins = new ArrayList<>();
            String[] protarray = prots.split(";");
            for (int i = 0; i < protarray.length; i++) {
                proteins.add(protarray[i]);
            }
            PepLib.put(peptideseq, proteins);
        }
        reader2.close();
    }

    public SortedPepCandidate ReadWholeIonLib() throws IOException, XmlPullParserException {

        ArrayList<PepIonCandidate> candidates = new ArrayList<>();
        File folder = new File(FilenameUtils.getFullPath(fasta) + FilenameUtils.getBaseName(fasta) + "/");
        for (final File fileEntry : folder.listFiles()) {
            if (fileEntry.isFile() && fileEntry.getName().endsWith(".ionlib")) {
                System.out.println("Loading ion library:" + fileEntry.getAbsolutePath());
                BufferedReader reader = new BufferedReader(new FileReader(fileEntry.getAbsolutePath()));
                String line = "";
                reader.readLine();
                while ((line = reader.readLine()) != null) {
                    float mzvalue = Float.parseFloat(line.split("\t")[0]);
                    //String peptideseq = line.split("\t")[1];
                    int charge = Integer.parseInt(line.split("\t")[2]);
                    PepIonCandidate candidate = new PepIonCandidate();
                    if (line.split("\t").length > 3) {
                        //candidate.ModString = line.split("\t")[3];
                    }
                    //candidate.Sequence=peptideseq;
                    candidate.Charge = charge;
                    candidate.SetMz(mzvalue);
                    candidates.add(candidate);
                    line = null;
                }
                reader.close();
                break;
            }
        }
        IonLib.clear();
        IonLib.addAll(candidates);
        IonLib.Finalize();
        return IonLib;
    }

    public SortedPepCandidate ReadIonLib(float mz) throws IOException, XmlPullParserException {

        if (LibMzRange != null && mz <= LibMzRange.getY()) {
            return IonLib;
        }

        IonLib.clear();
        File folder = new File(FilenameUtils.getFullPath(fasta) + FilenameUtils.getBaseName(fasta) + "/");
        for (final File fileEntry : folder.listFiles()) {
            if (fileEntry.isFile() && fileEntry.getName().endsWith(".ionlib")) {
                float lowmz = Float.parseFloat(FilenameUtils.getBaseName(fileEntry.getName()).split("_")[1]);
                if (mz >= lowmz && mz <= lowmz + 1200f) {
                    System.out.println("Loading ion library:" + fileEntry.getAbsolutePath());
                    BufferedReader reader = new BufferedReader(new FileReader(fileEntry.getAbsolutePath()));
                    String line = "";
                    reader.readLine();
                    while ((line = reader.readLine()) != null) {
                        float mzvalue = Float.parseFloat(line.split("\t")[0]);
                        String peptideseq = line.split("\t")[1];
                        int charge = Integer.parseInt(line.split("\t")[2]);
                        ArrayList<ModificationMatch> mods = new ArrayList<ModificationMatch>();
                        if (line.split("\t").length > 3) {
                            String[] modifications = line.split("\t")[3].split(";");
                            for (int i = 0; i < modifications.length; i++) {
                                String name = modifications[i].split("@")[0];
                                int site = Integer.parseInt(modifications[i].split("@")[1]);
                                mods.add(new ModificationMatch(name, true, site));
                            }
                        }
                        Peptide peptide = new Peptide(peptideseq, new ArrayList<String>(), mods);
                        peptide.setParentProteins(PepLib.get(peptide.getSequence()));
                        PepIonCandidate candidate = new PepIonCandidate();
                        candidate.peptide = peptide;
                        candidate.Charge = charge;
                        candidate.SetMz(mzvalue);
                        IonLib.add(candidate);
                    }
                    reader.close();
                    LibMzRange = new XYData(lowmz, lowmz + 1200);
                    break;
                }
            }
        }
        IonLib.Finalize();
        return IonLib;
    }

}
