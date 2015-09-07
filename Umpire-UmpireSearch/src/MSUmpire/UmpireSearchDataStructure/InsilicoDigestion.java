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

import MSUmpire.PSMDataStructure.ModificationInfo;
import MSUmpire.PSMDataStructure.PTMManager;
import com.compomics.util.experiment.biology.Enzyme;
import com.compomics.util.experiment.biology.EnzymeFactory;
import com.compomics.util.experiment.biology.PTM;
import com.compomics.util.experiment.biology.Peptide;
import com.compomics.util.experiment.identification.SequenceFactory;
import com.compomics.util.experiment.identification.matches.ModificationMatch;
import java.io.*;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import org.apache.commons.io.FilenameUtils;
import org.apache.commons.io.IOUtils;
import org.paukov.combinatorics.Factory;
import org.paukov.combinatorics.Generator;
import org.paukov.combinatorics.ICombinatoricsVector;
import org.xmlpull.v1.XmlPullParserException;

/**
 *
 * @author Chih-Chiang Tsou <chihchiang.tsou@gmail.com>
 */
public class InsilicoDigestion {

    String fasta;
    public int nMissedCleavages = 1;
    public int nMin = 4;
    public int nMax = 50;
    public int mzMin = 350;
    public int mzMax = 1200;
    public int StartCharge = 2;
    public int EndCharge = 4;
    public HashMap<String, ArrayList<String>> PepLib = new HashMap<>();
    //public SortedPepCandidate IonLib = new SortedPepCandidate();
    ArrayList<ModificationInfo> modifications = new ArrayList<>();
    public ArrayList<Integer> massrange;
    public ArrayList<FileWriter> MassPartitionWriter;

    public InsilicoDigestion(String fasta) throws IOException {
        this.fasta = fasta;
        massrange = new ArrayList<>();
        MassPartitionWriter = new ArrayList<>();
        int mz = mzMin;
        new File(FilenameUtils.getFullPath(fasta) + FilenameUtils.getBaseName(fasta) + "/").mkdirs();
        while (mz < mzMax) {
            FileWriter writer = new FileWriter(FilenameUtils.getFullPath(fasta) + FilenameUtils.getBaseName(fasta) + "/m_" + mz + ".ionlib");
            writer.write("mz\tseq\tcharge\tmod\n");
            massrange.add(mz);
            MassPartitionWriter.add(writer);
            mz += 1200;
        }
    }

    public void SetModifications(String site, float massdiff) throws IOException, XmlPullParserException {
        PTM ptm = PTMManager.GetInstance().GetPTM(site, massdiff);
        ModificationInfo modification = new ModificationInfo();
        modification.massdiff = massdiff;
        modification.site = site;
        modification.modification = ptm;
        modifications.add(modification);
    }

    public void AppendPeptideIntoLib(PepIonCandidate peptide) throws IOException {

        for (int i = 1; i < massrange.size(); i++) {
            if (massrange.get(i) > peptide.NeutralPrecursorMz()) {
                MassPartitionWriter.get(i - 1).write(peptide.NeutralPrecursorMz() + "\t" + peptide.peptide.getSequence() + "\t" + peptide.Charge + "\t" + peptide.GetModString() + "\n");
                //System.out.println(peptide.NeutralPrecursorMz() + "\t" + peptide.peptide.getSequence() + "\t" + peptide.Charge + "\t" + modstring + "\n");
                return;
            }
        }
        MassPartitionWriter.get(MassPartitionWriter.size() - 1).write(peptide.NeutralPrecursorMz() + "\t" + peptide.peptide.getSequence() + "\t" + peptide.Charge + "\t" + peptide.GetModString() + "\n");
    }

    public void GenerateModifications(Peptide peptide, int StartCharge, int EndCharge) throws IOException {
        String peptideseq = peptide.getSequence();

        ArrayList<ModificationMatch> mods = new ArrayList<>();
        for (ModificationInfo mod : modifications) {
            if ("C-term".equals(mod.site)) {
                mods.add(new ModificationMatch(mod.modification.getName(), true, peptideseq.length()));
            }
            if ("N-term".equals(mod.site)) {
                mods.add(new ModificationMatch(mod.modification.getName(), true, 1));
            }
            if (peptideseq.contains(mod.site)) {
                for (int i = 0; i < peptideseq.length(); i++) {
                    if (String.valueOf(peptideseq.charAt(i)).equals(mod.site)) {
                        mods.add(new ModificationMatch(mod.modification.getName(), true, i+1));
                    }
                }
            }
        }

        if (mods.size() > 0 && mods.size() < 6) {
            ICombinatoricsVector<ModificationMatch> initialSet = Factory.createVector(mods);
            Generator<ModificationMatch> gen = Factory.createSubSetGenerator(initialSet);
            for (ICombinatoricsVector<ModificationMatch> subSet : gen) {
                if (subSet.getSize() > 0) {
                    //System.out.println(peptide.getMass()+"\n");
                    PepIonCandidate modcandidate = new PepIonCandidate();
                    modcandidate.peptide = peptide;
                    peptide.clearModificationMatches();
                    for (ModificationMatch mod : (ArrayList<ModificationMatch>) subSet.getVector()) {
                        peptide.addModificationMatch(mod);
                    }

                    //Peptide modpeptide = new Peptide(peptideseq, new ArrayList<String>(), (ArrayList<ModificationMatch>) subSet.getVector());
                    peptide.estimateTheoreticMass();
                    //System.out.println(peptide.getMass()+"\n");  
                    for (int charge = StartCharge; charge <= EndCharge; charge++) {
                        modcandidate.Charge = charge;
                        //System.out.println(peptide.getMass()+"\n");
                        float mz = modcandidate.UpdateNeutralPrecursorMz();
                        if (mz >= mzMin && mz <= mzMax) {
                            AppendPeptideIntoLib(modcandidate);
                        }
                    }
                    //modcandidate.peptide=null;
                    //IonLib.add(modcandidate);
                }
            }
        }
    }

    private void SortAllIonLib() {
        ArrayList<PepIonCandidate> IonLib = new ArrayList<>();
        Collections.sort(IonLib, new Comparator<PepIonCandidate>() {
            @Override
            public int compare(PepIonCandidate o1, PepIonCandidate o2) {
                return -Float.compare(o1.NeutralPrecursorMz(), o2.NeutralPrecursorMz());
            }
        });
    }

    public void Perform() throws IOException, FileNotFoundException, ClassNotFoundException, XmlPullParserException, IOException, IOException, IOException, IOException, IOException, IOException, IOException, InterruptedException {

        SetModifications("M", 15.9949f);
        SetModifications("C", 57.021464f);
        SetModifications("E", -18.0106f);
        SetModifications("Q", -17.0265f);
        SetModifications("N-term", 42.0106f);

        File fastaFile = new File(fasta);
        SequenceFactory sequenceFactory = SequenceFactory.getInstance();
        sequenceFactory.loadFastaFile(fastaFile);

        EnzymeFactory enzymeFactory = EnzymeFactory.getInstance();
        InputStream is = InsilicoDigestion.class.getClassLoader().getResourceAsStream("Resource/enzymes.xml");

        String PREFIX = "stream2file";
        String SUFFIX = ".tmp";
        File enzymeFile = File.createTempFile(PREFIX, SUFFIX);
        enzymeFile.delete();
        FileOutputStream out = new FileOutputStream(enzymeFile);
        IOUtils.copy(is, out);

        enzymeFactory.importEnzymes(enzymeFile);
        Enzyme enzyme = enzymeFactory.getEnzyme("Trypsin");

        //int count=0;
        //int total=sequenceFactory.getNTargetSequences();
        for (String proteinKey : sequenceFactory.getAccessions()) {
            String sequence = sequenceFactory.getProtein(proteinKey).getSequence().replace("X", "").replace("U", "C");
            //count++;
            //System.out.print("\r processing protein:"+proteinKey+"("+sequence.length()+")--progress:"+count+"/"+total+":"+ Math.round(count/total)+"%");
            ArrayList<String> peptidelist = enzyme.digest(sequence, nMissedCleavages, nMin, nMax);
            for (String peptideseq : peptidelist) {
                if (!PepLib.containsKey(peptideseq)) {
                    Peptide peptide = new Peptide(peptideseq, new ArrayList<String>(), new ArrayList<ModificationMatch>());
                    peptide.clearModificationMatches();
                    peptide.estimateTheoreticMass();
                    float mass = peptide.getMass().floatValue();
                    if (mass >= mzMin * StartCharge && mass <= mzMax * EndCharge) {
                        PepIonCandidate candidate = new PepIonCandidate();
                        candidate.peptide = peptide;
                        for (int charge = StartCharge; charge <= EndCharge; charge++) {
                            candidate.Charge = charge;
                            float mz = candidate.UpdateNeutralPrecursorMz();
                            if (mz >= mzMin && mz <= mzMax) {
                                AppendPeptideIntoLib(candidate);
                            }
                        }
                        GenerateModifications(peptide, StartCharge, EndCharge);
                        //IonLib.add(candidate);
                        ArrayList<String> proteins = new ArrayList<>();
                        proteins.add(proteinKey);
                        PepLib.put(peptideseq, proteins);
                    }
                } else {
                    PepLib.get(peptideseq).add(proteinKey);
                }
            }
        }
        //System.out.println("No. of peptides:"+PepLib.size()+"\n");

        FileWriter writer2 = new FileWriter(FilenameUtils.getFullPath(fasta) + FilenameUtils.getBaseName(fasta) + ".pepprot");

        writer2.write("pep\tprot\n");
        for (String pep : PepLib.keySet()) {
            String prot = "";
            for (String Protein : PepLib.get(pep)) {
                prot += Protein + ";";
            }
            writer2.write(pep + "\t" + prot + "\n");
        }
        writer2.close();
        for (FileWriter writer : MassPartitionWriter) {
            writer.close();
        }
    }
}
