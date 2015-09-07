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
package Test;

import MSUmpire.UmpireSearchDataStructure.InsilicoDigestion;
import MSUmpire.UmpireSearchDataStructure.PepIonLib;
import com.compomics.util.experiment.biology.Ion;
import com.compomics.util.experiment.biology.IonFactory;
import com.compomics.util.experiment.biology.Peptide;
import com.compomics.util.experiment.biology.ions.PeptideFragmentIon;
import com.compomics.util.experiment.identification.matches.ModificationMatch;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import org.apache.commons.io.FilenameUtils;
import org.xmlpull.v1.XmlPullParserException;

/**
 *
 * @author Chih-Chiang Tsou <chihchiang.tsou@gmail.com>
 */
public class NewMain {

    /**
     * @param args the command line arguments
     * @throws java.lang.InterruptedException
     * @throws org.xmlpull.v1.XmlPullParserException
     * @throws java.io.IOException
     * @throws java.io.FileNotFoundException
     * @throws java.lang.ClassNotFoundException
     */
    public static void main(String[] args) throws InterruptedException, XmlPullParserException, IOException, FileNotFoundException, ClassNotFoundException {

        ArrayList<ModificationMatch> mods = new ArrayList<>();
        ArrayList<ModificationMatch> mods2 = new ArrayList<>();
        ArrayList<String> Prots = new ArrayList<>();
        mods.add(new ModificationMatch("oxidation of m", true, 5));
        Peptide pep = new Peptide("DPEQQMELNRESVR", Prots, mods);
        Peptide pep2 = new Peptide("DPEQQMELNRESVR", Prots, mods2);
        System.out.print(pep.getMass() + "\n");
        System.out.print(pep2.getMass() + "\n");
//        ArrayList<Ion> allfragment = IonFactory.getInstance().getFragmentIons(pep);
//        for (Ion frag : allfragment) {
//            if (frag.getType() == Ion.IonType.PEPTIDE_FRAGMENT_ION && (frag.getSubType() == PeptideFragmentIon.B_ION || frag.getSubType() == PeptideFragmentIon.Y_ION) && "".equals(frag.getNeutralLossesAsString())) {
//                System.out.println(frag.getName() + ((PeptideFragmentIon) frag).getNumber() + ";" + frag.getTheoreticMass());
//            }
//        }
//        System.out.print("-------------------------------------------------\n");
//        allfragment = IonFactory.getInstance().getFragmentIons(pep2);
//        for (Ion frag : allfragment) {
//            if (frag.getType() == Ion.IonType.PEPTIDE_FRAGMENT_ION && (frag.getSubType() == PeptideFragmentIon.B_ION || frag.getSubType() == PeptideFragmentIon.Y_ION) && "".equals(frag.getNeutralLossesAsString())) {
//                System.out.println(frag.getName() + ((PeptideFragmentIon) frag).getNumber() + ";" + frag.getTheoreticMass());
//            }
//        }
        //Factorial();
        String fasta = "F:/fasta/swissprot_Hs_plusREV.2013Jan09_2.fa";
        //String fasta = "C:/inetpub/wwwroot/ISB/data/Fdata/SWATH_IDA_Project94/New_Rep/0.01FDR/UPS_1_2_Ecoli_PlusRevTag2.fa";
        System.out.println("Loading fasta file:" + fasta + "\n");

        //insilico digestion
        if (!new File(FilenameUtils.getFullPath(fasta) + FilenameUtils.getBaseName(fasta) + ".ionidx").exists() || !new File(FilenameUtils.getFullPath(fasta) + FilenameUtils.getBaseName(fasta) + ".pepprot").exists()) {
            InsilicoDigestion digestion = new InsilicoDigestion(fasta);
            digestion.Perform();
        }

        PepIonLib lib = new PepIonLib(fasta);
    }

    private static void Factorial() throws IOException {
        FileWriter writer = new FileWriter("factorial.txt");
        double fact = 0;
        int x = 1000;
        writer.write("1\t" + fact + "\n");
        for (int i = 2; i <= x; i++) // loop
        {
            fact += Math.log(i); // shorthand for: fact = fact * i;
            writer.write(i + "\t" + fact + "\n");
        }
        writer.close();
    }
}
