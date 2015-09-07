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

import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.ArrayList;
import java.util.LinkedHashMap;
import org.biojava3.core.sequence.AccessionID;
import org.biojava3.core.sequence.ProteinSequence;
import org.biojava3.core.sequence.compound.AminoAcidCompound;
import org.biojava3.core.sequence.compound.AminoAcidCompoundSet;
import org.biojava3.core.sequence.io.FastaReader;
import org.biojava3.core.sequence.io.FastaWriterHelper;
import org.biojava3.core.sequence.io.GenericFastaHeaderParser;
import org.biojava3.core.sequence.io.ProteinSequenceCreator;

/**
 *
 * @author Chih-Chiang Tsou <chihchiang.tsou@gmail.com>
 */
public class UniprotExtraction {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) throws FileNotFoundException, IOException, ClassNotFoundException, Exception {
        String filename = "F:/fasta/uniprot_sprot.fasta";
        String resultname = "F:/fasta/sprot_Drosophila_plusREV.2014Jan24.fa";
        FileInputStream inStream = new FileInputStream(filename);
        FastaReader<ProteinSequence, AminoAcidCompound> fastaReader
                = new FastaReader<>(
                        inStream,
                        new GenericFastaHeaderParser<ProteinSequence, AminoAcidCompound>(),
                        new ProteinSequenceCreator(AminoAcidCompoundSet.getAminoAcidCompoundSet()));
        LinkedHashMap<String, ProteinSequence> b = fastaReader.process();
        ArrayList<ProteinSequence> output = new ArrayList<>();
        for (ProteinSequence seq : b.values()) {
            if (seq.getDescription().contains("OS=Drosophila")) {
                String rev = new StringBuilder(seq.getSequenceAsString()).reverse().toString();
                ProteinSequence revseq = new ProteinSequence(rev);
                revseq.setAccession(new AccessionID("rev_" + seq.getAccession().getID()));
                output.add(seq);
                output.add(revseq);
            }
        }
        FastaWriterHelper.writeProteinSequence(new File(resultname), output);

    }

}
