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

import MSUmpire.BaseDataStructure.InstrumentParameter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.concurrent.locks.ReentrantLock;
import org.apache.commons.io.FilenameUtils;
import org.xmlpull.v1.XmlPullParserException;

/**
 *
 * @author Chih-Chiang Tsou <chihchiang.tsou@gmail.com>
 */
public class PepIonLib {

    public HashMap<String, ArrayList<String>> PepSeqLib = new HashMap<>();
    public SortedPepCandidate IonMzLib;
    public String fasta;
    PepIonLibReader reader = new PepIonLibReader();
    ReentrantLock readerlock = new ReentrantLock(true);

    public PepIonLib(String fasta) throws IOException, FileNotFoundException, ClassNotFoundException, InterruptedException, XmlPullParserException {
        this.fasta = fasta;
        if (!new File(FilenameUtils.getFullPath(fasta) + FilenameUtils.getBaseName(fasta) + ".pepprot").exists()) {
            InsilicoDigestion digestion = new InsilicoDigestion(fasta);
            digestion.Perform();
        }
        reader.ReadPepLib(fasta);
        this.PepSeqLib = reader.PepLib;
    }

    public void GetWholeIonLib() throws IOException, XmlPullParserException {
        readerlock.lock();
        try {
            IonMzLib = reader.ReadWholeIonLib();
        } finally {
            readerlock.unlock();
        }
    }

    public SortedPepCandidate GetIonLib(float mz) throws IOException, XmlPullParserException {
        readerlock.lock();
        try {
            return reader.ReadIonLib(mz);
        } finally {
            readerlock.unlock();
        }
    }

    public ArrayList<PepIonCandidate> GetCandidates(float mz, int charge, float ppm) throws IOException, XmlPullParserException {
        SortedPepCandidate IonLib = GetIonLib(mz);
        ArrayList<PepIonCandidate> candidates = new ArrayList<>();
        float lowmz = InstrumentParameter.GetMzByPPM(mz, charge, ppm);
        int startidx = IonLib.BinarySearchLower(lowmz);

        for (int idx = startidx; idx < IonLib.size(); idx++) {
            PepIonCandidate candiate = IonLib.get(idx);
            if (candiate.Charge == charge) {
                if (InstrumentParameter.CalcPPM(mz, candiate.NeutralPrecursorMz()) <= ppm) {
                    candidates.add(candiate);
                } else if (candiate.NeutralPrecursorMz() > mz) {
                    break;
                }
            }
        }
        return candidates;
    }
}
