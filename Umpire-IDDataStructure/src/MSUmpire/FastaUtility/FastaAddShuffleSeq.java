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
import MSUmpire.SeqUtility.ShuffledSeqGen;
import Utility.UpdateProcess;
import jaligner.matrix.Matrix;
import jaligner.matrix.MatrixLoader;
import jaligner.matrix.MatrixLoaderException;
import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;
import org.apache.log4j.Logger;
import org.xmlpull.v1.XmlPullParserException;

/**
 *
 * @author Chih-Chiang Tsou <chihchiang.tsou@gmail.com>
 */
public class FastaAddShuffleSeq {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) throws IOException, FileNotFoundException, StringIndexOutOfBoundsException, ClassNotFoundException, IllegalArgumentException, InterruptedException, XmlPullParserException, MatrixLoaderException {
        String filename = "F:/fasta/uniprot_sp-human.fasta";
        String resultfile = "F:/fasta/uniprot_sp-human_sh.fasta";
        BufferedReader reader = new BufferedReader(new FileReader(filename));
        FileWriter resultfastaFile = new FileWriter(resultfile);
        FastaParser_V2 fastaparser=new FastaParser_V2(filename);                
        System.out.println("No. of proteins: "+fastaparser.ProteinList.size());
        StringBuilder outputBuilder=new StringBuilder();
        int missedcleave=1;
        int minlength=8;
        int maxlength=30;
        fastaparser.digestion(missedcleave, minlength, maxlength,"DECOY");
        //UpdateProcess update=new UpdateProcess();
        UpdateProcess update=null;
               
        Matrix blosum62=MatrixLoader.load("BLOSUM62");
        ExecutorService executorPool;
        executorPool = Executors.newFixedThreadPool(10);
        //update.SetTotal(fastaparser.PeptideList.size());
        //Thread thread = new Thread(update);
        //thread.start();
        ArrayList<ShuffledSeqGen> ResultList=new ArrayList<>();
        for (FastaParser_V2.PeptideEntry peptide : fastaparser.PeptideList.values()) {
            ShuffledSeqGen shuffledSeqGen = new ShuffledSeqGen(peptide.Sequence,blosum62,update);
            ResultList.add(shuffledSeqGen);
            executorPool.execute(shuffledSeqGen);
        }
        executorPool.shutdown();
//        while (!executorPool.isTerminated()) {
//        }
        try {
            executorPool.awaitTermination(Long.MAX_VALUE, TimeUnit.NANOSECONDS);
        } catch (InterruptedException e) {
            Logger.getRootLogger().info("interrupted..");
        }
        for(ShuffledSeqGen gen : ResultList){
            fastaparser.PeptideList.get(gen.seq).Decoy=gen.decoy;
        }
        
        for (FastaParser_V2.ProteinEntry protein : fastaparser.ProteinList.values()) {
            StringBuilder DecoyProtein=new StringBuilder();
            for(String pepseq: protein.Peptides){                
                DecoyProtein.append(fastaparser.PeptideList.get(pepseq).Decoy);
            }            
            outputBuilder.append(">" + protein.ACC+ " "+protein.Des+ "\n" + protein.Seq +"\n");
            outputBuilder.append(">DECOY_" + protein.ACC + "\n" + DecoyProtein.toString()+ "\n");
        }
        resultfastaFile.write(outputBuilder.toString());
        resultfastaFile.close();
        fastaparser.FasterSerialzationWrite(resultfile);
    }

}
