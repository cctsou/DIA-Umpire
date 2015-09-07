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
package MSUmpire.MSMSDBSearch;

import java.io.IOException;

import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;
import javax.xml.parsers.ParserConfigurationException;
import org.apache.commons.io.FilenameUtils;
import org.xml.sax.SAXException;

/**
 *
 * @author Chih-Chiang Tsou <chihchiang.tsou@gmail.com>
 */
public class TandemSearchFindBestThreshold {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) throws IOException, InterruptedException, ParserConfigurationException, SAXException {

        if (args.length < 10) {
            System.out.print("java -jar -d64 -Xms15G TandemSearchFindBestThreshold.jar mzXML fasta StartPrecursorPPM EndPrecursorPPM StartFragPPM EndFragPPM PPMInterval FilePath NoThread NoThreadTandem\n");
            return;
        }

        String Filepath = args[7];
        System.out.print("FilePath:" + Filepath + "\n");
        String mzxML = args[0];

        int StartPrecursorPPM = Integer.parseInt(args[2]);
        int EndPrecursorPPM = Integer.parseInt(args[3]);
        int StartFragPPM = Integer.parseInt(args[4]);
        int EndFragPPM = Integer.parseInt(args[5]);
        int PPMInterval = Integer.parseInt(args[6]);
        //String ExportCSV=args[7];

        ExecutorService executorPool = null;
        executorPool = Executors.newFixedThreadPool(Integer.parseInt(args[8]));

        for (int precursorppm = StartPrecursorPPM; precursorppm < EndPrecursorPPM; precursorppm += PPMInterval) {
            for (int fragppm = StartFragPPM; fragppm < EndFragPPM; fragppm += PPMInterval) {
                TandemParam para = new TandemParam(DBSearchParam.SearchInstrumentType.TOF5600);
                para.FastaPath = Filepath + args[1];
                para.NoCPUs = Integer.parseInt(args[9]);
                para.FragPPM = fragppm;
                para.PrecursorPPM = precursorppm;
                para.parameterPath = Filepath + FilenameUtils.getBaseName(mzxML) + "_" + precursorppm + "_" + fragppm + "_tandem.para";
                para.RawSearchResult = Filepath + FilenameUtils.getBaseName(mzxML) + "_" + precursorppm + "_" + fragppm + ".tandem";
                para.InteractPepXMLPath = Filepath + "interact-" + FilenameUtils.getBaseName(mzxML) + "_" + precursorppm + "_" + fragppm + ".pep.xml";
para.ProtXMLPath=para.InteractPepXMLPath.replace(".pep.xml", ".prot.xml");                
//para.OutputSeqPath=Filepath+FilenameUtils.getBaseName(mzxML)+"_"+precursorppm+"_"+fragppm+".seq";
                para.SpectrumPath = Filepath + mzxML;
                TandemSearch search = new TandemSearch(para);
                executorPool.execute(search);
            }
        }
        executorPool.shutdown();
//        while (!executorPool.isTerminated()) {
//        }
         try {
            executorPool.awaitTermination(Long.MAX_VALUE, TimeUnit.NANOSECONDS);
        } catch (InterruptedException e) {
             System.err.println("interrupted..");
        }
//        float errorrate=0.025f;
//        FileWriter writer=new FileWriter(Filepath+ExportCSV);
//        LCMSPepID pepID=new LCMSPepID(mzxML);
//         for (int precursorppm = StartPrecursorPPM; precursorppm < EndPrecursorPPM; precursorppm += PPMInterval) {
//            for (int fragppm = StartFragPPM; fragppm < EndFragPPM; fragppm += PPMInterval) {
//                pepID.ParseFromXML(Filepath+ "interact-"+FilenameUtils.getBaseName(mzxML)+"_"+precursorppm+"_"+fragppm+".pep.xml", 0f);
//                                
//                writer.write(precursorppm+","+fragppm+",");
//                TreeMap<Float, PSM> sortedPSMs=new TreeMap();
//                for(PSM psm :  pepID.PSMList.values())
//                {
//                    sortedPSMs.put(psm.Probability, psm);                    
//                }                
//                int total=0;
//                int decoy=0;
//                float prob=1f;
//                while(!sortedPSMs.isEmpty()){
//                    PSM topscore= sortedPSMs.pollLastEntry().getValue(); 
//                    total++;
//                    boolean isdecoy=true;
//                    for(String ProID :topscore.ParentProtIDs){
//                        if(!ProID.startsWith("rev_")){
//                            isdecoy=false;
//                            break;
//                        }
//                    }
//                    if(isdecoy){
//                        decoy++;
//                    }
//                    if(((float)decoy/(float)total)>errorrate){
//                        prob=topscore.Probability;
//                        break;
//                    }
//                }
//                writer.write(total+","+decoy+","+prob+"\n");
//            }
//        }
//         writer.close();
    }
}
