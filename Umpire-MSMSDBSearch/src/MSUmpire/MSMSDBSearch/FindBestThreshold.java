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

import MSUmpire.PSMDataStructure.LCMSID;
import MSUmpire.PSMDataStructure.PSM;
import MSUmpire.SearchResultParser.PepXMLParser;
import java.io.FileWriter;
import java.io.IOException;
import java.util.TreeMap;
import javax.xml.parsers.ParserConfigurationException;
import org.apache.commons.io.FilenameUtils;
import org.xml.sax.SAXException;
import org.xmlpull.v1.XmlPullParserException;

/**
 *
 * @author Chih-Chiang Tsou <chihchiang.tsou@gmail.com>
 */
public class FindBestThreshold {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) throws IOException, InterruptedException, ParserConfigurationException, SAXException, XmlPullParserException {
        String Filepath = "C:/inetpub/wwwroot/ISB/data/UPS1_Ecoli/";
        String mzxML1 = Filepath + "14343_UPS1_400fm_Ecolilysate_IDA_5600.mzXML";
        String mzxML2 = Filepath + "14343_UPS1_400fm_Ecolilysate_IDA_5600_ABConvert.mzXML";

        TandemParam para = new TandemParam(DBSearchParam.SearchInstrumentType.TOF5600);
        para.NoCPUs = 5;
        para.parameterPath = Filepath + "tandem.para";
        para.FastaPath = Filepath + "UPS_1_2_Ecoli.plusREV.fa";
//        TandemSearch search=new TandemSearch();
//       for (int precursorppm = 20; precursorppm < 40; precursorppm += 5) {
//            for (int fragppm = 5; fragppm < 40; fragppm += 5) {
//                para.FragPPM=fragppm;
//                para.PrecursorPPM=precursorppm;             
//                
//                para.OutputPath=  Filepath+FilenameUtils.getBaseName(mzxML1)+"_"+precursorppm+"_"+fragppm+".tandem";
//                para.InteractPepXMLPath= Filepath+ "interact-"+FilenameUtils.getBaseName(mzxML1)+"_"+precursorppm+"_"+fragppm+".pep.xml";
//                para.OutputSeqPath=Filepath+FilenameUtils.getBaseName(mzxML1)+"_"+precursorppm+"_"+fragppm+".seq";
//                para.SpectrumPath=mzxML1;
//                search.RunTandem(para);
//                para.OutputPath= Filepath+ FilenameUtils.getBaseName(mzxML2)+"_"+precursorppm+"_"+fragppm+".tandem";
//                para.InteractPepXMLPath=  Filepath+"interact-"+FilenameUtils.getBaseName(mzxML2)+"_"+precursorppm+"_"+fragppm+".pep.xml";
//                para.OutputSeqPath=Filepath+FilenameUtils.getBaseName(mzxML2)+"_"+precursorppm+"_"+fragppm+".seq";
//                para.SpectrumPath=mzxML2;
//                search.RunTandem(para);
//            }
//        }
        float errorrate = 0.025f;
        FileWriter writer = new FileWriter("UPS1_Ecoli_IDThreshold.csv");
        LCMSID pepID = new LCMSID(mzxML1,"rev","");
        LCMSID pepID2 = new LCMSID(mzxML2,"rev","");
        for (int precursorppm = 5; precursorppm < 40; precursorppm += 5) {
            for (int fragppm = 5; fragppm < 40; fragppm += 5) {
                PepXMLParser pepXMLParser = new PepXMLParser(pepID, Filepath + "interact-" + FilenameUtils.getBaseName(mzxML1) + "_" + precursorppm + "_" + fragppm + ".pep.xml", 0f);
                pepXMLParser = new PepXMLParser(pepID2, Filepath + "interact-" + FilenameUtils.getBaseName(mzxML2) + "_" + precursorppm + "_" + fragppm + ".pep.xml", 0f);

                writer.write(precursorppm + "," + fragppm + ",");
                TreeMap<Float, PSM> sortedPSMs = new TreeMap();
                TreeMap<Float, PSM> sortedPSMs2 = new TreeMap();
                for (PSM psm : pepID.PSMList.values()) {
                    sortedPSMs.put(psm.Probability, psm);
                }
                for (PSM psm : pepID2.PSMList.values()) {
                    sortedPSMs2.put(psm.Probability, psm);
                }
                int total = 0;
                int decoy = 0;
                float prob = 1f;
                while (!sortedPSMs.isEmpty()) {
                    PSM topscore = sortedPSMs.pollLastEntry().getValue();
                    total++;
                    boolean isdecoy = true;
                    for (String ProID : topscore.ParentProtIDs) {
                        if (!ProID.startsWith("rev_")) {
                            isdecoy = false;
                            break;
                        }
                    }
                    if (isdecoy) {
                        decoy++;
                    }
                    if (((float) decoy / (float) total) > errorrate) {
                        prob = topscore.Probability;
                        break;
                    }
                }
                writer.write(total + "," + decoy + "," + prob + ",");
                total = 0;
                decoy = 0;
                prob = 1f;
                while (!sortedPSMs2.isEmpty()) {
                    PSM topscore = sortedPSMs2.pollLastEntry().getValue();
                    total++;
                    boolean isdecoy = true;
                    for (String ProID : topscore.ParentProtIDs) {
                        if (!ProID.startsWith("rev_")) {
                            isdecoy = false;
                            break;
                        }
                    }
                    if (isdecoy) {
                        decoy++;
                    }
                    if (((float) decoy / (float) total) > errorrate) {
                        prob = topscore.Probability;
                        break;
                    }
                }
                writer.write(total + "," + decoy + "," + prob + "\n");
            }
        }
        writer.close();
    }
}
