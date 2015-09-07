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

import MSUmpire.PSMDataStructure.LCMSID;
import MSUmpire.PSMDataStructure.PepIonID;
import MSUmpire.SearchResultParser.PepXMLParser;
import java.io.*;
import java.sql.SQLException;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.concurrent.ExecutionException;
import java.util.zip.DataFormatException;
import javax.xml.parsers.ParserConfigurationException;
import org.apache.commons.io.FilenameUtils;
import org.xml.sax.SAXException;

/**
 *
 * @author Chih-Chiang Tsou <chihchiang.tsou@gmail.com>
 */
public class CombineTwoPepXML {

    /**
     * @param args the command line arguments
     *
     */
    public static void main(String[] args) throws InterruptedException, FileNotFoundException, ExecutionException, IOException, ParserConfigurationException, DataFormatException, SAXException, SQLException, Exception {

        String ResultFolder = "C:/Umich/My Box Files/Default Sync Folder/Prelim/";

        String Filename2 = "F:\\Data\\PepideAtlas\\PAe003776_Search_Results_7621_201211210846\\interact-prob.pep.xml";
        String Filename1 = "F:\\Data\\PepideAtlas\\PAe003736_Search_Results_7499_201211210941\\interact-prob.pep.xml";
        String Filename3 = "F:\\Data\\PepideAtlas\\PAe003765_Search_Results_7610_201211210622\\interact-prob.pep.xml";
        String Filename4 = "F:\\Data\\PepideAtlas\\PAe003773_Search_Results_7618_201211210423\\interact-prob.pep.xml";
        String mzxml1 = "846";
        String mzxml2 = "941";
        String mzxml3 = "622";
        String mzxml4 = "423";

        LCMSID IDsummary1 = new LCMSID(mzxml1,"rev","");
        LCMSID IDsummary2 = new LCMSID(mzxml2,"rev","");
        LCMSID IDsummary3 = new LCMSID(mzxml3,"rev","");
        LCMSID IDsummary4 = new LCMSID(mzxml4,"rev","");
        //TPPResult tppresult = new TPPResult();

        ArrayList<LCMSID> SummaryList = new ArrayList<>();
        SummaryList.add(IDsummary1);
        PepXMLParser pepxmlparser = new PepXMLParser(IDsummary1, Filename1, 0f);
        SummaryList.add(IDsummary2);
        pepxmlparser = new PepXMLParser(IDsummary2, Filename2, 0f);
        SummaryList.add(IDsummary3);
        pepxmlparser = new PepXMLParser(IDsummary3, Filename3, 0f);
        SummaryList.add(IDsummary4);
        pepxmlparser = new PepXMLParser(IDsummary4, Filename4, 0f);
        //tppresult.ReadSearchResult(IDsummary1, Filename1, Filename1.replace("pep.xml", "prot.xml"));
        //tppresult.ReadSearchResult(IDsummary2, Filename2, Filename2.replace("pep.xml", "prot.xml"));

        for (LCMSID IDsummary : SummaryList) {
            IDsummary.FilterByPepDecoyFDR("rev_", 0.01f);
            System.out.print("Protein No.:" + IDsummary.ProteinList.size() + "; Assigned Peptide No.:" + IDsummary.AssignedPepIonList.size() + "; All peptide No.:" + IDsummary.GetPepIonList().size() + "\n");
        }

        /////////////////////////////////////
        HashSet<String> IdentifiedPepMap = new HashSet<>();
        FileWriter writer = new FileWriter(ResultFolder + "ALLPeptideSummary.xls");
        FileWriter writer3 = new FileWriter(ResultFolder + "IDNoSummary.xls");

        writer3.write("File\tProteins\tNo. peptides\tNo. assigned peptides\n");
        //////////////////////////////////////////////////////////////////

        for (LCMSID IDsummary : SummaryList) {
            writer3.write(FilenameUtils.getBaseName(IDsummary.mzXMLFileName) + "\t" + IDsummary.ProteinList.size() + "\t" + IDsummary.GetPepIonList().size() + "\t" + IDsummary.AssignedPepIonList.size() + "\n");
            for (String key : IDsummary.GetPepIonList().keySet()) {
                if (!IdentifiedPepMap.contains(key)) {
                    IdentifiedPepMap.add(key);
                }
            }
        }
        writer3.close();

        writer.write("Peptide Key\tmz\tCharge\t");
        for (LCMSID IDSummary : SummaryList) {
            String file = FilenameUtils.getBaseName(IDSummary.mzXMLFileName);
            writer.write(file + "_Prob\t" + file + "_PSMs\t" + file + "_RT\t" + file + "_RTSD\t");
        }
        writer.write("\n");

        for (String key : IdentifiedPepMap) {
            writer.write(key + "\t");
            for (LCMSID IDSummary : SummaryList) {
                if (IDSummary.GetPepIonList().containsKey(key)) {
                    PepIonID peptide = IDSummary.GetPepIonList().get(key);
                    writer.write(peptide.NeutralPrecursorMz() + "\t" + peptide.Charge + "\t");
                    break;
                }
            }
            for (LCMSID IDSummary : SummaryList) {
                if (IDSummary.GetPepIonList().containsKey(key)) {
                    PepIonID peptide = IDSummary.GetPepIonList().get(key);
                    writer.write(peptide.MaxProbability + "\t" + peptide.GetPSMList().size() + "\t" + peptide.GetIDRT() + "\t" + peptide.GetRTSD() + "\t");
                } else {
                    writer.write("\t\t\t\t");
                }
            }
            writer.write("\n");
        }
        writer.close();
    }
}
