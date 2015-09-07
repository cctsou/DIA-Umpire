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
import MSUmpire.PSMDataStructure.ProtID;
import MSUmpire.SearchResultParser.PepXMLParser;
import MSUmpire.SearchResultParser.ProtXMLParser;
import MSUmpire.SearchResultParser.TPPResult;
import java.io.IOException;
import javax.xml.parsers.ParserConfigurationException;
import org.apache.log4j.ConsoleAppender;
import org.apache.log4j.Level;
import org.apache.log4j.Logger;
import org.apache.log4j.PatternLayout;
import org.xml.sax.SAXException;
import org.xmlpull.v1.XmlPullParserException;

/**
 *
 * @author Chih-Chiang Tsou <chihchiang.tsou@gmail.com>
 */
public class ParserTest {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) throws ParserConfigurationException, SAXException, IOException, XmlPullParserException, ClassNotFoundException, InterruptedException {

        Logger logger = Logger.getRootLogger();
        ConsoleAppender ca = new ConsoleAppender();
        ca.setThreshold(Level.INFO);
        ca.setName("ConsoleLogger_Info");
        ca.setLayout(new PatternLayout("%d %-5p [%c{1}] %m%n"));
        ca.activateOptions();

        logger.getLoggerRepository().resetConfiguration();
        logger.addAppender(ca);

        String protxml = "F:\\Data\\Fe_DIA\\Fusion\\DDA\\Hela1ug_QC_141226_01\\interact-Hela1ug_QC_141226_01.iproph.prot.xml";
        String pepxml = "F:\\Data\\Fe_DIA\\Fusion\\DDA\\Hela1ug_QC_141226_01\\interact-Hela1ug_QC_141226_01.iproph.pep.xml";
        
        String fasta = "F:\\Data\\Fe_DIA\\swissprot_Hs_plusREV.2013Jan09.fa";
        LCMSID protID = new LCMSID(protxml,"rev",fasta);
        protID.FastaPath = fasta;
        protID.DecoyTag="rev";
        protID.FDR=0.01f;
        
        TPPResult result=new TPPResult(0.01f, 0.01f, "ref");
        result.ReadSearchResult(protID, pepxml, protxml);
        
        int count=0;
        for(ProtID protein : protID.ProteinList.values()){
            if(protein.GetSpectralCount()>1){
                count++;
            }
        }
        logger.info("Protein No.:" + protID.ProteinList.size());        
    }
}
