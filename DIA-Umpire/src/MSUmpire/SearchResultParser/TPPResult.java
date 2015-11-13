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
package MSUmpire.SearchResultParser;

import MSUmpire.PSMDataStructure.LCMSID;
import MSUmpire.PSMDataStructure.PepIonID;
import java.io.IOException;
import java.util.ArrayList;
import javax.xml.parsers.ParserConfigurationException;
import org.apache.commons.io.FilenameUtils;
import org.apache.commons.lang.exception.ExceptionUtils;
import org.apache.log4j.Logger;
import org.xml.sax.SAXException;
import org.xmlpull.v1.XmlPullParserException;

/**
 *
 * @author Chih-Chiang Tsou <chihchiang.tsou@gmail.com>
 */
public class TPPResult implements Runnable {

    public TPPResult(float fdr, float probFDR,String DecoyTag) {
        this.FDR = fdr;
        this.ProtFDR = probFDR;
        this.DecoyTag=DecoyTag;
    }
    LCMSID lcmsid;
    String pepxml;
    String protxml;
    public boolean FilterIDBymzXMLname = false;

    public void SetData(LCMSID lcmsid, String pepxml, String protxml) {
        this.lcmsid = lcmsid;
        this.pepxml = pepxml;
        this.protxml = protxml;
    }

    String DecoyTag;
    float FDR = 0.01f;
    float ProtFDR = 0.01f;

    public void ReadSearchResult(LCMSID lcmsid, String pepxml, String protxml) throws ParserConfigurationException, IOException, SAXException, XmlPullParserException, ClassNotFoundException, InterruptedException {
        PepXMLParser pepxmlparser = new PepXMLParser(lcmsid, pepxml, 0f);
        pepxmlparser.FilteredID = FilterIDBymzXMLname;
        lcmsid.FilterByPepDecoyFDR(DecoyTag, FDR);
        ProtXMLParser protxmlparser = new ProtXMLParser(lcmsid, protxml, 0f);
        //lcmsid.FilterByProteinDecoyFDR(DecoyTag, ProtFDR);
        lcmsid.RemoveLowLocalPWProtein(0.5f);
        //IDsummary.RemoveLowMaxIniProbProtein(0.9f);
        lcmsid.FilterByProteinDecoyFDRUsingMaxIniProb(DecoyTag, ProtFDR);
        lcmsid.LoadSequence();
        lcmsid.ReMapProPep();
        lcmsid.CreateInstanceForAllPepIon();
        Logger.getRootLogger().info("Protein No.:" + lcmsid.ProteinList.size() + "; Assigned Peptide No.:" + lcmsid.AssignedPepIonList.size() + "; All peptide No.:" + lcmsid.GetPepIonList().size() + "; Spectrum level threshold: " + lcmsid.SpecProbThreshold + "; Peptide level threshold: " + lcmsid.PepProbThreshold + "; Protein level threshold: " + lcmsid.ProteinProbThreshold );
    }
   
    @Override
    public void run() {
        try {
            ReadSearchResult(lcmsid, pepxml, protxml);
        } catch (Exception ex) {
            Logger.getRootLogger().error(ExceptionUtils.getStackTrace(ex));
        }
    }
}
