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
    public void ReadSearchResultByRefPepID(LCMSID lcmsid, String pepxml, String protxml, LCMSID RefPepID) throws ParserConfigurationException, IOException, SAXException, XmlPullParserException, ClassNotFoundException, InterruptedException {
  
        LCMSID pepxmlid = new LCMSID(pepxml, lcmsid.DecoyTag, lcmsid.FastaPath);
        PepXMLParser pepxmlparser = new PepXMLParser(pepxmlid, pepxml, 0f);
        for (PepIonID pepion : pepxmlid.GetPepIonList().values()) {
            if (RefPepID.GetPepIonList().containsKey(pepion.GetKey())) {
                lcmsid.AddPeptideID(pepion);
            }
        }
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
    
    public void ReadSearchResultByRefIDProt(LCMSID lcmsid, String pepxml, LCMSID RefID) throws ParserConfigurationException, IOException, SAXException, XmlPullParserException, ClassNotFoundException, InterruptedException {
        PepXMLParser pepxmlparser = new PepXMLParser(lcmsid, pepxml, 0f);
        pepxmlparser.FilteredID = FilterIDBymzXMLname;
        lcmsid.FilterByPepDecoyFDR(DecoyTag, FDR);
        lcmsid.GenerateProteinByRefIDByPepSeq(RefID,false);
        lcmsid.LoadSequence();
        lcmsid.ReMapProPep();
        lcmsid.CreateInstanceForAllPepIon();
        Logger.getRootLogger().info("Protein No.:" + lcmsid.ProteinList.size() + "; Assigned Peptide No.:" + lcmsid.AssignedPepIonList.size() + "; All peptide No.:" + lcmsid.GetPepIonList().size() + "; Spectrum level threshold: " + lcmsid.SpecProbThreshold + "; Peptide level threshold: " + lcmsid.PepProbThreshold + "; Protein level threshold: " + lcmsid.ProteinProbThreshold );
    }
    public void ReadSearchResultByRefID(LCMSID lcmsid, String pepxml, String protxml, LCMSID RefID) throws ParserConfigurationException, IOException, SAXException, XmlPullParserException, ClassNotFoundException, InterruptedException {
        PepXMLParser pepxmlparser = new PepXMLParser(lcmsid, pepxml, 0f);
        pepxmlparser.FilteredID = FilterIDBymzXMLname;
        lcmsid.FilterByPepDecoyFDR(DecoyTag, FDR);
        ProtXMLParser protxmlparser = new ProtXMLParser(lcmsid, protxml, 0f);
        lcmsid.FilterProteinByRefID(RefID);
        lcmsid.LoadSequence();
        lcmsid.ReMapProPep();
        lcmsid.CreateInstanceForAllPepIon();
        Logger.getRootLogger().info("Protein No.:" + lcmsid.ProteinList.size() + "; Assigned Peptide No.:" + lcmsid.AssignedPepIonList.size() + "; All peptide No.:" + lcmsid.GetPepIonList().size() + "; Spectrum level threshold: " + lcmsid.SpecProbThreshold + "; Peptide level threshold: " + lcmsid.PepProbThreshold + "; Protein level threshold: " + lcmsid.ProteinProbThreshold );
    }

    public void ReadSearchResultAndFilterByProb(LCMSID lcmsid, String pepxml, String protxml, float prob) throws ParserConfigurationException, IOException, SAXException, XmlPullParserException, ClassNotFoundException, InterruptedException {
        PepXMLParser pepxmlparser = new PepXMLParser(lcmsid, pepxml, prob);
        pepxmlparser.FilteredID = FilterIDBymzXMLname;
        ProtXMLParser protxmlparser = new ProtXMLParser(lcmsid, protxml, prob);
        lcmsid.LoadSequence();
        lcmsid.ReMapProPep();
        lcmsid.CreateInstanceForAllPepIon();
        Logger.getRootLogger().info("Protein No.:" + lcmsid.ProteinList.size() + "; Assigned Peptide No.:" + lcmsid.AssignedPepIonList.size() + "; All peptide No.:" + lcmsid.GetPepIonList().size() + "; Spectrum level threshold: " + lcmsid.SpecProbThreshold + "; Peptide level threshold: " + lcmsid.PepProbThreshold + "; Protein level threshold: " + lcmsid.ProteinProbThreshold );
    }

    public void ReadSearchResultByRefPepID(LCMSID lcmsid, ArrayList<String> pepxmls, String protxml, final LCMSID RefPepID) throws ParserConfigurationException, IOException, SAXException, XmlPullParserException, ClassNotFoundException, InterruptedException {
        for (String pepxml : pepxmls) {
            LCMSID pepxmlid = new LCMSID(pepxml, lcmsid.DecoyTag, lcmsid.FastaPath);
            PepXMLParser pepxmlparser = new PepXMLParser(pepxmlid, pepxml, 0f);            
            for (PepIonID pepion : pepxmlid.GetPepIonList().values()) {
                if (RefPepID.GetPepIonList().containsKey(pepion.GetKey())) {
                    lcmsid.AddPeptideID(pepion);
                }
            }
        }
        ProtXMLParser protxmlparser = new ProtXMLParser(lcmsid, protxml, 0f);
        lcmsid.RemoveLowLocalPWProtein(0.5f);
        lcmsid.FilterByProteinDecoyFDRUsingMaxIniProb(DecoyTag, ProtFDR);
        lcmsid.LoadSequence();
        lcmsid.ReMapProPep();
        lcmsid.CreateInstanceForAllPepIon();
        Logger.getRootLogger().info("Protein No.:" + lcmsid.ProteinList.size() + "; Assigned Peptide No.:" + lcmsid.AssignedPepIonList.size() + "; All peptide No.:" + lcmsid.GetPepIonList().size() + "; Spectrum level threshold: " + lcmsid.SpecProbThreshold + "; Peptide level threshold: " + lcmsid.PepProbThreshold + "; Protein level threshold: " + lcmsid.ProteinProbThreshold );
    }
    public void ReadSearchResult(LCMSID lcmsid, ArrayList<String> pepxmls, String protxml) throws ParserConfigurationException, IOException, SAXException, XmlPullParserException, ClassNotFoundException, InterruptedException {
        for (String pepxml : pepxmls) {
            LCMSID pepxmlid = new LCMSID(FilenameUtils.getFullPath(pepxml) + FilenameUtils.getBaseName(pepxml),DecoyTag,lcmsid.FastaPath);
            PepXMLParser pepxmlparser = new PepXMLParser(pepxmlid, pepxml, 0f);
            pepxmlid.FilterByPepDecoyFDR(DecoyTag, FDR);
            Logger.getRootLogger().info( "peptide No.:" + pepxmlid.GetPepIonList().size() + "; Peptide level threshold: " + pepxmlid.PepProbThreshold);
            for (PepIonID pepID : pepxmlid.GetPepIonList().values()) {
                lcmsid.AddPeptideID(pepID);
            }
        }
        ProtXMLParser protxmlparser = new ProtXMLParser(lcmsid, protxml, 0f);
        lcmsid.RemoveLowLocalPWProtein(0.5f);
        //lcmsid.RemoveLowMaxIniProbProtein(0.9f);
        lcmsid.FilterByProteinDecoyFDRUsingMaxIniProb(DecoyTag, ProtFDR);
        lcmsid.LoadSequence();
        lcmsid.ReMapProPep();
        lcmsid.CreateInstanceForAllPepIon();
        Logger.getRootLogger().info("Protein No.:" + lcmsid.ProteinList.size() + "; Assigned Peptide No.:" + lcmsid.AssignedPepIonList.size() + "; All peptide No.:" + lcmsid.GetPepIonList().size() + "; Spectrum level threshold: " + lcmsid.SpecProbThreshold + "; Peptide level threshold: " + lcmsid.PepProbThreshold + "; Protein level threshold: " + lcmsid.ProteinProbThreshold);
    }
    
     public void ReadSearchResultByRefIDUseRefProtID(LCMSID lcmsid, ArrayList<String> pepxmls, LCMSID RefID,boolean UseMappIon) throws ParserConfigurationException, IOException, SAXException, XmlPullParserException, ClassNotFoundException, InterruptedException {
        for (String pepxml : pepxmls) {
            LCMSID pepxmlid = new LCMSID(FilenameUtils.getFullPath(lcmsid.mzXMLFileName) + FilenameUtils.getBaseName(lcmsid.mzXMLFileName),DecoyTag,lcmsid.FastaPath);
            PepXMLParser pepxmlparser = new PepXMLParser(pepxmlid, pepxml, 0f);
            pepxmlid.FilterByPepDecoyFDR(DecoyTag, FDR);
            Logger.getRootLogger().info( "peptide No.:" + pepxmlid.GetPepIonList().size() + "; Peptide level threshold: " + pepxmlid.PepProbThreshold);
            for (PepIonID pepID : pepxmlid.GetPepIonList().values()) {
                lcmsid.AddPeptideID(pepID);
            }           
        }        
        lcmsid.GenerateProteinByRefIDByPepSeq(RefID,UseMappIon);
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
