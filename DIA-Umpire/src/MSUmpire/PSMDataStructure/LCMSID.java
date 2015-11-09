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
package MSUmpire.PSMDataStructure;

import MSUmpire.BaseDataStructure.InstrumentParameter;
import MSUmpire.BaseDataStructure.ScanCollection;
import MSUmpire.BaseDataStructure.ScanData;
import MSUmpire.BaseDataStructure.XYData;
import MSUmpire.FastaParser.FastaParser;
import MSUmpire.PeakDataStructure.PeakCluster;
import MSUmpire.SortedListLib.SortedList;
import com.compomics.util.experiment.biology.AminoAcid;
import com.compomics.util.experiment.biology.Ion;
import com.compomics.util.experiment.biology.PTM;
import com.compomics.util.experiment.biology.ions.ElementaryIon;
import com.compomics.util.experiment.biology.ions.PeptideFragmentIon;
import org.nustaq.serialization.FSTObjectInput;
import org.nustaq.serialization.FSTObjectOutput;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;
import java.io.Serializable;
import java.sql.Connection;
import java.sql.ResultSet;
import java.sql.SQLException;
import java.sql.Statement;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.HashMap;
import java.util.logging.Level;
import org.apache.commons.io.FilenameUtils;
import org.apache.commons.lang.exception.ExceptionUtils;
import org.apache.log4j.Logger;
import org.xmlpull.v1.XmlPullParserException;

/**
 *
 * @author Chih-Chiang Tsou <chihchiang.tsou@gmail.com>
 */
public class LCMSID implements Serializable {

    private static final long serialVersionUID = 25494643464616L;

    public HashMap<String, PSM> PSMList;
    public HashMap<String, PSM> LowScorePSMByPepKey;
    public transient HashMap<String, PepIonID> LowScorePep;
    public transient HashMap<String, PSM> LowScorePSM;
    public transient HashMap<Integer, ArrayList<ProtID>> ProteinGroups;
    private HashMap<String, PepIonID> PepIonList;
    private HashMap<Integer, PepIonID> PepIonIndexList;
    public HashMap<String, HashMap<String, PepIonID>> PeptideList;
//    
    private HashMap<String, PepIonID> MappedPepIonList;
    private HashMap<Integer, PepIonID> MappedPepIonIndexList;
    public HashMap<String, HashMap<String, PepIonID>> MappedPeptideList;
//    
//    private HashMap<String, PepIonID> ExtLibPepIonList;
//    private HashMap<Integer, PepIonID> ExtLibPepIonIndexList;
//    public HashMap<String, HashMap<String, PepIonID>> ExtLibPeptideSeqList;

    
    
    public HashMap<String, PepIonID> AssignedPepIonList;
    public HashMap<String, PepIonID> ProtXMLPepIonList;
    
    public HashMap<String, ProtID> ProteinList;
    public HashMap<String, ProtID> IndisProteinIDList;
    
    public HashMap<String, ModificationInfo> ModificationList;
    public HashMap<String, ProtID> PepXMLProteinList;
    
    
    public String DataBase;
    public String SearchEngine;
    public String msModel;
    public String msManufacturer;
    public String msIonization;
    public String msMassAnalyzer;
    public String msDetector;
    public String mzXMLFileName;
    public float FDR = Float.POSITIVE_INFINITY;
    public float ProteinFDR = Float.POSITIVE_INFINITY;
    public float ExpectThreshold = Float.POSITIVE_INFINITY;
    public float SpecProbThreshold = 0f;
    public float PepProbThreshold = 0f;
    public float ProteinProbThreshold = 0f;
    public String DecoyTag = "rev_";
    public String FastaPath;
    private float NorFactor = 1f;
    //SequenceFactory sequenceFactory = null;
    FastaParser fastaParser;
    //transient FastaParser_V2 fastaParser;//removed 02032015
    public HashMap<String, String> LuciphorResult;
    public String Filename; //added 0828, needs to set as transient for older serialization

    private FastaParser GetFastaParser() {
        if (fastaParser == null) {
            fastaParser = new FastaParser(FastaPath);
        }
        return fastaParser;
    }
  
    public void WriteLCMSIDSerialization(String filepath) {
        WriteLCMSIDSerialization(filepath, "");
    }

    public void WriteLCMSIDSerialization(String filepath, String tag) {
        //JavaSerializationWrite(filepath);
        if (!FSWrite(filepath, tag)) {
            Logger.getRootLogger().debug("Writing LCMSID FS failed. writing standard serialization instead");
            JavaSerializationWrite(filepath, tag);
        }
    }

    public void RemoveLowWeightPep(float threshold) {
        for (ProtID protein : ProteinList.values()) {
            protein.RemoveLowWeightPepID(threshold);
        }
    }

    public LCMSID CreateEmptyLCMSID() {
        LCMSID newLcmsid = new LCMSID(mzXMLFileName, DecoyTag, FastaPath);
        newLcmsid.SearchEngine = SearchEngine;
        newLcmsid.msDetector = msDetector;

        newLcmsid.msIonization = msIonization;
        newLcmsid.msManufacturer = msManufacturer;
        newLcmsid.msMassAnalyzer = msMassAnalyzer;
        newLcmsid.msModel = msModel;
        newLcmsid.DataBase = DataBase;
        for (ModificationInfo mod : ModificationList.values()) {
            newLcmsid.ModificationList.put(mod.GetKey(), mod);
        }
        return newLcmsid;
    }
    
    public HashMap<String,LCMSID> GetLCMSIDFileMap(){
        
        HashMap<String,LCMSID> lcmsidmap=new HashMap<>();
        
        for(PSM psm : this.PSMList.values()){
            if (!lcmsidmap.containsKey(psm.GetRawNameString())) {
                LCMSID newLcmsid = CreateEmptyLCMSID();
                newLcmsid.mzXMLFileName = psm.GetRawNameString();
                lcmsidmap.put(psm.GetRawNameString(), newLcmsid);
            }            
            lcmsidmap.get(psm.GetRawNameString()).AddPSM(psm);            
        }        
        return lcmsidmap;
    }
    
    public LCMSID ExtractIDByFileName(String basename){
        LCMSID newLcmsid=CreateEmptyLCMSID();
        newLcmsid.mzXMLFileName=basename;
        for(PSM psm : this.PSMList.values()){
            if(psm.GetRawNameString().equals(basename)){
                newLcmsid.AddPSM(psm);
            }
        }        
        return newLcmsid;
    }
    
    private void JavaSerializationWrite(String filepath, String tag) {
        try {
            if (!tag.equals("")) {
                tag = "_" + tag;
            }
            Logger.getRootLogger().info("Writing ID results to file:" + FilenameUtils.getFullPath(filepath) + FilenameUtils.getBaseName(filepath) + tag + "_LCMSID.ser...");
            FileOutputStream fout = new FileOutputStream(FilenameUtils.getFullPath(filepath) + FilenameUtils.getBaseName(filepath) + tag + "_LCMSID.ser", false);
            ObjectOutputStream oos = new ObjectOutputStream(fout);
            oos.writeObject(this);
            oos.close();
            fout.close();
        } catch (Exception ex) {
            Logger.getRootLogger().error(ExceptionUtils.getStackTrace(ex));
        }
    }

    private boolean FSWrite(String filepath, String tag) {
        try {
            if (!tag.equals("")) {
                tag = "_" + tag;
            }
            Logger.getRootLogger().info("Writing ID results to file:" + FilenameUtils.getFullPath(filepath) + FilenameUtils.getBaseName(filepath) + tag + "_LCMSID.serFS...");
            FileOutputStream fout = new FileOutputStream(FilenameUtils.getFullPath(filepath) + FilenameUtils.getBaseName(filepath) + tag + "_LCMSID.serFS", false);
            FSTObjectOutput out = new FSTObjectOutput(fout);
            ReduceMemoryUsage();
            out.writeObject(this, LCMSID.class);
            out.close();
            fout.close();
        } catch (Exception ex) {
            Logger.getRootLogger().error(ExceptionUtils.getStackTrace(ex));
            return false;
        }
        return true;
    }

    private static LCMSID FS_Read(String filepath, String tag) throws Exception {
        if (!tag.equals("")) {
            tag = "_" + tag;
        }
        if (!new File(FilenameUtils.getFullPath(filepath) + FilenameUtils.getBaseName(filepath) + tag + "_LCMSID.serFS").exists()) {
            return null;
        }
        try {
            Logger.getRootLogger().info("Reading ID results from file:" + FilenameUtils.getFullPath(filepath) + FilenameUtils.getBaseName(filepath) + tag + "_LCMSID.serFS...");

            FileInputStream fileIn = new FileInputStream(FilenameUtils.getFullPath(filepath) + FilenameUtils.getBaseName(filepath) + tag + "_LCMSID.serFS");
            FSTObjectInput in = new FSTObjectInput(fileIn);
            LCMSID lcmsid = (LCMSID) in.readObject(LCMSID.class);
            in.close();
            fileIn.close();
            return lcmsid;

        } catch (Exception ex) {
            Logger.getRootLogger().info("Reading LCMSID FS results failed.");
            Logger.getRootLogger().error(ExceptionUtils.getStackTrace(ex));
            return null;
        }
    }

    private static LCMSID FS_Read_Old(String filepath, String tag) throws Exception {
        if (!tag.equals("")) {
            tag = "_" + tag;
        }
        if (!new File(FilenameUtils.getFullPath(filepath) + FilenameUtils.getBaseName(filepath) + tag + "_LCMSID.serFS").exists()) {
            return null;
        }
        try {
            Logger.getRootLogger().info("Reading ID results from file:" + FilenameUtils.getFullPath(filepath) + FilenameUtils.getBaseName(filepath) + tag + "_LCMSID.serFS...");

            FileInputStream fileIn = new FileInputStream(FilenameUtils.getFullPath(filepath) + FilenameUtils.getBaseName(filepath) + tag + "_LCMSID.serFS");
            org.nustaq_old.serialization.FSTObjectInput in = new org.nustaq_old.serialization.FSTObjectInput(fileIn);
            LCMSID lcmsid = (LCMSID) in.readObject(LCMSID.class);
            in.close();
            fileIn.close();
            return lcmsid;

        } catch (Exception ex) {
            Logger.getRootLogger().info("Reading LCMSID FS results failed.");
            Logger.getRootLogger().error(ExceptionUtils.getStackTrace(ex));
            return null;
        }
    }

    public static LCMSID ReadLCMSIDSerialization(String filepath) throws Exception {
        return ReadLCMSIDSerialization(filepath, "");
    }

    public static LCMSID ReadLCMSIDSerialization(String filepath, String tag) throws Exception {
        LCMSID lcmsid = FS_Read(filepath, tag);
        if (lcmsid == null) {
            //lcmsid = JavaSerializationRead(filepath);
            lcmsid = FS_Read_Old(filepath, tag);
            if (lcmsid != null) {
                lcmsid.WriteLCMSIDSerialization(filepath, tag);
            }
        }
        return lcmsid;
    }

    private static LCMSID JavaSerializationRead(String filepath) {
        if (!new File(FilenameUtils.getFullPath(filepath) + FilenameUtils.getBaseName(filepath) + "_LCMSID.ser").exists()) {
            return null;
        }
        try {
            Logger.getRootLogger().info("Reading ID results from file:" + FilenameUtils.getFullPath(filepath) + FilenameUtils.getBaseName(filepath) + "_LCMSID.ser...");

            FileInputStream fileIn = new FileInputStream(FilenameUtils.getFullPath(filepath) + FilenameUtils.getBaseName(filepath) + "_LCMSID.ser");
            ObjectInputStream in = new ObjectInputStream(fileIn);
            LCMSID lcmsid = (LCMSID) in.readObject();
            in.close();
            fileIn.close();
            return lcmsid;

        } catch (Exception ex) {
            Logger.getRootLogger().info("ID results don't exist, continue process.");
            Logger.getRootLogger().error(ExceptionUtils.getStackTrace(ex));
            return null;
        }
    }

    public float GetNorFactor() {
        if (NorFactor == 1f) {
            NorFactor = 0f;
            for (PepIonID pepIonID : GetPepIonList().values()) {
                NorFactor += pepIonID.PeakHeight[0];
            }
            NorFactor /= GetPepIonList().size();
        }
        return NorFactor;
    }

    public LCMSID(String mzXMLFileName, String Decoytag, String Fasta) {
        ModificationList = new HashMap<>();
        ProteinList = new HashMap<>();
        PSMList = new HashMap<>();
        PepIonList = new HashMap<>();
        AssignedPepIonList = new HashMap<>();
        ProtXMLPepIonList = new HashMap<>();
        IndisProteinIDList = new HashMap<>();
        ProteinGroups = new HashMap<>();
        PeptideList = new HashMap<>();
        MappedPepIonList = new HashMap<>();
        this.mzXMLFileName = mzXMLFileName;
        this.FastaPath = Fasta;
        this.DecoyTag = Decoytag;
    }

    public HashMap<String, PepIonID> GetPepIonList() {
        return PepIonList;
    }

    public HashMap<String, PepIonID> GetMappedPepIonList() {
        return MappedPepIonList;
    }

    public void ResetMappedPepProb() {
        for (PepIonID pepIonID : MappedPepIonList.values()) {
            pepIonID.MS1AlignmentLocalProbability = 0f;
            pepIonID.MS1AlignmentProbability = 0f;
            pepIonID.MS2AlignmentLocalProbability = 0f;
            pepIonID.MS2AlignmentProbability = 0f;
        }
    }

    public void SetMappedPepIonList(HashMap<String, PepIonID> list) {
        MappedPepIonList = list;
    }
//    public void ParseFromXML(String PepFile, float PepScore) throws ParserConfigurationException, SAXException, IOException {
//        PepXMLParser pepXMLParser = new PepXMLParser(this, PepFile, PepScore);
//    }
//    public void Parse(ArrayList<String> ProtFiles, ArrayList<String> PepFiles, float PepScore, float ProtScore) throws ParserConfigurationException, SAXException, IOException
//    {
//        for(String PepFile : PepFiles){
//         PepXMLParser pepXMLParser=new PepXMLParser(this, PepFile, PepScore);
//        }        
//         for(String ProtFile : ProtFiles){
//             ProtXMLParser protXMLParser=new ProtXMLParser(this, ProtFile, ProtScore);             
//        }
//    }

    public void GenerateFragmentPeakForPepIonByMSMS(ScanCollection scanCollection, float fragPPM) {

        double protonMass = ElementaryIon.proton.getTheoreticMass();
        for (PepIonID pepIonID : GetPepIonList().values()) {
            for (PSM psm : pepIonID.GetPSMList()) {
                ScanData scan = scanCollection.GetScan(psm.ScanNo);
                for (Ion frag : pepIonID.GetFragments()) {
                    XYData closetPeak = null;
                    float targetmz = (float) (frag.getTheoreticMass() + protonMass);
                    closetPeak = scan.GetHighestPeakInMzWindow(targetmz, fragPPM);
                    if (closetPeak != null) {
                        FragmentPeak fragmentpeak = new FragmentPeak();
                        fragmentpeak.ObservedMZ = closetPeak.getX();
                        fragmentpeak.FragMZ = targetmz;
                        fragmentpeak.intensity = closetPeak.getY();
                        fragmentpeak.Charge = 1;
                        fragmentpeak.ppm = InstrumentParameter.CalcSignedPPM(closetPeak.getX(), targetmz);
                        fragmentpeak.IonType = frag.getSubTypeAsString() + ((PeptideFragmentIon) frag).getNumber();
                        pepIonID.FragmentPeaks.add(fragmentpeak);
                    }

                    targetmz = (float) (frag.getTheoreticMass() + protonMass * 2) / 2;
                    closetPeak = scan.GetHighestPeakInMzWindow(targetmz, fragPPM);

                    if (closetPeak != null) {
                        FragmentPeak fragmentpeak = new FragmentPeak();
                        fragmentpeak.ObservedMZ = closetPeak.getX();
                        fragmentpeak.FragMZ = targetmz;
                        fragmentpeak.intensity = closetPeak.getY();
                        fragmentpeak.Charge = 2;
                        fragmentpeak.ppm = InstrumentParameter.CalcSignedPPM(closetPeak.getX(), targetmz);
                        fragmentpeak.IonType = frag.getSubTypeAsString() + ((PeptideFragmentIon) frag).getNumber();
                        pepIonID.FragmentPeaks.add(fragmentpeak);
                    }

                }
            }
        }
    }

    public void ReadFromSkylineExport(String SkylineExport) throws SQLException, FileNotFoundException, IOException {

        BufferedReader reader = new BufferedReader(new FileReader(SkylineExport));
        String line = "";
        while ((line = reader.readLine()) != null) {
            String[] Info = line.split(",");
            if (Info[0] == null ? FilenameUtils.getBaseName(mzXMLFileName) == null : FilenameUtils.getBaseName(Info[0]).equals(FilenameUtils.getBaseName(mzXMLFileName))) {

                if (!PepIonList.containsKey(Info[2] + "_" + Info[4])) {
                    PepIonID pepIonID = new PepIonID();
                    pepIonID.CreateQuantInstance(3);
                    pepIonID.Sequence = Info[1];
                    pepIonID.ModSequence = Info[2];
                    pepIonID.SetMz(Float.parseFloat(Info[6]));
                    PepIonList.put(Info[2] + "_" + Info[4], pepIonID);
                }
                PepIonID pepIon = PepIonList.get(Info[2] + "_" + Info[4]);

                if (Info[8].startsWith("precursor")) {
                    switch (Info[8]) {
                        case "precursor": {
                            pepIon.PeakArea[0] = Float.parseFloat(Info[9]);
                            pepIon.PeakHeight[0] = Float.parseFloat(Info[10]);
                            pepIon.SetRT(Float.parseFloat(Info[3]));
                            break;
                        }
                        case "precursor [M+1]": {
                            pepIon.PeakArea[1] = Float.parseFloat(Info[9]);
                            pepIon.PeakHeight[1] = Float.parseFloat(Info[10]);
                            break;
                        }
                        case "precursor [M+2]": {
                            pepIon.PeakArea[2] = Float.parseFloat(Info[9]);
                            pepIon.PeakHeight[2] = Float.parseFloat(Info[10]);
                            break;
                        }
                    }
                } else {
                    FragmentPeak fragment = new FragmentPeak();
                    fragment.IonType = Info[8];
                    fragment.Charge = Integer.parseInt(Info[5]);
                    fragment.intensity = Float.parseFloat(Info[10]);
                    pepIon.FragmentPeaks.add(fragment);
                }
            }
        }
        reader.close();
        for (PepIonID pepIonID : GetPepIonList().values()) {
            if (!PeptideList.containsKey(pepIonID.Sequence)) {
                PeptideList.put(pepIonID.Sequence, new HashMap<String, PepIonID>());
            }
            PeptideList.get(pepIonID.Sequence).put(pepIonID.GetKey(), pepIonID);
        }
    }

    public void ReadFromDBPepIon(Connection connection) throws SQLException {
        Logger.getRootLogger().info("Loading ID result from MySQL DB :" + FilenameUtils.getBaseName(mzXMLFileName) + "_PepIonIDs");
        PepIonIndexList = new HashMap<>();
        Statement state = connection.createStatement();
        ResultSet rsQuery = state.executeQuery("SELECT * FROM " + FilenameUtils.getBaseName(mzXMLFileName) + "_PepIonIDs");
        while (rsQuery.next()) {
            PepIonID pepIonID = new PepIonID();
            pepIonID.Index = rsQuery.getInt("PepIndex");
            pepIonID.Sequence = rsQuery.getString("Sequence");
            pepIonID.ModSequence = rsQuery.getString("ModSeq");
            try {
                //pepIonID.TPPModSeq = rsQuery.getString("TPPModSeq");
            } catch (Exception e) {
            }
            pepIonID.Is_NonDegenerate = (rsQuery.getInt("IsNonDegenerate") == 1 ? true : false);
            pepIonID.Charge = rsQuery.getInt("Charge");
            pepIonID.SetMz(rsQuery.getFloat("mz"));
            try {
                pepIonID.SetRT(rsQuery.getFloat("IDRT"));
            } catch (Exception e) {
                pepIonID.SetRT(rsQuery.getFloat("RT"));
            }
            try {
                pepIonID.PeakRT = rsQuery.getFloat("PeakRT");
            } catch (Exception e) {
            }
            pepIonID.CreateQuantInstance(3);
            pepIonID.PeakArea[0] = rsQuery.getFloat("PeakArea1");
            pepIonID.PeakArea[1] = rsQuery.getFloat("PeakArea2");
            pepIonID.PeakArea[2] = rsQuery.getFloat("PeakArea3");
            pepIonID.PeakHeight[0] = rsQuery.getFloat("PeakHeight1");
            pepIonID.PeakHeight[1] = rsQuery.getFloat("PeakHeight2");
            pepIonID.PeakHeight[2] = rsQuery.getFloat("PeakHeight3");

            pepIonID.MS1ClusIndex = rsQuery.getString("MS1ClusIndex");
            pepIonID.MS2ClusIndex = rsQuery.getString("MS2ClusIndex");

            try {
                pepIonID.PeakClusterScore = rsQuery.getFloat("PeakScore");
            } catch (SQLException e) {
            }

//            float score=Float.NEGATIVE_INFINITY;
//            for(PeakCluster cluster : pepIonID.MS1PeakClusters){
//                if(cluster.MS1Score>score){
//                    score=cluster.MS1Score;
//                }
//            }
//            pepIonID.PeakClusterScore=score;
            PepIonList.put(pepIonID.GetKey(), pepIonID);
            PepIonIndexList.put(pepIonID.Index, pepIonID);
        }
        GeneratePepSeqList();
    }

    public void DiscardMappedPepIonTable(Connection connection) throws SQLException {
        Statement state = connection.createStatement();
        try {
            state.execute("DROP TABLE " + FilenameUtils.getBaseName(mzXMLFileName) + "_MappedPepIonIDs");
        } catch (SQLException ex) {
            Logger.getRootLogger().error("(Discarding for " + FilenameUtils.getBaseName(mzXMLFileName) + "_MappedPepIonIDs failed.)");
            Logger.getRootLogger().error(ExceptionUtils.getStackTrace(ex));
        }
    }

    public void ReadFromDBMappedPepIon(Connection connection) throws SQLException {
        ReadFromDBMappedPepIon(connection, false, 0f);
    }

    public void ReadFromDBMappedPepIon(Connection connection, boolean filterbyProb, float MappedIonProbThreshold) throws SQLException {
        Logger.getRootLogger().info("Loading ID result from MySQL DB :" + FilenameUtils.getBaseName(mzXMLFileName) + "_MappedPepIonIDs");
        Statement state = connection.createStatement();
        MappedPepIonList = new HashMap<>();
        MappedPepIonIndexList = new HashMap<>();
        try {

            ResultSet rsQuery = state.executeQuery("SELECT * FROM " + FilenameUtils.getBaseName(mzXMLFileName) + "_MappedPepIonIDs");
            while (rsQuery.next()) {
                PepIonID pepIonID = new PepIonID();
                pepIonID.MS1ClusIndex = rsQuery.getString("MS1ClusIndex");
                pepIonID.MS2ClusIndex = rsQuery.getString("MS2ClusIndex");
                pepIonID.Sequence = rsQuery.getString("Sequence");
                pepIonID.ModSequence = rsQuery.getString("ModSeq");

                try {
                    //pepIonID.TPPModSeq = rsQuery.getString("TPPModSeq");
                    pepIonID.Index = rsQuery.getInt("PepIndex");
                    pepIonID.PeakRT = rsQuery.getFloat("PeakRT");
                    pepIonID.PeakClusterScore = rsQuery.getFloat("PeakScore");
                    pepIonID.Modifications = PTMManager.TranslateModificationString(rsQuery.getString("ModInfo"));
                } catch (Exception ex) {
                    Logger.getRootLogger().error(ExceptionUtils.getStackTrace(ex));
                }

                try {
                    pepIonID.MS1AlignmentProbability = rsQuery.getFloat("MS1AlignmentProb");
                    pepIonID.MS1AlignmentLocalProbability = rsQuery.getFloat("MS1AlignmentLProb");
                    pepIonID.MS2AlignmentProbability = rsQuery.getFloat("MS2AlignmentProb");
                    pepIonID.MS2AlignmentLocalProbability = rsQuery.getFloat("MS2AlignmentLProb");
                } catch (Exception e) {
                }

                pepIonID.Charge = rsQuery.getInt("Charge");
                pepIonID.SetMz(rsQuery.getFloat("mz"));

                String RTs = rsQuery.getString("PredictRT");

                for (String split : RTs.split(";")) {
                    if (!"".equals(split)) {
                        pepIonID.PredictRT.add(Float.parseFloat(split));
                    }
                }

                pepIonID.CreateQuantInstance(3);
                pepIonID.PeakArea[0] = rsQuery.getFloat("PeakArea1");
                pepIonID.PeakArea[1] = rsQuery.getFloat("PeakArea2");
                pepIonID.PeakArea[2] = rsQuery.getFloat("PeakArea3");
                pepIonID.PeakHeight[0] = rsQuery.getFloat("PeakHeight1");
                pepIonID.PeakHeight[1] = rsQuery.getFloat("PeakHeight2");
                pepIonID.PeakHeight[2] = rsQuery.getFloat("PeakHeight3");

                if (!filterbyProb || pepIonID.TargetedProbability() > MappedIonProbThreshold) {
                    MappedPepIonList.put(pepIonID.GetKey(), pepIonID);
                    MappedPepIonIndexList.put(pepIonID.Index, pepIonID);
                }
            }
        } catch (Exception ex) {
            Logger.getRootLogger().error("Loading " + FilenameUtils.getBaseName(mzXMLFileName) + "_MappedPepIonIDs failed");
            Logger.getRootLogger().error(ExceptionUtils.getStackTrace(ex));
        }
        GenerateMappedPepSeqList();
    }

    public void ReadFromCSVMappedPepIon() throws FileNotFoundException, IOException {
        Logger.getRootLogger().info("Loading ID result from CSV :" + FilenameUtils.getBaseName(mzXMLFileName) + "_MappedPepIonIDs");

        BufferedReader reader = new BufferedReader(new FileReader(FilenameUtils.getFullPath(mzXMLFileName) + FilenameUtils.getBaseName(mzXMLFileName) + "_MappedPepIonIDs.csv"));

        MappedPepIonList = new HashMap<>();
        MappedPepIonIndexList = new HashMap<>();
        String line = "";
        reader.readLine();
        while ((line = reader.readLine()) != null) {
            PepIonID pepIonID = new PepIonID();

            String[] info = line.split(",");
            pepIonID.Index = Integer.parseInt(info[0]);
            pepIonID.Sequence = info[1];
            pepIonID.ModSequence = info[2];
            pepIonID.TPPModSeq = info[3];
            pepIonID.Modifications = PTMManager.TranslateModificationString(info[4]);
            pepIonID.Charge = Integer.parseInt(info[5]);

            pepIonID.SetMz(Float.parseFloat(info[6]));

            String RTs = info[7];

            for (String split : RTs.split(";")) {
                if (!"".equals(split)) {
                    pepIonID.PredictRT.add(Float.parseFloat(split));
                }
            }

            pepIonID.PeakRT = Float.parseFloat(info[8]);

            MappedPepIonList.put(pepIonID.GetKey(), pepIonID);
            MappedPepIonIndexList.put(pepIonID.Index, pepIonID);
        }

        GenerateMappedPepSeqList();
    }

    public void ReadFromDBPSM(Connection connection) throws SQLException {
        Logger.getRootLogger().info("Loading ID result from MySQL DB :" + FilenameUtils.getBaseName(mzXMLFileName) + "_PSMs");
        Statement state = connection.createStatement();
        ResultSet rsQuery = state.executeQuery("SELECT * FROM " + FilenameUtils.getBaseName(mzXMLFileName) + "_PSMs");
        while (rsQuery.next()) {
            PSM psm = new PSM();
            psm.SpecNumber = rsQuery.getString("SpecID");
            psm.Sequence = rsQuery.getString("Sequence");
            psm.ModSeq = rsQuery.getString("ModSeq");
            try {
                psm.TPPModSeq = rsQuery.getString("TPPModSeq");
                psm.NeutralPepMass = rsQuery.getFloat("NeutralPepMass");
            } catch (Exception e) {
            }
            try {
                psm.LuciphorScore = rsQuery.getFloat("LuciphorScore");
                psm.LuciphorLFLR = rsQuery.getFloat("LuciphorLFLR");
                psm.LuciphorFLR = rsQuery.getFloat("LuciphorFLR");

            } catch (Exception ex) {
                Logger.getRootLogger().error(ExceptionUtils.getStackTrace(ex));
            }
            psm.Charge = rsQuery.getInt("Charge");
            psm.RetentionTime = rsQuery.getFloat("RT");
            psm.ObserPrecursorMass = rsQuery.getFloat("ObservedMass");

            psm.NeighborMaxRetentionTime = rsQuery.getFloat("AdjustedRT");
            psm.Rank = rsQuery.getInt("Rank");
            psm.ScanNo = rsQuery.getInt("ScanNo");
            psm.PreAA = rsQuery.getString("PreAA");
            psm.NextAA = rsQuery.getString("NextAA");
            psm.MissedCleavage = rsQuery.getInt("MissedCleavage");
            psm.expect = rsQuery.getFloat("ExpectValue");
            psm.Probability = rsQuery.getFloat("Prob");
            psm.MassError = rsQuery.getFloat("MassError");
            psm.RawDataName = rsQuery.getString("Rawname");
            String ModificationString = rsQuery.getString("Modification");
            psm.Modifications = PTMManager.TranslateModificationString(ModificationString);

            PepIonID pepIonID = PepIonList.get(psm.GetPepKey());

            if (!pepIonID.GetPSMList().contains(psm)) {
                pepIonID.AddPSM(psm);
            }
            PSMList.put(psm.SpecNumber, psm);
        }
    }

    public void ReadFromDBPepFragments(Connection connection) throws SQLException {
        Logger.getRootLogger().info("Loading ID result from MySQL DB :" + FilenameUtils.getBaseName(mzXMLFileName) + "_PepFragments");
        Statement state = connection.createStatement();
        ResultSet rsQuery = null;
        try {
            rsQuery = state.executeQuery("SELECT * FROM " + FilenameUtils.getBaseName(mzXMLFileName) + "_PepFragments");
        } catch (Exception ex) {
            Logger.getRootLogger().error("Querying table " + FilenameUtils.getBaseName(mzXMLFileName) + "_PepFragments failed.");
            Logger.getRootLogger().error(ExceptionUtils.getStackTrace(ex));
            return;
        }
        while (rsQuery.next()) {
            FragmentPeak fragment = new FragmentPeak();
            int pepindex = rsQuery.getInt("PepIndex");
            fragment.IonType = rsQuery.getString("IonType");
            try {
                fragment.FragMZ = rsQuery.getFloat("FragMZ");
            } catch (Exception e) {
                try {
                    fragment.FragMZ = rsQuery.getFloat("MZ");
                } catch (Exception e2) {
                    try {
                        fragment.FragMZ = rsQuery.getFloat("mz");
                    } catch (Exception ex) {
                        Logger.getRootLogger().error(ExceptionUtils.getStackTrace(ex));
                    }
                }
            }
            try {
                fragment.ObservedMZ = rsQuery.getFloat("ObservedMZ");
            } catch (Exception e) {
                fragment.ObservedMZ = fragment.FragMZ;
            }
            try {
                fragment.ppm = rsQuery.getFloat("PPM");
                fragment.ApexDelta = rsQuery.getFloat("ApexDelta");
                fragment.RTOverlapP = rsQuery.getFloat("RTOverlapP");
            } catch (Exception e) {
            }
            fragment.Charge = rsQuery.getInt("Charge");
            fragment.intensity = rsQuery.getFloat("Intensity");
            fragment.corr = rsQuery.getFloat("Correlation");

            PepIonID pepIonID = PepIonIndexList.get(pepindex);
            pepIonID.FragmentPeaks.add(fragment);
        }

    }

    public void ReadFromDBMappedPepFragments(Connection connection) throws SQLException {

        Statement state = connection.createStatement();

        try {
            ResultSet rsQuery = state.executeQuery("SELECT * FROM " + FilenameUtils.getBaseName(mzXMLFileName) + "_MappedPepFragments");
            while (rsQuery.next()) {
                FragmentPeak fragment = new FragmentPeak();
                int pepindex = rsQuery.getInt("PepIndex");
                fragment.IonType = rsQuery.getString("IonType");
                fragment.ObservedMZ = rsQuery.getFloat("ObservedMZ");
                fragment.FragMZ = rsQuery.getFloat("FragMZ");
                fragment.Charge = rsQuery.getInt("Charge");
                fragment.intensity = rsQuery.getFloat("Intensity");
                fragment.corr = rsQuery.getFloat("Correlation");
                fragment.ppm = rsQuery.getFloat("PPM");
                fragment.ApexDelta = rsQuery.getFloat("ApexDelta");
                fragment.RTOverlapP = rsQuery.getFloat("RTOverlapP");

                if (MappedPepIonIndexList.containsKey(pepindex)) {
                    PepIonID pepIonID = MappedPepIonIndexList.get(pepindex);
                    pepIonID.FragmentPeaks.add(fragment);
                }
            }

        } catch (Exception e) {
            Logger.getRootLogger().error("Loading " + FilenameUtils.getBaseName(mzXMLFileName) + "_MappedPepFragments failed");
        }
    }

    public void AssignProtForPepIon() {
        if (PeptideList != null) {
            for (String pepseq : PeptideList.keySet()) {                 
                for (ProtID protein : ProteinList.values()) {
                    if (protein.Sequence.contains(pepseq) || (protein.ProtPepSeq != null && protein.ProtPepSeq.contains(pepseq))) {
                        for (PepIonID pepIonID : PeptideList.get(pepseq).values()) {
                            protein.PeptideID.put(pepIonID.GetKey(), pepIonID);
                            if (!pepIonID.ParentProtID_ProtXML.contains(protein)) {
                                pepIonID.ParentProtID_ProtXML.add(protein);
                            }
                        }
                    }
                }
            }
        }
    }

    public void AssignProtForMappedIon() {
        if (MappedPeptideList != null) {
            for (String pepseq : MappedPeptideList.keySet()) {
                for (ProtID protein : ProteinList.values()) {
                    if (protein.Sequence.contains(pepseq) || (protein.ProtPepSeq != null && protein.ProtPepSeq.contains(pepseq))) {
                        for (PepIonID pepIonID : MappedPeptideList.get(pepseq).values()) {
                            protein.PeptideID.put(pepIonID.GetKey(), pepIonID);
                            if (!pepIonID.ParentProtID_ProtXML.contains(protein)) {
                                pepIonID.ParentProtID_ProtXML.add(protein);
                            }
                        }
                    }
                }
            }
        }
    }

    public void ReadFromDBProt(Connection connection) throws SQLException, IOException, XmlPullParserException, ClassNotFoundException, InterruptedException {
        Logger.getRootLogger().info("Loading ID result from MySQL DB :" + FilenameUtils.getBaseName(mzXMLFileName) + "_ProtIDs");
        Statement state = connection.createStatement();
        ResultSet rsQuery = state.executeQuery("SELECT * FROM " + FilenameUtils.getBaseName(mzXMLFileName) + "_ProtIDs");
        while (rsQuery.next()) {
            ProtID protein = new ProtID();
            protein.setAccNo(rsQuery.getString("AccNo"));
            protein.UniProtID = rsQuery.getString("UniProtID");
            protein.ProteinLength = rsQuery.getInt("ProteinLength");
            protein.ProteinGroup = rsQuery.getString("ProteinGroup");
            protein.Description = rsQuery.getString("Description");
            protein.Mass = rsQuery.getFloat("Mass");
            protein.Probability = rsQuery.getFloat("Score");
            //String pepString = rsQuery.getString("Peptides");
//                for (int i = 0; i < pepString.split(";").length; i++) {
//                    String key = pepString.split(";")[i];
//                    if (PepIonList.containsKey(key)) {
//                        protein.AddPeptideID(PepIonList.get(key));
//                        PepIonList.get(key).ParentProtID_ProtXML.add(protein);
//                    }
//                }
            try {
                String IndisProt = rsQuery.getString("IndisProt");
                for (int i = 0; i < IndisProt.split(";").length; i++) {
                    String key = IndisProt.split(";")[i];
                    protein.AddDisProteins(key);
                }
            } catch (Exception ex) {
                Logger.getRootLogger().error(ExceptionUtils.getStackTrace(ex));
            }
            protein.SetSequence(rsQuery.getString("Sequence"));
            AddProtID(protein);
        }
    }

    public void ExportMappedPepID() throws SQLException, IOException {
        ExportMappedPepIonCSV();
//        if (connectionManager != null) {
//            Connection connection = connectionManager.GetConnection();
//            CreateMappedPepIonTable(connection);
//            ExportTableDBBulkLoader("MappedPepIonIDs", connection, false);
//            connectionManager.CloseConnection();
//        }
    }

     public void ExportPepID() throws SQLException, IOException {
         ExportPepID(null);
     }
    public void ExportPepID(String folder) throws IOException {
        ExportPepIonCSV(folder);
        ExportPepPSMCSV(folder);
//        if (connectionManager != null) {
//            Connection connection = connectionManager.GetConnection();
//            CreatePepIonTable(connection);
//            ExportTableDBBulkLoader("PepIonIDs", connection, false);
//            CreatePepPSMTable(connection);
//            ExportTableDBBulkLoader("PSMs", connection, false);
//            connectionManager.CloseConnection();
//        }
    }

    public void ExportProtID() throws SQLException, IOException {
        ExportProtID(null);
    }
    public void ExportProtID(String folder) throws SQLException, IOException {
        ExportProtIDCSV(folder);
//        if (connectionManager != null) {
//            Connection connection = connectionManager.GetConnection();
//            CreateProtIDTable(connection);
//            ExportTableDBBulkLoader("ProtIDs", connection, false);
//            connectionManager.CloseConnection();
//        }
    }

    public void ExportPepFragmentPeak() throws SQLException, IOException {
        ExportPepFragmentCSV();
//        if (connectionManager != null) {
//            Connection connection = connectionManager.GetConnection();
//            CreatePepFragmentTable(connection);
//            ExportTableDBBulkLoader("PepFragments", connection, false);
//            connectionManager.CloseConnection();
//        }
    }

    public void ExportMappedPepFragmentPeak() throws SQLException, IOException {
        ExportMappedPepFragmentCSV();
//        if (connectionManager != null) {
//            Connection connection = connectionManager.GetConnection();
//            CreateMappedPepFragmentTable(connection);
//            ExportTableDBBulkLoader("MappedPepFragments", connection, false);
//            connectionManager.CloseConnection();
//        }
    }

    private void ExportPepPSMCSV(String folder) throws IOException {
        if(folder==null | "".equals(folder)){
            folder=FilenameUtils.getFullPath(mzXMLFileName);
        }        
        Logger.getRootLogger().info("Writing PSM result to file:" + folder + FilenameUtils.getBaseName(mzXMLFileName) + "_PSMs.csv...");
        FileWriter writer = new FileWriter(folder + FilenameUtils.getBaseName(mzXMLFileName) + "_PSMs.csv");
        writer.write("SpecID,Sequence,ModSeq,TPPModSeq,LuciphorScore,LuciphorLFLR,LuciphorFLR,Modification,Charge,mz,NeutralPepMass,ObservedMass,RT,AdjustedRT,Rank,ScanNo,PreAA,NextAA,MissedCleavage,ExpectValue,MassError,Prob,Rawname,ParentPepIndex,MS1Quant\n");
        for (PepIonID pepion : PepIonList.values()) {
            for (PSM psm : pepion.GetPSMList()) {
                writer.write(psm.SpecNumber + "," + psm.Sequence + "," + psm.ModSeq + "," + psm.TPPModSeq + "," + psm.LuciphorScore + "," + psm.LuciphorLFLR + "," + psm.LuciphorFLR + "," + psm.GetModificationString() + "," + psm.Charge + "," + psm.ObserPrecursorMz() + "," + psm.NeutralPepMass + "," + psm.ObserPrecursorMass + "," + psm.RetentionTime + "," + psm.NeighborMaxRetentionTime + "," + psm.Rank + "," + psm.ScanNo + "," + psm.PreAA + "," + psm.NextAA + "," + psm.MissedCleavage + "," + psm.expect + "," + psm.MassError + "," + psm.Probability + "," + psm.RawDataName + "," + pepion.Index + "," + pepion.GetMS1() + "\n");
            }
        }
        writer.close();
    }

    private void ExportMappedPepIonCSV() throws IOException {
        Logger.getRootLogger().info("Writing MappedPepIonIDs result to file:" + FilenameUtils.getFullPath(mzXMLFileName) + FilenameUtils.getBaseName(mzXMLFileName) + "_MappedPepIonIDs.csv...");
        FileWriter writer = new FileWriter(FilenameUtils.getFullPath(mzXMLFileName) + FilenameUtils.getBaseName(mzXMLFileName) + "_MappedPepIonIDs.csv");
        writer.write("PepIndex,Sequence,ModSeq,TPPModSeq,ModInfo,Charge,mz,PredictRT, PeakRT,MS1ClusIndex,MS2ClusIndex,PeakScore,PeakHeight1,PeakHeight2,PeakHeight3,PeakArea1,PeakArea2,PeakArea3,MS1AlignmentProb,MS1AlignmentLProb,MS2AlignmentProb,MS2AlignmentLProb\n");
        int index = 1;
        for (PepIonID pepion : MappedPepIonList.values()) {
            writer.write((index++) + "," + pepion.Sequence + "," + pepion.ModSequence + "," + pepion.TPPModSeq + "," + pepion.GetModificationString() + "," + pepion.Charge + "," + pepion.NeutralPrecursorMz() + "," + pepion.PredictRTString() + "," + pepion.PeakRT + "," + pepion.GetMS1ClusIndex() + "," + pepion.GetMS2ClusIndex() + "," + pepion.PeakClusterScore + "," + pepion.PeakHeight[0] + "," + pepion.PeakHeight[1] + "," + pepion.PeakHeight[2] + "," + pepion.PeakArea[0] + "," + pepion.PeakArea[1] + "," + pepion.PeakArea[2] + "," + pepion.MS1AlignmentProbability + "," + pepion.MS1AlignmentLocalProbability + "," + pepion.MS2AlignmentProbability + "," + pepion.MS2AlignmentLocalProbability + "\n");
        }
        writer.close();
    }

    private void ExportPepIonCSV(String folder) throws IOException {
        if(folder==null | "".equals(folder)){
            folder=FilenameUtils.getFullPath(mzXMLFileName);
        }
        
        Logger.getRootLogger().info("Writing PepIon result to file:" + folder + FilenameUtils.getBaseName(mzXMLFileName) + "_PepIonIDs.csv...");
        FileWriter writer = new FileWriter(folder + FilenameUtils.getBaseName(mzXMLFileName) + "_PepIonIDs.csv");
        writer.write("PepIndex,Sequence,ModSeq,TPPModSeq,LuciphorScore,LuciphorLFLR,LuciphorFLR, IsNonDegenerate,Charge,mz,IDRT,PeakRT,NoPSMs,MS1ClusIndex,MS2ClusIndex,PeakScore,PeakHeight1,PeakHeight2,PeakHeight3,PeakArea1,PeakArea2,PeakArea3\n");
        for (PepIonID pepion : PepIonList.values()) {
            writer.write(pepion.Index + "," + pepion.Sequence + "," + pepion.ModSequence + "," + pepion.TPPModSeq + "," + pepion.GetMaxLuciphorScore() + "," + pepion.GetMinLuciphorLFLR() + "," + pepion.GetMinLuciphorFLR() + "," + (pepion.Is_NonDegenerate ? 1 : 0) + "," + pepion.Charge + "," + pepion.NeutralPrecursorMz() + "," + pepion.GetIDRT() + "," + pepion.PeakRT + "," + pepion.GetSpectralCount() + "," + pepion.GetMS1ClusIndex() + "," + pepion.GetMS2ClusIndex() + "," + pepion.PeakClusterScore + "," + pepion.PeakHeight[0] + "," + pepion.PeakHeight[1] + "," + pepion.PeakHeight[2] + "," + pepion.PeakArea[0] + "," + pepion.PeakArea[1] + "," + pepion.PeakArea[2] + "\n");
        }
        writer.close();
    }

    private void ExportProtIDCSV(String folder) throws IOException {
        if(folder==null | "".equals(folder)){
            folder=FilenameUtils.getFullPath(mzXMLFileName);
        }
        Logger.getRootLogger().info("Writing ProteinID result to file:" + folder + FilenameUtils.getBaseName(mzXMLFileName) + "_ProtIDs.csv...");
        FileWriter writer = new FileWriter(folder + FilenameUtils.getBaseName(mzXMLFileName) + "_ProtIDs.csv");
        writer.write("AccNo,UniProtID,ProteinLength,ProteinGroup,IndisProt,Description,Mass,Score,Peptides,Sequence\n");
        for (ProtID protein : ProteinList.values()) {
            String pepstring = "";
            for (PepIonID pep : protein.PeptideID.values()) {
                pepstring += pep.GetKey() + ";";
            }
            String IndisProt = "";
            for (String indisprot : protein.IndisProteins) {
                IndisProt += indisprot + ";";
            }
            writer.write(protein.getAccNo() + "," + protein.UniProtID + "," + protein.ProteinLength + "," + protein.ProteinGroup + "," + IndisProt + "," + protein.Description + "," + protein.Mass + "," + protein.Probability + "," + pepstring + "," + protein.Sequence + "\n");
        }
        writer.close();
    }

    private void ExportMappedPepFragmentCSV() throws IOException {
        Logger.getRootLogger().info("Writing PepFragment result to file:" + FilenameUtils.getFullPath(mzXMLFileName) + FilenameUtils.getBaseName(mzXMLFileName) + "_MappedPepFragments.csv...");
        FileWriter writer = new FileWriter(FilenameUtils.getFullPath(mzXMLFileName) + FilenameUtils.getBaseName(mzXMLFileName) + "_MappedPepFragments.csv");
        writer.write("PepIndex,IonType,fragMZ,ObservedMZ,Charge,Intensity,Correlation,PPM,ApexDelta,RTOverlapP\n");
        for (PepIonID pepion : MappedPepIonList.values()) {
            for (FragmentPeak frag : pepion.FragmentPeaks) {
                writer.write(pepion.Index + "," + frag.IonType + "," + frag.FragMZ + "," + frag.ObservedMZ + "," + frag.Charge + "," + frag.intensity + "," + frag.corr + "," + frag.ppm + "," + frag.ApexDelta + "," + frag.RTOverlapP + "\n");
            }
        }
        writer.close();
    }

    private void ExportPepFragmentCSV() throws IOException {
        Logger.getRootLogger().info("Writing PepFragment result to file:" + FilenameUtils.getFullPath(mzXMLFileName) + FilenameUtils.getBaseName(mzXMLFileName) + "_PepFragments.csv...");
        FileWriter writer = new FileWriter(FilenameUtils.getFullPath(mzXMLFileName) + FilenameUtils.getBaseName(mzXMLFileName) + "_PepFragments.csv");
        writer.write("PepIndex,IonType,fragMZ,ObservedMZ,Charge,Intensity,Correlation,PPM,ApexDelta,RTOverlapP\n");
        for (PepIonID pepion : PepIonList.values()) {
            for (FragmentPeak frag : pepion.FragmentPeaks) {
                writer.write(pepion.Index + "," + frag.IonType + "," + frag.FragMZ + "," + frag.ObservedMZ + "," + frag.Charge + "," + frag.intensity + "," + frag.corr + "," + frag.ppm + "," + frag.ApexDelta + "," + frag.RTOverlapP + "\n");
            }
        }
        writer.close();
    }

    protected void CreatePepFragmentTable(Connection connection) throws SQLException {
        connection.createStatement().execute("DROP TABLE IF EXISTS " + FilenameUtils.getBaseName(mzXMLFileName) + "_PepFragments;");
        connection.createStatement().execute("CREATE TABLE " + FilenameUtils.getBaseName(mzXMLFileName) + "_PepFragments (PepIndex INT NOT NULL, IonType VARCHAR(4) NOT NULL, fragMZ DOUBLE NOT NULL,ObservedMZ DOUBLE NOT NULL, Charge INT NOT NULL, Intensity DOUBLE NOT NULL, Correlation DOUBLE NOT NULL,PPM DOUBLE NOT NULL, ApexDelta DOUBLE NOT NULL,RTOverlapP DOUBLE NOT NULL, PRIMARY KEY (PepIndex,IonType,Charge));");
    }

    protected void CreateMappedPepFragmentTable(Connection connection) throws SQLException {
        connection.createStatement().execute("DROP TABLE IF EXISTS " + FilenameUtils.getBaseName(mzXMLFileName) + "_MappedPepFragments;");
        connection.createStatement().execute("CREATE TABLE " + FilenameUtils.getBaseName(mzXMLFileName) + "_MappedPepFragments (PepIndex INT NOT NULL, IonType VARCHAR(4) NOT NULL, fragMZ DOUBLE NOT NULL,ObservedMZ DOUBLE NOT NULL, Charge INT NOT NULL, Intensity DOUBLE NOT NULL, Correlation DOUBLE NOT NULL,PPM DOUBLE NOT NULL,ApexDelta DOUBLE NOT NULL,RTOverlapP DOUBLE NOT NULL, PRIMARY KEY (PepIndex,IonType,Charge));");
    }

    protected void CreateMappedPepIonTable(Connection connection) throws SQLException {
        connection.createStatement().execute("DROP TABLE IF EXISTS " + FilenameUtils.getBaseName(mzXMLFileName) + "_MappedPepIonIDs;");
        connection.createStatement().execute("CREATE TABLE " + FilenameUtils.getBaseName(mzXMLFileName) + "_MappedPepIonIDs (PepIndex INT NOT NULL, Sequence VARCHAR(100) NOT NULL, ModSeq VARCHAR(200) NOT NULL, TPPModSeq VARCHAR(200) NOT NULL, ModInfo TEXT, Charge INT NOT NULL, mz DOUBLE NOT NULL, PredictRT VARCHAR(200), PeakRT DOUBLE NOT NULL, MS1ClusIndex VARCHAR(100),MS2ClusIndex VARCHAR(100),PeakScore DOUBLE,PeakHeight1 DOUBLE, PeakHeight2 DOUBLE, PeakHeight3 DOUBLE, PeakArea1 DOUBLE, PeakArea2 DOUBLE, PeakArea3 DOUBLE, MS1AlignmentProb DOUBLE,MS1AlignmentLProb DOUBLE,MS2AlignmentProb DOUBLE,MS2AlignmentLProb DOUBLE,PRIMARY KEY (ModSeq, Charge));");
    }

    protected void CreatePepIonTable(Connection connection) throws SQLException {
        connection.createStatement().execute("DROP TABLE IF EXISTS " + FilenameUtils.getBaseName(mzXMLFileName) + "_PepIonIDs;");
        connection.createStatement().execute("CREATE TABLE " + FilenameUtils.getBaseName(mzXMLFileName) + "_PepIonIDs (PepIndex INT NOT NULL, Sequence VARCHAR(100) NOT NULL, ModSeq VARCHAR(200) NOT NULL,TPPModSeq VARCHAR(200) NOT NULL, LuciphorScore DOUBLE NOT NULL,LuciphorLFLR DOUBLE NOT NULL,LuciphorFLR DOUBLE NOT NULL,  IsNonDegenerate TINYINT(1) NOT NULL,Charge INT NOT NULL, mz DOUBLE NOT NULL, IDRT DOUBLE NOT NULL,PeakRT DOUBLE NOT NULL, NoPSMs INT NOT NULL, MS1ClusIndex VARCHAR(100),MS2ClusIndex VARCHAR(100), PeakScore DOUBLE, PeakHeight1 DOUBLE, PeakHeight2 DOUBLE, PeakHeight3 DOUBLE, PeakArea1 DOUBLE, PeakArea2 DOUBLE, PeakArea3 DOUBLE, PRIMARY KEY (ModSeq, Charge));");
    }

    protected void CreateProtIDTable(Connection connection) throws SQLException {
        connection.createStatement().execute("DROP TABLE IF EXISTS " + FilenameUtils.getBaseName(mzXMLFileName) + "_ProtIDs;");
        connection.createStatement().execute("CREATE TABLE " + FilenameUtils.getBaseName(mzXMLFileName) + "_ProtIDs (AccNo VARCHAR(50) NOT NULL,UniProtID VARCHAR(20),ProteinLength INT NOT NULL,ProteinGroup TEXT,IndisProt TEXT, Description TEXT,Mass DOUBLE NOT NULL,Score DOUBLE NOT NULL,Peptides TEXT NOT NULL,Sequence TEXT, PRIMARY KEY (AccNo));");
    }

    protected void CreatePepPSMTable(Connection connection) throws SQLException {
        connection.createStatement().execute("DROP TABLE IF EXISTS " + FilenameUtils.getBaseName(mzXMLFileName) + "_PSMs;");
        connection.createStatement().execute("CREATE TABLE " + FilenameUtils.getBaseName(mzXMLFileName) + "_PSMs (SpecID VARCHAR(100) NOT NULL, Sequence VARCHAR(100) NOT NULL, ModSeq VARCHAR(200) NOT NULL, TPPModSeq VARCHAR(200) NOT NULL, LuciphorScore DOUBLE NOT NULL,LuciphorLFLR DOUBLE NOT NULL,LuciphorFLR DOUBLE NOT NULL,Modification  TEXT, Charge INT NOT NULL, mz DOUBLE NOT NULL, NeutralPepMass DOUBLE NOT NULL, ObservedMass DOUBLE NOT NULL, RT DOUBLE NOT NULL, AdjustedRT DOUBLE, Rank INT NOT NULL, ScanNo INT NOT NULL, PreAA VARCHAR(2), NextAA VARCHAR(2), MissedCleavage INT NOT NULL, ExpectValue DOUBLE NOT NULL, MassError DOUBLE NOT NULL, Prob DOUBLE NOT NULL, Rawname VARCHAR(50) NOT NULL, ParentPepIndex INT NOT NULL, PRIMARY KEY (SpecID, Charge));");
    }

    protected void ExportTableDBBulkLoader(String TableSuffix, Connection connection, boolean filedelete) throws SQLException {

        Logger.getRootLogger().info("Writing " + TableSuffix + " result to MySQL DB:" + FilenameUtils.getBaseName(mzXMLFileName) + "_" + TableSuffix + "...");

        Statement state = connection.createStatement();
        state.executeUpdate("LOAD DATA LOCAL INFILE '" + FilenameUtils.separatorsToUnix(FilenameUtils.getFullPath(mzXMLFileName)) + FilenameUtils.getBaseName(mzXMLFileName) + "_" + TableSuffix + ".csv'" + " INTO TABLE " + FilenameUtils.getBaseName(mzXMLFileName) + "_" + TableSuffix + " FIELDS TERMINATED BY ','" + " LINES TERMINATED BY '\\n' IGNORE 1 LINES");
        state.close();
        File file = new File(FilenameUtils.separatorsToUnix(FilenameUtils.getFullPath(mzXMLFileName)) + FilenameUtils.getBaseName(mzXMLFileName) + "_" + TableSuffix + ".csv");
        if (filedelete) {
            file.delete();
        }
    }

    public PepIonID GetPepID(PSM psm) {
        return PepIonList.get(psm.GetPepKey());
    }

    public void AddPSM(PSM psm) {
        PSMList.put(psm.SpecNumber, psm);
        PepIonID pepIonID = GetPepID(psm);
        if (pepIonID == null) {
            pepIonID = new PepIonID();
            pepIonID.SetInfobyPSM(psm);
            AddPeptideID(pepIonID);
        }
        pepIonID.AddPSM(psm);
    }

    public void LoadSequence() throws IOException, XmlPullParserException {
        for (ProtID protID : ProteinList.values()) {
            if (protID.Sequence == null || "".equals(protID.Sequence) || "null".equals(protID.Sequence)) {
                if (!protID.IsDecoy(DecoyTag)) {
                    if (FastaPath != null && !"".equals(FastaPath)) {
                        //String Sequence = GetSequenceFactory().getProtein(protID.UniProtID).getSequence();                    
                        //String Sequence = GetFastaParser().ProteinList.get(protID.getAccNo()).Seq;
                        try {
                            String Sequence = GetFastaParser().ProteinList.get(protID.getAccNo())[0];
                            if (Sequence != null) {
                                protID.SetSequence(Sequence);
                            } else {
                                Logger.getRootLogger().error("Can't find sequence in fasta file for protein:" + protID.getAccNo());
                            }
                        } catch (Exception e) {
                            Logger.getRootLogger().error("Can't find sequence in fasta file for protein:" + protID.getAccNo());
                        }
                    }
                }
            }
        }
    }

    public void AddProtID(ProtID protID) throws ClassNotFoundException, InterruptedException, IOException, XmlPullParserException {
        if (!ProteinList.containsKey(protID.getAccNo())) {
            ProteinList.put(protID.getAccNo(), protID);
        }
    }

    public ArrayList<PSM> FindPsmsBymzRT(float mz, float RT, int charge, float PPM, float RTtol) {
        ArrayList<PSM> psms = new ArrayList<>();
        for (PSM psm : PSMList.values()) {
            if (psm.Charge == charge && InstrumentParameter.CalcPPM(mz, psm.ObserPrecursorMz()) < PPM && Math.abs(psm.RetentionTime - RT) < RTtol) {
                psms.add(psm);
            }
        }
        return psms;
    }

    private void FindSpecExpectThresholdByFDR() {

        SortedPSMListEvalue sortedlist = new SortedPSMListEvalue();
        sortedlist.addAll(PSMList.values());

        int positive = 0;
        int negative = 0;
        for (int i = 0; i < sortedlist.size(); i++) {
            PSM psm = sortedlist.get(i);
            if (psm.IsDecoy(DecoyTag)) {
                negative++;
            } else {
                positive++;
            }
            if ((float) negative / (float) (positive + negative) >= FDR) {
                ExpectThreshold = psm.expect;
                return;
            }
        }
    }

    public void ROCProtByMaxIniProb(String decoytag) throws IOException {
        if (ProteinList.isEmpty()) {
            return;
        }
        SortedProteinListMaxIniProb sortedlist = new SortedProteinListMaxIniProb();
        sortedlist.addAll(ProteinList.values());

        FileWriter writer = new FileWriter(FilenameUtils.getFullPath(mzXMLFileName) + "/" + FilenameUtils.getBaseName(mzXMLFileName) + "_ROCProt.txt");
        int positive = 0;
        int negative = 0;
        ProtID protein = sortedlist.get(0);
        if (protein.IsDecoy(decoytag)) {
            negative++;
        } else {
            positive++;
        }
        for (int i = 1; i < sortedlist.size(); i++) {
            protein = sortedlist.get(i);
            if (protein.IsDecoy(decoytag)) {
                negative++;
                //System.out.println(protein.getAccNo());
            } else {
                positive++;
            }
            if (protein.MaxIniProb < sortedlist.get(i - 1).MaxIniProb) {
                writer.write(positive + "\t" + negative + "\t" + protein.MaxIniProb + "\t" + ((float) negative / (float) (positive)) + "\n");
                if ((float) negative / (float) (positive) >= 0.5f) {
                    break;
                }
            }
        }
        writer.close();
    }

    public void ROCPepByMaxIniProb(String decoytag) throws IOException {

        if (PepIonList.isEmpty()) {
            return;
        }
        SortedPepListProb sortedlist = new SortedPepListProb();
        sortedlist.addAll(PepIonList.values());

        int positive = 0;
        int negative = 0;

        FileWriter writer = new FileWriter(FilenameUtils.getFullPath(mzXMLFileName) + FilenameUtils.getBaseName(mzXMLFileName) + "_ROCPep.txt");
        PepIonID pep = sortedlist.get(0);
        if (pep.IsDecoy(DecoyTag)) {
            negative++;
        } else {
            positive++;
        }
        for (int i = 1; i < sortedlist.size(); i++) {
            pep = sortedlist.get(i);
            if (pep.MaxProbability < 0.1f) {
                break;
            }
            if (pep.IsDecoy(DecoyTag)) {
                negative++;
            } else {
                positive++;
            }

            if (pep.MaxProbability < sortedlist.get(i - 1).MaxProbability) {
                writer.write(positive + "\t" + negative + "\t" + pep.MaxProbability + "\t" + ((float) negative / (float) (positive)) + "\n");
                if ((float) negative / (float) (positive) >= 0.5f) {
                    break;
                }
            }
        }
        writer.close();
    }

    private void FindMaxIniProbThresholdByFDR() {
        if (ProteinList.isEmpty()) {
            return;
        }
        SortedProteinListMaxIniProb sortedlist = new SortedProteinListMaxIniProb();
        sortedlist.addAll(ProteinList.values());

        //writer = new FileWriter(FilenameUtils.getFullPath(mzXMLFileName)+"/" + FilenameUtils.getBaseName(mzXMLFileName)+"_Pro.txt");
        int positive = 0;
        int negative = 0;
        ProtID protein = sortedlist.get(0);
        if (protein.IsDecoy(DecoyTag)) {
            negative++;
        } else {
            positive++;
        }
        for (int i = 1; i < sortedlist.size(); i++) {
            protein = sortedlist.get(i);
            if (protein.IsDecoy(DecoyTag)) {
                negative++;
                //System.out.println(protein.getAccNo()+"-"+protein.ProteinGroup+"-Decoy");
            } else {
                positive++;
                //System.out.println(protein.getAccNo()+"-"+ protein.ProteinGroup);
            }
            if (protein.MaxIniProb < sortedlist.get(i - 1).MaxIniProb && (float) negative / (float) (positive) >= ProteinFDR) {
                ProteinProbThreshold = protein.MaxIniProb;
                Logger.getRootLogger().info("Protein maxiniprob threshold=" + ProteinProbThreshold + " Estimated raw protein FDR:" + (float) negative / (float) (positive) + "(Target/Decoy)=(" + positive + "/" + negative + ")");
                return;
            }
        }
        Logger.getRootLogger().info("Protein maxiniprob threshold=" + ProteinProbThreshold + " Estimated raw protein FDR:" + (float) negative / (float) (positive) + "(Target/Decoy)=(" + positive + "/" + negative + ")");
    }

    private void FindLocalPWThresholdByFDR() {
        //FileWriter writer = null;
        //try {
        if (ProteinList.isEmpty()) {
            return;
        }
        SortedProteinListProb sortedlist = new SortedProteinListProb();
        sortedlist.addAll(ProteinList.values());

        //writer = new FileWriter(FilenameUtils.getFullPath(mzXMLFileName)+"/" + FilenameUtils.getBaseName(mzXMLFileName)+"_Pro.txt");
        int positive = 0;
        int negative = 0;
        ProtID protein = sortedlist.get(0);
        if (protein.IsDecoy(DecoyTag)) {
            negative++;
        } else {
            positive++;
        }
        for (int i = 1; i < sortedlist.size(); i++) {
            protein = sortedlist.get(i);
            if (protein.IsDecoy(DecoyTag)) {
                negative++;
                //System.out.println(protein.getAccNo()+"-"+protein.ProteinGroup+"-Decoy");
            } else {
                positive++;
                //System.out.println(protein.getAccNo()+"-"+ protein.ProteinGroup);
            }
            if (protein.Probability < sortedlist.get(i - 1).Probability && (float) negative / (float) (positive) >= ProteinFDR) {
                ProteinProbThreshold = protein.Probability;
                Logger.getRootLogger().info("Protein probability threshold=" + ProteinProbThreshold + " Estimated raw protein FDR:" + (float) negative / (float) (positive) + "(Target/Decoy)=(" + positive + "/" + negative + ")");
                return;
            }
        }
    }

    private void FindProteinMaxLocalPWThresholdByFDR() {
        if (ProteinList.isEmpty()) {
            return;
        }
        SortedProteinListMaxLocalPW sortedlist = new SortedProteinListMaxLocalPW();
        sortedlist.addAll(ProteinList.values());

        //writer = new FileWriter(FilenameUtils.getFullPath(mzXMLFileName)+"/" + FilenameUtils.getBaseName(mzXMLFileName)+"_Pro.txt");
        int positive = 0;
        int negative = 0;
        ProtID protein = sortedlist.get(0);
        if (protein.IsDecoy(DecoyTag)) {
            negative++;
        } else {
            positive++;
        }
        for (int i = 1; i < sortedlist.size(); i++) {
            protein = sortedlist.get(i);
            if (protein.IsDecoy(DecoyTag)) {
                negative++;
            } else {
                positive++;
                //System.out.println(protein.getAccNo()+"-"+ protein.ProteinGroup);
            }
            if (protein.MaxLocalPW < sortedlist.get(i - 1).MaxLocalPW && (float) negative / (float) (positive) >= ProteinFDR) {
                ProteinProbThreshold = protein.MaxLocalPW;
                Logger.getRootLogger().info("Protein probability threshold=" + ProteinProbThreshold + " Estimated raw protein FDR:" + (float) negative / (float) (positive) + "(Target/Decoy)=(" + positive + "/" + negative + ")");
                return;
            }
        }
    }

    private void FindSpecProbThresholdByFDR() {
        if (PSMList.isEmpty()) {
            return;
        }
        SortedPSMListProb sortedlist = new SortedPSMListProb();
        sortedlist.addAll(PSMList.values());

        //writer = new FileWriter(FilenameUtils.getFullPath(mzXMLFileName)+"/" + FilenameUtils.getBaseName(mzXMLFileName)+"_Pep.txt");
        int positive = 0;
        int negative = 0;
        PSM psm = sortedlist.get(0);
        if (psm.IsDecoy(DecoyTag)) {
            negative++;
        } else {
            positive++;
        }
        for (int i = 1; i < sortedlist.size(); i++) {
            psm = sortedlist.get(i);
            if (psm.Probability < 0.1f) {
                break;
            }
            if (psm.IsDecoy(DecoyTag)) {
                negative++;
            } else {
                positive++;
            }
            if (psm.Probability < sortedlist.get(i - 1).Probability && ((float) negative / (float) (positive) >= FDR)) {
                SpecProbThreshold = psm.Probability;
                Logger.getRootLogger().info("Spectrum probability threshold=" + SpecProbThreshold + " Estimated Spectrum FDR:" + (float) negative / (float) (positive) + "(Target/Decoy)=(" + positive + "/" + negative + ")");
                return;
            }
        }
    }

    public void ClearPeakData() {
        for (PepIonID pepIonID : PepIonList.values()) {
            for (PeakCluster peak : pepIonID.MS1PeakClusters) {
                peak.IsoPeaksCurves = null;
                peak.MonoIsotopePeak = null;
            }
            for (PeakCluster peak : pepIonID.MS2UnfragPeakClusters) {
                peak.IsoPeaksCurves = null;
                peak.MonoIsotopePeak = null;
            }
        }
    }

    public void ClearAssignPeakCluster() {
        for (PepIonID pepIonID : PepIonList.values()) {
            pepIonID.MS1PeakClusters.clear();
            pepIonID.MS2UnfragPeakClusters.clear();
        }
    }

    public void ClearAssignPeakClusterMappedion() {
        for (PepIonID pepIonID : MappedPepIonList.values()) {
            pepIonID.MS1PeakClusters.clear();
            pepIonID.MS2UnfragPeakClusters.clear();
        }
    }

    public void ReduceMemoryUsage() {
        fastaParser = null;
        for (PepIonID pepIonID : PepIonList.values()) {
            pepIonID.ClearPepFragFactory();
            for (PeakCluster peakCluster : pepIonID.MS1PeakClusters) {
                peakCluster.GroupedFragmentPeaks.clear();
                peakCluster.MonoIsotopePeak = null;
            }
            for (PeakCluster peakCluster : pepIonID.MS2UnfragPeakClusters) {
                peakCluster.GroupedFragmentPeaks.clear();
                peakCluster.MonoIsotopePeak = null;
            }
        }
        for (PepIonID pepIonID : MappedPepIonList.values()) {
            pepIonID.ClearPepFragFactory();
            for (PeakCluster peakCluster : pepIonID.MS1PeakClusters) {
                peakCluster.GroupedFragmentPeaks.clear();
                peakCluster.MonoIsotopePeak = null;
            }
            for (PeakCluster peakCluster : pepIonID.MS2UnfragPeakClusters) {
                peakCluster.GroupedFragmentPeaks.clear();
                peakCluster.MonoIsotopePeak = null;
            }
        }
    }

    public void ClearPSMs() {
        PSMList.clear();
        for (PepIonID pepIonID : PepIonList.values()) {
            pepIonID.CalcRT();
            pepIonID.GetSpectralCount();
            pepIonID.GetPSMList().clear();
        }
    }

    public void ReMapProPep() {
        ClearProPeplist();
        GeneratePepSeqList();
        AssignProtForPepIon();
        GenerateMappedPepSeqList();
        AssignProtForMappedIon();
        GenearteAssignIonList();
    }

    public void RemoveRedundantFrag() {
        for (PepIonID pep : GetPepIonList().values()) {
            pep.RemoveRedundantFrag();
        }
        for (PepIonID pep : GetMappedPepIonList().values()) {
            pep.RemoveRedundantFrag();
        }
    }

    public void SetFilterByGroupWeight() {
        for (PepIonID pep : GetPepIonList().values()) {
            pep.FilteringWeight = pep.GroupWeight;
        }
        for (PepIonID pep : GetMappedPepIonList().values()) {
            pep.FilteringWeight = pep.GroupWeight;
        }
    }

    public void SetFilterByWeight() {
        for (PepIonID pep : GetPepIonList().values()) {
            pep.FilteringWeight = pep.Weight;
        }
        for (PepIonID pep : GetMappedPepIonList().values()) {
            pep.FilteringWeight = pep.Weight;
        }
    }

    public void CheckRT(String step) {

        for (PSM psm : PSMList.values()) {
            if (psm.RetentionTime == -1f) {
                System.out.println("CheckRTfailed " + step);
                return;
            }
        }
        for (PepIonID pepIonID : GetPepIonList().values()) {
            for (PSM psm : pepIonID.GetPSMList()) {
                if (psm.RetentionTime == -1f) {
                    System.out.println("CheckRTfailed " + step);
                    return;
                }
            }
        }
        System.out.println("done");
    }

    public void RemoveLowFragPep(int i) {
        ArrayList<PepIonID> removelist = new ArrayList<>();
        for (PepIonID pep : PepIonList.values()) {
            if (pep.FragmentPeaks.size()<i) {
                removelist.add(pep);
            }
        }
        for (PepIonID pep : removelist) {
            PepIonList.remove(pep.GetKey());
        }
        removelist = new ArrayList<>();
        for (PepIonID pep : MappedPepIonList.values()) {
            if (pep.FragmentPeaks.size()<i) {
                removelist.add(pep);
            }
        }
        for (PepIonID pep : removelist) {
            MappedPepIonList.remove(pep.GetKey());
        }
    }

    public class SeqDecoyObj {

        public String Seq;
        public boolean IsDecoy = false;
        public float prob;

        public SeqDecoyObj(String Seq, boolean IsDecoy, float prob) {
            this.Seq = Seq;
            this.IsDecoy = IsDecoy;
            this.prob = prob;
        }
    }

    public class SortedSeqDecoyObj extends SortedList<SeqDecoyObj> {

        public SortedSeqDecoyObj() {
            super(new Comparator<SeqDecoyObj>() {
                @Override
                public int compare(SeqDecoyObj x, SeqDecoyObj y) {
                    return Float.compare(y.prob, x.prob);
                }
            });
        }

    }

    public void FindPepProbThresholdByFDRAtPepSeq() {

        SortedSeqDecoyObj PepSeqScore = new SortedSeqDecoyObj();
        for (String pepseq : PeptideList.keySet()) {
            float maxscore = 0f;
            boolean isdecoy = true;
            for (PepIonID pepion : PeptideList.get(pepseq).values()) {
                if (maxscore < pepion.MaxProbability) {
                    maxscore = pepion.MaxProbability;
                }
                if (!pepion.IsDecoy(DecoyTag)) {
                    isdecoy = false;
                }
            }
            PepSeqScore.add(new SeqDecoyObj(pepseq, isdecoy, maxscore));
        }

        int positive = 0;
        int negative = 0;
        SeqDecoyObj pep = PepSeqScore.get(0);
        if (pep.IsDecoy) {
            negative++;
        } else {
            positive++;
        }
        for (int i = 1; i < PepSeqScore.size(); i++) {
            pep = PepSeqScore.get(i);
            if (pep.IsDecoy) {
                negative++;
            } else {
                positive++;
            }
            if (pep.prob < PepSeqScore.get(i - 1).prob && ((float) negative / (float) (positive) >= FDR)) {
                PepProbThreshold = pep.prob;
                Logger.getRootLogger().info("Probability threshold=" + PepProbThreshold + " Estimated FDR:" + (float) negative / (float) (positive) + "(Target/Decoy)=(" + positive + "/" + negative + ")");
                return;
            }
        }
    }

    public void FindPepProbThresholdByFDR() {
        if (PepIonList.isEmpty()) {
            return;
        }
        SortedPepListProb sortedlist = new SortedPepListProb();
        sortedlist.addAll(PepIonList.values());

        int positive = 0;
        int negative = 0;

        PepIonID pep = sortedlist.get(0);
        if (pep.IsDecoy(DecoyTag)) {
            negative++;
        } else {
            positive++;
        }
        for (int i = 1; i < sortedlist.size(); i++) {
            pep = sortedlist.get(i);
            if (pep.MaxProbability < 0.1f) {
                break;
            }
            if (pep.IsDecoy(DecoyTag)) {
                negative++;
            } else {
                positive++;
            }
            if (pep.MaxProbability < sortedlist.get(i - 1).MaxProbability && ((float) negative / (float) (positive) >= FDR)) {
                PepProbThreshold = pep.MaxProbability;
                Logger.getRootLogger().info("Probability threshold=" + PepProbThreshold + " Estimated FDR:" + (float) negative / (float) (positive) + "(Target/Decoy)=(" + positive + "/" + negative + ")");
                return;
            }
        }
        Logger.getRootLogger().info("Probability threshold=" + PepProbThreshold + " Estimated FDR:" + (float) negative / (float) (positive) + "(Target/Decoy)=(" + positive + "/" + negative + ")");
    }

    private void RemoveHighExpPSM() {

        ArrayList<PSM> removelist = new ArrayList<>();
        for (PSM psm : PSMList.values()) {
            if (psm.expect > ExpectThreshold) {
                psm.pepIonID.GetPSMList().remove(psm);
                removelist.add(psm);
            }
        }
        for (PSM psm : removelist) {
            PSMList.remove(psm.SpecNumber);
        }
        ArrayList<PepIonID> emptypeplist = new ArrayList<>();
        for (PepIonID pep : PepIonList.values()) {
            if (pep.GetPSMList().isEmpty()) {
                emptypeplist.add(pep);
            }
        }
        for (PepIonID pep : emptypeplist) {
            PepIonList.remove(pep.GetKey());
        }
    }

    public void GenearteAssignIonList() {
        AssignedPepIonList.clear();
        for (ProtID protein : ProteinList.values()) {
            for (PepIonID pepIonID : protein.PeptideID.values()) {
                AssignedPepIonList.put(pepIonID.GetKey(), pepIonID);
            }
        }
    }

    public void RemoveLowLocalPWProtein(float LocalPW) {
        ArrayList<ProtID> removelist = new ArrayList<>();
        for (ProtID protein : ProteinList.values()) {
            if (protein.Probability < LocalPW) {
                removelist.add(protein);
            }
        }
        for (ProtID protein : removelist) {
            ProteinList.remove(protein.getAccNo());
        }
        GenearteAssignIonList();
    }

    private void RemoveLowMaxLocaPWProteinDecoy() {

        ArrayList<ProtID> removelist = new ArrayList<>();
        for (ProtID protein : ProteinList.values()) {
            if (protein.MaxLocalPW < ProteinProbThreshold || protein.IsDecoy(DecoyTag)) {
                removelist.add(protein);
            }
        }
        for (ProtID protein : removelist) {
            ProteinList.remove(protein.getAccNo());
        }
        GenearteAssignIonList();
    }

    public void RemoveLowMaxIniProbProteinDecoy() {

        ArrayList<ProtID> removelist = new ArrayList<>();
        for (ProtID protein : ProteinList.values()) {
            if (protein.MaxIniProb < ProteinProbThreshold || protein.IsDecoy(DecoyTag)) {
                removelist.add(protein);
            }
        }
        for (ProtID protein : removelist) {
            ProteinList.remove(protein.getAccNo());
        }
        GenearteAssignIonList();
    }

    public void RemoveLowMaxIniProbProtein(float maxiniprob) {

        ArrayList<ProtID> removelist = new ArrayList<>();
        for (ProtID protein : ProteinList.values()) {
            if (protein.MaxIniProb < maxiniprob) {
                removelist.add(protein);
            }
        }
        for (ProtID protein : removelist) {
            ProteinList.remove(protein.getAccNo());
        }
    }

    private void RemoveLowProbProteinDecoy() {

        ArrayList<ProtID> removelist = new ArrayList<>();
        for (ProtID protein : ProteinList.values()) {
            if (protein.Probability < ProteinProbThreshold || protein.IsDecoy(DecoyTag)) {
                removelist.add(protein);
            }
        }
        for (ProtID protein : removelist) {
            ProteinList.remove(protein.getAccNo());
        }
        GenearteAssignIonList();
    }

    public void DetermineAssignIonListByProtPepSeq() {
        for (ProtID proid : ProteinList.values()) {
            for (String seq : proid.ProtPepSeq) {
                if (PeptideList.containsKey(seq)) {
                    for (PepIonID pep : PeptideList.get(seq).values()) {
                        proid.AddPeptideID(pep);
                        if (!pep.ParentProtID_ProtXML.contains(proid)) {
                            pep.ParentProtID_ProtXML.add(proid);
                        }
                        if (!AssignedPepIonList.containsKey(pep.GetKey())) {
                            AssignedPepIonList.put(pep.GetKey(), pep);
                        }
                    }
                }
            }
        }
    }

    public void RemoveLowProbPepSeq() {

        // = new HashMap<>();
        ArrayList<PepIonID> removelist = new ArrayList<>();
        for (PepIonID pep : GetPepIonList().values()) {
            if (pep.MaxProbability < PepProbThreshold) {
                removelist.add(pep);
            }
        }
        for (PepIonID pep : removelist) {
            GetPepIonList().remove(pep.GetKey());
            //LowScorePep.put(pep.GetKey(), pep);
        }
    }

    public void RemoveLowProbPep() {

        // = new HashMap<>();
        ArrayList<PepIonID> removelist = new ArrayList<>();
        for (PepIonID pep : GetPepIonList().values()) {
            if (pep.MaxProbability < PepProbThreshold) {
                removelist.add(pep);
            }
        }
        for (PepIonID pep : removelist) {
            GetPepIonList().remove(pep.GetKey());
            //LowScorePep.put(pep.GetKey(), pep);
        }
    }

    private void RemoveLowProbPSM() {

        LowScorePSMByPepKey = new HashMap<>();
        LowScorePSM = new HashMap<>();
        ArrayList<PSM> removelist = new ArrayList<>();
        for (PSM psm : PSMList.values()) {
            if (psm.Probability < SpecProbThreshold) {
                psm.pepIonID.GetPSMList().remove(psm);
                removelist.add(psm);
            }
        }
        for (PSM psm : removelist) {
            PSMList.remove(psm.SpecNumber);
            LowScorePSM.put(psm.SpecNumber, psm);
            if (LowScorePSMByPepKey.containsKey(psm.GetPepKey())) {
                if (LowScorePSMByPepKey.get(psm.GetPepKey()).Probability < psm.Probability) {
                    LowScorePSMByPepKey.put(psm.GetPepKey(), psm);
                }
            } else {
                LowScorePSMByPepKey.put(psm.GetPepKey(), psm);
            }
        }
        ArrayList<PepIonID> emptypeplist = new ArrayList<>();
        for (PepIonID pep : PepIonList.values()) {
            if (pep.GetPSMList().isEmpty()) {
                emptypeplist.add(pep);
            }
        }
        for (PepIonID pep : emptypeplist) {
            PepIonList.remove(pep.GetKey());
        }
    }

    public void RemoveDecoyPep() {
        ArrayList<PepIonID> removelist = new ArrayList<>();
        for (PepIonID pep : PepIonList.values()) {
            if (pep.IsDecoy(DecoyTag)) {
                removelist.add(pep);
            }
        }
        for (PepIonID pep : removelist) {
            PepIonList.remove(pep.GetKey());
        }
    }

    public void RemoveDecoyProtein() {
        ArrayList<ProtID> removelist = new ArrayList<>();
        for (ProtID protein : ProteinList.values()) {
            if (protein.IsDecoy(DecoyTag)) {
                removelist.add(protein);
            }
        }
        for (ProtID protein : removelist) {
            ProteinList.remove(protein.getAccNo());
        }
    }

    public void FilterByPepDecoyFDR(String DecoyTag, float fdr) {
        this.DecoyTag = DecoyTag;
        this.FDR = fdr;
        try {
            ROCPepByMaxIniProb(DecoyTag);
        } catch (IOException ex) {
            java.util.logging.Logger.getLogger(LCMSID.class.getName()).log(Level.SEVERE, null, ex);
        }
        FindPepProbThresholdByFDR();
        RemoveLowProbPep();
        RemoveDecoyPep();
        GeneratePepSeqList();
    }

    private void GeneratePepSeqList() {
        PeptideList = new HashMap();
        for (PepIonID pepID : GetPepIonList().values()) {
            if (!PeptideList.containsKey(pepID.Sequence)) {
                PeptideList.put(pepID.Sequence, new HashMap<String, PepIonID>());
            }
            if (!PeptideList.get(pepID.Sequence).containsKey(pepID.GetKey())) {
                PeptideList.get(pepID.Sequence).put(pepID.GetKey(), pepID);
            }
        }
    }

    public void GenerateMappedPepSeqList() {
        MappedPeptideList = new HashMap<>();
        for (PepIonID pepID : GetMappedPepIonList().values()) {
            if (!MappedPeptideList.containsKey(pepID.Sequence)) {
                MappedPeptideList.put(pepID.Sequence, new HashMap<String, PepIonID>());
            }
            if (!MappedPeptideList.get(pepID.Sequence).containsKey(pepID.GetKey())) {
                MappedPeptideList.get(pepID.Sequence).put(pepID.GetKey(), pepID);
            }
        }
    }

    public void GenerateProteinByRefIDByProbPepSeq(LCMSID RefID, boolean UseMappedIon) {
        ProteinList.clear();
        AssignedPepIonList.clear();
        for (PepIonID pepIonID : GetPepIonList().values()) {
            pepIonID.Weight = 0f;
            pepIonID.GroupWeight = 0f;
        }
        for (ProtID protein : RefID.ProteinList.values()) {
            if (!protein.IsDecoy(DecoyTag)) {
                ProtID newprot = protein.CloneProtein();
                newprot.Probability = protein.Probability;
                for (String Pepseq : protein.ProtPepSeq) {
                    float Maxweight = 0f;
                    float MaxGweight = 0f;
                    for (PepIonID pep : protein.ProtPeptideID.values()) {
                        if (pep.Sequence.equals(Pepseq)) {
                            if (pep.Weight > Maxweight) {
                                Maxweight = pep.Weight;
                            }
                            if (pep.GroupWeight > MaxGweight) {
                                MaxGweight = pep.GroupWeight;
                            }
                        }
                    }
                    if (PeptideList.containsKey(Pepseq)) {
                        for (PepIonID pepIonID : PeptideList.get(Pepseq).values()) {
                            pepIonID.Weight = Maxweight;
                            pepIonID.GroupWeight = MaxGweight;
                            newprot.AddPeptideID(pepIonID);
                            newprot.IDByDBSearch = true;
                            if (!pepIonID.ParentProtID_ProtXML.contains(newprot)) {
                                pepIonID.ParentProtID_ProtXML.add(newprot);
                            }
                        }
                    }
                    if (UseMappedIon) {
                        if (MappedPeptideList.containsKey(Pepseq)) {
                            for (PepIonID pepIonID : MappedPeptideList.get(Pepseq).values()) {
                                pepIonID.Weight = Maxweight;
                                pepIonID.GroupWeight = MaxGweight;
                                newprot.AddPeptideID(pepIonID);
                                if (!pepIonID.ParentProtID_ProtXML.contains(newprot)) {
                                    pepIonID.ParentProtID_ProtXML.add(newprot);
                                }
                            }
                        }
                    }
                }
                if (!newprot.PeptideID.isEmpty()) {
                    ProteinList.put(newprot.getAccNo(), newprot);
                }
            }
        }
    }

    public void GenerateProteinByRefIDByTrypticPep(LCMSID RefID, boolean UseMappedIon) {
        ProteinList.clear();
        AssignedPepIonList.clear();
        for (PepIonID pepIonID : GetPepIonList().values()) {
            pepIonID.Weight = 0f;
            pepIonID.GroupWeight = 0f;
        }
        for (ProtID protein : RefID.ProteinList.values()) {
            if (!protein.IsDecoy(DecoyTag)) {
                ProtID newprot = protein.CloneProtein();
                newprot.Probability = protein.Probability;
                for (String Pepseq : protein.TheoPeptides) {
                    float Maxweight = 0f;
                    float MaxGweight = 0f;
                    for (PepIonID pep : protein.ProtPeptideID.values()) {
                        if (pep.Sequence.equals(Pepseq)) {
                            if (pep.Weight > Maxweight) {
                                Maxweight = pep.Weight;
                            }
                            if (pep.GroupWeight > MaxGweight) {
                                MaxGweight = pep.GroupWeight;
                            }
                        }
                    }
                    if (PeptideList.containsKey(Pepseq)) {
                        for (PepIonID pepIonID : PeptideList.get(Pepseq).values()) {
                            pepIonID.Weight = Maxweight;
                            pepIonID.GroupWeight = MaxGweight;
                            newprot.AddPeptideID(pepIonID);
                            newprot.IDByDBSearch = true;
                            if (!pepIonID.ParentProtID_ProtXML.contains(newprot)) {
                                pepIonID.ParentProtID_ProtXML.add(newprot);
                            }
                        }
                    }
                    if (UseMappedIon) {
                        if (MappedPeptideList.containsKey(Pepseq)) {
                            for (PepIonID pepIonID : MappedPeptideList.get(Pepseq).values()) {
                                pepIonID.Weight = Maxweight;
                                pepIonID.GroupWeight = MaxGweight;
                                newprot.AddPeptideID(pepIonID);
                                if (!pepIonID.ParentProtID_ProtXML.contains(newprot)) {
                                    pepIonID.ParentProtID_ProtXML.add(newprot);
                                }
                            }
                        }
                    }
                }
                if (!newprot.PeptideID.isEmpty()) {
                    ProteinList.put(newprot.getAccNo(), newprot);
                }
            }
        }
    }

    public void GenerateProteinByRefIDByPepSeq(LCMSID RefID, boolean UseMappedIon) {
        ProteinList.clear();
        AssignedPepIonList.clear();
        for (PepIonID pepIonID : GetPepIonList().values()) {
            pepIonID.Weight = 0f;
            pepIonID.GroupWeight = 0f;
        }
        for (ProtID protein : RefID.ProteinList.values()) {            
            if (!protein.IsDecoy(DecoyTag)) {
                ProtID newprot = protein.CloneProtein();
                newprot.Probability = protein.Probability;
                newprot.ProtPepSeq = (ArrayList<String>) protein.ProtPepSeq.clone();

                for (PepIonID pep : protein.ProtPeptideID.values()) {                        
                    if (PeptideList.containsKey(pep.Sequence)) {
                        for (PepIonID pepIonID : PeptideList.get(pep.Sequence).values()) {
                            pepIonID.Weight = pep.Weight;
                            pepIonID.GroupWeight = pep.GroupWeight;
                            newprot.AddPeptideID(pepIonID);
                            newprot.IDByDBSearch = true;
                            if (!pepIonID.ParentProtID_ProtXML.contains(newprot)) {
                                pepIonID.ParentProtID_ProtXML.add(newprot);
                            }
                        }
                    }
                    if (UseMappedIon) {
                        if (MappedPeptideList.containsKey(pep.Sequence)) {
                            for (PepIonID pepIonID : MappedPeptideList.get(pep.Sequence).values()) {
                                pepIonID.Weight = pep.Weight;
                                pepIonID.GroupWeight = pep.GroupWeight;
                                newprot.AddPeptideID(pepIonID);
                                if (!pepIonID.ParentProtID_ProtXML.contains(newprot)) {
                                    pepIonID.ParentProtID_ProtXML.add(newprot);
                                }
                            }
                        }
                    }
                }
                if (!newprot.PeptideID.isEmpty()) {
                    ProteinList.put(newprot.getAccNo(), newprot);
                }
            }
        }
    }

    public void FilterProteinByRefID(LCMSID RefID) {

        ArrayList<ProtID> Premovelist = new ArrayList<>();

        for (ProtID protein : ProteinList.values()) {
            boolean found = false;
            if (!protein.IsDecoy(DecoyTag)) {
                for (String Indisprot : protein.IndisProteins) {
                    if (RefID.IndisProteinIDList.containsKey(Indisprot)) {
                        found = true;
                        break;
                    }
                }
            }
            if (!found) {
                Premovelist.add(protein);
            }
        }
        for (ProtID protein : Premovelist) {
            ProteinList.remove(protein.getAccNo());
        }
    }

    public void FilterPepIonByRefID(LCMSID RefID) {
        ArrayList<PepIonID> removelist = new ArrayList<>();
        for (PepIonID pep : GetPepIonList().values()) {
            if (!RefID.ProtXMLPepIonList.containsKey(pep.GetKey())) {
                removelist.add(pep);
            }
        }
        for (PepIonID pep : removelist) {
            GetPepIonList().remove(pep.GetKey());
        }
    }

    public void FilterBySpecDecoyFDR(String DecoyTag, float fdr) {
        this.DecoyTag = DecoyTag;
        this.FDR = fdr;
        //FindExpectThresholdByFDR();
        //RemoveHighExpPSM();
        FindSpecProbThresholdByFDR();
        RemoveLowProbPSM();
        //RemoveDecoy();
    }

    public void GenerateIndisProtMap() {
        IndisProteinIDList.clear();
        for (ProtID protID : ProteinList.values()) {
            for (String id : protID.IndisProteins) {
                if (!id.startsWith(DecoyTag)) {
                    IndisProteinIDList.put(id, protID);
                }
            }
        }
    }

    public void FilterByProteinDecoyFDRUsingMaxLocalPW(ArrayList<LCMSID> ProtIDList, String DecoyTag, float fdr) {
        this.DecoyTag = DecoyTag;
        this.ProteinFDR = fdr;
        for (ProtID protein : ProteinList.values()) {

            for (LCMSID protid : ProtIDList) {
                if (protid.IndisProteinIDList.containsKey(protein.getAccNo())) {
                    if (protein.MaxLocalPW < protid.IndisProteinIDList.get(protein.getAccNo()).Probability) {
                        protein.MaxLocalPW = protid.IndisProteinIDList.get(protein.getAccNo()).Probability;
                    }
                }
            }
        }
        FindProteinMaxLocalPWThresholdByFDR();
        RemoveLowMaxLocaPWProteinDecoy();
    }

    public void UpdateProteinKey() throws ClassNotFoundException, InterruptedException, IOException, XmlPullParserException {
        ArrayList<ProtID> temp = new ArrayList<>();
        for (ProtID protein : ProteinList.values()) {
            temp.add(protein);
        }
        ProteinList.clear();
        for (ProtID protein : temp) {
            AddProtID(protein);
        }
    }

    public void FixProteinWithDecoyHead() throws ClassNotFoundException, InterruptedException, IOException, XmlPullParserException {
        for (ProtID protein : ProteinList.values()) {
            if (protein.IsDecoy(DecoyTag)) {
                for (int i = 0; i < protein.IndisProteins.size(); i++) {
                    if (!protein.IndisProteins.get(i).startsWith(DecoyTag)) {
                        protein.setAccNo(protein.IndisProteins.get(i));
                        protein.SetDescription(protein.IndisProtDes.get(i));
                        break;
                    }
                }
            }
        }
        UpdateProteinKey();
    }

    public void UpdateDecoyMaxIniProb() {
        for (ProtID protein : ProteinList.values()) {
            if (protein.IsDecoy(DecoyTag)) {
                protein.MaxIniProb = 0f;
                for (PepIonID pep : protein.ProtPeptideID.values()) {
                    boolean include = true;
                    for (String prot : pep.ParentProtString_ProtXML) {
                        if (!prot.startsWith(DecoyTag)) {
                            include = false;
                        }
                    }
                    if (include && protein.MaxIniProb < pep.MaxProbability) {
                        protein.MaxIniProb = pep.MaxProbability;
                    }
                }
            }
        }
    }

    public void SetGroupProbForNonDecoyGroupHead() {
        for (ArrayList<ProtID> group : ProteinGroups.values()) {
            boolean allzero = true;
            for (ProtID protein : group) {
                if (protein.Probability > 0) {
                    allzero = false;
                    break;
                }
            }
            if (allzero) {
                for (ProtID protein : group) {
                    if (!protein.IsDecoy(DecoyTag)) {
                        protein.Probability = protein.GroupProb;
                        break;
                    }
                }
            }
        }
    }

    public void FilterByProteinDecoyFDRUsingMaxIniProb(String DecoyTag, float fdr) {
        this.DecoyTag = DecoyTag;
        this.ProteinFDR = fdr;
        FindMaxIniProbThresholdByFDR();
        RemoveLowMaxIniProbProteinDecoy();
    }

    public void FilterByProteinDecoyFDRUsingLocalPW(String DecoyTag, float fdr) {
        this.DecoyTag = DecoyTag;
        this.ProteinFDR = fdr;
        FindLocalPWThresholdByFDR();
        RemoveLowProbProteinDecoy();
    }

    public void GeneratePepXMLProtList() {
        PepXMLProteinList = new HashMap<>();
        for (PepIonID pepIonID : PepIonList.values()) {
            InsertProteinPepXML(pepIonID);
        }
    }

    private void InsertProteinPepXML(PepIonID pepIonID) {
        for (String prot : pepIonID.ParentProtID_PepXML) {
            if (!PepXMLProteinList.containsKey(prot)) {
                ProtID protein = new ProtID();
                protein.setAccNo(prot);
                PepXMLProteinList.put(prot, protein);
            }

            if (!PepXMLProteinList.get(prot).PeptideID.containsKey(pepIonID.GetKey())) {
                PepXMLProteinList.get(prot).PeptideID.put(pepIonID.GetKey(), pepIonID);
            }
        }
    }

    public void AddPeptideID(PepIonID pepID) {
        if (!PepIonList.containsKey(pepID.GetKey())) {
            pepID.Index = PepIonList.size();
            PepIonList.put(pepID.GetKey(), pepID);
        }
        for (PSM psm : pepID.GetPSMList()) {
            PSMList.put(psm.SpecNumber, psm);
        }
    }

    public PSM GetPSM(int SpecNum) {
        return PSMList.get(SpecNum);
    }

    public void AddModification(PTM ptm, String site) {
        ModificationInfo modification = new ModificationInfo();
        modification.site = site;
        modification.massdiff = (float) ptm.getMass();
        modification.mass = (float) (ptm.getMass() + AminoAcid.getAminoAcid(site).monoisotopicMass);
        modification.modification = ptm;
        if (!ModificationList.containsKey(ptm.getName() + "_" + site)) {
            ModificationList.put(ptm.getName() + "_" + site, modification);
        }
    }

    public void CreateProtByRefID(LCMSID RefProtID) {
        for (ProtID protein : RefProtID.ProteinList.values()) {
            ProtID newprot = protein.CloneProtein();
            newprot.Probability = protein.Probability;
            for (String pepseq : protein.ProtPepSeq) {
                if (PeptideList.containsKey(pepseq)) {
                    for (PepIonID pepion : PeptideList.get(pepseq).values()) {
                        newprot.PeptideID.put(pepion.GetKey(), pepion);
                    }
                }
            }
            ProteinList.put(newprot.getAccNo(), newprot);
        }
    }

    public void CreateInstanceForAllPepIon() {
        for (PepIonID pepIonID : GetPepIonList().values()) {
            pepIonID.CreateQuantInstance(4);
        }
    }

    public void ClearMappedPep() {
        MappedPepIonList.clear();
        MappedPeptideList.clear();
        if (MappedPepIonIndexList != null) {
            MappedPepIonIndexList.clear();
        }
    }

    public void ClearProPeplist() {
        for (ProtID protein : ProteinList.values()) {
            protein.PeptideID.clear();
        }
        for (PepIonID pep : GetPepIonList().values()) {
            pep.ParentProtID_ProtXML.clear();
        }
        for (PepIonID pep : GetMappedPepIonList().values()) {
            pep.ParentProtID_ProtXML.clear();
        }
    }

    public void RemoveLowProbMappedIon(float ProbThreshold) {
        //LowScorePep = new HashMap<>();
        ArrayList<PepIonID> removelist = new ArrayList<>();
        for (PepIonID pep : GetMappedPepIonList().values()) {
            if (pep.TargetedProbability() < ProbThreshold) {
                removelist.add(pep);
            }
        }
        for (PepIonID pep : removelist) {
            GetMappedPepIonList().remove(pep.GetKey());
            //LowScorePep.put(pep.GetKey(), pep);
        }
        GenerateMappedPepSeqList();
    }

    public void ReleaseIDs() {
        PSMList = null;
        LowScorePSMByPepKey = null;
        LowScorePep = null;
        LowScorePSM = null;
        PepIonList = null;
        PepIonIndexList = null;
        MappedPepIonIndexList = null;
        AssignedPepIonList = null;
        ProtXMLPepIonList = null;
        ProteinList = null;
        IndisProteinIDList = null;
        MappedPepIonList = null;
        PepXMLProteinList = null;
        PeptideList = null;
        MappedPeptideList = null;
    }

    public float GetRFactor(float ProbThreshold) {
        int forward = 0;
        int decoy = 0;
        for (ProtID protID : ProteinList.values()) {
            if (protID.Probability < ProbThreshold) {
                if (protID.IsDecoy(DecoyTag)) {
                    decoy++;
                } else {
                    forward++;
                }
            }
        }
        float rf = (float) forward / decoy;
        Logger.getRootLogger().info("Caculating R factor: probability threshold =" + ProbThreshold);
        Logger.getRootLogger().info("R factor=" + rf + " (forward/decoy=" + forward + "/" + decoy + ")");
        return rf;
    }
}
