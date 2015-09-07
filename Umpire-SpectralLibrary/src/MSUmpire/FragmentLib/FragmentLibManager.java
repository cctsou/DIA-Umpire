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
package MSUmpire.FragmentLib;

import MSUmpire.MySQLTool.ConnectionManager;
import MSUmpire.PSMDataStructure.FragmentPeak;
import MSUmpire.PSMDataStructure.FragmentPeakGroup;
import MSUmpire.PSMDataStructure.FragmentSelection;
import MSUmpire.PSMDataStructure.LCMSID;
import MSUmpire.PSMDataStructure.ModStringConvert;
import MSUmpire.PSMDataStructure.ModificationInfo;
import MSUmpire.PSMDataStructure.PTMManager;
import MSUmpire.PSMDataStructure.PepFragmentLib;
import MSUmpire.PSMDataStructure.PepIonID;
import MSUmpire.SeqUtility.ShuffledSeqGen;
import com.compomics.util.experiment.biology.AminoAcid;
import com.compomics.util.experiment.biology.Ion;
import com.compomics.util.experiment.biology.IonFactory;
import com.compomics.util.experiment.biology.Peptide;
import com.compomics.util.experiment.biology.ions.PeptideFragmentIon;
import com.compomics.util.experiment.identification.matches.ModificationMatch;
import jaligner.matrix.Matrix;
import jaligner.matrix.MatrixLoader;
import jaligner.matrix.MatrixLoaderException;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
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
import java.util.HashMap;
import java.util.Random;
import org.apache.commons.io.FilenameUtils;
import org.apache.commons.lang.exception.ExceptionUtils;
import org.apache.log4j.Logger;
import org.hupo.psi.ms.traml.CvParamType;
import org.hupo.psi.ms.traml.ModificationType;
import org.hupo.psi.ms.traml.PeptideType;
import org.hupo.psi.ms.traml.RetentionTimeType;
import org.hupo.psi.ms.traml.TransitionType;
import org.nustaq.serialization.FSTObjectInput;
import org.nustaq.serialization.FSTObjectOutput;
import org.systemsbiology.apps.tramlparser.TraMLParser;
import org.xmlpull.v1.XmlPullParserException;

/**
 *
 * @author Chih-Chiang Tsou <chihchiang.tsou@gmail.com>
 */
public class FragmentLibManager implements Serializable {
    private static final long serialVersionUID = -1866384504919716378L;

    transient ConnectionManager connectionManager;
    public HashMap<String, PepFragmentLib> PeptideFragmentLib = new HashMap<>();
    HashMap<String, PepFragmentLib> PeptideDecoyFragmentLib = new HashMap<>();
        
    public String LibID = "Test";
    transient FragmentSelection fragselection;
    
    public void ReduceMemoryUsage(){
        for(PepFragmentLib pep : PeptideFragmentLib.values()){
            for(FragmentPeakGroup frag: pep.FragmentGroups.values()){
                frag.ClearGroups();
            }
        }
        for(PepFragmentLib pep : PeptideDecoyFragmentLib.values()){
            for(FragmentPeakGroup frag: pep.FragmentGroups.values()){
                frag.ClearGroups();
            }
        }
    }

    public FragmentLibManager(String LibID) {
        this.LibID = LibID;
    }

    public FragmentLibManager(String LibID, ConnectionManager connectionManager) {
        this.LibID = LibID;
        this.connectionManager = connectionManager;
    }

    public void WriteFragmentLibSerialization(String path) {
        //JavaSerializationFragmentLibWrite(path, LibID);
        FSFragmentLibWrite(path, LibID);
    }

    private void JavaSerializationFragmentLibWrite(String path, String LibID1) {
        try {
            Logger.getRootLogger().info("Writing FragmentLib to file:" + path + LibID1 + ".ser...");
            FileOutputStream fout = new FileOutputStream(path + LibID1 + ".ser", false);
            ObjectOutputStream oos = new ObjectOutputStream(fout);
            oos.writeObject(this);
            oos.close();
            fout.close();
        } catch (Exception ex) {
            Logger.getRootLogger().error(ExceptionUtils.getStackTrace(ex));
        }
    }

    private void FSFragmentLibWrite(String path, String LibID1) {
        try {
            Logger.getRootLogger().info("Writing FragmentLib to file:" + path + LibID1 + ".serFS...");
            FileOutputStream fout = new FileOutputStream(path + LibID1 + ".serFS", false);
            FSTObjectOutput oos = new FSTObjectOutput(fout);
            oos.writeObject(this);
            oos.close();
            fout.close();
        } catch (Exception ex) {
            Logger.getRootLogger().error(ExceptionUtils.getStackTrace(ex));
        }
    }

    public static FragmentLibManager ReadFragmentLibSerialization(String path, String LibID) {
        FragmentLibManager lib = FSFragmentLibRead(path, LibID);
        if (lib == null) {
            //lib = JavaSerializationFragmentLibRead(path, LibID);
            lib=FSFragmentLibRead_Old(path, LibID);
            if (lib != null) {
                lib.WriteFragmentLibSerialization(path);
            }
        }
        if(lib!=null){
            lib.LibID=LibID;
        }
        return lib;
    }

    private static FragmentLibManager FSFragmentLibRead(String path, String LibID1) {
        if (!new File(path + LibID1 + ".serFS").exists()) {
            Logger.getRootLogger().debug(path + LibID1 + ".serFS does not exsit.");
            return null;
        }
        try {
            Logger.getRootLogger().info("Reading spectral library from file:" + path + LibID1 + ".serFS...");
            FileInputStream fileIn = new FileInputStream(path + LibID1 + ".serFS");
            FSTObjectInput in = new FSTObjectInput(fileIn);
            FragmentLibManager FragLib = (FragmentLibManager) in.readObject();
            in.close();
            fileIn.close();            
            return FragLib;
        } catch (Exception ex) {
            Logger.getRootLogger().error(ExceptionUtils.getStackTrace(ex));
            return null;
        }
    }
    
    private static FragmentLibManager FSFragmentLibRead_Old(String path, String LibID1) {
        if (!new File(path + LibID1 + ".serFS").exists()) {
            Logger.getRootLogger().debug(path + LibID1 + ".serFS does not exsit.");
            return null;
        }
        try {
            Logger.getRootLogger().info("Reading internal spectral library from file:" + path + LibID1 + ".serFS...");
            FileInputStream fileIn = new FileInputStream(path + LibID1 + ".serFS");
            org.nustaq_old.serialization.FSTObjectInput in = new org.nustaq_old.serialization.FSTObjectInput(fileIn);
            FragmentLibManager FragLib = (FragmentLibManager) in.readObject();
            in.close();
            fileIn.close();
            return FragLib;
        } catch (Exception ex) {
            Logger.getRootLogger().error(ExceptionUtils.getStackTrace(ex));
            return null;
        }
    }

    private static FragmentLibManager JavaSerializationFragmentLibRead(String path, String LibID1) {
        if (!new File("Reading FragmentLib from file:" + path + LibID1 + ".ser...").exists()) {
            return null;
        }
        try {
            Logger.getRootLogger().info("Reading FragmentLib from file:" + path + LibID1 + ".ser...");
            FileInputStream fileIn = new FileInputStream(path + LibID1 + ".ser");
            ObjectInputStream in = new ObjectInputStream(fileIn);
            FragmentLibManager FragLib = (FragmentLibManager) in.readObject();
            in.close();
            fileIn.close();
            return FragLib;
        } catch (Exception ex) {
            Logger.getRootLogger().error(ExceptionUtils.getStackTrace(ex));
            return null;
        }
    }

    

    public PepFragmentLib GetFragmentLib(String PepKey) {
        return PeptideFragmentLib.get(PepKey);
    }

    public PepFragmentLib GetDecoyFragmentLib(String PepKey) {
        return PeptideDecoyFragmentLib.get("decoy_" + PepKey);
    }

    private PepFragmentLib GenerateDecoy(PepFragmentLib target, String recoyseq) {
        PepFragmentLib fragmentLibdecoy = new PepFragmentLib();
        fragmentLibdecoy.Sequence = recoyseq;
        fragmentLibdecoy.ModSequence = "decoy_" + target.ModSequence;
        //fragmentLibdecoy.ModificationString="decoy_" +target.ModificationString;
        fragmentLibdecoy.Charge = target.Charge;
        fragmentLibdecoy.PrecursorMz = target.PrecursorMz;
        ArrayList<ModificationMatch> Modifications = PTMManager.TranslateModificationString(target.ModificationString);

        for (ModificationMatch mod : Modifications) {
            if (mod.getModificationSite() > 0 && mod.getModificationSite() <= target.Sequence.length()) {
                char site = target.Sequence.charAt(mod.getModificationSite() - 1);
                ArrayList<Integer> Position = new ArrayList<>();
                for (int i = 0; i < recoyseq.length(); i++) {
                    if (recoyseq.charAt(i) == site) {
                        Position.add(i + 1);
                    }
                }
                int rnd = new Random().nextInt(Position.size());
                mod.setModificationSite(Position.get(rnd));
            }
        }

        Peptide peptide = peptide = new Peptide(recoyseq, Modifications);
        HashMap<Integer, HashMap<Integer, ArrayList<Ion>>> allfragment = IonFactory.getInstance().getFragmentIons(peptide);        
        fragmentLibdecoy.FragmentGroups = target.CloneFragmentGroup();
        ArrayList<Ion> fragments=new ArrayList<>();
        fragments.addAll(allfragment.get(Ion.IonType.PEPTIDE_FRAGMENT_ION.index).get(PeptideFragmentIon.B_ION));
        fragments.addAll(allfragment.get(Ion.IonType.PEPTIDE_FRAGMENT_ION.index).get(PeptideFragmentIon.Y_ION));
        
        for (Ion frag : fragments) {
            float targetmz = (float) frag.getTheoreticMz(1);
            String IonString = frag.getSubTypeAsString() + ((PeptideFragmentIon) frag).getNumber() + "_1";
            if (fragmentLibdecoy.FragmentGroups.containsKey(IonString)) {
                FragmentPeakGroup fragmentPeakGroup = fragmentLibdecoy.FragmentGroups.get(IonString);
                fragmentPeakGroup.FragMZ = targetmz;
            }
            targetmz = (float) frag.getTheoreticMz(2);
            IonString = frag.getSubTypeAsString() + ((PeptideFragmentIon) frag).getNumber() + "_2";
            if (fragmentLibdecoy.FragmentGroups.containsKey(IonString)) {
                FragmentPeakGroup fragmentPeakGroup = fragmentLibdecoy.FragmentGroups.get(IonString);
                fragmentPeakGroup.FragMZ = targetmz;
            }
        }
        return fragmentLibdecoy;
    }

    public void GenerateDecoyLib() throws MatrixLoaderException {
        Logger.getRootLogger().info("generating decoy spectra");       
        PeptideDecoyFragmentLib = new HashMap<>();
        Matrix blosum62=MatrixLoader.load("BLOSUM62");                
        for (PepFragmentLib fragmentLib : PeptideFragmentLib.values()) {            
            ShuffledSeqGen shufen=new ShuffledSeqGen(fragmentLib.Sequence,blosum62,null);
            shufen.Generate();
            String decoyseq=shufen.decoy;
            PepFragmentLib fragmentLibdecoy = GenerateDecoy(fragmentLib, decoyseq);
            PeptideDecoyFragmentLib.put(fragmentLibdecoy.GetKey(), fragmentLibdecoy);            
        }
    }

    public void ImportDB(ConnectionManager connectionManager, String Path) throws IOException, SQLException {
        ImportLib(connectionManager, Path);
        ImportDecoyLib(connectionManager, Path);
    }

    private void ImportLib(ConnectionManager connectionManager, String Path) throws IOException, SQLException {
        ExportFragLibCSV(Path);
        this.connectionManager = connectionManager;
        CreateFragLibTable();
        ExportTableDBBulkLoader(Path, LibID + "_PepInfo.csv", false);
        ExportTableDBBulkLoader(Path, LibID + "_FragLib.csv", false);
    }

    private void ImportDecoyLib(ConnectionManager connectionManager, String Path) throws IOException, SQLException {
        ExportdecoyFragLibCSV(Path);
        this.connectionManager = connectionManager;
        CreateDecoyFragLibTable();
        ExportTableDBBulkLoader(Path, LibID + "_DecoyPepInfo.csv", false);
        ExportTableDBBulkLoader(Path, LibID + "_DecoyFragLib.csv", false);
    }

    private void ExportdecoyFragLibCSV(String Path) throws IOException {
        FileWriter writer = new FileWriter(Path + LibID + "_DecoyPepInfo.csv");
        FileWriter writer2 = new FileWriter(Path + LibID + "_DecoyFragLib.csv");
        writer.write("Sequence,ModSeq,ModificationString,Charge,PrecursorMz,MaxProbability,MS1Score\n");
        writer2.write("PepKey,IonType,FragMz,Corr,PPM,Intensity,ApexDelta,RTOverlapP\n");
        for (PepFragmentLib PepFrag : PeptideDecoyFragmentLib.values()) {
            writer.write(PepFrag.Sequence + "," + PepFrag.ModSequence + "," + PepFrag.ModificationString + "," + PepFrag.Charge + "," + PepFrag.PrecursorMz + "," + PepFrag.MaxProbability + "," + PepFrag.MS1Score + "\n");
            for (FragmentPeakGroup frag : PepFrag.FragmentGroups.values()) {
                writer2.write(PepFrag.GetKey() + "," + frag.GetFragKey() + "," + frag.FragMZ + "," + frag.GetCorrString() + "," + frag.GetPPMString() + "," + frag.GetIntString() + "," + frag.GetApexDeltaString() + "," + frag.GetRTOverlapString() + "\n");
            }
        }
        writer.close();
        writer2.close();
    }

    private void ExportFragLibCSV(String Path) throws IOException {
        FileWriter writer = new FileWriter(Path + LibID + "_PepInfo.csv");
        FileWriter writer2 = new FileWriter(Path + LibID + "_FragLib.csv");
        writer.write("Sequence,ModSeq,ModificationString,Charge,PrecursorMz,MaxProbability,MS1Score\n");
        writer2.write("PepKey,IonType,FragMz,Corr,PPM,Intensity,ApexDelta,RTOverlapP\n");
        for (PepFragmentLib PepFrag : PeptideFragmentLib.values()) {
            writer.write(PepFrag.Sequence + "," + PepFrag.ModSequence + "," + PepFrag.ModificationString + "," + PepFrag.Charge + "," + PepFrag.PrecursorMz + "," + PepFrag.MaxProbability + "," + PepFrag.MS1Score + "\n");
            for (FragmentPeakGroup frag : PepFrag.FragmentGroups.values()) {
                writer2.write(PepFrag.GetKey() + "," + frag.GetFragKey() + "," + frag.FragMZ + "," + frag.GetCorrString() + "," + frag.GetPPMString() + "," + frag.GetIntString() + "," + frag.GetApexDeltaString() + "," + frag.GetRTOverlapString() + "\n");
            }
        }
        writer.close();
        writer2.close();
    }

    private void ReadDecoyLibFromDB() throws SQLException {
        PeptideDecoyFragmentLib = new HashMap<>();
        Connection connection = connectionManager.GetConnection();
        try (Statement state = connection.createStatement()) {
            ResultSet rsQuery = state.executeQuery("SELECT * FROM " + LibID + "_DecoyPepInfo");
            while (rsQuery.next()) {
                PepFragmentLib PepFrag = new PepFragmentLib();
                PepFrag.Sequence = rsQuery.getString("Sequence");
                PepFrag.ModSequence = rsQuery.getString("ModSeq");
                PepFrag.ModificationString = rsQuery.getString("ModificationString");
                PepFrag.Charge = rsQuery.getInt("Charge");
                PepFrag.PrecursorMz = rsQuery.getFloat("PrecursorMz");
                PepFrag.MS1Score = rsQuery.getFloat("MS1Score");
                PepFrag.MaxProbability = rsQuery.getFloat("MaxProbability");
                PeptideDecoyFragmentLib.put(PepFrag.GetKey(), PepFrag);
            }

            rsQuery = state.executeQuery("SELECT * FROM " + LibID + "_DecoyFragLib");
            while (rsQuery.next()) {
                PepFragmentLib PepFrag = PeptideDecoyFragmentLib.get(rsQuery.getString("PepKey"));
                FragmentPeakGroup fragmentgroup = new FragmentPeakGroup();
                String IonString = rsQuery.getString("IonType");
                fragmentgroup.IonType = IonString.split("_")[0];
                fragmentgroup.Charge = Integer.parseInt(IonString.split("_")[1]);
                fragmentgroup.FragMZ = rsQuery.getFloat("FragMz");
                String[] info = rsQuery.getString("Corr").split(";");
                for (int i = 0; i < info.length; i++) {
                    fragmentgroup.CorrGroup.add(Float.parseFloat(info[i]));
                }
                info = rsQuery.getString("PPM").split(";");
                for (int i = 0; i < info.length; i++) {
                    fragmentgroup.PPMGroup.add(Float.parseFloat(info[i]));
                }
                info = rsQuery.getString("Intensity").split(";");
                for (int i = 0; i < info.length; i++) {
                    fragmentgroup.IntensityGroup.add(Float.parseFloat(info[i]));
                }
                info = rsQuery.getString("ApexDelta").split(";");
                for (int i = 0; i < info.length; i++) {
                    fragmentgroup.ApexDeltaGroup.add(Float.parseFloat(info[i]));
                }
                info = rsQuery.getString("RTOverlapP").split(";");
                for (int i = 0; i < info.length; i++) {
                    fragmentgroup.RTOverlapPGroup.add(Float.parseFloat(info[i]));
                }
                PepFrag.FragmentGroups.put(fragmentgroup.GetFragKey(), fragmentgroup);
            }
        }
        connection.close();
    }

    public void ReadFromDB() throws SQLException {
        ReadFragLibFromDB();
        ReadDecoyLibFromDB();
    }

    private void ReadFragLibFromDB() throws SQLException {
        Connection connection = connectionManager.GetConnection();
        try (Statement state = connection.createStatement()) {
            ResultSet rsQuery = state.executeQuery("SELECT * FROM " + LibID + "_PepInfo");
            while (rsQuery.next()) {
                PepFragmentLib PepFrag = new PepFragmentLib();
                PepFrag.Sequence = rsQuery.getString("Sequence");
                PepFrag.ModSequence = rsQuery.getString("ModSeq");
                PepFrag.MS1Score = rsQuery.getFloat("MS1Score");
                PepFrag.ModificationString = rsQuery.getString("ModificationString");
                PepFrag.Charge = rsQuery.getInt("Charge");
                PepFrag.PrecursorMz = rsQuery.getFloat("PrecursorMz");
                PepFrag.MaxProbability = rsQuery.getFloat("MaxProbability");
                PeptideFragmentLib.put(PepFrag.GetKey(), PepFrag);
            }

            rsQuery = state.executeQuery("SELECT * FROM " + LibID + "_FragLib");
            while (rsQuery.next()) {
                PepFragmentLib PepFrag = PeptideFragmentLib.get(rsQuery.getString("PepKey"));
                FragmentPeakGroup fragmentgroup = new FragmentPeakGroup();
                String IonString = rsQuery.getString("IonType");
                fragmentgroup.IonType = IonString.split("_")[0];
                fragmentgroup.Charge = Integer.parseInt(IonString.split("_")[1]);
                fragmentgroup.FragMZ = rsQuery.getFloat("FragMz");
                String[] info = rsQuery.getString("Corr").split(";");
                for (int i = 0; i < info.length; i++) {
                    fragmentgroup.CorrGroup.add(Float.parseFloat(info[i]));
                }
                info = rsQuery.getString("PPM").split(";");
                for (int i = 0; i < info.length; i++) {
                    fragmentgroup.PPMGroup.add(Float.parseFloat(info[i]));
                }
                info = rsQuery.getString("Intensity").split(";");
                for (int i = 0; i < info.length; i++) {
                    fragmentgroup.IntensityGroup.add(Float.parseFloat(info[i]));
                }
                info = rsQuery.getString("ApexDelta").split(";");
                for (int i = 0; i < info.length; i++) {
                    fragmentgroup.ApexDeltaGroup.add(Float.parseFloat(info[i]));
                }
                info = rsQuery.getString("RTOverlapP").split(";");
                for (int i = 0; i < info.length; i++) {
                    fragmentgroup.RTOverlapPGroup.add(Float.parseFloat(info[i]));
                }
                PepFrag.FragmentGroups.put(fragmentgroup.GetFragKey(), fragmentgroup);
            }
        }
        connection.close();
    }

    protected void CreateFragLibTable() throws SQLException {
        Connection connection = this.connectionManager.GetConnection();
        connection.createStatement().execute("DROP TABLE IF EXISTS " + LibID + "_PepInfo;");
        connection.createStatement().execute("DROP TABLE IF EXISTS " + LibID + "_FragLib;");
        connection.createStatement().execute("CREATE TABLE " + LibID + "_PepInfo (Sequence VARCHAR(100), ModSeq VARCHAR(200), ModificationString VARCHAR(200), Charge INT NOT NULL, PrecursorMz DOUBLE NOT NULL, MaxProbability DOUBLE, MS1Score DOUBLE, PRIMARY KEY (ModSeq, Charge));");
        connection.createStatement().execute("CREATE TABLE " + LibID + "_FragLib (PepKey VARCHAR(200), IonType VARCHAR(10), FragMz DOUBLE NOT NULL, Corr TEXT, PPM TEXT, Intensity TEXT,ApexDelta TEXT,RTOverlapP TEXT, PRIMARY KEY (PepKey, IonType));");
        connection.close();
    }

    protected void CreateDecoyFragLibTable() throws SQLException {
        Connection connection = this.connectionManager.GetConnection();
        connection.createStatement().execute("DROP TABLE IF EXISTS " + LibID + "_DecoyPepInfo;");
        connection.createStatement().execute("DROP TABLE IF EXISTS " + LibID + "_DecoyFragLib;");
        connection.createStatement().execute("CREATE TABLE " + LibID + "_DecoyPepInfo (Sequence VARCHAR(100), ModSeq VARCHAR(200), ModificationString VARCHAR(200), Charge INT NOT NULL, PrecursorMz DOUBLE NOT NULL, MaxProbability DOUBLE, MS1Score DOUBLE, PRIMARY KEY (ModSeq, Charge));");
        connection.createStatement().execute("CREATE TABLE " + LibID + "_DecoyFragLib (PepKey VARCHAR(200), IonType VARCHAR(10), FragMz DOUBLE NOT NULL, Corr TEXT, PPM TEXT, Intensity TEXT,ApexDelta TEXT,RTOverlapP TEXT, PRIMARY KEY (PepKey, IonType));");
        connection.close();
    }

    protected void ExportTableDBBulkLoader(String Path, String Filename, boolean filedelete) throws SQLException {
        Connection connection = this.connectionManager.GetConnection();
        Statement state = connection.createStatement();
        state.executeUpdate("LOAD DATA LOCAL INFILE '" + Path + Filename + "' INTO TABLE " + FilenameUtils.getBaseName(Filename) + " FIELDS TERMINATED BY ','" + " LINES TERMINATED BY '\\n' IGNORE 1 LINES");
        state.close();
        connection.close();
        File file = new File(Path + Filename);
        if (filedelete) {
            file.delete();
        }
    }
  
    public void ImportFragLibByTraML(String tramlpath, String DecoyPrefix) throws IOException, XmlPullParserException, Exception {
        Logger.getRootLogger().info("Parsing " + tramlpath);
        try {
            TraMLParser traMLParser = new TraMLParser();
            traMLParser.parse_file(tramlpath, Logger.getRootLogger());
            PTMManager.GetInstance();
            HashMap<String, PepFragmentLib> TraMLMap = new HashMap<>();
            HashMap<String, PepFragmentLib> Decoys=new HashMap<>();

            for (PeptideType peptide : traMLParser.getTraML().getCompoundList().getPeptide()) {
                PepFragmentLib fraglib = new PepFragmentLib();
                fraglib.Sequence = peptide.getSequence();
                fraglib.ModSequence = peptide.getSequence();
                for (RetentionTimeType rt : peptide.getRetentionTimeList().getRetentionTime()) {
                    fraglib.RetentionTime.add(Float.parseFloat(rt.getCvParam().get(0).getValue()));
                }
                for (CvParamType param : peptide.getCvParam()) {
                    if (param.getName().equals("charge state")) {
                        fraglib.Charge = Integer.parseInt(param.getValue());
                    }
                }

                if (peptide.getModification() != null) {
                    for (ModificationType mod : peptide.getModification()) {
                        ModificationInfo modinfo = new ModificationInfo();
                        int idx=mod.getLocation();
                        if (idx== 0) {
                            modinfo.site = "N-term";
                            idx=1;
                        } else if (idx== peptide.getSequence().length() + 1) {
                            modinfo.site = "C-term";                            
                            idx=peptide.getSequence().length();
                            if (mod.getCvParam().get(0).getAccession().equals("UNIMOD:35")) {
                                modinfo.site = "M";
                            }
                        } else {
                            modinfo.site = String.valueOf(peptide.getSequence().charAt(idx - 1));
                        }
                        modinfo.modification = PTMManager.GetInstance().GetPTM(modinfo.site, mod.getMonoisotopicMassDelta().floatValue());
                        if(modinfo.modification==null){
                            Logger.getRootLogger().error("Modification was not found in the library: site:"+modinfo.site+", massdiff="+mod.getMonoisotopicMassDelta().floatValue());
                        }
                        modinfo.massdiff = (float) modinfo.modification.getMass();
                        modinfo.mass = (float) (modinfo.modification.getMass() + AminoAcid.getAminoAcid(modinfo.site).monoisotopicMass);
                        ModificationMatch modmatch = new ModificationMatch(modinfo.modification.getName(), true, idx - 1);
                        fraglib.Modifications.add(modmatch);
                        fraglib.ModSequence = ModStringConvert.AddModIntoSeqBeforeSite(fraglib.ModSequence, modinfo.GetKey(), idx - 1);
                    }
                }
                if (peptide.getId().startsWith(DecoyPrefix)) {
                    //PeptideDecoyFragmentLib.put("decoy_" +fraglib.GetKey(), fraglib);
                    fraglib.ModSequence ="decoy_" + fraglib.ModSequence;
                    Decoys.put(peptide.getId(), fraglib);
                } else {
                    PeptideFragmentLib.put(fraglib.GetKey(), fraglib);
                }
                TraMLMap.put(peptide.getId(), fraglib);
            }
            for (String key : Decoys.keySet()) {
                String pepkey = key.replace(DecoyPrefix + "_", "");
                Decoys.get(key).ModSequence = "decoy_" + TraMLMap.get(pepkey).ModSequence;
                PeptideDecoyFragmentLib.put(Decoys.get(key).GetKey(), Decoys.get(key));
            }

            HashMap<String, ArrayList<FragmentPeak>> TransitionList = new HashMap<>();
            HashMap<String, Float> PrecursorMZList = new HashMap<>();

            for (TransitionType trans : traMLParser.getTraML().getTransitionList().getTransition()) {
                String pepid = ((PeptideType) trans.getPeptideRef()).getId();
                PrecursorMZList.put(pepid, Float.parseFloat(trans.getPrecursor().getCvParam().get(0).getValue()));
                FragmentPeak fragment = new FragmentPeak();
                if (!trans.getUserParam().isEmpty()) {
                    fragment.IonType = trans.getUserParam().get(0).getValue().split("/")[0];
                }
                else{
                    fragment.IonType = trans.getId().split("_")[1]+"_"+trans.getId().split("_")[0];
                    fragment.IonType = fragment.IonType.replace("_noloss", "");
                }
                for (CvParamType cv : trans.getProduct().getCvParam()) {
                    if (cv.getName().equals("charge state")) {
                        fragment.Charge = Integer.parseInt(cv.getValue());
                    }
                    if (cv.getName().equals("isolation window target m/z")) {
                        fragment.FragMZ = Float.parseFloat(cv.getValue());
                    }
                }

                for (CvParamType cv : trans.getCvParam()) {
                    if (cv.getName().equals("product ion intensity")) {
                        fragment.intensity = Float.parseFloat(cv.getValue());
                    }
                }
                
                if(!TransitionList.containsKey(pepid)){
                    TransitionList.put(pepid, new ArrayList<FragmentPeak>());
                }
                TransitionList.get(pepid).add(fragment);
            }
            for (String pepid : TransitionList.keySet()) {
                PepFragmentLib fraglib = TraMLMap.get(pepid);
                fraglib.PrecursorMz = PrecursorMZList.get(pepid);
                fraglib.AddFragments(TransitionList.get(pepid));
            }
            Logger.getRootLogger().info("No. of peptide ions in the imported library:"+PeptideFragmentLib.size());
        } catch (MatrixLoaderException ex) {
            Logger.getRootLogger().error(ExceptionUtils.getStackTrace(ex));
        }
    }

    public void ImportFragLibBySPTXT(String sptxtpath) throws IOException, XmlPullParserException {
        Logger.getRootLogger().info("Parsing " + sptxtpath);
        try {
            BufferedReader reader=new BufferedReader(new FileReader(sptxtpath));
            
            PTMManager.GetInstance();

            String line="";
            boolean Header=false;
            boolean Peak=false;
            PepFragmentLib fraglib =null;
            ArrayList<FragmentPeak> FragmentPeaks=null;
            while ((line = reader.readLine()) != null) {
                if (line.startsWith("#")) {
                    continue;
                } else if (line.startsWith("Name:")) {
                    Header = true;
                    fraglib = new PepFragmentLib();
                    FragmentPeaks = new ArrayList<>();
                } else if (line.startsWith("FullName:") && Header) {
                    String temp = line.substring(line.indexOf(".") + 1);
                    String modstring = temp.substring(0, temp.indexOf("."));
                    fraglib.Sequence = modstring.replace("n[", "[").replace("c[", "").replaceAll("[\\[0-9\\]]", "");
                    fraglib.ModSequence = ModStringConvert.ConvertTPPModString(modstring, fraglib.Modifications);
                    fraglib.Charge = Integer.parseInt(line.split("/")[1].subSequence(0, 1).toString());
                } else if (line.startsWith("Comment:") && Header) {
                    //String RTs = line.split("RetentionTime=")[1].split(" ")[0];
                    String RTs = line.split("iRT=")[1].split(" ")[0];
                    for (String rt : RTs.split(",")) {
                        //fraglib.RetentionTime.add(Float.parseFloat(rt)/60f);
                        fraglib.RetentionTime.add(Float.parseFloat(rt));
                    }
                } else if (line.startsWith("PrecursorMZ:") && Header) {
                    fraglib.PrecursorMz = Float.parseFloat(line.split(":")[1].trim());
                } else if (line.startsWith("NumPeaks:") && Header) {
                    Peak = true;
                } else if (!"".equals(line) && Peak) {
                    String mz = line.split("\t")[0];
                    String intensity = line.split("\t")[1];
                    String type = line.split("\t")[2];
                                        
                    FragmentPeak fragment = new FragmentPeak();
                    fragment.IonType = "?";
                    float delta=Float.MAX_VALUE;
                    String iontype="";
                    for (String ion : type.split(",")) {
                        if (ion.startsWith("b") || ion.startsWith("y")) {
                            String temp = ion.split("/")[0].split("-")[0];
                            if (temp.contains("^")) {
                                fragment.IonType = temp.substring(0, temp.indexOf("^"));
                            } else {
                                fragment.IonType = temp.replace("i", "");
                            }
                            break;
                        }
                        fragment.FragMZ = Float.parseFloat(mz.trim());
                        fragment.intensity = Float.parseFloat(intensity.trim());
                        fragment.Charge = 1;
                        FragmentPeaks.add(fragment);
                    }
                } else if ("".equals(line) && Peak) {
                    Header = false;
                    Peak = false;
                    fraglib.AddFragments(FragmentPeaks);
                    PeptideFragmentLib.put(fraglib.GetKey(), fraglib);
                }
            }
            Logger.getRootLogger().info("No. of peptide ions in the imported library:"+PeptideFragmentLib.size());       
            GenerateDecoyLib();
        } catch (MatrixLoaderException ex) {
            Logger.getRootLogger().error(ExceptionUtils.getStackTrace(ex));
        }
    }

    public void FragmentSelection(ArrayList<LCMSID> LCMSIDList, float Freq, int TopNFrag) {
        fragselection = new FragmentSelection(LCMSIDList);
        fragselection.freqPercent = Freq;
        fragselection.GeneratePepFragScoreMap();
        fragselection.FillMissingFragScoreMap();
        fragselection.GenerateTopFragMap(TopNFrag);
    }
    
    public void ImportFragLibTopFrag(ArrayList<LCMSID> LCMSIDList, float Freq, int topNFrag) {
        FragmentSelection(LCMSIDList, Freq, topNFrag);
        for (LCMSID lcmsid : LCMSIDList) {
            for (PepIonID pepIonID : lcmsid.GetPepIonList().values()) {
                if (!PeptideFragmentLib.containsKey(pepIonID.GetKey())) {
                    PepFragmentLib fraglib = new PepFragmentLib();
                    fraglib.Sequence = pepIonID.Sequence;
                    fraglib.ModificationString = pepIonID.GetModificationString();
                    fraglib.Charge = pepIonID.Charge;
                    fraglib.ModSequence = pepIonID.ModSequence;
                    fraglib.PrecursorMz = pepIonID.NeutralPrecursorMz();
                    fraglib.MS1Score = pepIonID.PeakClusterScore;
                    fraglib.RetentionTime.add(pepIonID.PeakRT);
                    if (pepIonID.MaxProbability > fraglib.MaxProbability) {
                        fraglib.MaxProbability = pepIonID.MaxProbability;
                    }
                    if (pepIonID.PeakClusterScore > fraglib.MS1Score) {
                        fraglib.MS1Score = pepIonID.PeakClusterScore;
                    }
                    PeptideFragmentLib.put(pepIonID.GetKey(), fraglib);
                }
                                
                if (pepIonID.FragmentPeaks != null && !pepIonID.FragmentPeaks.isEmpty()) {                                        
                    //PeptideFragmentLib.get(pepIonID.GetKey()).AddFragments(pepIonID.FragmentPeaks);
                    ArrayList<FragmentPeak> frags=new ArrayList<>();
                    for(FragmentPeak fra : pepIonID.FragmentPeaks){
                        if(fragselection.TopFrags.get(pepIonID.GetKey()).contains(fra.GetFragKey())){
                            frags.add(fra);
                        }
                    }
                    if(!frags.isEmpty()){
                        PeptideFragmentLib.get(pepIonID.GetKey()).AddFragments(frags);
                    }
                    else{
                        Logger.getRootLogger().warn("Skipped peptide ion: " + pepIonID.GetKey() + " because it does not have enough matched fragments from file: " + lcmsid.mzXMLFileName);
                    }                    
                } else {
                    Logger.getRootLogger().warn("Skipped peptide ion: " + pepIonID.GetKey() + " because it does not have any matched fragment from file: " + lcmsid.mzXMLFileName);
                }
            }
        }
        try {
            GenerateDecoyLib();
        } catch (MatrixLoaderException ex) {
            Logger.getRootLogger().error(ExceptionUtils.getStackTrace(ex));
        }
    }
    
    public void ImportFragLib(ArrayList<LCMSID> LCMSIDList) {
        for (LCMSID lcmsid : LCMSIDList) {
            for (PepIonID pepIonID : lcmsid.GetPepIonList().values()) {
                if (!PeptideFragmentLib.containsKey(pepIonID.GetKey())) {
                    PepFragmentLib fraglib = new PepFragmentLib();
                    fraglib.Sequence = pepIonID.Sequence;
                    fraglib.ModificationString = pepIonID.GetModificationString();
                    fraglib.Charge = pepIonID.Charge;
                    fraglib.ModSequence = pepIonID.ModSequence;
                    fraglib.PrecursorMz = pepIonID.NeutralPrecursorMz();
                    fraglib.MS1Score = pepIonID.PeakClusterScore;
                    fraglib.RetentionTime.add(pepIonID.PeakRT);
                    if (pepIonID.MaxProbability > fraglib.MaxProbability) {
                        fraglib.MaxProbability = pepIonID.MaxProbability;
                    }
                    if (pepIonID.PeakClusterScore > fraglib.MS1Score) {
                        fraglib.MS1Score = pepIonID.PeakClusterScore;
                    }
                    PeptideFragmentLib.put(pepIonID.GetKey(), fraglib);
                }
                                
                if (pepIonID.FragmentPeaks != null && !pepIonID.FragmentPeaks.isEmpty()) {                                        
                    PeptideFragmentLib.get(pepIonID.GetKey()).AddFragments(pepIonID.FragmentPeaks);                    
                } else {
                    Logger.getRootLogger().warn("Skipped peptide ion: " + pepIonID.GetKey() + " because it does not have any matched fragment from file: " + lcmsid.mzXMLFileName);
                }
            }
        }
        try {
            GenerateDecoyLib();
        } catch (MatrixLoaderException ex) {
            Logger.getRootLogger().error(ExceptionUtils.getStackTrace(ex));
        }
    }
}
