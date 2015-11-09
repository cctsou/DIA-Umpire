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
package MSUmpire.LCMSBaseStructure;

import MSUmpire.BaseDataStructure.InstrumentParameter;
import MSUmpire.BaseDataStructure.SpectralDataType;
import MSUmpire.PSMDataStructure.PSM;
import MSUmpire.PSMDataStructure.PepIonID;
import MSUmpire.PeakDataStructure.PeakCluster;
import MSUmpire.PeakDataStructure.PeakCurve;
import MSUmpire.PeakDataStructure.SortedClusterCollectionClassApexRT;
import MSUmpire.PeakDataStructure.SortedCurveCollectionApexRT;
import MSUmpire.SpectrumParser.mzXMLParser;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;
import java.sql.SQLException;
import java.util.*;
import java.util.logging.Level;
import net.sf.javaml.core.kdtree.KDTree;
import net.sf.javaml.core.kdtree.KeyDuplicateException;
import net.sf.javaml.core.kdtree.KeySizeException;
import org.apache.commons.io.FilenameUtils;
import org.apache.commons.lang.exception.ExceptionUtils;
import org.apache.commons.math3.analysis.polynomials.PolynomialSplineFunction;
import org.apache.log4j.Logger;
import org.nustaq.serialization.FSTObjectInput;
import org.nustaq.serialization.FSTObjectOutput;


/**
 *
 * @author Chih-Chiang Tsou <chihchiang.tsou@gmail.com>
 */
public class LCMSPeakBase {

    public ArrayList<PeakCluster> PeakClusters = new ArrayList<>(1000);    
    //public SortedClusterCollectionClassMZ MZSortedClusters;
    private SortedClusterCollectionClassApexRT PeakClusterListRT=null;
    private KDTree PeakClusterMassSearchTree=null;
    //private KDTree PeakClusterMZSearchTree=null;
    private KDTree PeakCurveSearchTree=null;
    public String ScanCollectionName;
    public String ParentmzXMLName;
    protected mzXMLParser mzxml;
    public int MaxNoPeakCluster;
    public int MinNoPeakCluster;
    public int StartCharge;
    public int EndCharge;
    public float MiniIntensity;
    public float SNR;
    //public ConnectionManager connectionManager;
    //public SortedCurveCollectionMZ PeakCurveListMZ = new SortedCurveCollectionMZ();
    private SortedCurveCollectionApexRT PeakCurveListRT = null;
    
    //KDTree PeakClusterSearchTree =null;
    public ArrayList<PeakCurve> UnSortedPeakCurves;    
    //public ArrayList<PeakCurve> IsolatedMS1PeakCurves;
    public InstrumentParameter parameter;
    public SpectralDataType.DataType datattype;
    public boolean Resume = true;
    public boolean ExportPeakClusterTable=true;
    public boolean SaveSerializationFile=true;
    public boolean ExportPeakCurveTable=false;
    int NoCPUs = Runtime.getRuntime().availableProcessors() - 2;
    public PolynomialSplineFunction Masscalibrationfunction;

    public void ClearAllPeaks(){
        BaseClearAllPeaks();
    }
    public void BaseClearAllPeaks() {        
        PeakClusters = null;
        PeakClusterListRT = null;
        PeakCurveListRT = null;
        UnSortedPeakCurves=null;
        //IsolatedMS1PeakCurves = null;
        PeakCurveSearchTree=null;
        PeakClusterMassSearchTree=null;
        //PeakClusterMZSearchTree=null;
    }

    public void ClearRawPeaks() {
        for (PeakCurve peakCurve : UnSortedPeakCurves) {
            peakCurve.GetPeakList().clear();
        }
    }
     
//    public void SetMySQLConnection(ConnectionManager connectionManager){
//        this.connectionManager = connectionManager;
//    }
//    public void GenerateIsolatedPeakCurve() {
//        ArrayList<Integer> IncludedIndex = new ArrayList<>();
//        IsolatedMS1PeakCurves = new ArrayList<>();
//        for (PeakCluster cluster : PeakClusters) {
//            for (PeakCurve peak : cluster.IsoPeaksCurves) {
//                if (peak != null) {
//                    if (!IncludedIndex.contains(peak.Index)) {
//                        IncludedIndex.add(peak.Index);
//                    }
//                }
//            }
//        }
//
//        for (PeakCurve peak : UnSortedPeakCurves) {
//            if (!IncludedIndex.contains(peak.Index)) {
//                IsolatedMS1PeakCurves.add(peak);
//            }
//        }
//        Logger.getRootLogger().info("No. of isolated peak curve: " + IsolatedMS1PeakCurves.size());
//    }

//    public void ReadPeakCurveFromDB(boolean ReadPeak) throws SQLException {
//        Connection connection = connectionManager.GetConnection();
//        Logger.getRootLogger().info("Loading peak curve table....");
//        Statement state = connection.createStatement();
//        ResultSet rsCurve = state.executeQuery("SELECT * FROM " + FilenameUtils.getBaseName(ScanCollectionName) + "_PeakCurve Order by Curve_index");
//
//        UnSortedPeakCurves = new ArrayList<>();
//        while (rsCurve.next()) {
//            PeakCurve peakCurve = new PeakCurve(parameter);
//            peakCurve.TargetMz = rsCurve.getFloat("mz");
//
//            peakCurve.Index = rsCurve.getInt("Curve_index");
//            peakCurve.SetRTs(rsCurve.getFloat("StartRT"), rsCurve.getFloat("EndRT"));
//            peakCurve.ApexRT = rsCurve.getFloat("ApexRT");
//            peakCurve.ApexInt = rsCurve.getFloat("ApexInt");
//
//            try {
//                peakCurve.MzVar = rsCurve.getFloat("mzVar");
//            } catch (Exception e) {
//            }
//            try {
//                String Ridges = rsCurve.getString("NoRidges");
//                for (String ridges : Ridges.split("#")) {
//                    if (!"".equals(ridges)) {
//                        peakCurve.RegionRidge.add(Float.parseFloat(ridges));
//                    }
//                }
//
//            } catch (Exception e) {
//            }
//            //peakCurve.SetRTs(rsCurve.getFloat("StartRT"),rsCurve.getFloat("EndRT"));
//            if (ReadPeak) {
//                peakCurve.ReadPeakResultMySQL(rsCurve);
//            }
//            UnSortedPeakCurves.add(peakCurve);            
//        }
////        PeakCurveListMZ.addAll(peakCurves);
////        PeakCurveListRT.addAll(peakCurves);
////        PeakCurveListMZ.Finalize();
////        PeakCurveListRT.Finalize();
//                
//        state.close();
//        connectionManager.CloseConnection();
//    }

//    public void GenerateRTSoretedClusterList(boolean reindex) {
//        ApexRTSortedClusters = new SortedClusterCollectionClassApexRT();
//        ApexRTSortedClusters.addAll(PeakClusters);
//
//        if (reindex) {
//            for (int i = 0; i < ApexRTSortedClusters.size(); i++) {
//                ApexRTSortedClusters.get(i).Index = i + 1;
//            }
//        }
//    }

    public ArrayList<PeakCluster> FindAllPeakClustersForPepByPSM(PepIonID pep) {
        ArrayList<PeakCluster> allclusterList = new ArrayList<>();
        for (PSM psm : pep.GetPSMList()) {
            ArrayList<PeakCluster> clusterList = FindPeakClustersByMassRTTol(psm.ObserPrecursorMass, psm.Charge, psm.RetentionTime, parameter.MaxCurveRTRange / 2);
            for (PeakCluster cluster : clusterList) {
                if (!allclusterList.contains(cluster)) {
                    allclusterList.add(cluster);
                }
            }
        }
        return allclusterList;
    }

    private float GetMassError(float RT) {

        if (RT > Masscalibrationfunction.getKnots()[Masscalibrationfunction.getN()]) {
            RT = (float) Masscalibrationfunction.getKnots()[Masscalibrationfunction.getN()];
        }
        if (RT < Masscalibrationfunction.getKnots()[Masscalibrationfunction.getN()]) {
            RT = (float) Masscalibrationfunction.getKnots()[0];
        }
        return (float) Masscalibrationfunction.value(RT);
    }

    public ArrayList<PeakCluster> FindAllPeakClustersForMappedPep(PepIonID pep) {
        ArrayList<PeakCluster> allclusterList = new ArrayList<>();        
        float idrt=pep.GetRT();
        if (idrt != -1) {
            float calibratedmass = InstrumentParameter.GetMzByPPM(pep.CalcNeutralPepMass(), 1, -GetMassError(idrt));
            ArrayList<PeakCluster> clusterList = FindPeakClustersByMassIDTime(calibratedmass, pep.Charge, idrt, parameter.RTtol);
            return clusterList;
        }
        for (float rt : pep.PredictRT) {
            float calibratedmass = InstrumentParameter.GetMzByPPM(pep.CalcNeutralPepMass(), 1, -GetMassError(rt));
            float rtrange=parameter.RT_window_Targeted;            
            if (rtrange == -1) {
                rtrange = Math.min(5, Math.max(pep.GetRTSD() * 2, parameter.RTtol));
            }
            ArrayList<PeakCluster> clusterList = FindPeakClustersByMassRTTol(calibratedmass, pep.Charge, rt, rtrange);
            for (PeakCluster cluster : clusterList) {
                if (!cluster.Identified && !allclusterList.contains(cluster)) {
                    allclusterList.add(cluster);
                }
            }
        }
        return allclusterList;
    }

    public ArrayList<PeakCurve> FindPeakCurveByMZRT(float mz, float RT, float ppm) {
        ArrayList<PeakCurve> ReturnList = new ArrayList<>();
        float lowrt = RT - parameter.MaxCurveRTRange;
        float highrt = RT + parameter.MaxCurveRTRange;
        float lowmz = InstrumentParameter.GetMzByPPM(mz, 1, ppm);
        float highmz = InstrumentParameter.GetMzByPPM(mz, 1, -ppm);
       
        Object[] found=null;
        try {
            found = GetPeakCurveSearchTree().range(new double[]{lowrt,lowmz}, new double[]{highrt,highmz});
        } catch (KeySizeException ex) {
            
        }
        if(found==null || found.length==0){
            return ReturnList;
        }
       for(Object obj : found){         
            PeakCurve peakcurve = (PeakCurve)obj;
            if (InstrumentParameter.CalcPPM(peakcurve.TargetMz, mz) < ppm) {
                if (peakcurve.StartRT() <= RT && peakcurve.EndRT() >= RT) {
                    ReturnList.add(peakcurve);
                }
            } else {
                if (peakcurve.TargetMz > mz) {
                    return ReturnList;
                }
            }
        }
        return ReturnList;
    }
    
    public ArrayList<PeakCluster> FindPeakClustersByMassIDTime(float mass, int charge, float RT, float RTtol) {
        ArrayList<PeakCluster> ReturnList = new ArrayList<>();
        float lowrt = RT - RTtol;
        float highrt = RT + RTtol;
        float lowmass = InstrumentParameter.GetMzByPPM(mass, 1, parameter.MS1PPM);
        float highmass = InstrumentParameter.GetMzByPPM(mass, 1, -parameter.MS1PPM);
       
        Object[] found=null;
        try {
            found = GetPeakClusterMassSearchTree().range(new double[]{lowrt,lowmass}, new double[]{highrt,highmass});
        } catch (KeySizeException ex) {
            
        }
        if(found==null || found.length==0){
            return ReturnList;
        }
        for (Object obj : found) {
            PeakCluster candidateCluster = (PeakCluster) obj;            
            if (candidateCluster.Charge == charge && Math.abs(candidateCluster.PeakHeightRT[0] - RT) <= RTtol) {                
                ReturnList.add(candidateCluster);
            }
        }
        return ReturnList;
    }

    public ArrayList<PeakCluster> FindPeakClustersByMassRTTol(float mass, int charge, float RT, float RTtol) {
        
        ArrayList<PeakCluster> Clusters = new ArrayList<>();
        ArrayList<PeakCluster> TightRangeClusters = new ArrayList<>();
        
        float lowrt = RT - RTtol;
        float highrt = RT + RTtol;
        float lowmass = InstrumentParameter.GetMzByPPM(mass, 1, parameter.MS1PPM);
        float highmass = InstrumentParameter.GetMzByPPM(mass, 1, -parameter.MS1PPM);
       
        Object[] found=null;
        try {
            found = GetPeakClusterMassSearchTree().range(new double[]{lowrt,lowmass}, new double[]{highrt,highmass});
        } catch (KeySizeException ex) {
            
        }
        if(found==null || found.length==0){
            return Clusters;
        }
        for (Object obj : found) {
            PeakCluster candidateCluster = (PeakCluster) obj;
            if (candidateCluster.Charge == charge && candidateCluster.startRT - RTtol <= RT && candidateCluster.endRT + RTtol >= RT) {
                Clusters.add(candidateCluster);
            }
            if (candidateCluster.Charge == charge && candidateCluster.startRT <= RT && candidateCluster.endRT >= RT) {
                TightRangeClusters.add(candidateCluster);
            }
        }
        if (!TightRangeClusters.isEmpty()) {
            return TightRangeClusters;
        }
        return Clusters;
    }

    public ArrayList<PeakCluster> FindPeakClustersByMassRTRange(float mass, int charge, float startrt, float endrt) {
        ArrayList<PeakCluster> Clusters = new ArrayList<>();
        float lowrt = startrt;
        float highrt = endrt;
        float lowmass = InstrumentParameter.GetMzByPPM(mass, 1, parameter.MS1PPM);
        float highmass = InstrumentParameter.GetMzByPPM(mass, 1, -parameter.MS1PPM);

        Object[] found = null;
        try {
            found = GetPeakClusterMassSearchTree().range(new double[]{lowrt, lowmass}, new double[]{highrt, highmass});
        } catch (KeySizeException ex) {

        }
        if (found == null || found.length == 0) {
            return Clusters;
        }
        for (Object obj : found) {
            PeakCluster candidateCluster = (PeakCluster) obj;

            if (candidateCluster.Charge == charge && candidateCluster.PeakHeightRT[0] >= startrt && candidateCluster.PeakHeightRT[0] <= endrt) {
                Clusters.add(candidateCluster);
            }

        }
        return Clusters;
    }


//    public void GenerateMZSortedClusterList(boolean reindex) {
//        MZSortedClusters = new SortedClusterCollectionClassMZ();
//        MZSortedClusters.addAll(PeakClusters);
//        if (reindex) {
//            for (int i = 0; i < MZSortedClusters.size(); i++) {
//                MZSortedClusters.get(i).Index = i + 1;
//            }
//        }
//    }

    public void ExportPeakClusterResultCSV() throws IOException {

        Logger.getRootLogger().info("Writing PeakCluster CSV:" + FilenameUtils.getFullPath(ScanCollectionName) + FilenameUtils.getBaseName(ScanCollectionName) + "_PeakCluster.csv...");
        FileWriter writer = new FileWriter(FilenameUtils.getFullPath(ScanCollectionName) + FilenameUtils.getBaseName(ScanCollectionName) + "_PeakCluster.csv");

        String mzstring = "";
        String Idxstring = "";
        String CorrString = "";
        String SNRString = "";
        String PeakheightString = "";
        String PeakheightRTString = "";
        String PeakAreaString = "";
        String IdentifiedString = "0";
        //String IDIsoPatternString="";
        //String MapIsoPatternString="";

        for (int i = 0; i < MaxNoPeakCluster; i++) {
            mzstring += ",mz" + (i + 1);
            Idxstring += ",PeakIdx" + (i + 1);
            if (i > 0) {
                CorrString += ",Corr" + (i + 1);
                //IDIsoPatternString += ",IsoPatternID" + (i + 1);
                //MapIsoPatternString += ",IsoPatternMap" + (i + 1);
            }
            SNRString += ",SNR" + (i + 1);
            PeakheightString += ",PeakHeight" + (i + 1);
            PeakheightRTString += ",PeakHeightRT" + (i + 1);
            PeakAreaString += ",PeakArea" + (i + 1);
        }

        writer.write("Cluster_Index,StartRT,EndRT,StartScan,EndScan,Identified,Charge" + mzstring + Idxstring + CorrString + SNRString + PeakheightString + PeakheightRTString + PeakAreaString /*+ IDIsoPatternString+MapIsoPatternString*/ + ",IDIsoPatternProb,IsoMapProb,ConflictCorr,LeftInt,RightInt,NoRidges,MS1Score,MS1Prob,MS1LProb\n");

        for (PeakCluster cluster : PeakClusters) {
            IdentifiedString = "0";
            if (cluster.Identified) {
                IdentifiedString = "1";
            }

            String statementString = cluster.Index + "," + cluster.startRT + "," + cluster.endRT + "," + cluster.StartScan + "," + cluster.EndScan + "," + IdentifiedString + "," + cluster.Charge + ",";

            mzstring = "";
            Idxstring = "";
            CorrString = "";
            SNRString = "";
            PeakheightString = "";
            PeakheightRTString = "";
            PeakAreaString = "";
            //IDIsoPatternString="";
            //MapIsoPatternString="";

            for (int i = 0; i < MaxNoPeakCluster; i++) {
                mzstring += cluster.mz[i] + ",";
                Idxstring += cluster.IsoPeakIndex[i] + ",";
                if (i > 0) {
                    CorrString += cluster.Corrs[i - 1] + ",";
                    //IDIsoPatternString += ",";
                    //MapIsoPatternString += ",";
                }
                SNRString += cluster.GetSNR(i) + ",";
                PeakheightString += cluster.PeakHeight[i] + ",";
                PeakheightRTString += cluster.PeakHeightRT[i] + ",";
                PeakAreaString += cluster.PeakArea[i] + ",";
            }
            statementString += mzstring + Idxstring + CorrString + SNRString + PeakheightString + PeakheightRTString + PeakAreaString /*+ IDIsoPatternString+MapIsoPatternString*/ + cluster.IDIsoPatternProb + "," + cluster.IsoMapProb + "," + cluster.GetConflictCorr() + "," + cluster.LeftInt + "," + cluster.RightInt + "," + cluster.NoRidges + "," + cluster.MS1Score + "," + cluster.MS1ScoreProbability + "," + cluster.MS1ScoreLocalProb + "\n";
            writer.write(statementString);
        }
        writer.close();
        //System.out.print("Finished multithreading\n");
    }

    
    public void ExportPeakCluster() throws IOException {
        WritePeakClusterSerialization();      
//        if (ExportPeakClusterTable) {
//            ExportPeakClusterResultCSV();
//        }
//        if (connectionManager != null) {
//            ExportPeakClusterResultCSV();
//            CreatePeakClusterTable();
//            ExportTableDBBulkLoader("PeakCluster",false);
//        }
    }

     public void CreatePeakFolder(){
        new File(FilenameUtils.getFullPath(ParentmzXMLName)+FilenameUtils.getBaseName(ParentmzXMLName) + "_Peak/").mkdir();        
    }
    
    public void WritePeakClusterSerialization() {
        //JavaSerializationPeakClusterWrite();
        FS_PeakClusterWrite();
    }

    private void FS_PeakClusterWrite() {
        try {
            Logger.getRootLogger().info("Writing PeakCluster serialization to file:" +  FilenameUtils.getBaseName(ScanCollectionName) + "_PeakCluster.serFS...");
            FileOutputStream fout = new FileOutputStream(FilenameUtils.getFullPath(ParentmzXMLName)+ FilenameUtils.getBaseName(ParentmzXMLName)+"_Peak/" + FilenameUtils.getBaseName(ScanCollectionName) + "_PeakCluster.serFS", false);
            FSTObjectOutput out = new FSTObjectOutput(fout);
            out.writeObject(PeakClusters);
            out.close();
            fout.close();
        } catch (Exception ex) {
            Logger.getRootLogger().error(ExceptionUtils.getStackTrace(ex));
            JavaSerializationPeakClusterWrite();
        }
    }
    
    private void JavaSerializationPeakClusterWrite() {
        try {
            Logger.getRootLogger().info("Writing PeakCluster serialization to file:" +  FilenameUtils.getBaseName(ScanCollectionName) + "_PeakCluster.ser...");
            FileOutputStream fout = new FileOutputStream(FilenameUtils.getFullPath(ParentmzXMLName)+ FilenameUtils.getBaseName(ParentmzXMLName)+"_Peak/" + FilenameUtils.getBaseName(ScanCollectionName) + "_PeakCluster.ser", false);
            ObjectOutputStream oos = new ObjectOutputStream(fout);
            oos.writeObject(PeakClusters);
            oos.close();
            fout.close();
        } catch (Exception ex) {
            Logger.getRootLogger().error(ExceptionUtils.getStackTrace(ex));
        }
    }

     public boolean ReadPeakCluster() {
        return ReadPeakClusterSerialization();
    }
    private boolean ReadPeakClusterSerialization() {
        //return JavaSerializationPeakClusterRead();
        if(!FS_PeakClusterRead()){
            if (JavaSerializationPeakClusterRead()) {
                //FS_PeakClusterWrite();
                return true;
            }
            return false;
        }
        for(PeakCluster cluster : PeakClusters){
            cluster.CreateLock();
        }
        return true;        
    }

    private boolean JavaSerializationPeakClusterRead() {
        if (!new File(FilenameUtils.getFullPath(ParentmzXMLName)+ FilenameUtils.getBaseName(ParentmzXMLName)+"_Peak/" + FilenameUtils.getBaseName(ScanCollectionName) + "_PeakCluster.ser").exists()) {
            return false;
        }
        try {
            Logger.getRootLogger().info("Reading PeakCluster serialization from file:" + FilenameUtils.getBaseName(ScanCollectionName) + "_PeakCluster.ser...");

            FileInputStream fileIn = new FileInputStream(FilenameUtils.getFullPath(ParentmzXMLName)+ FilenameUtils.getBaseName(ParentmzXMLName)+"_Peak/" + FilenameUtils.getBaseName(ScanCollectionName) + "_PeakCluster.ser");
            ObjectInputStream in = new ObjectInputStream(fileIn);
            PeakClusters = (ArrayList<PeakCluster>) in.readObject();
            in.close();
            fileIn.close();
        } catch (Exception ex) {
            Logger.getRootLogger().error(ExceptionUtils.getStackTrace(ex));
            return false;
        }
        return true;
    }

    private boolean FS_PeakClusterRead() {
        if (!new File(FilenameUtils.getFullPath(ParentmzXMLName)+ FilenameUtils.getBaseName(ParentmzXMLName)+"_Peak/" + FilenameUtils.getBaseName(ScanCollectionName) + "_PeakCluster.serFS").exists()) {
            return false;
        }
        try {
            Logger.getRootLogger().info("Reading PeakCluster serialization from file:" + FilenameUtils.getBaseName(ScanCollectionName) + "_PeakCluster.serFS...");

            FileInputStream fileIn = new FileInputStream(FilenameUtils.getFullPath(ParentmzXMLName)+ FilenameUtils.getBaseName(ParentmzXMLName)+"_Peak/" + FilenameUtils.getBaseName(ScanCollectionName)+ "_PeakCluster.serFS");
            FSTObjectInput in = new FSTObjectInput(fileIn);
            PeakClusters = (ArrayList<PeakCluster>) in.readObject();
            in.close();
            fileIn.close();            
        } catch (Exception ex) {            
            Logger.getRootLogger().error(ExceptionUtils.getStackTrace(ex));
            if(FS_PeakClusterRead_Old()){
                WritePeakClusterSerialization();                
                return true;
            }
            Logger.getRootLogger().error(ExceptionUtils.getStackTrace(ex));
            return false;
        } 
        return true;
    }
    
    private boolean FS_PeakClusterRead_Old() {
        if (!new File(FilenameUtils.getFullPath(ParentmzXMLName)+ FilenameUtils.getBaseName(ParentmzXMLName)+"_Peak/" + FilenameUtils.getBaseName(ScanCollectionName) + "_PeakCluster.serFS").exists()) {
            return false;
        }
        try {
            Logger.getRootLogger().info("Old PeakCluster serialization from file:" + FilenameUtils.getBaseName(ScanCollectionName) + "_PeakCluster.serFS...");

            FileInputStream fileIn = new FileInputStream(FilenameUtils.getFullPath(ParentmzXMLName)+ FilenameUtils.getBaseName(ParentmzXMLName)+"_Peak/" + FilenameUtils.getBaseName(ScanCollectionName) + "_PeakCluster.serFS");
            org.nustaq_old.serialization.FSTObjectInput in = new org.nustaq_old.serialization.FSTObjectInput(fileIn);
            PeakClusters = (ArrayList<PeakCluster>) in.readObject();
            in.close();
            fileIn.close();            
        } catch (Exception ex) {
            Logger.getRootLogger().error("Old version reader still failed.");
            Logger.getRootLogger().error(ExceptionUtils.getStackTrace(ex));
            return false;
        } 
        return true;
    }
    
    private void WritePeakCurveSerialization() {
        //JavaSerializationPeakCurveWrite();
        FSPeakCurveWrite();
    }

     private void FSPeakCurveWrite() {
        try {
            Logger.getRootLogger().info("Writing PeakCurve serialization to file:" + FilenameUtils.getBaseName(ScanCollectionName) + "_PeakCurve.serFS...");
            FileOutputStream fout = new FileOutputStream(FilenameUtils.getFullPath(ParentmzXMLName)+ FilenameUtils.getBaseName(ParentmzXMLName)+"_Peak/" + FilenameUtils.getBaseName(ScanCollectionName) + "_PeakCurve.serFS", false);
            FSTObjectOutput out = new FSTObjectOutput(fout);
            out.writeObject(UnSortedPeakCurves);
            out.close();            
            fout.close();
        } catch (Exception ex) {
            Logger.getRootLogger().error(ExceptionUtils.getStackTrace(ex));
        }
    }
    
    private void JavaSerializationPeakCurveWrite() {
        try {
            Logger.getRootLogger().info("Writing PeakCurve serialization to file:" + FilenameUtils.getBaseName(ScanCollectionName) + "_PeakCurve.ser...");
            FileOutputStream fout = new FileOutputStream(FilenameUtils.getFullPath(ParentmzXMLName)+ FilenameUtils.getBaseName(ParentmzXMLName)+"_Peak/" + FilenameUtils.getBaseName(ScanCollectionName) + "_PeakCurve.ser", false);
            ObjectOutputStream oos = new ObjectOutputStream(fout);
            oos.writeObject(UnSortedPeakCurves);
            oos.close();
            fout.close();
        } catch (Exception ex) {
            Logger.getRootLogger().error(ExceptionUtils.getStackTrace(ex));
        }
    }

    private boolean ReadPeakCurveSerialization() {
        if (!FSPeakCurveRead()) {
            if (JavaSerializationPeakCurveRead()) {
                FSPeakCurveWrite();
                return true;
            }
            return false;
        }
        return true;
    }

    private boolean FSPeakCurveRead() {
        if (!new File(FilenameUtils.getFullPath(ParentmzXMLName)+ FilenameUtils.getBaseName(ParentmzXMLName)+"_Peak/" + FilenameUtils.getBaseName(ScanCollectionName) + "_PeakCurve.serFS").exists()) {
            return false;
        }
        try {
            Logger.getRootLogger().info("Reading PeakCurve serialization from file:" + FilenameUtils.getBaseName(ScanCollectionName) + "_PeakCurve.serFS...");
            FileInputStream fileIn = new FileInputStream(FilenameUtils.getFullPath(ParentmzXMLName)+ FilenameUtils.getBaseName(ParentmzXMLName)+"_Peak/" + FilenameUtils.getBaseName(ScanCollectionName) + "_PeakCurve.serFS");
            FSTObjectInput in = new FSTObjectInput(fileIn);
            UnSortedPeakCurves = (ArrayList<PeakCurve>) in.readObject();
            in.close();
            fileIn.close();         
        } catch (Exception ex) {
            Logger.getRootLogger().error(ExceptionUtils.getStackTrace(ex));
            return false;
        }
        return true;
    }
    
    private boolean JavaSerializationPeakCurveRead() {
        if (!new File(FilenameUtils.getFullPath(ParentmzXMLName)+ FilenameUtils.getBaseName(ParentmzXMLName)+"_Peak/" + FilenameUtils.getBaseName(ScanCollectionName) + "_PeakCurve.ser").exists()) {
            return false;
        }
        try {
            Logger.getRootLogger().info("Reading PeakCurve serialization from file:" + FilenameUtils.getBaseName(ScanCollectionName) + "_PeakCurve.ser...");
            FileInputStream fileIn = new FileInputStream(FilenameUtils.getFullPath(ParentmzXMLName)+ FilenameUtils.getBaseName(ParentmzXMLName)+"_Peak/" + FilenameUtils.getBaseName(ScanCollectionName) + "_PeakCurve.ser");
            ObjectInputStream in = new ObjectInputStream(fileIn);
            UnSortedPeakCurves = (ArrayList<PeakCurve>) in.readObject();
            in.close();
            fileIn.close();
        } catch (Exception ex) {
            Logger.getRootLogger().error(ExceptionUtils.getStackTrace(ex));
            return false;
        } 
        return true;
    }

    public void ExportPeakCurveResult() throws IOException, SQLException {
        WritePeakCurveSerialization();
//        if (ExportPeakCurveTable) {
//            ExportPeakCurveResultCSV_V2();
////            if (connectionManager != null) {
////                CreatePeakCurveTable_V2();
////                ExportTableDBBulkLoader("PeakCurve", true);
////            }
//        }
    }
//<editor-fold defaultstate="collapsed" desc="DB process">
    
//    protected void ExportTableDBBulkLoader(String TableSuffix, boolean delete) throws SQLException {
//        if (connectionManager == null) {
//            return;
//        }
//        Logger.getRootLogger().info("Writing " + TableSuffix + " result to MySQL DB:" + FilenameUtils.getBaseName(ScanCollectionName) + "_" + TableSuffix + "...");
//        Connection connection = connectionManager.GetConnection();
//        Statement state = connection.createStatement();
//        state.executeUpdate("LOAD DATA LOCAL INFILE '" + FilenameUtils.separatorsToUnix(FilenameUtils.getFullPath(ScanCollectionName)) + FilenameUtils.getBaseName(ScanCollectionName) + "_" + TableSuffix + ".csv'" + " INTO TABLE " + FilenameUtils.getBaseName(ScanCollectionName) + "_" + TableSuffix + " FIELDS TERMINATED BY ','" + " LINES TERMINATED BY '\\n' IGNORE 1 LINES");
//        state.close();
//        connectionManager.CloseConnection();
//        if (delete) {
//            File file = new File(FilenameUtils.separatorsToUnix(FilenameUtils.getFullPath(ScanCollectionName)) + FilenameUtils.getBaseName(ScanCollectionName) + "_" + TableSuffix + ".csv");
//            file.delete();
//        }
//    }
    
//    protected void CreatePeakClusterTable() throws SQLException {
//        if (connectionManager == null) {
//            return;
//        }
//        Connection connection = connectionManager.GetConnection();
//        connection.createStatement().execute("DROP TABLE IF EXISTS " + FilenameUtils.getBaseName(ScanCollectionName) + "_PeakCluster;");
//        String statementString = "CREATE TABLE " + FilenameUtils.getBaseName(ScanCollectionName) + "_PeakCluster (Cluster_Index INT NOT NULL, StartRT FLOAT NOT NULL, EndRT FLOAT NOT NULL, Identified TINYINT, Charge INT";
//
//        String mzstring = "";
//        String Idxstring = "";
//        String CorrString = "";
//        String SNRString = "";
//        String PeakheightString = "";
//        String PeakheightRTString = "";
//        String PeakAreaString = "";
//        String IdentifiedString = "0";
//        //String IDIsoPatternString = "";
//        //String MapIsoPatternString = "";
//
//        for (int i = 0; i < MaxNoPeakCluster; i++) {
//
//            mzstring += ", mz" + (i + 1) + " DOUBLE";
//            Idxstring += ", PeakIdx" + (i + 1) + " INT";
//            if (i > 0) {
//                CorrString += ", Corr" + (i + 1) + " DOUBLE";
//                //IDIsoPatternString += ", IsoPatternID" + (i + 1) + " DOUBLE";
//                //MapIsoPatternString += ", IsoPatternMap" + (i + 1) + " DOUBLE";
//            }
//            SNRString += ", SNR" + (i + 1) + " DOUBLE";
//            PeakheightString += ", PeakHeight" + (i + 1) + " DOUBLE";
//            PeakheightRTString += ", PeakHeightRT" + (i + 1) + " DOUBLE";
//            PeakAreaString += ", PeakArea" + (i + 1) + " DOUBLE";
//
//        }
//        statementString += mzstring + Idxstring + CorrString + SNRString + PeakheightString + PeakheightRTString + PeakAreaString /*+ IDIsoPatternString + MapIsoPatternString*/ + ", IDIsoPatternProb DOUBLE, IsoMapProb DOUBLE, ConflictCorr DOUBLE, LeftInt DOUBLE, RightInt DOUBLE, NoRidges INT, MS1Score DOUBLE, MS1Prob DOUBLE, MS1LProb DOUBLE, PRIMARY KEY (Cluster_Index));";
//
//        connection.createStatement().execute(statementString);
//        connectionManager.CloseConnection();
//    }
//
//    protected void CreatePeakCurveTable_V2() throws SQLException {
//        if (connectionManager == null) {
//            return;
//        }
//        Connection connection = connectionManager.GetConnection();
//        connection.createStatement().execute("DROP TABLE IF EXISTS " + FilenameUtils.getBaseName(ScanCollectionName) + "_PeakCurve;");
//        connection.createStatement().execute("CREATE TABLE " + FilenameUtils.getBaseName(ScanCollectionName) + "_PeakCurve (Curve_index INT NOT NULL, StartRT DOUBLE NOT NULL, EndRT DOUBLE NOT NULL, mz DOUBLE NOT NULL,mzVar DOUBLE NOT NULL, ApexRT DOUBLE NOT NULL, ApexInt DOUBLE NOT NULL, Ridges TEXT NOT NULL, Smoothed_Peak LONGTEXT NOT NULL, PRIMARY KEY (Curve_index));");
//        connectionManager.CloseConnection();
//    }
//
//    protected void CreatePeakCurveTable() throws SQLException {
//        if (connectionManager == null) {
//            return;
//        }
//        Connection connection = connectionManager.GetConnection();
//        connection.createStatement().execute("DROP TABLE IF EXISTS " + FilenameUtils.getBaseName(ScanCollectionName) + "_PeakCurve;");
//        connection.createStatement().execute("CREATE TABLE " + FilenameUtils.getBaseName(ScanCollectionName) + "_PeakCurve (Curve_index INT NOT NULL, StartRT DOUBLE NOT NULL, EndRT DOUBLE NOT NULL, StartScan INT NOT NULL, EndScan INT NOT NULL, mz DOUBLE NOT NULL, RAW_Peak TEXT NOT NULL, Smoothed_Peak LONGTEXT NOT NULL , PeakRegion TEXT, PSMs TEXT, PRIMARY KEY (Curve_index));");
//        connectionManager.CloseConnection();
//    }
//
//    protected void CreateCorrMatrixTable() throws SQLException {
//        if (connectionManager == null) {
//            return;
//        }
//        Connection connection = connectionManager.GetConnection();
//        connection.createStatement().execute("DROP TABLE IF EXISTS " + FilenameUtils.getBaseName(ScanCollectionName) + "_PeakCorr;");
//        connection.createStatement().execute("CREATE TABLE " + FilenameUtils.getBaseName(ScanCollectionName) + "_PeakCorr (ID INT NOT NULL AUTO_INCREMENT , IndexA INT NOT NULL, IndexB INT NOT NULL, AStart INT NOT NULL, Aend INT NOT NULL, Bstart INT NOT NULL, Bend INT NOT NULL,Corr DOUBLE NOT NULL, PRIMARY KEY (ID));");
//        connectionManager.CloseConnection();
//    }
//
//    protected void ExportCorrMatrix_V2(ArrayList<PeakCurveClusteringCorrV2Unit> ResultArrayList) throws IOException, SQLException {
//        ExportCorrMatrixCSV_V2(ResultArrayList);
//        CreateCorrMatrixTable();
//        ExportTableDBBulkLoader("PeakCorr",false);
//    }
//
//    protected void ExportCorrMatrixCSV_V2(ArrayList<PeakCurveClusteringCorrV2Unit> ResultArrayList) throws IOException {
//        Logger.getRootLogger().info("Writing PeakCorr result to file:" + ScanCollectionName + "_PeakCorr.csv...");
//        FileWriter writer = new FileWriter(FilenameUtils.getFullPath(ScanCollectionName) + FilenameUtils.getBaseName(ScanCollectionName) + "_PeakCorr.csv");
//        writer.write("ID,IndexA,IndexB,Corr\n");
//        int id = 0;
//        for (PeakCurveClusteringCorrV2Unit result : ResultArrayList) {
////            for (PeakOverlapRegion region : result.overlapRegions) {
////                writer.write((id++) + "," + region.PeakCurveIndexA + "," + region.PeakCurveIndexB + "," + region.Astart + "," + region.Aend + "," + region.Bstart + "," + region.Bend + "," + region.Correlation + "\n");
////            }
//        }
//        writer.close();
//    }
    
    //    public void ReadPeakClusterFromDB() throws SQLException {
//        Logger.getRootLogger().info("Read peak cluster table........"+FilenameUtils.getBaseName(ScanCollectionName) + "_PeakCluster");
//        Connection connection = connectionManager.GetConnection();
//        Statement state = connection.createStatement();
//        //UnIDpeakClusters = new HashMap<>();
//        PeakClusters = new ArrayList<>();
//        //RTSortedClusters=new SortedClusterCollectionClassRT();
//
//        ResultSet rsCluster = state.executeQuery("SELECT * FROM " + FilenameUtils.getBaseName(ScanCollectionName) + "_PeakCluster Order by Cluster_Index");
//        while (rsCluster.next()) {
//            int Charge = rsCluster.getInt("Charge");
//            PeakCluster cluster = new PeakCluster(MaxNoPeakCluster, Charge);
//            cluster.Index = rsCluster.getInt("Cluster_Index");
//            cluster.startRT = rsCluster.getFloat("StartRT");
//            cluster.endRT = rsCluster.getFloat("EndRT");
//            cluster.Identified = rsCluster.getInt("Identified") == 1;
//
//            for (int i = 0; i < MaxNoPeakCluster; i++) {
//                cluster.PeakHeight[i] = rsCluster.getFloat("PeakHeight" + (i + 1));
//                cluster.PeakArea[i] = rsCluster.getFloat("PeakArea" + (i + 1));
//                cluster.SetSNR(i, rsCluster.getFloat("SNR" + (i + 1)));
//                cluster.SetMz(i, rsCluster.getFloat("mz" + (i + 1)));
//                cluster.IsoPeakIndex[i] = rsCluster.getInt("PeakIdx" + (i + 1));
//                cluster.PeakHeightRT[i] = rsCluster.getFloat("PeakHeightRT" + (i + 1));
//                if (i > 0) {
//                    cluster.Corrs[i - 1] = rsCluster.getFloat("Corr" + (i + 1));
//                    try {
//                        //cluster.IsoPatternErrorID[i-1]=rsCluster.getFloat("IsoPatternID" + (i + 1));
//                        //cluster.IsoPatternErrorMap[i]=rsCluster.getFloat("IsoPatternMap" + (i + 1));
//                    } catch (Exception e) {
//                    }
//                }
//            }
//
//            cluster.IDIsoPatternProb = rsCluster.getFloat("IDIsoPatternProb");
//            cluster.IsoMapProb = rsCluster.getFloat("IsoMapProb");
//            cluster.SetConflictCorr(rsCluster.getFloat("ConflictCorr"));
//
//            try {
//                cluster.LeftInt = rsCluster.getFloat("LeftInt");
//                cluster.RightInt = rsCluster.getFloat("RightInt");
//                cluster.NoRidges = rsCluster.getInt("NoRidges");
//            } catch (Exception e) {
//            }
//
//            try {
//                cluster.MS1Score = rsCluster.getFloat("LDA");
//            } catch (Exception e) {
//            }
//            try {
//                cluster.MS1Score = rsCluster.getFloat("MS1Score");
//                cluster.MS1ScoreProbability = rsCluster.getFloat("MS1Prob");
//                cluster.MS1ScoreLocalProb = rsCluster.getFloat("MS1LProb");
//            } catch (Exception e) {
//            }
//
//            if (UnSortedPeakCurves != null && UnSortedPeakCurves.size() > 0) {
//
//                for (int i = 0; i < MaxNoPeakCluster; i++) {
//                    if (cluster.IsoPeakIndex[i] != 0) {
//                        cluster.IsoPeaksCurves[i] = UnSortedPeakCurves.get(cluster.IsoPeakIndex[i] - 1);
//                    }
//                }
//            }
//
////            cluster.IsoPeaksCurves[0] = PeakCurveList.get(rsCluster.getInt("FirstPeakIdx"));
////            cluster.IsoPeaksCurves[1] = PeakCurveList.get(rsCluster.getInt("SecondPeakIdx"));
////            cluster.IsoPeaksCurves[2] = PeakCurveList.get(rsCluster.getInt("ThirdPeakIdx"));
//            //cluster.IsoPeaksCurves[0] = curveDBReader.Query(ScanCollectionName, rsCluster.getInt("FirstPeakIdx"), parameter);
//            //cluster.IsoPeaksCurves[1] = curveDBReader.Query(ScanCollectionName, rsCluster.getInt("SecondPeakIdx"), parameter);
//            //cluster.IsoPeaksCurves[2] = curveDBReader.Query(ScanCollectionName, rsCluster.getInt("ThirdPeakIdx"), parameter);
//            //PeakClusters.put(GetPeakClusterHashKey(cluster), cluster);
//            PeakClusters.add(cluster);
//            //RTSortedClusters.add(cluster);
//        }
//        state.close();
//        connectionManager.CloseConnection();
//        //GenerateMZSortedClusterList(false);
//        //GenerateRTSoretedClusterList(false);
//        //System.out.print("done\n");
//    }
//
//</editor-fold>
    
    public void ExportPeakClusterRegionTXT() throws IOException {
        FileWriter writer = new FileWriter(FilenameUtils.getFullPath(ScanCollectionName) + FilenameUtils.getBaseName(ScanCollectionName) + "_PeakClusterRange.txt");
        for (PeakCluster cluster : PeakClusters) {
            writer.write(cluster.startRT + "\t" + cluster.endRT + "\t" + (cluster.mz[0] + 1f) + "\t" + (cluster.mz[0]) + "\n");
        }
        writer.close();
    }

    public void ExportPeakCurveTXT() throws IOException {
        FileWriter writer = new FileWriter(FilenameUtils.getFullPath(ScanCollectionName) + FilenameUtils.getBaseName(ScanCollectionName) + "_PeakCurveText.txt");
        for (PeakCurve curve : UnSortedPeakCurves) {
            writer.write(curve.TargetMz + "\t" + curve.ApexRT + "\t" + curve.ApexInt + "\n");
        }
        writer.close();
    }
   
    public void ClearMonoisotopicPeakOfCluster() {
        for (PeakCluster peak : PeakClusters) {
            peak.IsoPeaksCurves = null;
            peak.MonoIsotopePeak = null;
        }
    }

//    public KDTree GetPeakClusterMZSearchTree() {
//        if (PeakClusterMZSearchTree == null) {
//            PeakClusterMZSearchTree = new KDTree(2);
//            for (int i = 0; i < PeakClusters.size(); i++) {
//                PeakCluster peakCluster = PeakClusters.get(i);
//                try {
//                    PeakClusterMZSearchTree.insert(new double[]{peakCluster.PeakHeightRT[0], peakCluster.TargetMz()}, peakCluster);
//                } catch (KeyDuplicateException | KeySizeException ex) {
//                    Logger.getRootLogger().error(ex.getMessage());
//                }
//            }
//        }
//        return PeakClusterMZSearchTree;
//    }
    
    public KDTree GetPeakClusterMassSearchTree() {
        if (PeakClusterMassSearchTree == null) {
            Logger.getRootLogger().info("Building PeakCluster Mass-RT KD tree");
            PeakClusterMassSearchTree = new KDTree(2);
            for (int i = 0; i < PeakClusters.size(); i++) {
                PeakCluster peakCluster = PeakClusters.get(i);
                try {
                    PeakClusterMassSearchTree.insert(new double[]{peakCluster.PeakHeightRT[0], peakCluster.NeutralMass()}, peakCluster);
                } catch (KeyDuplicateException | KeySizeException ex) {
                    Logger.getRootLogger().error(ex.getMessage());
                }
            }
        }
        return PeakClusterMassSearchTree;
    }
       
    public KDTree GetPeakCurveSearchTree(){
        if(PeakCurveSearchTree==null){
            Logger.getRootLogger().info("Building PeakCurve Mass-RT KD tree");
            PeakCurveSearchTree=new KDTree(2);
            for(PeakCurve peakCurve : UnSortedPeakCurves){
                try {
                    PeakCurveSearchTree.insert(new double[]{peakCurve.ApexRT,peakCurve.TargetMz}, peakCurve);
                } catch (KeySizeException ex) {
                    java.util.logging.Logger.getLogger(LCMSPeakBase.class.getName()).log(Level.SEVERE, null, ex);
                } catch (KeyDuplicateException ex) {
                    java.util.logging.Logger.getLogger(LCMSPeakBase.class.getName()).log(Level.SEVERE, null, ex);
                }
            }
        }
        return PeakCurveSearchTree;
    }

    public SortedCurveCollectionApexRT GetPeakCurveListRT() {
        if (PeakCurveListRT == null) {
            PeakCurveListRT=new SortedCurveCollectionApexRT();
            PeakCurveListRT.addAll(UnSortedPeakCurves);
        }
        return PeakCurveListRT;
    }
     public SortedClusterCollectionClassApexRT GetPeakClusterListRT() {
        if (PeakClusterListRT == null) {
            PeakClusterListRT=new SortedClusterCollectionClassApexRT();
            PeakClusterListRT.addAll(PeakClusters);
        }
        return PeakClusterListRT;
    }

}
