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
import MSUmpire.BaseDataStructure.ScanCollection;
import MSUmpire.BaseDataStructure.ScanData;
import MSUmpire.BaseDataStructure.SpectralDataType;
import MSUmpire.BaseDataStructure.XYData;
import MSUmpire.MSMSDBSearch.DBSearchParam;
import MSUmpire.MSMSDBSearch.MSMSDBSearch;
import MSUmpire.MSMSDBSearch.TandemParam;
import MSUmpire.MSMSDBSearch.TandemSearch;
import MSUmpire.MSMSDBSearch.iProphet;
import MSUmpire.MathPackage.ChiSquareGOF;
import MSUmpire.MathPackage.KMeans;
import MSUmpire.PSMDataStructure.LCMSID;
import MSUmpire.PSMDataStructure.PSM;
import MSUmpire.PSMDataStructure.PepIonID;
import MSUmpire.PeakDataStructure.MS1ScoreModel;
import MSUmpire.PeakDataStructure.PeakCluster;
import MSUmpire.PeptidePeakClusterDetection.PDHandlerMS1;
import MSUmpire.SearchResultParser.PepXMLParser;
import MSUmpire.SearchResultParser.TPPResult;
import MSUmpire.SortedListLib.SortedList;
import MSUmpire.SpectraST.SpectraSTCreateLib;
import MSUmpire.SpectraST.SpectraSTSearch;
import MSUmpire.UmpireSearchDataStructure.PepIonCandidate;
import MSUmpire.UmpireSearchDataStructure.PepIonLib;
import MSUmpire.spectrumparser.mzXMLParser;
import jMEF.ExpectationMaximization1D;
import jMEF.MixtureModel;
import jMEF.PVector;
import java.awt.BasicStroke;
import java.awt.Color;
import java.awt.geom.Ellipse2D;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.sql.SQLException;
import java.text.DecimalFormat;
import java.util.*;
import java.util.concurrent.ExecutionException;
import java.util.zip.DataFormatException;
import javastat.multivariate.DiscriminantAnalysis;
import javax.xml.parsers.ParserConfigurationException;
import org.apache.commons.io.FilenameUtils;
import org.apache.commons.lang.exception.ExceptionUtils;
import org.apache.commons.math3.analysis.interpolation.LoessInterpolator;
import org.apache.log4j.Logger;
import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartUtilities;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.chart.plot.XYPlot;
import org.jfree.chart.renderer.xy.XYItemRenderer;
import org.jfree.data.statistics.HistogramDataset;
import org.jfree.data.statistics.HistogramType;
import org.jfree.data.xy.XYDataItem;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;
import org.xml.sax.SAXException;
import org.xmlpull.v1.XmlPullParserException;

/**
 *
 * @author Chih-Chiang Tsou <chihchiang.tsou@gmail.com>
 */
public class LCMSPeakMS1 extends LCMSPeakBase {

    public LCMSID IDsummary;
    //public ArrayList<LCMSPeakDIAMS2> SwathWindows;        
    public HashMap<String, MS1ScoreModel> ScoreModels;
    public boolean UseMappedIon = false;
    public TreeMap<XYData, ArrayList<Integer>> MS1Windows;

//    public LCMSPeakMS1(String Filename, ConnectionManager connectionManager, InstrumentParameter parameter, int NoCPUs) {
//        this.ScanCollectionName = Filename;
//        this.ParentmzXMLName = Filename;
//        this.datattype = SpectralDataType.DataType.DDA;
//        this.NoCPUs = NoCPUs;
//        SetParameter(parameter);
//        //SetMySQLConnection(connectionManager);
//    }
     public LCMSPeakMS1(String Filename, InstrumentParameter parameter, int NoCPUs) {
        this.ScanCollectionName = Filename;
        this.ParentmzXMLName = Filename;
        this.datattype = SpectralDataType.DataType.DDA;
        this.NoCPUs = NoCPUs;
        SetParameter(parameter);
    }
    
    public LCMSPeakMS1(String Filename, int NoCPUs) {
        this.ScanCollectionName = Filename;
        this.ParentmzXMLName=Filename;
        this.datattype=SpectralDataType.DataType.DDA;        
        this.NoCPUs = NoCPUs;
    }
    
    public void SetMS1Windows(TreeMap<XYData, ArrayList<Integer>> MS1Windows){
        this.MS1Windows=MS1Windows;
    }

    public void SetParameter(InstrumentParameter parameter){
        this.parameter = parameter;
        this.MaxNoPeakCluster = parameter.MaxNoPeakCluster;
        this.MinNoPeakCluster = parameter.MinNoPeakCluster;
        ChiSquareGOF.GetInstance(this.MaxNoPeakCluster);
        this.StartCharge = parameter.StartCharge;
        this.EndCharge = parameter.EndCharge;
        this.MiniIntensity = parameter.MinMSIntensity;
        this.SNR = parameter.SNThreshold;
    }
    
    public void DoTandem(TandemParam para) throws IOException, InterruptedException, ParserConfigurationException, SAXException, XmlPullParserException, ClassNotFoundException {
        para.NoCPUs = NoCPUs;
        para.parameterPath = FilenameUtils.separatorsToUnix(FilenameUtils.getFullPath(ScanCollectionName) + FilenameUtils.getBaseName(ScanCollectionName) + "_tandem.param");
        para.RawSearchResult = FilenameUtils.separatorsToUnix(FilenameUtils.getFullPath(ScanCollectionName) + FilenameUtils.getBaseName(ScanCollectionName) + ".tandem");
        para.InteractPepXMLPath = FilenameUtils.separatorsToUnix(FilenameUtils.getFullPath(ScanCollectionName) + "interact-" + FilenameUtils.getBaseName(ScanCollectionName) + ".pep.xml");
        para.ProtXMLPath=para.InteractPepXMLPath.replace(".pep.xml", ".prot.xml");
        para.SpectrumPath = ScanCollectionName;
        TandemSearch search = new TandemSearch(para);
        Logger.getRootLogger().info("Identifying protein by X!Tandem.....");
        search.RunTandem();
    }
    
    public void DBSearch(MSMSDBSearch dbsearch) throws IOException, InterruptedException, ParserConfigurationException, SAXException, XmlPullParserException, ClassNotFoundException {
        dbsearch.SetResultFilePath(ParentmzXMLName);
        dbsearch.DBSearch();        
    }

    public void GeneratePeptideFragmentByMSMS(TandemParam para) throws InterruptedException, ExecutionException, IOException, SQLException {
        GetmzXML().GetAllScanCollectionMS2Only(true, true);
        IDsummary.GenerateFragmentPeakForPepIonByMSMS(GetmzXML().scanCollection, para.FragPPM);
        //IDsummary.ExportPepFragmentPeak(connectionManager);
    }

    public void WriteLCMSIDSerialization(){
        WriteLCMSIDSerialization("");
    }
    public void WriteLCMSIDSerialization(String tag){
        this.IDsummary.WriteLCMSIDSerialization(ParentmzXMLName,tag);
    }
    public boolean ReadSerializedLCMSID() throws Exception{
        return ReadSerializedLCMSID("");
    }
    public boolean ReadSerializedLCMSID(String tag) throws Exception{
        this.IDsummary=LCMSID.ReadLCMSIDSerialization(ParentmzXMLName,tag);            
        if(this.IDsummary==null){
            return false;
        }
        this.IDsummary.mzXMLFileName=ParentmzXMLName;
        this.IDsummary.Filename=ParentmzXMLName;
        return true;
    }
    
    public void ConstructSpecLib(TandemParam para) throws IOException, InterruptedException, ParserConfigurationException, SAXException, XmlPullParserException, ClassNotFoundException {
        para.NoCPUs = NoCPUs;
        para.parameterPath = FilenameUtils.separatorsToUnix(FilenameUtils.getFullPath(ScanCollectionName) + FilenameUtils.getBaseName(ScanCollectionName) + "_tandem.param");
        para.RawSearchResult = FilenameUtils.separatorsToUnix(FilenameUtils.getFullPath(ScanCollectionName) + FilenameUtils.getBaseName(para.InteractPepXMLPath) + ".tandem");
        para.InteractPepXMLPath = FilenameUtils.separatorsToUnix(FilenameUtils.getFullPath(ScanCollectionName) + "interact-" + FilenameUtils.getBaseName(ScanCollectionName)+ ".pep.xml");
        para.ProtXMLPath = para.InteractPepXMLPath.replace(".pep.xml", ".prot.xml");
        para.SpectrumPath = ScanCollectionName;
        SpectraSTCreateLib spectraSTCreateLib = new SpectraSTCreateLib();
        spectraSTCreateLib.CreateLib(para.PepXMLPath);
    }
   
    public String GetDefaultPepXML() {
        return FilenameUtils.separatorsToUnix(FilenameUtils.getFullPath(ScanCollectionName) + "interact-" + FilenameUtils.getBaseName(ScanCollectionName) + ".pep.xml");
    }

    public String GetDefaultProtXML() {
        return GetDefaultPepXML().replace(".pep.xml",".prot.xml");
    }

    public void SetDefaultPath(DBSearchParam para) {
        para.PepXMLPath = GetDefaultPepXML();
        para.ProtXMLPath = GetDefaultProtXML();
    }

    public void ParseSearchEngineResult(MSMSDBSearch dbsearch) throws ParserConfigurationException, SAXException, IOException, XmlPullParserException, ClassNotFoundException, InterruptedException {   
        dbsearch.SetResultFilePath(ParentmzXMLName);        
        this.IDsummary = new LCMSID(ScanCollectionName,dbsearch.GetParameter().DecoyPrefix,dbsearch.GetParameter().FastaPath);        
        TPPResult tppresult = new TPPResult(dbsearch.GetParameter().PepFDR, dbsearch.GetParameter().ProtFDR, dbsearch.GetParameter().DecoyPrefix);
        ParseLuciphor(IDsummary);
        tppresult.ReadSearchResult(IDsummary, dbsearch.GetParameter().InteractPepXMLPath, dbsearch.GetParameter().ProtXMLPath);        
    }
    
    public void ParseTPP(DBSearchParam para, boolean UserDefinedPath) throws ParserConfigurationException, IOException, InterruptedException, ClassNotFoundException, XmlPullParserException, SAXException {
        this.IDsummary = new LCMSID(ScanCollectionName,para.DecoyPrefix,para.FastaPath);        
        if (!UserDefinedPath) {
            SetDefaultPath(para);
        }
        TPPResult tppresult = new TPPResult(para.PepFDR, para.ProtFDR, para.DecoyPrefix);
        ParseLuciphor(IDsummary);
        tppresult.ReadSearchResult(IDsummary, para.PepXMLPath, para.ProtXMLPath);
    }

    public void ParseTPP(DBSearchParam para) throws ParserConfigurationException, IOException, InterruptedException, ClassNotFoundException, XmlPullParserException, SAXException {
        ParseTPP(para, false);
    }

    public void ParseTPPByRefID(DBSearchParam para, LCMSID ReferenceID) throws ParserConfigurationException, IOException, InterruptedException, XmlPullParserException, SAXException, ClassNotFoundException {
        SetDefaultPath(para);
        //para.PepXMLPath = FilenameUtils.separatorsToUnix(FilenameUtils.getFullPath(ScanCollectionName)  + FilenameUtils.getBaseName(ScanCollectionName) + ".Comet.pep.xml");
        this.IDsummary = new LCMSID(ScanCollectionName,para.DecoyPrefix,para.FastaPath);        
        TPPResult tppresult = new TPPResult(para.PepFDR, para.ProtFDR,para.DecoyPrefix);
        ParseLuciphor(IDsummary);
        tppresult.ReadSearchResultByRefID(IDsummary, para.PepXMLPath, para.ProtXMLPath, ReferenceID);
    }

    public void ParseLuciphor(LCMSID IDsummary) throws FileNotFoundException, IOException {
        if (!new File(ScanCollectionName.replace(".mzXML", ".luciphor.tsv")).exists()) {
            return;
        }
        IDsummary.LuciphorResult = new HashMap<>();
        BufferedReader reader = new BufferedReader(new FileReader(ScanCollectionName.replace(".mzXML", ".luciphor.tsv")));
        String line = "";
        reader.readLine();
        while ((line = reader.readLine()) != null) {
            IDsummary.LuciphorResult.put(line.split("\t")[0], line);
        }
        reader.close();
    }

    public void ParseTPPByProbThreshold(DBSearchParam para, float prob) throws  InterruptedException, ClassNotFoundException, XmlPullParserException, SAXException, ParserConfigurationException, IOException {
        para.PepXMLPath = FilenameUtils.separatorsToUnix(FilenameUtils.getFullPath(ScanCollectionName) + "interact-" + FilenameUtils.getBaseName(ScanCollectionName) + ".pep.xml");
        para.ProtXMLPath = para.PepXMLPath.replace(".pep.xml", ".prot.xml");
//para.PepXMLPath = FilenameUtils.separatorsToUnix(FilenameUtils.getFullPath(ScanCollectionName)  + FilenameUtils.getBaseName(ScanCollectionName) + ".Comet.pep.xml");
        this.IDsummary = new LCMSID(ScanCollectionName,para.DecoyPrefix,para.FastaPath);
        IDsummary.FastaPath = para.FastaPath;
        TPPResult tppresult = new TPPResult(-1f, -1f, para.DecoyPrefix);
        tppresult.ReadSearchResultAndFilterByProb(IDsummary, para.PepXMLPath, para.ProtXMLPath, prob);
    }

    public void ParseAssignedTPP(DBSearchParam para, float fdr) throws ParserConfigurationException, IOException, InterruptedException, ClassNotFoundException, XmlPullParserException, SAXException {
        this.IDsummary = new LCMSID(ScanCollectionName, para.DecoyPrefix,para.FastaPath);        
       TPPResult tppresult = new TPPResult(fdr, fdr, para.DecoyPrefix);
        tppresult.FilterIDBymzXMLname = true;
        tppresult.ReadSearchResult(IDsummary, para.PepXMLPath, para.ProtXMLPath);
    }

    public void UmpireSearch(PepIonLib ionlib) throws InterruptedException, ExecutionException, IOException {
        GetmzXML().GetAllScanCollectionMS2Only(true, true);
        ScanCollection Scans = GetmzXML().scanCollection;
        for (ScanData scan : Scans.ScanHashMap.values()) {
            float lowmz = InstrumentParameter.GetMzByPPM(scan.PrecursorMz, scan.PrecursorCharge, parameter.MS1PPM);
            float highmz = InstrumentParameter.GetMzByPPM(scan.PrecursorMz, scan.PrecursorCharge, -parameter.MS1PPM);
            ArrayList<PepIonCandidate> candidates = ionlib.IonMzLib.GetCandidate(lowmz, highmz, scan.PrecursorCharge);
        }
    }

//    public void DiscardMappedIonTable() throws SQLException {
//        if (connectionManager != null) {
//            Connection connection = connectionManager.GetConnection();
////         this.IDsummary.DiscardMappedPepIonTable(connection);
//            try (Statement state = connection.createStatement()) {
//                state.execute("DROP TABLE " + FilenameUtils.getBaseName(ScanCollectionName) + "_MappedPepIonIDs");
//
//            } catch (Exception ex) {
//                Logger.getRootLogger().error("(Discarding table: " + FilenameUtils.getBaseName(ScanCollectionName) + "_MappedPepIonIDs failed.)");
//                Logger.getRootLogger().error(ExceptionUtils.getStackTrace(ex));
//            }
//            connectionManager.CloseConnection();
//        }
//    }

     public void AssignMS1Cluster() {
        for (PepIonID pepIonID : IDsummary.GetPepIonList().values()) {
            for (String indexint : pepIonID.MS1ClusIndex.split(";")) {
                if (!"".equals(indexint)) {
                    pepIonID.MS1PeakClusters.add(PeakClusters.get(Integer.parseInt(indexint) - 1));
                }
            }
        }
        if (IDsummary.GetMappedPepIonList() != null) {
            for (PepIonID pepIonID : IDsummary.GetMappedPepIonList().values()) {
                for (String indexint : pepIonID.MS1ClusIndex.split(";")) {
                    if (!"".equals(indexint)) {
                        pepIonID.MS1PeakClusters.add(PeakClusters.get(Integer.parseInt(indexint) - 1));
                    }
                }
            }
        }
    }
//    public void ReadIDFromDB(String mzxmlname, String FastaFile) throws SQLException, IOException, XmlPullParserException, ClassNotFoundException, InterruptedException {
//
//        Logger.getRootLogger().info("Loading ID result from MySQL DB :" + FilenameUtils.getBaseName(mzxmlname) + "...");
//        this.IDsummary = new LCMSID(FilenameUtils.getFullPath(ScanCollectionName) + FilenameUtils.getBaseName(mzxmlname),"rev",FastaFile);        
//        Connection connection = connectionManager.GetConnection();
//        this.IDsummary.ReadFromDBPepIon(connection);
//        this.IDsummary.ReadFromDBPepFragments(connection);
//        if (UseMappedIon) {
//            this.IDsummary.ReadFromDBMappedPepIon(connection, false,0f);
//        }        
//        this.IDsummary.ReadFromDBPSM(connection);
//        this.IDsummary.ReadFromDBProt(connection);
//        this.IDsummary.AssignProtForPepIon();
//        if (UseMappedIon) {
//            this.IDsummary.AssignProtForMappedIon();
//        }
//        this.IDsummary.GenearteAssignIonList();
//        this.IDsummary.GenerateIndisProtMap();
//        //connectionManager.CloseConnection();
//        Logger.getRootLogger().info("Protein No.:" + IDsummary.ProteinList.size() + "; All assigned peptide ions:" + IDsummary.AssignedPepIonList.size() + "; Identified peptide ions:" + IDsummary.GetPepIonList().size());
//    }

    public void AssignIDResult(LCMSID summary) {
        this.IDsummary = summary;
    }

    public void SetmzXML(mzXMLParser mzxml) {
        this.mzxml = mzxml;
    }

    public mzXMLParser GetmzXML() {
        if (mzxml == null) {
            try {
                mzxml = new mzXMLParser(ScanCollectionName, parameter, datattype, null, NoCPUs);
            } catch (Exception ex) {
                Logger.getRootLogger().error(ExceptionUtils.getStackTrace(ex));
            }
        }
        return mzxml;
    }

    
    
    public void PeakClusterDetection() throws FileNotFoundException, IOException, InterruptedException, ExecutionException, ParserConfigurationException, SAXException, DataFormatException, XmlPullParserException {

        if (Resume && ReadIfProcessed()) {
            return;
        }
        CreatePeakFolder();
        ArrayList<ScanCollection> scanCollections = new ArrayList<>();
        parameter.NoPeakPerMin = (int) (5f / GetmzXML().GetMS1CycleTime());
        Logger.getRootLogger().info("MS1 average cycle time : "+GetmzXML().GetMS1CycleTime()*60+ " seconds");
        if (MS1Windows == null || MS1Windows.isEmpty()) {
            GetmzXML().GetAllScanCollectionByMSLabel(true, true, true, false,parameter.startRT, parameter.endRT);
            scanCollections.add(GetmzXML().scanCollection);
            if (parameter.DetermineBGByID) {
                DetermineBGByLowIntID(GetmzXML().scanCollection);
            }
        }
        else{
            for( XYData window : MS1Windows.keySet()){
                scanCollections.add(GetmzXML().GetScanCollectionMS1Window(window, true));
            }
        }
        
        PDHandlerMS1 detection = new PDHandlerMS1(this, NoCPUs, parameter.MS1PPM);        
        if(parameter.TargetIDOnly){
            detection.SetTargetedDetectionOnly();
        }
        if (IDsummary != null) {
            for (PepIonID pepIonID : IDsummary.GetPepIonList().values()) {
                for (int i = 0; i < parameter.MaxNoPeakCluster; i++) {
                    detection.AddToInclusionList(pepIonID.GetPeakMz(i), pepIonID.GetRT());
                }
            }
        }
        detection.DetectPeakClusters(scanCollections);
        MapScanNoForPeakClusters();
        
        if (ExportPeakClusterTable) {
            //ExportPeakCurveResult();
            ExportPeakClusterResultCSV();
        }
        //GenerateMZSortedClusterList(false);
        
        ExportPeakCluster();
    }

    public void MapScanNoForPeakClusters() {
        for (PeakCluster peakCluster : PeakClusters) {
            peakCluster.StartScan = GetmzXML().GetScanNoByRT(peakCluster.startRT);
            peakCluster.EndScan = GetmzXML().GetScanNoByRT(peakCluster.endRT);
        }
    }
    
    public void DetermineBGByLowIntID(ScanCollection scanCollection) {
        if(IDsummary==null){
            return;
        }
        Logger.getRootLogger().info("Estimating background by median of 10 lowest intensity precursors of identified peptide ions...");
        
        ArrayList<Float> LowList=new ArrayList<>();
        float lowest = Float.MAX_VALUE;
        
        int ScanNo=0;
        for (PSM psm : IDsummary.PSMList.values()) {
            ScanData MS2Scan = scanCollection.GetScan(psm.ScanNo);
            if (MS2Scan.PrecursorIntensity > 0f) {
                LowList.add(MS2Scan.PrecursorIntensity);
                if (lowest > MS2Scan.PrecursorIntensity) {
                    lowest = MS2Scan.PrecursorIntensity;
                    ScanNo = psm.ScanNo;
                }
            } else {
                ScanData MS1Scan = scanCollection.GetParentMSScan(psm.ScanNo);                
                float lowmz = InstrumentParameter.GetMzByPPM(psm.ObserPrecursorMz(), psm.Charge, parameter.MS1PPM);
                float highmz = InstrumentParameter.GetMzByPPM(psm.ObserPrecursorMz(), psm.Charge, -parameter.MS1PPM);
                int startindx = MS1Scan.GetLowerIndexOfX(lowmz);
                int endindx = MS1Scan.GetHigherIndexOfX(highmz);
                float maxint = -1f;
                for (int idx = startindx; idx <= endindx; idx++) {
                    if (InstrumentParameter.CalcPPM(psm.ObserPrecursorMz(), MS1Scan.Data.get(idx).getX()) < parameter.MS1PPM) {
                        if (MS1Scan.Data.get(idx).getY() > maxint) {
                            maxint = MS1Scan.Data.get(idx).getY();
                        }
                    }
                }
                if (maxint != -1) {
                    LowList.add(maxint);
                    if (lowest > maxint) {
                        lowest = maxint;
                        ScanNo = psm.ScanNo;
                    }
                }
            }
        }
        SortedList<Float> LowQueue=new SortedList(new Comparator<Float>() {

            @Override
            public int compare(Float o1, Float o2) {
                return Float.compare(o1, o2);
            }
        });
        LowQueue.addAll(LowList);
        
        parameter.MinMSIntensity=LowQueue.get(5)*0.3f;
        MiniIntensity = parameter.MinMSIntensity;
        scanCollection.RemoveBackground(1, MiniIntensity);
        Logger.getRootLogger().info("Lowest identified precursor intensity: "+lowest+" at ScanNo:"+ScanNo+ ", median intensity of lowest 10 precursor is "+LowQueue.get(5));
        Logger.getRootLogger().info("Use 30% of the median as background intensity");
    }

    public void DetermineMS1PeakClusterScore() throws IOException {
        if (IDsummary == null) {
            return;
        }
        GenerateScoreLDAModelsByQ();
        MixtureModeling();
        EstimateMS1Prob();
        AssignProb();
        GenerateCharts();
    }

    public void AssignProb() {
        int NoPoints = 500;
        for (PeakCluster cluster : PeakClusters) {
            MS1ScoreModel model = ScoreModels.get(String.valueOf(cluster.GetQualityCategory()));
            if (model.IDList.isEmpty()) {
                continue;
            }

            if (model.valid) {
                for (int i = 0; i < NoPoints; i++) {
                    if (cluster.MS1Score >= model.Probability[i][0]) {
                        cluster.MS1ScoreProbability = model.Probability[i][1];
                        cluster.MS1ScoreLocalProb = model.Probability[i][2];
                        break;
                    }
                }
            } else {
                cluster.MS1ScoreProbability = Float.NaN;
                cluster.MS1ScoreLocalProb = Float.NaN;
            }
        }
    }

    public void EstimateMS1Prob() {
        for (String key : ScoreModels.keySet()) {
            MS1ScoreModel model = ScoreModels.get(key);
            if (model.IDList.isEmpty()) {
                continue;
            }
            if (((PVector) model.mmc.param[0]).array[0] < ((PVector) model.mmc.param[1]).array[0]) {
                model.valid = false;
                continue;
            }
            String pngfile = FilenameUtils.getFullPath(ScanCollectionName) + "/" + FilenameUtils.getBaseName(ScanCollectionName) + "_Model_" + key + ".png";
            XYSeries model1 = new XYSeries("Positive model");
            XYSeries model2 = new XYSeries("Negative model");

            int NoPoints = 500;
            float intv = (model.max - model.min) / NoPoints;
            PVector point = new PVector(2);
            for (int i = 0; i < NoPoints; i++) {
                point.array[0] = model.max - i * intv;
                point.array[1] = model.mmc.EF.density(point, model.mmc.param[0]) * model.mmc.weight[0];

                model1.add(point.array[0], point.array[1]);
                point.array[1] = model.mmc.EF.density(point, model.mmc.param[1]) * model.mmc.weight[1];
                model2.add(point.array[0], point.array[1]);
            }

            model.Probability = new Float[NoPoints][3];
            float positiveaccu = 0f;
            float negativeaccu = 0f;

            for (int i = 0; i < NoPoints; i++) {
                float positiveNumber = model1.getY(NoPoints - 1 - i).floatValue();
                float negativeNumber = model2.getY(NoPoints - 1 - i).floatValue();
                model.Probability[i][0] = model1.getX(NoPoints - 1 - i).floatValue();
                positiveaccu += positiveNumber;
                negativeaccu += negativeNumber;
                model.Probability[i][2] = positiveNumber / (negativeNumber + positiveNumber);
                model.Probability[i][1] = positiveaccu / (negativeaccu + positiveaccu);
            }

            XYSeriesCollection dataset = new XYSeriesCollection();

            dataset.addSeries(model1);
            dataset.addSeries(model2);

            JFreeChart chart = ChartFactory.createXYLineChart(FilenameUtils.getBaseName(pngfile), "Score", "Frequency", dataset, PlotOrientation.VERTICAL, true, true, false);
            XYPlot plot = chart.getXYPlot();
            plot.setBackgroundPaint(Color.white);
            plot.setDomainGridlinePaint(Color.white);
            plot.setRangeGridlinePaint(Color.white);
            plot.setForegroundAlpha(0.8f);
            chart.setBackgroundPaint(Color.white);
            try {
                ChartUtilities.saveChartAsPNG(new File(pngfile), chart, 1000, 600);
            } catch (IOException e) {
            }
        }
    }

    public void MixtureModeling() {
        for (String key : ScoreModels.keySet()) {
            MS1ScoreModel model = ScoreModels.get(key);

            if (model.IDList.isEmpty()) {
                continue;
            }

            PVector[] points = new PVector[model.UnIDList.size()];
            PVector[] centroids = new PVector[2];

            double[] IDScore = new double[model.IDList.size()];
            double[] UnIDScore = new double[model.UnIDList.size()];
            double IDmean = 0d;
            double UnIDmean = 0d;
            for (int i = 0; i < model.IDList.size(); i++) {
                IDScore[i] = model.IDList.get(i).MS1Score;
                IDmean += model.IDList.get(i).MS1Score;
            }
            for (int i = 0; i < model.UnIDList.size(); i++) {
                UnIDScore[i] = model.UnIDList.get(i).MS1Score;
                UnIDmean += model.UnIDList.get(i).MS1Score;
            }
            Arrays.sort(IDScore);

            IDmean /= model.IDList.size();
            UnIDmean /= model.UnIDList.size();
            //model.ttest = new TwoSampMeansTTest(0, "equal", IDScore, UnIDScore);

            ArrayList<Double> lowerQueue = new ArrayList<>();

            for (int i = 0; i < model.UnIDList.size(); i++) {
                points[i] = new PVector(1);
                points[i].array[0] = UnIDScore[i];

                if (UnIDScore[i] < IDmean) {
                    lowerQueue.add(UnIDScore[i]);
                }
            }
            Collections.sort(lowerQueue);
            centroids[0] = new PVector(1);
            centroids[0].array[0] = IDmean;
            centroids[1] = new PVector(1);
            centroids[1].array[0] = lowerQueue.get(lowerQueue.size() / 2);
            Vector<PVector>[] clusters = KMeans.run(points, 2, centroids);
            MixtureModel mm = ExpectationMaximization1D.initialize(clusters);
            model.mmc = ExpectationMaximization1D.run(points, mm);
            DecimalFormat df = new DecimalFormat("#.####");
            Logger.getRootLogger().debug("----------------------------------------------------------------------------------------");
            Logger.getRootLogger().debug("Mixture modeling for clusters of quality category " + key );
            Logger.getRootLogger().debug("Mean of identified clusters=" + df.format(IDmean) );
            Logger.getRootLogger().debug("Mean of unidentified clusters=" + df.format(UnIDmean) );
            //.getRootLogger().debug("T-test: p-value=" + df.format(model.ttest.pValue).toString());
            Logger.getRootLogger().debug("Model 1 mean=" + df.format(((PVector) model.mmc.param[0]).array[0]) + " sd=" + df.format(((PVector) model.mmc.param[0]).array[1]) + " weight=" + df.format(model.mmc.weight[0]) );
            Logger.getRootLogger().debug("Model 2 mean=" + df.format(((PVector) model.mmc.param[1]).array[0]) + " sd=" + df.format(((PVector) model.mmc.param[1]).array[1]) + " weight=" + df.format(model.mmc.weight[1]) );
        }
    }

    public void AssignQuant() throws SQLException, IOException {
        AssignQuant(true);
    }

     public void DoTPPWithName(TandemParam para, String FastaFile, String SpecifiedFilename) throws IOException, InterruptedException, ParserConfigurationException, SAXException, XmlPullParserException, ClassNotFoundException {
        para.FastaPath = FilenameUtils.separatorsToUnix(FilenameUtils.getFullPath(FastaFile) + FilenameUtils.getName(FastaFile));
        para.NoCPUs = NoCPUs;
        para.parameterPath = FilenameUtils.separatorsToUnix(FilenameUtils.getFullPath(ScanCollectionName) + FilenameUtils.getBaseName(SpecifiedFilename) + "_tandem.param");
        para.RawSearchResult = FilenameUtils.separatorsToUnix(FilenameUtils.getFullPath(ScanCollectionName) + FilenameUtils.getBaseName(SpecifiedFilename) + ".tandem");
        para.PepXMLPath = FilenameUtils.separatorsToUnix(FilenameUtils.getFullPath(ScanCollectionName) + "interact-" + FilenameUtils.getBaseName(SpecifiedFilename) + ".pep.xml");
        para.ProtXMLPath = para.PepXMLPath.replace(".pep.xml", ".prot.xml");
        para.SpectrumPath = ScanCollectionName;
        TandemSearch search = new TandemSearch(para);
        Logger.getRootLogger().info("Processing TPP.....");
        search.RunTandem();
    }
     public void DoSpecLibSearch(TandemParam para) throws IOException, InterruptedException, ParserConfigurationException, SAXException, XmlPullParserException, ClassNotFoundException {
        para.NoCPUs = NoCPUs;
        para.parameterPath = FilenameUtils.separatorsToUnix(FilenameUtils.getFullPath(ScanCollectionName) + FilenameUtils.getBaseName(ScanCollectionName) + "_tandem.param");
        para.RawSearchResult = FilenameUtils.separatorsToUnix(FilenameUtils.getFullPath(ScanCollectionName) + FilenameUtils.getBaseName(ScanCollectionName));
        para.InteractPepXMLPath = FilenameUtils.separatorsToUnix(FilenameUtils.getFullPath(ScanCollectionName) + "interact-" + FilenameUtils.getBaseName(ScanCollectionName) + ".pep.xml");
        para.ProtXMLPath = para.InteractPepXMLPath.replace(".pep.xml", ".prot.xml");
        para.SpectrumPath = ScanCollectionName;
        SpectraSTSearch spsearch = new SpectraSTSearch();
        spsearch.RunSpectraST(ScanCollectionName);
    }
     
    public void AssignQuant(boolean export) throws SQLException, IOException {
        if (IDsummary == null) {
            return;
        }
        if (!FilenameUtils.getBaseName(IDsummary.mzXMLFileName).equals(FilenameUtils.getBaseName(ScanCollectionName))) {
            return;
        }
        Logger.getRootLogger().info("Assigning peak cluster to peptide IDs......");
        Logger.getRootLogger().info("Total peptide ions to be quantified..:"+IDsummary.GetPepIonList().size());
        int count=0;
        for (PepIonID pepIonID : IDsummary.GetPepIonList().values()) {
            pepIonID.CreateQuantInstance(MaxNoPeakCluster);
            pepIonID.MS1PeakClusters.clear();
            ArrayList<PeakCluster> clusterList = FindAllPeakClustersForPepByPSM(pepIonID);
            if (!clusterList.isEmpty()) {
                count++;
                PeakCluster targetCluster = clusterList.get(0);
                float Intensity = targetCluster.PeakHeight[0];

                for (int i = 1; i < clusterList.size(); i++) {
                    if (clusterList.get(i).PeakHeight[0] > Intensity) {
                        targetCluster = clusterList.get(i);
                        Intensity = clusterList.get(i).PeakHeight[0];
                    }
                }
                pepIonID.PeakArea = targetCluster.PeakArea;
                pepIonID.PeakHeight = targetCluster.PeakHeight;
                if (!pepIonID.MS1PeakClusters.contains(targetCluster)) {
                    pepIonID.MS1PeakClusters.add(targetCluster);
                }
                pepIonID.PeakClusterScore = targetCluster.MS1Score;
                pepIonID.PeakRT = targetCluster.PeakHeightRT[0];
                pepIonID.ObservedMz = targetCluster.mz[0];
                targetCluster.AssignedPepIon = pepIonID.GetKey();
                targetCluster.Identified = true;
            }
            else {
                Logger.getRootLogger().warn("Feature for "+pepIonID.ModSequence+ " not found: Charge="+pepIonID.Charge+", mz="+pepIonID.ObservedMz+", RT="+pepIonID.GetRT());
            }
        }
        DecimalFormat df = new DecimalFormat("###.##");
              
        Logger.getRootLogger().info("No. of peptide ions quantified..:"+count+"("+ df.format((float)count/IDsummary.GetPepIonList().size()*100)+"%)");
        
        //DetermineMS1PeakClusterScore();        
        //ExportPeakCluster();
        if (export) {
            ExportID();
        }
        //System.out.print("done\n");
    }

    public void AssignMappedPepQuant() throws SQLException, IOException {
        if (IDsummary == null || IDsummary.GetMappedPepIonList().isEmpty()) {
            return;
        }
        if (!FilenameUtils.getBaseName(IDsummary.mzXMLFileName).equals(FilenameUtils.getBaseName(ScanCollectionName))) {
            return;
        }
        Logger.getRootLogger().info("Assigning peak cluster to mapped peptide IDs......");
        for (PepIonID pepIonID : IDsummary.GetMappedPepIonList().values()) {
//            if(pepIonID.ModSequence.equals("TTLLH[15.994522(M)]MLKDDR")){
//                System.out.print("");
//            }
            pepIonID.CreateQuantInstance(MaxNoPeakCluster);

            ArrayList<PeakCluster> clusterList = FindAllPeakClustersForMappedPep(pepIonID);
//            if(clusterList.isEmpty()){
//                for(PSM psm : pepIonID.GetPSMList()){
//                    System.out.print(pepIonID.ModSequence+" mz:"+psm.ObserPrecursorMz()+" RT:"+psm.RetentionTime+"\n");
//                }
//                System.out.print("");
//            }

            if (!clusterList.isEmpty()) {
                PeakCluster targetCluster = null;
                float Score = 0f;

                for (int i = 0; i < clusterList.size(); i++) {
                    PeakCluster cluster = clusterList.get(i);
                    if ("".equals(cluster.AssignedPepIon)) {
                        Score = cluster.PeakHeight[0] * cluster.MS1Score;
                        if (targetCluster == null || clusterList.get(i).MS1Score > Score) {
                            targetCluster = cluster;
                            Score = cluster.MS1Score;
                        }
                    }
                }
                if (targetCluster != null) {
                    pepIonID.PeakArea = targetCluster.PeakArea;
                    pepIonID.PeakHeight = targetCluster.PeakHeight;
                    pepIonID.MS1PeakClusters.add(targetCluster);
                    pepIonID.PeakClusterScore = targetCluster.MS1Score;
                    pepIonID.PeakRT = targetCluster.PeakHeightRT[0];
                    pepIonID.ObservedMz = targetCluster.mz[0];
                    targetCluster.AssignedPepIon = pepIonID.GetKey();
                }
            }
        }
        //ExportMappedIDQuant();
        //System.out.print("done\n");
    }

    public boolean ReadIfProcessed(){
        return ReadPeakCluster();        
//        Connection connection = connectionManager.GetConnection();
//        Statement state = connection.createStatement();
//        ResultSet rsCluster = state.executeQuery("SHOW TABLES LIKE '" + FilenameUtils.getBaseName(ScanCollectionName) + "_PeakCluster'");
//        if (rsCluster.next()) {
//            ReadCurvePeakFromDB(true);
//            ReadPeakCluster();
//            return true;
//        }
//        return false;
    }

    private void GenerateCharts() {

        for (String key : ScoreModels.keySet()) {
            MS1ScoreModel model = ScoreModels.get(key);
            if (model.IDList.isEmpty()) {
                continue;
            }

            double[] IDCorr = new double[model.IDList.size()];
            double[] UnIDCorr = new double[model.UnIDList.size()];
            double[] IDIsoPattern = new double[model.IDList.size()];
            double[] UnIDIsoPattern = new double[model.UnIDList.size()];
            double[] IDConflictCorr = new double[model.IDList.size()];
            double[] UnIDConflictCorr = new double[model.UnIDList.size()];
            double[] IDPeakApexVar = new double[model.IDList.size()];
            double[] UnIDPeakApexVar = new double[model.UnIDList.size()];
            double[] IDNoRidges = new double[model.IDList.size()];
            double[] UnIDNoRidges = new double[model.UnIDList.size()];
            double[] IDSymScore = new double[model.IDList.size()];
            double[] UnIDSymScore = new double[model.UnIDList.size()];
            double[] IDSNR = new double[model.IDList.size()];
            double[] UnIDSNR = new double[model.UnIDList.size()];
            double[] IDLDA = new double[model.IDList.size()];
            double[] UnIDLDA = new double[model.UnIDList.size()];
            int idx = 0;
            for (PeakCluster cluster : model.IDList) {
                if (cluster.GetQualityCategory() == 2) {
                    IDCorr[idx] = cluster.Corrs[0];
                } else {
                    IDCorr[idx] = cluster.Corrs[0] + cluster.Corrs[1];
                }
                IDSNR[idx] = Math.log(cluster.GetSNR(0));
                IDIsoPattern[idx] = cluster.IsoMapProb;
                IDConflictCorr[idx] = cluster.GetConflictCorr();
                IDPeakApexVar[idx] = cluster.GetApexVar();
                IDNoRidges[idx] = cluster.NoRidges;
                IDSymScore[idx] = cluster.GetSymScore();
                IDLDA[idx] = cluster.MS1Score;
                idx++;
            }
            idx = 0;
            for (PeakCluster cluster : model.UnIDList) {
                if (cluster.GetQualityCategory() == 2) {
                    UnIDCorr[idx] = cluster.Corrs[0];
                } else {
                    UnIDCorr[idx] = cluster.Corrs[0] + cluster.Corrs[1];
                }
                UnIDSNR[idx] = Math.log(cluster.GetSNR(0));
                UnIDIsoPattern[idx] = cluster.IsoMapProb;
                UnIDConflictCorr[idx] = cluster.GetConflictCorr();
                UnIDPeakApexVar[idx] = cluster.GetApexVar();
                UnIDNoRidges[idx] = cluster.NoRidges;
                UnIDSymScore[idx] = cluster.GetSymScore();
                UnIDLDA[idx] = cluster.MS1Score;
                idx++;
            }

            GenerateHistogram(FilenameUtils.getFullPath(ScanCollectionName) + "/" + FilenameUtils.getBaseName(ScanCollectionName) + "_SNR_" + key + ".png", IDSNR, UnIDSNR);
            GenerateHistogram(FilenameUtils.getFullPath(ScanCollectionName) + "/" + FilenameUtils.getBaseName(ScanCollectionName) + "_Corr_" + key + ".png", IDCorr, UnIDCorr);
            GenerateHistogram(FilenameUtils.getFullPath(ScanCollectionName) + "/" + FilenameUtils.getBaseName(ScanCollectionName) + "_IsoPattern_" + key + ".png", IDIsoPattern, UnIDIsoPattern);
            GenerateHistogram(FilenameUtils.getFullPath(ScanCollectionName) + "/" + FilenameUtils.getBaseName(ScanCollectionName) + "_ConflictCorr_" + key + ".png", IDConflictCorr, UnIDConflictCorr);
            GenerateHistogram(FilenameUtils.getFullPath(ScanCollectionName) + "/" + FilenameUtils.getBaseName(ScanCollectionName) + "_ApexVar_" + key + ".png", IDPeakApexVar, UnIDPeakApexVar);
            GenerateHistogram(FilenameUtils.getFullPath(ScanCollectionName) + "/" + FilenameUtils.getBaseName(ScanCollectionName) + "_NoRidges_" + key + ".png", IDNoRidges, UnIDNoRidges);
            GenerateHistogram(FilenameUtils.getFullPath(ScanCollectionName) + "/" + FilenameUtils.getBaseName(ScanCollectionName) + "_SymScore_" + key + ".png", IDSymScore, UnIDSymScore);
            GenerateHistogram(FilenameUtils.getFullPath(ScanCollectionName) + "/" + FilenameUtils.getBaseName(ScanCollectionName) + "_LDA_" + key + ".png", IDLDA, UnIDLDA);

        }
    }

    private void GenerateScoreLDAModelsByQ() throws IOException {
        ScoreModels = new HashMap<>();
        if (IDsummary != null) {
            for (int Q = 1; Q <= 2; Q++) {
                ScoreModels.put(String.valueOf(Q), new MS1ScoreModel());
            }

            for (PeakCluster cluster : PeakClusters) {
                MS1ScoreModel model = ScoreModels.get(String.valueOf(cluster.GetQualityCategory()));
                if (cluster.Identified) {
                    model.IDList.add(cluster);
                } else {
                    model.UnIDList.add(cluster);
                }
            }
            ArrayList<String> emptyList = new ArrayList<>();

            for (String key : ScoreModels.keySet()) {
                MS1ScoreModel model = ScoreModels.get(key);
                if (model.IDList.isEmpty()) {
                    Logger.getRootLogger().info("No. of identified cluster of Q " + key + " is zero");
                    emptyList.add(key);
                }
            }
            for (String key : emptyList) {
                ScoreModels.remove(key);
            }

            for (String key : ScoreModels.keySet()) {
                MS1ScoreModel model = ScoreModels.get(key);
                double[][] testdata = new double[7][model.IDList.size() + model.UnIDList.size()];
                double[] testgroup = new double[model.IDList.size() + model.UnIDList.size()];

                int idx = 0;
                for (PeakCluster cluster : model.IDList) {
                    if (cluster.GetQualityCategory() == 2) {
                        testdata[0][idx] = cluster.Corrs[0];
                    } else {
                        testdata[0][idx] = cluster.Corrs[0] + cluster.Corrs[1];
                    }
                    testdata[1][idx] = Math.log(cluster.GetSNR(0));
                    testdata[2][idx] = cluster.IsoMapProb;
                    testdata[3][idx] = cluster.GetConflictCorr();
                    testdata[4][idx] = cluster.GetApexVar();
                    testdata[5][idx] = cluster.NoRidges;
                    testdata[6][idx] = cluster.GetSymScore();
                    testgroup[idx] = 1;
                    idx++;
                }

                for (PeakCluster cluster : model.UnIDList) {
                    if (cluster.GetQualityCategory() == 2) {
                        testdata[0][idx] = cluster.Corrs[0];
                    } else {
                        testdata[0][idx] = cluster.Corrs[0] + cluster.Corrs[1];
                    }
                    testdata[1][idx] = Math.log(cluster.GetSNR(0));
                    testdata[2][idx] = cluster.IsoMapProb;
                    testdata[3][idx] = cluster.GetConflictCorr();
                    testdata[4][idx] = cluster.GetApexVar();
                    testdata[5][idx] = cluster.NoRidges;
                    testdata[6][idx] = cluster.GetSymScore();
                    testgroup[idx] = 0;
                    idx++;
                }
                DiscriminantAnalysis LDA = new DiscriminantAnalysis();
                int[] group = LDA.predictedGroup(testgroup, testdata, testdata);

                model.LDACoeff[0] = -LDA.linearDiscriminants[0][0];
                model.LDACoeff[1] = -LDA.linearDiscriminants[0][1];
                model.LDACoeff[2] = -LDA.linearDiscriminants[0][2];
                model.LDACoeff[3] = -LDA.linearDiscriminants[0][3];
                model.LDACoeff[4] = -LDA.linearDiscriminants[0][4];
                model.LDACoeff[5] = -LDA.linearDiscriminants[0][5];
                model.LDACoeff[6] = -LDA.linearDiscriminants[0][6];

                for (PeakCluster cluster : model.IDList) {
                    if (cluster.GetQualityCategory() == 2) {
                        cluster.MS1Score = (float) ((model.LDACoeff[0] * cluster.Corrs[0]) + (model.LDACoeff[1] * Math.log(cluster.GetSNR(0))) + (model.LDACoeff[2] * cluster.IsoMapProb) + (model.LDACoeff[3] * cluster.GetConflictCorr()) + (model.LDACoeff[4] * cluster.GetApexVar()) + (model.LDACoeff[5] * cluster.NoRidges) + (model.LDACoeff[6] * cluster.GetSymScore()));
                    } else {
                        cluster.MS1Score = (float) ((model.LDACoeff[0] * (cluster.Corrs[0] + cluster.Corrs[1])) + (model.LDACoeff[1] * Math.log(cluster.GetSNR(0))) + (model.LDACoeff[2] * cluster.IsoMapProb) + (model.LDACoeff[3] * cluster.GetConflictCorr()) + (model.LDACoeff[4] * cluster.GetApexVar()) + (model.LDACoeff[5] * cluster.NoRidges) + (model.LDACoeff[6] * cluster.GetSymScore()));
                    }
                    if (cluster.MS1Score > model.max) {
                        model.max = cluster.MS1Score;
                    }
                    if (cluster.MS1Score < model.min) {
                        model.min = cluster.MS1Score;
                    }
                }
                for (PeakCluster cluster : model.UnIDList) {
                    if (cluster.GetQualityCategory() == 2) {
                        cluster.MS1Score = (float) ((model.LDACoeff[0] * cluster.Corrs[0]) + (model.LDACoeff[1] * Math.log(cluster.GetSNR(0))) + (model.LDACoeff[2] * cluster.IsoMapProb) + (model.LDACoeff[3] * cluster.GetConflictCorr()) + (model.LDACoeff[4] * cluster.GetApexVar()) + (model.LDACoeff[5] * cluster.NoRidges) + (model.LDACoeff[6] * cluster.GetSymScore()));
                    } else {
                        cluster.MS1Score = (float) ((model.LDACoeff[0] * (cluster.Corrs[0] + cluster.Corrs[1])) + (model.LDACoeff[1] * Math.log(cluster.GetSNR(0))) + (model.LDACoeff[2] * cluster.IsoMapProb) + (model.LDACoeff[3] * cluster.GetConflictCorr()) + (model.LDACoeff[4] * cluster.GetApexVar()) + (model.LDACoeff[5] * cluster.NoRidges) + (model.LDACoeff[6] * cluster.GetSymScore()));
                    }
                    if (cluster.MS1Score > model.max) {
                        model.max = cluster.MS1Score;
                    }
                    if (cluster.MS1Score < model.min) {
                        model.min = cluster.MS1Score;
                    }
                }

                DecimalFormat df = new DecimalFormat("#.####");

                Logger.getRootLogger().debug("----------------------------------------------------------------------------------------");
                Logger.getRootLogger().debug("Linear discriminate analysis for clusters of quality category " + key );
                Logger.getRootLogger().debug("No. of identified peak clusters:" + model.IDList.size());
                Logger.getRootLogger().debug("No. of unidentified peak clusters:" + model.UnIDList.size());
                Logger.getRootLogger().debug("Trained weights:");
                Logger.getRootLogger().debug("Corr=" + df.format(model.LDACoeff[0]).toString() );                
                Logger.getRootLogger().debug("SNR=" + df.format(model.LDACoeff[1]).toString() );
                Logger.getRootLogger().debug("Isotopic pattern=" + df.format(model.LDACoeff[2]).toString() );
                Logger.getRootLogger().debug("ConflictCorr=" + df.format(model.LDACoeff[3]).toString() );
                Logger.getRootLogger().debug("ApexVar=" + df.format(model.LDACoeff[4]).toString() );
                Logger.getRootLogger().debug("NoRidges=" + df.format(model.LDACoeff[5]).toString() );
                Logger.getRootLogger().debug("SymScore=" + df.format(model.LDACoeff[6]).toString() );
            }
        }
    }

    private void GenerateScoreLDAModelsByChargeAndQ() throws IOException {

        ScoreModels = new HashMap<>();
        if (IDsummary != null) {
            for (int Q = 1; Q <= 2; Q++) {
                for (int charge = parameter.StartCharge; charge <= parameter.EndCharge; charge++) {
                    ScoreModels.put(Q + "_" + charge, new MS1ScoreModel());
                }
            }

            for (PeakCluster cluster : PeakClusters) {
                MS1ScoreModel model = ScoreModels.get(cluster.GetQualityCategory() + "_" + cluster.Charge);
                if (cluster.Identified) {
                    model.IDList.add(cluster);
                } else {
                    model.UnIDList.add(cluster);
                }
            }

            for (String key : ScoreModels.keySet()) {
                MS1ScoreModel model = ScoreModels.get(key);
                double[][] testdata = new double[7][model.IDList.size() + model.UnIDList.size()];
                double[] testgroup = new double[model.IDList.size() + model.UnIDList.size()];

                int idx = 0;
                for (PeakCluster cluster : model.IDList) {
                    if (cluster.GetQualityCategory() == 2) {
                        testdata[0][idx] = cluster.Corrs[0];
                    } else {
                        testdata[0][idx] = cluster.Corrs[0] + cluster.Corrs[1];
                    }
                    testdata[1][idx] = Math.log(cluster.GetSNR(0));
                    testdata[2][idx] = cluster.IsoMapProb;
                    testdata[3][idx] = cluster.GetConflictCorr();
                    testdata[4][idx] = cluster.GetApexVar();
                    testdata[5][idx] = cluster.NoRidges;
                    testdata[6][idx] = cluster.MS1Score;
                    testgroup[idx] = 1;
                    idx++;
                }
                for (PeakCluster cluster : model.UnIDList) {
                    if (cluster.GetQualityCategory() == 2) {
                        testdata[0][idx] = cluster.Corrs[0];
                    } else {
                        testdata[0][idx] = cluster.Corrs[0] + cluster.Corrs[1];
                    }
                    testdata[1][idx] = Math.log(cluster.GetSNR(0));
                    testdata[2][idx] = cluster.IsoMapProb;
                    testdata[3][idx] = cluster.GetConflictCorr();
                    testdata[4][idx] = cluster.GetApexVar();
                    testdata[5][idx] = cluster.NoRidges;
                    testdata[6][idx] = Math.log(cluster.GetSymScore());
                    testgroup[idx] = 0;
                    idx++;
                }
                DiscriminantAnalysis LDA = new DiscriminantAnalysis();
                int[] group = LDA.predictedGroup(testgroup, testdata, testdata);

                model.LDACoeff[0] = -LDA.linearDiscriminants[0][0];
                model.LDACoeff[1] = -LDA.linearDiscriminants[0][1];
                model.LDACoeff[2] = -LDA.linearDiscriminants[0][2];
                model.LDACoeff[3] = -LDA.linearDiscriminants[0][3];
                model.LDACoeff[4] = -LDA.linearDiscriminants[0][4];
                model.LDACoeff[5] = -LDA.linearDiscriminants[0][5];
                model.LDACoeff[6] = -LDA.linearDiscriminants[0][6];

                for (PeakCluster cluster : model.IDList) {
                    if (cluster.GetQualityCategory() == 2) {
                        cluster.MS1Score = (float) ((model.LDACoeff[0] * cluster.Corrs[0]) + (model.LDACoeff[1] * Math.log(cluster.GetSNR(0))) + (model.LDACoeff[2] * cluster.IsoMapProb) + (model.LDACoeff[3] * cluster.GetConflictCorr()) + (model.LDACoeff[4] * cluster.GetApexVar()) + (model.LDACoeff[5] * cluster.NoRidges) + (model.LDACoeff[6] * Math.log(cluster.GetSymScore())));
                    } else {
                        cluster.MS1Score = (float) ((model.LDACoeff[0] * (cluster.Corrs[0] + cluster.Corrs[1])) + (model.LDACoeff[1] * Math.log(cluster.GetSNR(0))) + (model.LDACoeff[2] * cluster.IsoMapProb) + (model.LDACoeff[3] * cluster.GetConflictCorr()) + (model.LDACoeff[4] * cluster.GetApexVar()) + (model.LDACoeff[5] * cluster.NoRidges) + (model.LDACoeff[6] * Math.log(cluster.GetSymScore())));
                    }
                }
                for (PeakCluster cluster : model.UnIDList) {
                    if (cluster.GetQualityCategory() == 2) {
                        cluster.MS1Score = (float) ((model.LDACoeff[0] * cluster.Corrs[0]) + (model.LDACoeff[1] * Math.log(cluster.GetSNR(0))) + (model.LDACoeff[2] * cluster.IsoMapProb) + (model.LDACoeff[3] * cluster.GetConflictCorr()) + (model.LDACoeff[4] * cluster.GetApexVar()) + (model.LDACoeff[5] * cluster.NoRidges) + (model.LDACoeff[6] * Math.log(cluster.GetSymScore())));
                    } else {
                        cluster.MS1Score = (float) ((model.LDACoeff[0] * (cluster.Corrs[0] + cluster.Corrs[1])) + (model.LDACoeff[1] * Math.log(cluster.GetSNR(0))) + (model.LDACoeff[2] * cluster.IsoMapProb) + (model.LDACoeff[3] * cluster.GetConflictCorr()) + (model.LDACoeff[4] * cluster.GetApexVar()) + (model.LDACoeff[5] * cluster.NoRidges) + (model.LDACoeff[6] * Math.log(cluster.GetSymScore())));
                    }
                }

                DecimalFormat df = new DecimalFormat("#.####");
                Logger.getRootLogger().debug("LDA weights:" + key );
                Logger.getRootLogger().debug("Corr(" + df.format(model.LDACoeff[0]).toString() + ")");
                //Logger.getRootLogger().debug("Corr2(" + df.format(LDACoeff[1]).toString() + "), ");
                Logger.getRootLogger().debug("SNR(" + df.format(model.LDACoeff[1]).toString() + ")");
                Logger.getRootLogger().debug("Isotopic pattern(" + df.format(model.LDACoeff[2]).toString() + ")");
                Logger.getRootLogger().debug("ConflictCorr (" + df.format(model.LDACoeff[3]).toString() + ")");
                Logger.getRootLogger().debug("RTVar(" + df.format(model.LDACoeff[4]).toString() + ")");
                Logger.getRootLogger().info("NoRidges(" + df.format(model.LDACoeff[5]).toString() + ")");
                Logger.getRootLogger().info("SymScore(" + df.format(model.LDACoeff[6]).toString() + ")");
            }
        }
    }

    public void ExportMappedIDQuant() throws SQLException, IOException {
        if (IDsummary == null) {
            return;
        }
        if (!FilenameUtils.getBaseName(IDsummary.mzXMLFileName).equals(FilenameUtils.getBaseName(ScanCollectionName))) {
            return;
        }
        IDsummary.ExportMappedPepID();
    }

    public void ExportID()throws SQLException, IOException{
        ExportID("");
    }
            
    public void ExportID(String tag) throws SQLException, IOException {
        if (IDsummary == null) {
            return;
        }
        if (!FilenameUtils.getBaseName(IDsummary.mzXMLFileName).equals(FilenameUtils.getBaseName(ScanCollectionName))) {
            return;
        }
        IDsummary.WriteLCMSIDSerialization(ParentmzXMLName,tag);
        //IDsummary.ExportPepID(connectionManager);
        //IDsummary.ExportProtID(connectionManager);
        //ExportMappedIDQuant();
    }

    private void GenerateHistogram(String pngfile, double[] IDvalues, double[] UnIDvalues) {

        //System.out.print("Generating histogram :" + pngfile + "\n");
        HistogramDataset histogramDataset = new HistogramDataset();
        histogramDataset.setType(HistogramType.RELATIVE_FREQUENCY);
        histogramDataset.addSeries("Identified clusters", IDvalues, 100);
        histogramDataset.addSeries("Unidentified clusters", UnIDvalues, 100);

        JFreeChart chart = ChartFactory.createHistogram(FilenameUtils.getBaseName(pngfile), "Score", "Frequency", histogramDataset, PlotOrientation.VERTICAL, true, false, false);
        XYPlot plot = chart.getXYPlot();
        plot.setBackgroundPaint(Color.white);
        plot.setDomainGridlinePaint(Color.white);
        plot.setRangeGridlinePaint(Color.white);
        plot.setForegroundAlpha(0.8f);
        chart.setBackgroundPaint(Color.white);
        try {
            ChartUtilities.saveChartAsPNG(new File(pngfile), chart, 1000, 600);
        } catch (IOException e) {
        }
    }

    public void GenerateMassCalibrationMWMap() throws IOException {
        //FileWriter writer = new FileWriter(FilenameUtils.getFullPath(ScanCollectionName) + "/" + FilenameUtils.getBaseName(ScanCollectionName) + "_masscali.txt");
        String pngfile = FilenameUtils.getFullPath(ScanCollectionName) + "/" + FilenameUtils.getBaseName(ScanCollectionName) + "_masscaliMW.png";
        XYSeries series = new XYSeries("PSM");
        XYSeriesCollection xySeriesCollection = new XYSeriesCollection();
        for (PSM psm : this.IDsummary.PSMList.values()) {
            float ppm = InstrumentParameter.CalcSignedPPM(psm.ObserPrecursorMass, psm.NeutralPepMass);
            series.add(new XYDataItem(psm.NeutralPepMass, ppm));
        }
//        writer.write("RT\tMW\tIntensity\tppm\n");
//        for (PepIonID pep : this.IDsummary.GetPepIonList().values()) {
//            float ppm = InstrumentParameter.CalcSignedPPM(pep.ObservedMass(), pep.CalcNeutralPepMass());
//            
//            writer.write(pep.GetRT() + "\t" + pep.CalcNeutralPepMass() + "\t" + pep.PeakHeight[0] + "\t" + ppm + "\n");
//        }
//        writer.close();
        xySeriesCollection.addSeries(series);
        JFreeChart chart = ChartFactory.createScatterPlot("Mass calibration", "Mass", "Mass error (ppm)", xySeriesCollection,
                PlotOrientation.VERTICAL, true, true, false);
        XYPlot xyPlot = (XYPlot) chart.getPlot();
        xyPlot.setDomainCrosshairVisible(true);
        xyPlot.setRangeCrosshairVisible(true);
        XYItemRenderer renderer = xyPlot.getRenderer();
        renderer.setSeriesPaint(0, Color.blue);
        renderer.setSeriesShape(0, new Ellipse2D.Double(0, 0, 3, 3));
        renderer.setSeriesStroke(0, new BasicStroke(1.0f));
        xyPlot.setBackgroundPaint(Color.white);
        ChartUtilities.saveChartAsPNG(new File(pngfile), chart, 1000, 600);
    }

    public void GenerateMassCalibrationRTMap() throws IOException {
        String pngfile = FilenameUtils.getFullPath(ScanCollectionName) + "/" + FilenameUtils.getBaseName(ScanCollectionName) + "_masscaliRT.png";
        XYSeries series = new XYSeries("PSM");
        XYSeriesCollection xySeriesCollection = new XYSeriesCollection();
        LoessInterpolator loessInterpolator = new LoessInterpolator(
                0.75,//bandwidth,
                2//robustnessIters
        );

        for (PSM psm : this.IDsummary.PSMList.values()) {
            float ppm = InstrumentParameter.CalcSignedPPM(psm.ObserPrecursorMass, psm.NeutralPepMass);            
            series.add(new XYDataItem(psm.RetentionTime, ppm));
        }
        double x[] = new double[IDsummary.PSMList.size()];
        double y[] = new double[x.length];
        double currentmin = 0f;
        for (int i = 0; i < series.getItemCount(); i++) {
            x[i] = (double) series.getX(i);
            if (x[i] <= currentmin) {
                x[i] = currentmin + 0.0001f;
            }
            currentmin = x[i];
            y[i] = (double) series.getY(i);
        }

        Masscalibrationfunction = loessInterpolator.interpolate(x, y);
        XYSeries smoothline = new XYSeries("Loess Regression");

        double xvalue = series.getMinX();

        while (xvalue < series.getMaxX()) {
            smoothline.add(xvalue, Masscalibrationfunction.value(xvalue));
            xvalue += 0.05d;
        }
        xySeriesCollection.addSeries(smoothline);
        xySeriesCollection.addSeries(series);

        JFreeChart chart = ChartFactory.createScatterPlot("Mass calibration", "RT", "Mass error (ppm)", xySeriesCollection,
                PlotOrientation.VERTICAL, true, true, false);
        XYPlot xyPlot = (XYPlot) chart.getPlot();
        xyPlot.setDomainCrosshairVisible(true);
        xyPlot.setRangeCrosshairVisible(true);

        XYItemRenderer renderer = xyPlot.getRenderer();
        renderer.setSeriesPaint(1, Color.blue);
        renderer.setSeriesPaint(0, Color.BLACK);
        renderer.setSeriesShape(1, new Ellipse2D.Double(0, 0, 3, 3));
        renderer.setSeriesStroke(1, new BasicStroke(1.0f));
        xyPlot.setBackgroundPaint(Color.white);
        ChartUtilities.saveChartAsPNG(new File(pngfile), chart, 1000, 600);
    }
   
    public void ReplaceProtByRefID(LCMSID protID) {
        IDsummary.GenerateProteinByRefIDByPepSeq(protID, UseMappedIon);
        //IDsummary.GenearteAssignIonList();        
    }

    public void iProphet(ArrayList<MSMSDBSearch> searches) throws InterruptedException, IOException {
        iProphet ipro = new iProphet(GetIPROPHETPepXML(), searches);
        ipro.DoiProphetPepXML();
        ipro.DoiProphetProteinXML();
    }
    public String GetIPROPHETPepXML(){
        return FilenameUtils.getFullPath(ParentmzXMLName) + "interact-"+FilenameUtils.getBaseName(ParentmzXMLName) + ".iproph.pep.xml";
    }
    
    public String GetIPROPHETProtXML(){
        return FilenameUtils.getFullPath(ParentmzXMLName) + "interact-"+FilenameUtils.getBaseName(ParentmzXMLName) + ".iproph.prot.xml";
    }

    public void ParseiProphet(DBSearchParam param) throws ParserConfigurationException, IOException, SAXException, XmlPullParserException, ClassNotFoundException, InterruptedException {
        ParseiProphet(param,null);
    }
    
    public void ParseiProphet(DBSearchParam param,LCMSID RefPepID) throws ParserConfigurationException, IOException, SAXException, XmlPullParserException, ClassNotFoundException, InterruptedException {
        IDsummary = new LCMSID(FilenameUtils.getFullPath(ParentmzXMLName) + FilenameUtils.getBaseName(ParentmzXMLName),param.DecoyPrefix,param.FastaPath);                
        TPPResult tppresult = new TPPResult(param.PepFDR,param.ProtFDR, param.DecoyPrefix);
        if (RefPepID == null) {
            tppresult.ReadSearchResult(IDsummary, GetIPROPHETPepXML(), GetIPROPHETProtXML());
        } else {
            tppresult.ReadSearchResultByRefPepID(IDsummary, GetIPROPHETPepXML(), GetIPROPHETProtXML(), RefPepID);
        }
        CheckPSMRT();
    }
    
    public void ParsePepXML(DBSearchParam param, float prob) throws ParserConfigurationException, IOException, SAXException, XmlPullParserException, ClassNotFoundException, InterruptedException {
        IDsummary = new LCMSID(FilenameUtils.getFullPath(ParentmzXMLName) + FilenameUtils.getBaseName(ParentmzXMLName),param.DecoyPrefix,param.FastaPath);                
        PepXMLParser pepxmlparser = new PepXMLParser(IDsummary, param.InteractPepXMLPath, prob);
        CheckPSMRT();
    }
    public void ParseiProphetPepXML(DBSearchParam param, float prob) throws ParserConfigurationException, IOException, SAXException, XmlPullParserException, ClassNotFoundException, InterruptedException {
        IDsummary = new LCMSID(FilenameUtils.getFullPath(ParentmzXMLName) + FilenameUtils.getBaseName(ParentmzXMLName),param.DecoyPrefix,param.FastaPath);                
        PepXMLParser pepxmlparser = new PepXMLParser(IDsummary, GetIPROPHETPepXML(), prob);
        CheckPSMRT();
    }

    private void CheckPSMRT() {
        for(PSM psm : IDsummary.PSMList.values()){
            if(psm.RetentionTime==-1f){
                psm.RetentionTime=GetmzXML().GetScanElutionTimeMap().get(psm.ScanNo);
                psm.NeighborMaxRetentionTime=psm.RetentionTime;
            }
        }
    }

    public void RemoveContaminantPeaks(float proportion) {
        Logger.getRootLogger().info("Removing peak clusters whose m/z apprear more than " +proportion*100+ "% chromatography. No. of peak clusters : "+PeakClusters.size());
        float minmz=Float.MAX_VALUE;
        float maxmz=Float.MIN_VALUE;
        float minrt=Float.MAX_VALUE;
        float maxrt=Float.MIN_VALUE;
        for(PeakCluster peak : PeakClusters){
            if(peak.TargetMz()>maxmz){
                maxmz=peak.TargetMz();
            }
            if(peak.TargetMz()<minmz){
                minmz=peak.TargetMz();
            }
            if(peak.endRT>maxrt){
                maxrt=peak.endRT;
            }
            if(peak.startRT<minrt){
                minrt=peak.startRT;
            }
        }
        HashMap<Integer,ArrayList<PeakCluster>> map=new HashMap<>();
        
        float[] MzBin=new float[(int)Math.ceil((maxmz-minmz)*10)+1];       
        for(PeakCluster peak : PeakClusters){
            int binkey=(int)Math.ceil((peak.TargetMz()-minmz)*10);
            MzBin[binkey]+=peak.endRT-peak.startRT;
            if(!map.containsKey(binkey)){
                map.put(binkey, new ArrayList<PeakCluster>());
            }
            map.get(binkey).add(peak);
        }
        float threshold=proportion*(maxrt-minrt);
        for(int i=0;i<MzBin.length;i++){
            if(MzBin[i]>threshold){
                for (PeakCluster peakCluster : map.get(i)){ 
                    PeakClusters.remove(peakCluster);
                    Logger.getRootLogger().debug("Removing the cluster m/z: "+ peakCluster.TargetMz()+", StartRT: "+ peakCluster.startRT +", EndRT: "+peakCluster.endRT);
                }
            }
        }
        Logger.getRootLogger().info("Remaining peak clusters : "+PeakClusters.size());
    }
}
