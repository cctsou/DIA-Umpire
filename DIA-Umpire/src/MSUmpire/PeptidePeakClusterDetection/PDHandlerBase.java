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
package MSUmpire.PeptidePeakClusterDetection;

import MSUmpire.BaseDataStructure.InstrumentParameter;
import MSUmpire.BaseDataStructure.ScanCollection;
import MSUmpire.BaseDataStructure.ScanData;
import MSUmpire.BaseDataStructure.XYData;
import MSUmpire.BaseDataStructure.XYPointCollection;
import MSUmpire.BaseDataStructure.XYZData;
import MSUmpire.LCMSBaseStructure.LCMSPeakBase;
import MSUmpire.PeakDataStructure.PeakCluster;
import MSUmpire.PeakDataStructure.PeakCurve;
import java.io.*;
import java.util.*;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;
import net.sf.javaml.core.kdtree.KDTree;
import net.sf.javaml.core.kdtree.KeySizeException;
import org.apache.commons.lang.exception.ExceptionUtils;
import org.apache.log4j.Logger;

/**
 * Peak detection processing parent class
 * @author Chih-Chiang Tsou <chihchiang.tsou@gmail.com>
 */
public class PDHandlerBase {
    
    protected HashSet<String> IncludedHashMap;
    protected HashMap<XYData,Boolean> InclusionFound=new HashMap();
    protected XYPointCollection InclusionRT=new XYPointCollection();
    protected KDTree InclusionRange=new KDTree(2);    
    protected int NoCPUs = 4;
    public float minSNR;
    public TreeMap<Float, XYData>[] IsotopePatternMap;
    public TreeMap<Float, XYData>[] IsotopePatternFragMap;
    protected LCMSPeakBase LCMSPeakBase;
    protected InstrumentParameter parameter;
    protected boolean ReleaseScans = true;
    protected boolean TargetedOnly = false;
    protected float PPM;
    public int MSlevel=1;

    public PDHandlerBase() {
    }

    public void SetTargetedDetectionOnly(){
        TargetedOnly=true;
    }
    
    protected String InclusionCheckInfo(){
        int count=0;
        
        for(XYData point : InclusionFound.keySet()){
            boolean value = InclusionFound.get(point);
            if(value){
                count++;
            }
            else{
                //Logger.getRootLogger().warn("Missing signals: mz="+point.getX()+", RT="+point.getY());
            }
        }
        return count+"/"+InclusionFound.size();                
    }
    
    //Add mz and RT coordinate to inclusion list
     public void AddToInclusionList(float mz, float rt){
         XYData point=new XYData(mz, rt);
        InclusionFound.put(point, false);
        InclusionRT.AddPoint(rt,mz);        
        try {
            InclusionRange.insert(new double[]{rt,mz}, point);
        } catch (Exception ex) {
            Logger.getRootLogger().error(ExceptionUtils.getStackTrace(ex));
        }
    }
          
    public void ClearAllPeaks() {
        LCMSPeakBase.BaseClearAllPeaks();
    }

    public void ClearRawPeaks() {
        for (PeakCurve peakCurve : LCMSPeakBase.UnSortedPeakCurves) {
            peakCurve.CalculateMzVar();
            peakCurve.StartRT();
            peakCurve.EndRT();
            peakCurve.ReleaseRawPeak();
        }
    }

    //Detect all m/z trace peak curves given a list of ScanCollection
    public void FindAllMzTracePeakCurvesForScanCollections(ArrayList<ScanCollection> scanCollections) throws IOException {
        //Read peptide isotope pattern table
        ReadPepIsoMS1PatternMap();
        LCMSPeakBase.UnSortedPeakCurves = new ArrayList<>();
        
        for (ScanCollection scanCollection : scanCollections) {
            //Detect mz trace peak curves for each ScanCollection
            FindAllMzTracePeakCurves(scanCollection);
        }
        Logger.getRootLogger().info("Inclusion mz values found: "+InclusionCheckInfo());        
        //Perform peak smoothing for each detected peak curve
        PeakCurveSmoothing();      
        ClearRawPeaks();
    }
    
    //Detect all m/z trace / peak curves
    protected void FindAllMzTracePeakCurves(ScanCollection scanCollection) throws IOException {

        IncludedHashMap = new HashSet<>();
        Logger.getRootLogger().info("Processing all scans to detect possible m/z peak curves....");
        
        float preRT = 0f;
        
        //Loop for each scan in the ScanCollection
        for (int idx = 0; idx < scanCollection.GetScanNoArray(MSlevel).size(); idx++) {
            int scanNO = scanCollection.GetScanNoArray(MSlevel).get(idx);
            ScanData scanData = scanCollection.GetScan(scanNO);
            
            //If we are doing targeted peak detection and the RT of current scan is not in the range of targeted list, jump to the next scan 
            if(TargetedOnly && !FoundInInclusionRTList(scanData.RetentionTime)){
                continue;
            }
            if (idx == 0) {
                preRT = scanData.RetentionTime - 0.01f;
            }
            for (int i = 0; i < scanData.PointCount(); i++) {
                XYData peak = scanData.Data.get(i);
                //If we are doing targeted peak detection and the RT and m/z of current peak is not in the range of targeted list, jump to the next peak 
                if (TargetedOnly && !FoundInInclusionMZList(scanData.RetentionTime,peak.getX())) {
                    continue;
                }
                
                if(peak.getX()<parameter.MinMZ){
                    continue;
                }
                
                //Check if the current peak has been included in previously developed peak curves
                if (!IncludedHashMap.contains(scanNO + "_" + peak.getX())) {//The peak hasn't been included
                   
                    //The current peak will be the starting peak of a new peak curve
                    //Add it to the hash table
                    IncludedHashMap.add(scanNO + "_" + peak.getX());

                    float startmz = peak.getX();
                    float startint = peak.getY();
                    
                   //Find the maximum peak within PPM window as the starting peak
                    for (int j = i + 1; j < scanData.PointCount(); j++) {
                        XYData currentpeak = scanData.Data.get(j);
                        if (!IncludedHashMap.contains(scanNO + "_" + currentpeak.getX())) {
                            if (InstrumentParameter.CalcPPM(currentpeak.getX(), startmz) <= PPM) {
                                IncludedHashMap.add(scanNO + "_" + currentpeak.getX());

                                if (currentpeak.getY() >= startint) {
                                    startmz = currentpeak.getX();
                                    startint = currentpeak.getY();
                                }
                            } else {
                                break;
                            }
                        }
                    }

                    //Initialize a new peak curve
                    PeakCurve Peakcurve = new PeakCurve(parameter);
                    //Add a background peak
                    Peakcurve.AddPeak(new XYZData(preRT, startmz, scanData.background));
                    //Add the starting peak
                    Peakcurve.AddPeak(new XYZData(scanData.RetentionTime, startmz, startint));
                    Peakcurve.StartScan = scanNO;

                    int missedScan = 0;
                    float endrt=scanData.RetentionTime;
                    int endScan=scanData.ScanNum;
                    float bk=0f;

                     //Starting from the next scan, find the following peaks given the starting peak
                    for (int idx2 = idx + 1; idx2 < scanCollection.GetScanNoArray(MSlevel).size() && (missedScan < parameter.NoMissedScan /*|| (TargetedOnly && Peakcurve.RTWidth()<parameter.MaxCurveRTRange)*/); idx2++) {
                        int scanNO2 = scanCollection.GetScanNoArray(MSlevel).get(idx2);
                        ScanData scanData2 = scanCollection.GetScan(scanNO2);

                        endrt=scanData2.RetentionTime;
                        endScan=scanData2.ScanNum;
                        bk=scanData2.background;
                        float currentmz = 0f;
                        float currentint = 0f;

                        //If the scan is empty
                        if (scanData2.PointCount() == 0) {                            
                            if (parameter.FillGapByBK) {
                                Peakcurve.AddPeak(new XYZData(scanData2.RetentionTime, Peakcurve.TargetMz, scanData2.background));
                            }
                            missedScan++;
                            continue;
                        }

                        //Find the m/z index 
                        int mzidx = scanData2.GetLowerIndexOfX(Peakcurve.TargetMz);
                        for (int pkidx = mzidx; pkidx < scanData2.Data.size(); pkidx++) {
                            XYData currentpeak = scanData2.Data.get(pkidx);
                            if (currentpeak.getX() < parameter.MinMZ) {
                                continue;
                            }
                            //Check if the peak has been included or not
                            if (!IncludedHashMap.contains(scanNO2 + "_" + currentpeak.getX())) {
                                if (InstrumentParameter.CalcPPM(currentpeak.getX(), Peakcurve.TargetMz) > PPM) {
                                    if (currentpeak.getX() > Peakcurve.TargetMz) {
                                        break;
                                    }
                                } else {
                                    //////////The peak is in the ppm window, select the highest peak
                                    IncludedHashMap.add(scanNO2 + "_" + currentpeak.getX());
                                    if (currentint < currentpeak.getY()) {
                                        currentmz = currentpeak.getX();
                                        currentint = currentpeak.getY();
                                    }
                                }
                            }
                        }
                        
                        //No peak in the PPM window has been found
                        if (currentmz == 0f) {
                            if (parameter.FillGapByBK) {
                                Peakcurve.AddPeak(new XYZData(scanData2.RetentionTime, Peakcurve.TargetMz, scanData2.background));
                            }
                            missedScan++;
                        } else {
                            missedScan = 0;
                            Peakcurve.AddPeak(new XYZData(scanData2.RetentionTime, currentmz, currentint));
                        }
                    }
                    Peakcurve.AddPeak(new XYZData(endrt, Peakcurve.TargetMz, bk));
                    Peakcurve.EndScan=endScan;

                    //First check if the peak curve is in targeted list
                    if (FoundInInclusionList(Peakcurve.TargetMz, Peakcurve.StartRT(), Peakcurve.EndRT())) {
                        LCMSPeakBase.UnSortedPeakCurves.add(Peakcurve);
                    //Then check if the peak curve passes the criteria
                    } else if (Peakcurve.GetRawSNR() > LCMSPeakBase.SNR && Peakcurve.GetPeakList().size() >= parameter.MinPeakPerPeakCurve + 2) {
                        LCMSPeakBase.UnSortedPeakCurves.add(Peakcurve);
                    } else {
                        Peakcurve = null;
                    }
                }
            }
            preRT = scanData.RetentionTime;
            if (ReleaseScans) {
                scanData.dispose();
            }
        }

        //System.out.print("PSM removed (PeakCurve generation):" + PSMRemoved );         
        IncludedHashMap.clear();
        IncludedHashMap = null;

        int i = 1;
        //Assign peak curve index
        for (PeakCurve peakCurve : LCMSPeakBase.UnSortedPeakCurves) {
            peakCurve.Index = i++;
        }
        
        System.gc();
        Logger.getRootLogger().info(LCMSPeakBase.UnSortedPeakCurves.size() + " Peak curves found (Memory usage:" + Math.round((Runtime.getRuntime().totalMemory() - Runtime.getRuntime().freeMemory()) / 1048576) + "MB)");
    }
        
    private boolean FoundInInclusionRTList(float rt){              
        return Math.abs(InclusionRT.Data.get(InclusionRT.GetClosetIndexOfX(rt)).getX()-rt)<parameter.MaxCurveRTRange;
    }
    
    private boolean FoundInInclusionMZList(float rt, float mz) {
        if(InclusionRT.PointCount()==0){
            return false;
        }
        float lowrt = rt - parameter.MaxCurveRTRange;
        float highrt = rt + parameter.MaxCurveRTRange;
        float lowmz = InstrumentParameter.GetMzByPPM(mz, 1, PPM);
        float highmz = InstrumentParameter.GetMzByPPM(mz, 1, -PPM);
       
        Object[] found=null;
        try {
            found = InclusionRange.range(new double[]{lowrt,lowmz}, new double[]{highrt,highmz});
        } catch (KeySizeException ex) {
            Logger.getRootLogger().error(ExceptionUtils.getStackTrace(ex));
        }
        if(found!=null && found.length>0){
            return true;
        }
        return false;        
    }
    
    private boolean FoundInInclusionList(float mz, float startrt, float endrt){       
        if(InclusionRT.PointCount()==0){
            return false;
        }
        float lowmz = InstrumentParameter.GetMzByPPM(mz, 1, PPM);
        float highmz = InstrumentParameter.GetMzByPPM(mz, 1, -PPM);
        float lowrt=startrt-parameter.RTtol ;
        float highrt=endrt+parameter.RTtol;
                
        Object[] found=null;
        try {
            found = InclusionRange.range(new double[]{lowrt,lowmz}, new double[]{highrt,highmz});
        } catch (KeySizeException ex) {
            Logger.getRootLogger().error(ExceptionUtils.getStackTrace(ex));
        }
        if(found!=null && found.length>0){
            for(Object point: found){
                InclusionFound.put((XYData) point,true);
            }
            return true;
        }
        return false;     
    }
    
    //Signal smoothing for each detected peak curve
    protected void PeakCurveSmoothing() {
        //System.out.print("Using multithreading now: " + NoCPUs + " processors");
        Logger.getRootLogger().info("Smoothing detected signals......");
        
        //Threading pool
        ExecutorService executorPool = null;
        executorPool = Executors.newFixedThreadPool(NoCPUs);
        ArrayList<PeakCurveSmoothingUnit> ResultList = new ArrayList<>();
        for (PeakCurve Peakcurve : LCMSPeakBase.UnSortedPeakCurves) {
            PeakCurveSmoothingUnit unit = new PeakCurveSmoothingUnit(Peakcurve, parameter);
            ResultList.add(unit);

            executorPool.execute(unit);
        }
        executorPool.shutdown();
        try {
            executorPool.awaitTermination(Long.MAX_VALUE, TimeUnit.NANOSECONDS);
        } catch (InterruptedException e) {
            Logger.getRootLogger().info("interrupted..");
        }
        executorPool = null;
        
        LCMSPeakBase.UnSortedPeakCurves.clear();
        for (PeakCurveSmoothingUnit result : ResultList) {
            LCMSPeakBase.UnSortedPeakCurves.addAll(result.ResultCurves);
        }
        
        int i = 1;
        for (PeakCurve peakCurve : LCMSPeakBase.UnSortedPeakCurves) {            
            peakCurve.Index = i++;
        }        
    }

    //Load pre-built peptide isotope pattern table
    protected void ReadPepIsoMS1PatternMap() throws FileNotFoundException, IOException {

        InputStream is = this.getClass().getClassLoader().getResourceAsStream("resource/IsotopicPatternRange.csv");
        BufferedReader reader = new BufferedReader(new InputStreamReader(is));
        IsotopePatternMap = new TreeMap[LCMSPeakBase.MaxNoPeakCluster - 1];
        for (int i = 0; i < IsotopePatternMap.length; i++) {
            IsotopePatternMap[i] = new TreeMap<>();
        }
        String line = "";
        while ((line = reader.readLine()) != null) {
            float MW = Float.parseFloat(line.split(",")[0]);

            for (int i = 0; i < LCMSPeakBase.MaxNoPeakCluster - 1; i++) {
                float Mean = Float.parseFloat(line.split(",")[1 + (i * 2)]);
                float SD = Float.parseFloat(line.split(",")[2 + (i * 2)]);

                if (!Float.isNaN(Mean)) {
                    IsotopePatternMap[i].put(MW, new XYData(Mean + 3.3f * SD, Mean - 3.3f * SD));
                }
            }
        }
        reader.close();
    }
    
    //Group peak curves based on peak profile correlation of isotope peaks
    protected void PeakCurveCorrClustering(XYData mzRange) throws IOException {
        Logger.getRootLogger().info("Grouping isotopic peak curves........");

        LCMSPeakBase.PeakClusters = new ArrayList<>();
        
        //Thread pool
        ExecutorService executorPool = null;
        executorPool = Executors.newFixedThreadPool(NoCPUs);        
        ArrayList<PeakCurveClusteringCorrKDtree> ResultList = new ArrayList<>();

        //For each peak curve
        for (PeakCurve Peakcurve : LCMSPeakBase.UnSortedPeakCurves) {
            if (Peakcurve.TargetMz >= mzRange.getX() && Peakcurve.TargetMz <= mzRange.getY()) {
                //Create a thread unit for doing isotope clustering given a peak curve as the monoisotope peak
                PeakCurveClusteringCorrKDtree unit = new PeakCurveClusteringCorrKDtree(Peakcurve, LCMSPeakBase.GetPeakCurveSearchTree(), parameter, IsotopePatternMap, LCMSPeakBase.StartCharge, LCMSPeakBase.EndCharge, LCMSPeakBase.MaxNoPeakCluster, LCMSPeakBase.MinNoPeakCluster);
                ResultList.add(unit);
                executorPool.execute(unit);                
            }
        }

        executorPool.shutdown();
        try {
            executorPool.awaitTermination(Long.MAX_VALUE, TimeUnit.NANOSECONDS);
        } catch (InterruptedException e) {
            Logger.getRootLogger().info("interrupted..");
        }
        for (PeakCurveClusteringCorrKDtree unit : ResultList) {
            for (PeakCluster peakCluster : unit.ResultClusters) {
                //Check if the monoistope peak of cluster has been grouped in other isotope cluster, if yes, remove the peak cluster
                if (!parameter.RemoveGroupedPeaks || !peakCluster.MonoIsotopePeak.ChargeGrouped.contains(peakCluster.Charge)) {
                    peakCluster.Index = LCMSPeakBase.PeakClusters.size() + 1;
                    peakCluster.GetConflictCorr();
                    LCMSPeakBase.PeakClusters.add(peakCluster);
                }               
            }
        }
        ResultList.clear();
        ResultList = null;
        System.gc();
        Logger.getRootLogger().info("No of ion clusters:" + LCMSPeakBase.PeakClusters.size() + " (Memory usage:" + Math.round((Runtime.getRuntime().totalMemory() - Runtime.getRuntime().freeMemory()) / 1048576) + "MB)");
    }
}
