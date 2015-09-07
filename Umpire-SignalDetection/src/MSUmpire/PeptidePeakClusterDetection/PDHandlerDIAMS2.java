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

import MSUmpire.BaseDataStructure.ScanCollection;
import MSUmpire.BaseDataStructure.ScanData;
import MSUmpire.BaseDataStructure.SortedListInteger;
import MSUmpire.BaseDataStructure.XYData;
import MSUmpire.PeakDataStructure.PeakCluster;
import MSUmpire.LCMSBaseStructure.LCMSPeakMS1;
import MSUmpire.LCMSBaseStructure.LCMSPeakDIAMS2;
import MSUmpire.DIA.CorrCalcCluster2CurveUnit;
import MSUmpire.PeakDataStructure.PeakCurve;
import MSUmpire.PeakDataStructure.SortedCurveCollectionIntensity;
import MSUmpire.PeakDataStructure.SortedCurveCollectionApexRT;
import MSUmpire.DIA.CorrCalcCluster2ClusterUnit;
import MSUmpire.DIA.FragDirectedGroupingUnit;
import MSUmpire.PeakDataStructure.PrecursorFragmentPairEdge;
import Utility.UpdateProcess;
import java.io.*;
import java.sql.SQLException;
import java.util.*;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;
import org.apache.log4j.Logger;

/**
 *
 * @author Chih-Chiang Tsou <chihchiang.tsou@gmail.com>
 */
public class PDHandlerDIAMS2 extends PDHandlerBase {

    private LCMSPeakMS1 ms1lcms;
    private XYData DIAWindowMz;

    public PDHandlerDIAMS2(LCMSPeakDIAMS2 swathWindow, int NoCPUs, LCMSPeakMS1 ms1lcms, float PPM){
        this.ms1lcms = ms1lcms;
        this.PPM = PPM;
        this.NoCPUs = NoCPUs;
        this.DIAWindowMz = swathWindow.DIA_MZ_Range;
        this.LCMSPeakBase = swathWindow;
        this.parameter = swathWindow.parameter;
        //this.connectionManager = swathWindow.connectionManager;
    }

    public void pSMARTGrouping(ScanCollection scanCollection) throws FileNotFoundException, IOException {
        ((LCMSPeakDIAMS2) LCMSPeakBase).FragmentsClu2Cur = new HashMap<>();
        for (PeakCluster peakCluster : ms1lcms.PeakClusters) {
            if (peakCluster.GetMaxMz()>= DIAWindowMz.getX() && peakCluster.TargetMz() <= DIAWindowMz.getY()) {
                ScanCollection SearchScans = scanCollection.GetSubCollectionByElutionTimeAndMZ(peakCluster.startRT, peakCluster.endRT, -1, -1, 2, false);
                ArrayList<PrecursorFragmentPairEdge> ResultList = new ArrayList<>();
                for (ScanData scan : SearchScans.ScanHashMap.values()) {
                    for (int i = 0; i < scan.PointCount(); i++) {
                        XYData peak = scan.Data.get(i);
                        PrecursorFragmentPairEdge PrecursorFragmentPair = new PrecursorFragmentPairEdge();
                        PrecursorFragmentPair.PeakCurveIndexA = peakCluster.Index;
                        PrecursorFragmentPair.PeakCurveIndexB = scan.Num;
                        PrecursorFragmentPair.FragmentMz = peak.getX();
                        PrecursorFragmentPair.Intensity = peak.getY();
                        PrecursorFragmentPair.RTOverlapP = 1f;
                        PrecursorFragmentPair.ApexDelta = Math.abs(peakCluster.MonoIsotopePeak.ApexRT - scan.RetentionTime);
                        float rtrange = peakCluster.endRT - peakCluster.startRT;
                        PrecursorFragmentPair.Correlation = (rtrange - PrecursorFragmentPair.ApexDelta) / rtrange;
                        //FragmentPeaks.put(peakCurve.Index, peakCurve);
                        ResultList.add(PrecursorFragmentPair);
                    }
                }
                ((LCMSPeakDIAMS2) LCMSPeakBase).FragmentsClu2Cur.put(peakCluster.Index, ResultList);
            }
        }
        
        ((LCMSPeakDIAMS2) LCMSPeakBase).ExportCluster2CurveCorr();
    }
    
    public void DetectPeakCurves(ScanCollection scanCollection) throws FileNotFoundException, IOException {        
        //ReadFragIsoPatternMap();
        LCMSPeakBase.UnSortedPeakCurves = new ArrayList<>();
        FindAllPeakCurve(scanCollection);
        //PeaksCheckingStage1AND2(1);
        WaveletDetectMax();
        //CreateSortedPeakCurveList();
        ClearRawPeaks();
        ReadPepIsoMS1PatternMap();
        PeakCurveCorrClustering_V2(DIAWindowMz);
    }
        
    public void FragmentGrouping() throws SQLException, IOException {
        PrecursorFragmentPairBuildingForMS1();
        PrecursorFragmentPairBuildingForUnfragmentedIon();
        //PrecursorFragmentPairBuildingForIsolatedPeakCurve();
    }

    public void FragmentDirectedGrouping() {

        SortedCurveCollectionIntensity fragmentCurveCollectionIntensity = new SortedCurveCollectionIntensity();
        fragmentCurveCollectionIntensity.addAll(LCMSPeakBase.UnSortedPeakCurves);

        ArrayList<FragDirectedGroupingUnit> ResultList = new ArrayList<>();
        ExecutorService executorPool = null;
        executorPool = Executors.newFixedThreadPool(NoCPUs);
        float RTtol = 0.5f;
        //UpdateProcess progress = new UpdateProcess();

        while (!fragmentCurveCollectionIntensity.isEmpty()) {
            SortedListInteger RTList = new SortedListInteger();
            SortedCurveCollectionApexRT RTCurveList = new SortedCurveCollectionApexRT();
            RTCurveList.addAll(LCMSPeakBase.UnSortedPeakCurves);

            for (PeakCurve fragment : fragmentCurveCollectionIntensity) {
                if (Math.abs(RTList.get(RTList.BinarySearchClosest(fragment.ApexRT)) - fragment.ApexRT) > 2 * RTtol) {
                    FragDirectedGroupingUnit groupunit = new FragDirectedGroupingUnit(fragment, RTCurveList, parameter, null);
                    ResultList.add(groupunit);
                    executorPool.execute(groupunit);
                }
            }

            //Thread thread = new Thread(progress);
            //thread.start();

            executorPool.shutdown();
//            while (!executorPool.isTerminated()) {
//            }
            try {
                executorPool.awaitTermination(Long.MAX_VALUE, TimeUnit.NANOSECONDS);
            } catch (InterruptedException e) {
                Logger.getRootLogger().info("interrupted..");
            }

            //thread = null;
            //progress.ClearMSG();
            //progress = null;
        }
    }

    protected void PeakCurveCorrClustering_IsolationOnly() throws SQLException, IOException {
        Logger.getRootLogger().info("Grouping unfragmented precursor ion peaks in SWATH isolation window........");

        if (LCMSPeakBase.PeakClusters == null) {
            LCMSPeakBase.PeakClusters = new ArrayList<>();
        }

        ExecutorService executorPool = null;
        executorPool = Executors.newFixedThreadPool(NoCPUs);
        //executorPool = Executors.newFixedThreadPool(1);
        ArrayList<PeakCurveClusteringCorrKDtree> ResultList = new ArrayList<>();

        //UpdateProcess progress = new UpdateProcess();
        UpdateProcess progress = null;

        //progress.SetTotal(LCMSPeakBase.PeakCurveListMZ.size());
        //Thread thread = new Thread(progress);
        //thread.start();

        for (PeakCurve peakCurve : LCMSPeakBase.UnSortedPeakCurves) {
            if (peakCurve.TargetMz >= DIAWindowMz.getX() && peakCurve.TargetMz <= DIAWindowMz.getY()) {
                //PeakCurveClusteringCorrV2Unit unit = new PeakCurveClusteringCorrV2Unit(peakCurve, LCMSPeakBase.PeakCurveListMZ, LCMSPeakBase.PeakCurveListRT, parameter, IsotopePatternMap, LCMSPeakBase.StartCharge, LCMSPeakBase.EndCharge, LCMSPeakBase.MaxNoPeakCluster, LCMSPeakBase.MinNoPeakCluster, progress);
                PeakCurveClusteringCorrKDtree unit = new PeakCurveClusteringCorrKDtree(peakCurve, LCMSPeakBase.GetPeakCurveSearchTree(), parameter, IsotopePatternMap, LCMSPeakBase.StartCharge, LCMSPeakBase.EndCharge, LCMSPeakBase.MaxNoPeakCluster, LCMSPeakBase.MinNoPeakCluster, progress);
                ResultList.add(unit);
                executorPool.execute(unit);
            }
        }

        executorPool.shutdown();
//        while (!executorPool.isTerminated()) {
//        }
        try {
            executorPool.awaitTermination(Long.MAX_VALUE, TimeUnit.NANOSECONDS);
        } catch (InterruptedException e) {
            Logger.getRootLogger().info("interrupted..");
        }

        //thread = null;
        //progress.ClearMSG();
        //progress = null;

        //ExportCorrMatrixToDB(GroupedFragmentList);
        for (PeakCurveClusteringCorrKDtree unit : ResultList) {
            for (PeakCluster peakCluster : unit.ResultClusters) {
//                if (!LCMSPeakBase.PeakClusters.containsKey(LCMSPeakBase.GetPeakClusterHashKey(peakCluster))) {
//                    LCMSPeakBase.PeakClusters.put(LCMSPeakBase.GetPeakClusterHashKey(peakCluster), peakCluster);
                peakCluster.Index = LCMSPeakBase.PeakClusters.size() + 1;
                LCMSPeakBase.PeakClusters.add(peakCluster);
                //}
            }
        }
        Logger.getRootLogger().info("No of peak clusters:" + LCMSPeakBase.PeakClusters.size() + "\n");

        //ExportCorrMatrix_V2(GroupedFragmentList);
        ResultList.clear();
        ResultList = null;
        //System.out.print("Finished multithreading\n");    
    }

    private void PrecursorFragmentPairBuildingForUnfragmentedIon() throws SQLException, IOException {
        //System.out.print("Using multithreading now: " + NoCPUs + " processors\n");
        Logger.getRootLogger().info("Building precursor-fragment pairs for unfragmented ions....");
        ExecutorService executorPool = null;
        ArrayList<CorrCalcCluster2CurveUnit> UnfragmentedIonPairList = new ArrayList<>();
        executorPool = Executors.newFixedThreadPool(NoCPUs);
        //UpdateProcess progress = new UpdateProcess();
        UpdateProcess progress = null;
        for (PeakCluster peakCluster : LCMSPeakBase.PeakClusters) {
            if (peakCluster.Charge >= parameter.StartCharge && peakCluster.Charge <= parameter.EndCharge && peakCluster.TargetMz() >= DIAWindowMz.getX() && peakCluster.TargetMz() <= DIAWindowMz.getY()) {
                CorrCalcCluster2CurveUnit unit = new CorrCalcCluster2CurveUnit(peakCluster, LCMSPeakBase.GetPeakCurveListRT(), parameter, progress);
                UnfragmentedIonPairList.add(unit);
            }
        }

        //progress.SetTotal(UnfragmentedIonPairList.size());
        //Thread thread = new Thread(progress);
        //thread.start();
        for (CorrCalcCluster2CurveUnit unit : UnfragmentedIonPairList) {
            executorPool.execute(unit);
        }
        executorPool.shutdown();
//        while (!executorPool.isTerminated()) {
//        }
        try {
            executorPool.awaitTermination(Long.MAX_VALUE, TimeUnit.NANOSECONDS);
        } catch (InterruptedException e) {
            Logger.getRootLogger().info("interrupted..");
        }

        ((LCMSPeakDIAMS2) LCMSPeakBase).UnFragIonClu2Cur = new HashMap<>();
        for (CorrCalcCluster2CurveUnit unit : UnfragmentedIonPairList) {
            if (!unit.GroupedFragmentList.isEmpty()) {
                ((LCMSPeakDIAMS2) LCMSPeakBase).UnFragIonClu2Cur.put(unit.MS1PeakCluster.Index, unit.GroupedFragmentList);
            }
        }
        //thread.stop();
        //thread = null;
        //progress.ClearMSG();
        //progress = null;
        //System.out.print("done\n");
        executorPool = null;
        
        ((LCMSPeakDIAMS2) LCMSPeakBase).BuildFragmentUnfragranking();
        ((LCMSPeakDIAMS2) LCMSPeakBase).FilterByCriteriaUnfrag();
        ((LCMSPeakDIAMS2) LCMSPeakBase).ExportUnfragmentedClusterCurve();
        //ExportParentClusterCurveCorrToDB(ResultArrayList);
        //GenerateMGF(ResultArrayList);
        UnfragmentedIonPairList.clear();
        UnfragmentedIonPairList = null;
        executorPool = null;

        //System.out.print("Finished multithreading\n");
    }

    public void CalcCorrByUnfragmentedCluster2Cluster() throws SQLException, IOException {
        //System.out.print("Using multithreading now: " + NoCPUs + " processors\n");
        Logger.getRootLogger().info("Clustering fragments by unfragmented clusters....");
        ExecutorService executorPool = null;
        ArrayList<CorrCalcCluster2ClusterUnit> ResultArrayList = new ArrayList<>();
        executorPool = Executors.newFixedThreadPool(NoCPUs);
        //UpdateProcess progress = new UpdateProcess();
        UpdateProcess progress = null;
        for (PeakCluster peakCluster : LCMSPeakBase.PeakClusters) {
            if (peakCluster.TargetMz() >= DIAWindowMz.getX() && peakCluster.TargetMz() <= DIAWindowMz.getY()) {
                CorrCalcCluster2ClusterUnit unit = new CorrCalcCluster2ClusterUnit(peakCluster, LCMSPeakBase.GetPeakClusterListRT(), parameter, progress);
                ResultArrayList.add(unit);
            }
        }

        //progress.SetTotal(ResultArrayList.size());
        //Thread thread = new Thread(progress);
        //thread.start();
        for (CorrCalcCluster2ClusterUnit unit : ResultArrayList) {
            executorPool.execute(unit);
        }
        executorPool.shutdown();
//        while (!executorPool.isTerminated()) {
//        }
        try {
            executorPool.awaitTermination(Long.MAX_VALUE, TimeUnit.NANOSECONDS);
        } catch (InterruptedException e) {
            Logger.getRootLogger().info("interrupted..");
        }

        ((LCMSPeakDIAMS2) LCMSPeakBase).UnFragIonClu2Clu = new HashMap<>();
        for (CorrCalcCluster2ClusterUnit unit : ResultArrayList) {
            ((LCMSPeakDIAMS2) LCMSPeakBase).UnFragIonClu2Clu.put(unit.MS1PeakCluster.Index, unit.ResultList);
        }
        //thread.stop();
        //thread = null;
        //progress.ClearMSG();
        //progress = null;
        //System.out.print("done\n");

        //ExportUnfragmentedCluster2Cluster(ResultArrayList);
        //ExportParentClusterCurveCorrToDB(ResultArrayList);
        //GenerateMGF(ResultArrayList);
        ResultArrayList.clear();
        ResultArrayList = null;
        executorPool = null;

        //System.out.print("Finished multithreading\n");
    }

    private void FragmentGroupingCluster2Cluster() throws SQLException, IOException {
        //System.out.print("Using multithreading now: " + NoCPUs + " processors\n");
        Logger.getRootLogger().info("Clustering fragments by MS1 curve....");
        ExecutorService executorPool = null;
        ArrayList<CorrCalcCluster2ClusterUnit> ResultArrayList = new ArrayList<>();
        executorPool = Executors.newFixedThreadPool(NoCPUs);
        //UpdateProcess progress = new UpdateProcess();
        UpdateProcess progress = null;
        for (PeakCluster peakCluster : ms1lcms.PeakClusters) {
            if (peakCluster.TargetMz() >= DIAWindowMz.getX() && peakCluster.TargetMz() <= DIAWindowMz.getY()) {
                CorrCalcCluster2ClusterUnit unit = new CorrCalcCluster2ClusterUnit(peakCluster, LCMSPeakBase.GetPeakClusterListRT(), parameter, progress);
                ResultArrayList.add(unit);
            }
        }

        //progress.SetTotal(ResultArrayList.size());
        //Thread thread = new Thread(progress);
        //thread.start();
        for (CorrCalcCluster2ClusterUnit unit : ResultArrayList) {
            executorPool.execute(unit);
        }
        executorPool.shutdown();
//        while (!executorPool.isTerminated()) {
//        }
        try {
            executorPool.awaitTermination(Long.MAX_VALUE, TimeUnit.NANOSECONDS);
        } catch (InterruptedException e) {
            Logger.getRootLogger().info("interrupted..");
        }

        ((LCMSPeakDIAMS2) LCMSPeakBase).FragmentsClu2Clu = new HashMap<>();
        for (CorrCalcCluster2ClusterUnit unit : ResultArrayList) {
            ((LCMSPeakDIAMS2) LCMSPeakBase).FragmentsClu2Clu.put(unit.MS1PeakCluster.Index, unit.ResultList);
        }
        //thread.stop();
        //thread = null;
        //progress.ClearMSG();
        //progress = null;
        //System.out.print("done\n");
        //ExportCluster2ClusterCorr(ResultArrayList);
        //ExportParentClusterCurveCorrToDB(ResultArrayList);
        //GenerateMGF(ResultArrayList);
        ResultArrayList.clear();
        ResultArrayList = null;
        executorPool = null;

        //System.out.print("Finished multithreading\n");
    }

    private void PrecursorFragmentPairBuildingForMS1() throws SQLException, IOException {
        //System.out.print("Using multithreading now: " + NoCPUs + " processors\n");
        Logger.getRootLogger().info("Building precursor-fragment pairs for MS1 features....");
        ExecutorService executorPool = null;
        ArrayList<CorrCalcCluster2CurveUnit> PrecursorPairList = new ArrayList<>();
        executorPool = Executors.newFixedThreadPool(NoCPUs);
        //UpdateProcess progress = new UpdateProcess();
        UpdateProcess progress = null;
        for (PeakCluster peakCluster : ms1lcms.PeakClusters) {
            if (peakCluster.GetMaxMz()>= DIAWindowMz.getX() && peakCluster.TargetMz() <= DIAWindowMz.getY()) {
                CorrCalcCluster2CurveUnit unit = new CorrCalcCluster2CurveUnit(peakCluster, LCMSPeakBase.GetPeakCurveListRT(), parameter, progress);
                PrecursorPairList.add(unit);
            }
        }

        //progress.SetTotal(PrecursorPairList.size());
        //Thread thread = new Thread(progress);
        //thread.start();
        for (CorrCalcCluster2CurveUnit unit : PrecursorPairList) {
            executorPool.execute(unit);
        }
        executorPool.shutdown();
//        while (!executorPool.isTerminated()) {
//        }
        try {
            executorPool.awaitTermination(Long.MAX_VALUE, TimeUnit.NANOSECONDS);
        } catch (InterruptedException e) {
            Logger.getRootLogger().info("interrupted..");
        }

        ((LCMSPeakDIAMS2) LCMSPeakBase).FragmentsClu2Cur = new HashMap<>();
        for (CorrCalcCluster2CurveUnit unit : PrecursorPairList) {
            if (!unit.GroupedFragmentList.isEmpty()) {
                ((LCMSPeakDIAMS2) LCMSPeakBase).FragmentsClu2Cur.put(unit.MS1PeakCluster.Index, unit.GroupedFragmentList);
            }
        }
        //thread.stop();
        //thread = null;
        //progress.ClearMSG();
        //progress = null;
        ((LCMSPeakDIAMS2) LCMSPeakBase).BuildFragmentMS1ranking();
        ((LCMSPeakDIAMS2) LCMSPeakBase).FilterByCriteria();
        ((LCMSPeakDIAMS2) LCMSPeakBase).ExportCluster2CurveCorr();
        //ExportParentClusterCurveCorrToDB(ResultArrayList);
        //GenerateMGF(ResultArrayList);
        PrecursorPairList.clear();
        PrecursorPairList = null;
        executorPool = null;

        //System.out.print("Finished multithreading\n");
    }
    
//    //<editor-fold defaultstate="collapsed" desc="Grouping isolated signal">
//    private void PrecursorFragmentPairBuildingForIsolatedPeakCurve() throws SQLException, IOException {
//        //System.out.print("Using multithreading now: " + NoCPUs + " processors\n");
//        if(ms1lcms.IsolatedMS1PeakCurves==null || ms1lcms.IsolatedMS1PeakCurves.isEmpty()){
//            return;
//        }
//        Logger.getRootLogger().info("Building precursor-fragment pairs for isolated MS1 curve....");
//        ExecutorService executorPool = null;
//        ArrayList<IsolatedCurveCorrCalcUnit> ResultArrayList = new ArrayList<>();
//        executorPool = Executors.newFixedThreadPool(NoCPUs);
//        UpdateProcess progress = new UpdateProcess();
//        for (PeakCurve peak : ms1lcms.IsolatedMS1PeakCurves) {
//            if (peak.TargetMz >= DIAWindowMz.getX() && peak.TargetMz <= DIAWindowMz.getY()) {
//                IsolatedCurveCorrCalcUnit unit = new IsolatedCurveCorrCalcUnit(peak, LCMSPeakBase.PeakCurveListRT, parameter, progress);
//                ResultArrayList.add(unit);
//            }
//        }
//
//        progress.SetTotal(ResultArrayList.size());
//        Thread thread = new Thread(progress);
//        thread.start();
//        for (IsolatedCurveCorrCalcUnit unit : ResultArrayList) {
//            executorPool.execute(unit);
//        }
//        executorPool.shutdown();
////        while (!executorPool.isTerminated()) {
////        }
//        try {
//            executorPool.awaitTermination(Long.MAX_VALUE, TimeUnit.NANOSECONDS);
//        } catch (InterruptedException e) {
//            Logger.getRootLogger().info("interrupted..");
//        }
//
//        ((LCMSPeakDIAMS2) LCMSPeakBase).IsolatedMS1IonClu2Cur = new HashMap<>();
//        for (IsolatedCurveCorrCalcUnit unit : ResultArrayList) {
//            ((LCMSPeakDIAMS2) LCMSPeakBase).IsolatedMS1IonClu2Cur.put(unit.targetMS1Curve.Index, unit.ResultList);
//        }
//        thread.stop();
//        thread = null;
//        progress.ClearMSG();
//        progress = null;
//        System.out.print("done\n");
//        ((LCMSPeakDIAMS2) LCMSPeakBase).ExportIsolatedClusterCurve(ResultArrayList);
//        //ExportParentClusterCurveCorrToDB(ResultArrayList);
//        //GenerateMGF(ResultArrayList);
//        ResultArrayList.clear();
//        ResultArrayList = null;
//        executorPool = null;
//
//    }
//</editor-fold>
}
