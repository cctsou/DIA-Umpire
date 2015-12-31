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
import MSUmpire.BaseDataStructure.XYData;
import MSUmpire.PeakDataStructure.PeakCluster;
import MSUmpire.LCMSPeakStructure.LCMSPeakMS1;
import MSUmpire.LCMSPeakStructure.LCMSPeakDIAMS2;
import MSUmpire.DIA.CorrCalcCluster2CurveUnit;
import MSUmpire.PeakDataStructure.PrecursorFragmentPairEdge;
import java.io.*;
import java.sql.SQLException;
import java.util.*;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;
import org.apache.log4j.Logger;

/**
 * Peak detection processing class for DIA MS2 peak
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
                        PrecursorFragmentPair.PeakCurveIndexB = scan.ScanNum;
                        PrecursorFragmentPair.FragmentMz = peak.getX();
                        PrecursorFragmentPair.Intensity = peak.getY();
                        PrecursorFragmentPair.RTOverlapP = 1f;
                        PrecursorFragmentPair.ApexDelta = Math.abs(peakCluster.MonoIsotopePeak.ApexRT - scan.RetentionTime);
                        float rtrange = peakCluster.endRT - peakCluster.startRT;
                        PrecursorFragmentPair.Correlation = (rtrange - PrecursorFragmentPair.ApexDelta) / rtrange;                        
                        ResultList.add(PrecursorFragmentPair);
                    }
                }
                ((LCMSPeakDIAMS2) LCMSPeakBase).FragmentsClu2Cur.put(peakCluster.Index, ResultList);
            }
        }
        
        ((LCMSPeakDIAMS2) LCMSPeakBase).ExportCluster2CurveCorr();
    }
    
    public void DetectPeakCurves(ScanCollection scanCollection) throws FileNotFoundException, IOException {        
        
        LCMSPeakBase.UnSortedPeakCurves = new ArrayList<>();
        FindAllMzTracePeakCurves(scanCollection);
        PeakCurveSmoothing();
        ClearRawPeaks();
        ReadPepIsoMS1PatternMap();
        PeakCurveCorrClustering(DIAWindowMz);
    }
        
    public void FragmentGrouping() throws SQLException, IOException {
        PrecursorFragmentPairBuildingForMS1();
        PrecursorFragmentPairBuildingForUnfragmentedIon();
    }


    private void PrecursorFragmentPairBuildingForUnfragmentedIon() throws SQLException, IOException {
        //System.out.print("Using multithreading now: " + NoCPUs + " processors\n");
        Logger.getRootLogger().info("Building precursor-fragment pairs for unfragmented ions....");
        ExecutorService executorPool = null;
        ArrayList<CorrCalcCluster2CurveUnit> UnfragmentedIonPairList = new ArrayList<>();
        executorPool = Executors.newFixedThreadPool(NoCPUs);
        for (PeakCluster peakCluster : LCMSPeakBase.PeakClusters) {
            if (peakCluster.Charge >= parameter.StartCharge && peakCluster.Charge <= parameter.EndCharge && peakCluster.TargetMz() >= DIAWindowMz.getX() && peakCluster.TargetMz() <= DIAWindowMz.getY()) {
                CorrCalcCluster2CurveUnit unit = new CorrCalcCluster2CurveUnit(peakCluster, LCMSPeakBase.GetPeakCurveListRT(), parameter);
                UnfragmentedIonPairList.add(unit);
            }
        }

        for (CorrCalcCluster2CurveUnit unit : UnfragmentedIonPairList) {
            executorPool.execute(unit);
        }
        executorPool.shutdown();

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
        executorPool = null;
        
        ((LCMSPeakDIAMS2) LCMSPeakBase).BuildFragmentUnfragranking();
        ((LCMSPeakDIAMS2) LCMSPeakBase).FilterByCriteriaUnfrag();
        ((LCMSPeakDIAMS2) LCMSPeakBase).ExportUnfragmentedClusterCurve();
        
        UnfragmentedIonPairList.clear();
        UnfragmentedIonPairList = null;
        executorPool = null;

        //System.out.print("Finished multithreading\n");
    }

    private void PrecursorFragmentPairBuildingForMS1() throws SQLException, IOException {
        //System.out.print("Using multithreading now: " + NoCPUs + " processors\n");
        Logger.getRootLogger().info("Building precursor-fragment pairs for MS1 features....");
        ExecutorService executorPool = null;
        ArrayList<CorrCalcCluster2CurveUnit> PrecursorPairList = new ArrayList<>();
        executorPool = Executors.newFixedThreadPool(NoCPUs);
        for (PeakCluster peakCluster : ms1lcms.PeakClusters) {
            if (peakCluster.GetMaxMz()>= DIAWindowMz.getX() && peakCluster.TargetMz() <= DIAWindowMz.getY()) {
                CorrCalcCluster2CurveUnit unit = new CorrCalcCluster2CurveUnit(peakCluster, LCMSPeakBase.GetPeakCurveListRT(), parameter);
                PrecursorPairList.add(unit);
            }
        }

        for (CorrCalcCluster2CurveUnit unit : PrecursorPairList) {
            executorPool.execute(unit);
        }
        executorPool.shutdown();
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
        ((LCMSPeakDIAMS2) LCMSPeakBase).BuildFragmentMS1ranking();
        ((LCMSPeakDIAMS2) LCMSPeakBase).FilterByCriteria();
        ((LCMSPeakDIAMS2) LCMSPeakBase).ExportCluster2CurveCorr();
        PrecursorPairList.clear();
        PrecursorPairList = null;
        executorPool = null;

        //System.out.print("Finished multithreading\n");
    }
    
}
