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
package MSUmpire.PeakDataStructure;

import MSUmpire.BaseDataStructure.InstrumentParameter;
import MSUmpire.BaseDataStructure.XYData;
import MSUmpire.BaseDataStructure.XYPointCollection;
import MSUmpire.BaseDataStructure.XYZData;
import java.io.FileWriter;
import java.io.IOException;
import java.io.Serializable;
import java.sql.Connection;
import java.sql.ResultSet;
import java.sql.SQLException;
import java.sql.Statement;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.PriorityQueue;
import org.apache.commons.io.FilenameUtils;
import org.jfree.data.xy.XYSeries;

/**
 *
 * @author Chih-Chiang Tsou <chihchiang.tsou@gmail.com>
 */
public class PeakCurve implements Serializable  {
    private static final long serialVersionUID = 6498163564821L;

    private ArrayList<XYZData> PeakList;
    private XYPointCollection SmoothData;
    private float startint = 0f;
    public int Index;
    private float endrt = -1f;
    private float startrt = -1f;
    public int StartScan=-1;
    public int EndScan=-1;
    private float TotalIntMzF;
    private float TotalIntF;
    public float TargetMz;
    public float ApexInt;
    public float minIntF = Float.POSITIVE_INFINITY;
    public float ApexRT;
    public float MaxCorr = 0f;
    public boolean CheckState = false;
    public float ConflictCorr = 0f;
    public boolean Grouped = false;
    public transient HashSet<Integer> ChargeGrouped=new HashSet<>();
    public float MzVar = -1f;
    //private ArrayList<XYPoint> FitData = null;
    //public ArrayList<PSM> PSMs = new ArrayList<>();
    public SortedRidgeCollectionClass PeakRidgeList;
    public transient WaveletMassDetector waveletMassDetector;
    private ArrayList<XYZData> PeakRegionList;
    public ArrayList<Float> RegionRidge;
    private ArrayList<ArrayList<Float>> NoRidgeRegion;
    //private CurveFitter fitter = null;
    //private float fitness=-1f;
    public InstrumentParameter parameter;

    //using B-spline to generate smoothed peak signals
    public void DoBspline() {
        for (XYZData point : PeakList) {
            XYData pt = new XYData(point.getX(), point.getZ());
            SmoothData.AddPoint(pt);
        }
        Bspline bspline = new Bspline();
        SmoothData = bspline.Run(SmoothData, (int) Math.max((RTWidth() * parameter.NoPeakPerMin), PeakList.size()), 2);
        bspline = null;
    }

    public void DoInterpolation() {
        for (XYZData point : PeakList) {
            XYData pt = new XYData(point.getX(), point.getZ());
            SmoothData.AddPoint(pt);
        }
        LinearInterpolation interpo = new LinearInterpolation();
        SmoothData = interpo.Run(SmoothData, (int) Math.max((RTWidth() * parameter.NoPeakPerMin), PeakList.size()));
        interpo = null;
    }

    public void AddConflictScore(float corr) {
        synchronized (this) {
            ConflictCorr += corr;
        }
    }

    public float GetRawSNR() {
        return ApexInt / minIntF;
    }

    public void SetRTs(float StartRT, float EndRT) {
        startrt = StartRT;
        endrt = EndRT;
    }

//    public void CheckRidge(SortedRidgeCollectionClass list){
//        for(int i=0;i<list.size()-1;i++){
//            if(list.get(i+1).RT<list.get(i).RT){
//                System.out.print("");
//            }
//            if(list.get(i+1).equals(list.get(i))){
//                System.out.print("");
//            }
//        }
//    }
    public void DetectPeakRegion() {
        ArrayList<XYData> PeakArrayList = new ArrayList<>();
        PeakRidgeList = new SortedRidgeCollectionClass();
        PeakRegionList = new ArrayList<>();
        NoRidgeRegion = new ArrayList<>();
        if (RTWidth() * parameter.NoPeakPerMin < 1) {
            return;
        }
        //CheckRidge(PeakRidgeList);
//        for (int i = 0; i < PeakList.size(); i++) {
//            PeakArrayList.add(new XYPoint(PeakList.get(i).X, PeakList.get(i).Z));
//        }

        for (int i = 0; i < SmoothData.PointCount(); i++) {
            PeakArrayList.add(new XYData(SmoothData.Data.get(i).getX(), SmoothData.Data.get(i).getY()));
        }
        waveletMassDetector = new WaveletMassDetector(parameter, PeakArrayList, (int) (RTWidth() * parameter.NoPeakPerMin));
        waveletMassDetector.Run();

        int maxScale = waveletMassDetector.PeakRidge.length - 1;

        //trace peak ridge from maximum wavelet scale to minimum scale
        for (int i = maxScale; i >= 0; i--) {
            ArrayList<XYData> PeakRidgeArray = waveletMassDetector.PeakRidge[i];

            if (PeakRidgeArray == null) {
                maxScale = i;
                continue;
            }
            if (PeakRidgeArray.isEmpty()) {
                continue;
            }

            float[][] DisMatrixF = new float[PeakRidgeList.size()][PeakRidgeArray.size()];

            for (int k = 0; k < PeakRidgeList.size(); k++) {///For each existing peak ridge line
                for (int l = 0; l < PeakRidgeArray.size(); l++) {
                    DisMatrixF[k][l] = Math.abs(PeakRidgeList.get(k).RT - PeakRidgeArray.get(l).getX());
                }
            }

            boolean conti = true;
            ArrayList<XYData> RemovedRidgeList = new ArrayList<>();
            while (conti) {
                float closest = Float.MAX_VALUE;
                int ExistingRideIdx = -1;
                int PeakRidgeInx = -1;
                for (int k = 0; k < PeakRidgeList.size(); k++) {
                    for (int l = 0; l < PeakRidgeArray.size(); l++) {
                        {
                            if (DisMatrixF[k][l] < closest) {
                                closest = DisMatrixF[k][l];
                                ExistingRideIdx = k;
                                PeakRidgeInx = l;
                            }
                        }
                    }
                }

                if (closest < Float.MAX_VALUE && closest <= parameter.MinRTRange) {
                    PeakRidge ridge = PeakRidgeList.remove(ExistingRideIdx);
                    ridge.lowScale = i;
                    ridge.ContinueousLevel++;
                    XYData nearestRidge = PeakRidgeArray.get(PeakRidgeInx);
                    ridge.RT = nearestRidge.getX();
                    PeakRidgeList.add(ridge);
                    RemovedRidgeList.add(nearestRidge);
                    for (int k = 0; k < PeakRidgeList.size(); k++) {
                        DisMatrixF[k][PeakRidgeInx] = Float.MAX_VALUE;
                    }
                    for (int l = 0; l < PeakRidgeArray.size(); l++) {
                        DisMatrixF[ExistingRideIdx][l] = Float.MAX_VALUE;
                    }
                } else {
                    conti = false;
                }
            }

            for (XYData removeridge : RemovedRidgeList) {
                PeakRidgeArray.remove(removeridge);
            }
            RemovedRidgeList.clear();
            RemovedRidgeList = null;
            ArrayList<PeakRidge> removelist = new ArrayList<>();
            for (int k = 0; k < PeakRidgeList.size(); k++) {
                PeakRidge existridge = PeakRidgeList.get(k);
                if (existridge.lowScale - i > 2 && existridge.ContinueousLevel < maxScale / 2) {
                    removelist.add(existridge);
                }
            }
            for (int k = 0; k < removelist.size(); k++) {
                PeakRidgeList.remove(removelist.get(k));
            }
            removelist.clear();
            removelist = null;
            if (i > maxScale / 2) {
                for (XYData ridge : PeakRidgeArray) {
                    PeakRidge newRidge = new PeakRidge();
                    newRidge.RT = ridge.getX();
                    newRidge.lowScale = i;
                    newRidge.ContinueousLevel++;
                    newRidge.intensity = SmoothData.GetPoinByXCloset(newRidge.RT).getY();
                    PeakRidgeList.add(newRidge);
                }
            }
            PeakRidgeArray.clear();
            PeakRidgeArray = null;
        }

        if (PeakRidgeList.size() <= 1) {
            PeakRegionList.add(new XYZData(SmoothData.Data.get(0).getX(), ApexRT, SmoothData.Data.get(SmoothData.PointCount() - 1).getX()));
            ArrayList<Float> RidgeRTs = new ArrayList<>();
            RidgeRTs.add(ApexRT);
            NoRidgeRegion.add(RidgeRTs);
        }
        if (PeakRidgeList.size() > 1) {
            XYData[] ValleyPoints = new XYData[PeakRidgeList.size() + 1];
            ValleyPoints[0] = SmoothData.Data.get(0).cloneXYData();
            PeakRidge currentridge = PeakRidgeList.get(0);
            XYData localmin = new XYData(-1f, Float.MAX_VALUE);
            int startidx = SmoothData.GetLowerIndexOfX(currentridge.RT);

            for (int j = 1; j < PeakRidgeList.size(); j++) {
                PeakRidge nextridge = PeakRidgeList.get(j);
                for (int i = startidx; i < SmoothData.Data.size(); i++) {
                    XYData point = SmoothData.Data.get(i);
                    if (point.getX() > currentridge.RT && point.getX() < nextridge.RT) {
                        if (localmin.getY() > point.getY()) {
                            localmin = point.cloneXYData();
                        }
                    }
                    if (point.getX() >= nextridge.RT) {
                        startidx = i;
                        break;
                    }
                }
                ValleyPoints[j] = localmin;
                localmin = new XYData(-1f, Float.MAX_VALUE);
                currentridge = nextridge;
            }
            ValleyPoints[PeakRidgeList.size()] = SmoothData.Data.get(SmoothData.PointCount() - 1).cloneXYData();

            //Correct ridge rt and intensity
            startidx = 0;
            for (int i = 0; i < PeakRidgeList.size(); i++) {
                PeakRidge ridge = PeakRidgeList.get(i);
                for (int j = startidx; j < SmoothData.Data.size(); j++) {
                    XYData point = SmoothData.Data.get(j);
                    if (point.getX() < ValleyPoints[i + 1].getX()) {
                        if (ridge.intensity < point.getY()) {
                            ridge.intensity = point.getY();
                            ridge.RT = point.getX();
                        }
                    } else {
                        startidx = j;
                        break;
                    }
                }
            }

            //Find split points to generate peak regions
            boolean[] Splitpoints = new boolean[PeakRidgeList.size() - 1];
            int left = 0;
            int right = PeakRidgeList.size() - 1;
            FindSplitPoint(left, right, ValleyPoints, Splitpoints);

            ArrayList<Float> RidgeRTs = new ArrayList<>();
            startidx = 0;
            PeakRidge maxridge = PeakRidgeList.get(0);

            for (int i = 0; i < PeakRidgeList.size() - 1; i++) {
                RidgeRTs.add(PeakRidgeList.get(i).RT);
                if (PeakRidgeList.get(i).intensity > maxridge.intensity) {
                    maxridge = PeakRidgeList.get(i);
                }
                if (Splitpoints[i]) {
                    PeakRegionList.add(new XYZData(ValleyPoints[startidx].getX(), maxridge.RT, ValleyPoints[i + 1].getX()));
                    NoRidgeRegion.add(RidgeRTs);

                    maxridge = PeakRidgeList.get(i + 1);
                    RidgeRTs = new ArrayList<>();
                    startidx = i + 1;
                }
            }
            RidgeRTs.add(PeakRidgeList.get(PeakRidgeList.size() - 1).RT);
            if (PeakRidgeList.get(PeakRidgeList.size() - 1).intensity > maxridge.intensity) {
                maxridge = PeakRidgeList.get(PeakRidgeList.size() - 1);
            }
            PeakRegionList.add(new XYZData(ValleyPoints[startidx].getX(), maxridge.RT, ValleyPoints[PeakRidgeList.size()].getX()));

            NoRidgeRegion.add(RidgeRTs);
        }
        waveletMassDetector = null;
        PeakArrayList.clear();
        PeakArrayList = null;
        PeakRidgeList.clear();
        PeakRidgeList = null;
    }

    private void FindSplitPoint(int left, int right, XYData[] ValleyPoints, boolean[] splitpoints) {
        for (int i = left; i < right; i++) {
            if (ValidSplitPoint(left, right, i, ValleyPoints)) {
                splitpoints[i] = true;
                FindSplitPoint(left, i, ValleyPoints, splitpoints);
                FindSplitPoint(i + 1, right, ValleyPoints, splitpoints);
                break;
            }
        }
    }

    private boolean ValidSplitPoint(int left, int right, int cut, XYData[] ValleyPoints) {

        PeakRidge leftridge = PeakRidgeList.get(left);
        PeakRidge rightridge = PeakRidgeList.get(cut + 1);

        for (int i = left; i <= cut; i++) {
            if (PeakRidgeList.get(i).intensity > leftridge.intensity) {
                leftridge = PeakRidgeList.get(i);
            }
        }
        for (int i = cut + 1; i <= right; i++) {
            if (PeakRidgeList.get(i).intensity > rightridge.intensity) {
                rightridge = PeakRidgeList.get(i);
            }
        }
        return (Math.abs(ValleyPoints[left].getY() - ValleyPoints[cut + 1].getY()) / leftridge.intensity < parameter.SymThreshold && Math.abs(ValleyPoints[cut + 1].getY() - ValleyPoints[right + 1].getY()) / rightridge.intensity < parameter.SymThreshold);
    }

    public ArrayList<PeakCurve> SeparatePeakByRegion(float SN) {

        ArrayList<PeakCurve> tempArrayList = new ArrayList<>();
        ArrayList<PeakCurve> returnArrayList = new ArrayList<>();

//        if(TargetMz>466.73 && TargetMz<466.8 && startrt<20 && endrt>14){
//            System.out.print("");
//        }
        for (int i = 0; i < GetPeakRegionList().size(); i++) {
            PeakCurve peakCurve = new PeakCurve(parameter);
            peakCurve.RegionRidge = NoRidgeRegion.get(i);
            tempArrayList.add(peakCurve);
            XYZData region = GetPeakRegionList().get(i);
            if (region.getZ() - region.getX() > parameter.MaxCurveRTRange) {
                int leftidx = GetSmoothedList().GetLowerIndexOfX(region.getX());
                int rightidx = GetSmoothedList().GetHigherIndexOfX(region.getZ());
                XYData left = GetSmoothedList().Data.get(leftidx);
                XYData right = GetSmoothedList().Data.get(rightidx);
                while ((right.getX() - left.getX()) > parameter.MaxCurveRTRange) {
                    if (right.getX() - region.getY() <= parameter.MaxCurveRTRange / 4f) {
                        leftidx++;
                    } else if (region.getY() - left.getX() <= parameter.MaxCurveRTRange / 4f) {
                        rightidx--;
                    } else if (left.getY() < right.getY()) {
                        leftidx++;
                    } else {
                        rightidx--;
                    }
                    left = GetSmoothedList().Data.get(leftidx);
                    right = GetSmoothedList().Data.get(rightidx);
                }
                region.setX(left.getX());
                region.setZ(right.getX());
            }
        }

        for (int i = 0; i < GetPeakList().size(); i++) {
            XYZData peak = GetPeakList().get(i);
            for (int j = 0; j < GetPeakRegionList().size(); j++) {
                XYZData region = GetPeakRegionList().get(j);
                if (peak.getX() >= region.getX() && peak.getX() <= region.getZ()) {
                    tempArrayList.get(j).AddPeak(peak);
                    break;
                }
            }
        }

        for (int i = 0; i < GetSmoothedList().Data.size(); i++) {
            XYData peak = GetSmoothedList().Data.get(i);
            for (int j = 0; j < GetPeakRegionList().size(); j++) {
                XYZData region = GetPeakRegionList().get(j);
                if (peak.getX() >= region.getX() && peak.getX() <= region.getZ()) {
                    tempArrayList.get(j).GetSmoothedList().Data.add(peak);
                    break;
                }
            }
        }

        for (PeakCurve peak : tempArrayList) {
            if (peak.PeakList.size() > 2) {
                peak.GetSmoothedList().Data.Finalize();                
                returnArrayList.add(peak);
            }
        }        
        return returnArrayList;
    }

    public PeakCurve(InstrumentParameter parameter) {
        this.parameter = parameter;
        SmoothData = new XYPointCollection();
        PeakList = new ArrayList<>();
        PeakRegionList = new ArrayList<>();
    }

    public float StartInt() {
        if (startint == 0f) {
            startint = PeakList.get(0).getZ();
        }
        return startint;
    }

    public float StartRT() {
        if (startrt == -1f) {
            if (SmoothData != null && SmoothData.Data.size() > 0) {
                startrt = SmoothData.Data.get(0).getX();
            } else {
                startrt = PeakList.get(1).getX();
            }
        }
        return startrt;
    }
    float _snr = -1f;

    public float GetSNR() {
        if (_snr == -1) {
            //_snr = (ApexInt - GetBaseLine()) / GetNoiseLevel();
            //_snr = (float) Math.log(ApexInt);
            _snr = ApexInt;
        }
        return _snr;
    }

    public float GetMaxIntensityByRegionRange(float StartRT, float EndRT) {
        float max = 0f;
        for (int j = 0; j < GetSmoothedList().PointCount(); j++) {
            XYData pt = GetSmoothedList().Data.get(j);
            if (pt.getX() >= StartRT && pt.getX() <= EndRT && pt.getY() > max) {
                max = pt.getY();
            }
        }
        return max;
    }

    private void CalculateBaseLine() {
        _baseLine = 0f;
        PriorityQueue<Float> IntensityQueue = new PriorityQueue<>();
        for (XYData point : SmoothData.Data) {
            IntensityQueue.add(point.getY());
        }

        if (IntensityQueue.size() > 10) {
            for (int i = 0; i < IntensityQueue.size() / 10; i++) {
                _baseLine += IntensityQueue.poll();
            }
            _baseLine /= (IntensityQueue.size() / 10);
        } else {
            _baseLine = IntensityQueue.poll();
        }
//        _baseLine /= SmoothData.PointCount();
//        _noiseLevel = 0f;
//
//        for (XYPoint point : SmoothData.Data) {
//            _noiseLevel += (point.Y - _baseLine) * (point.Y - _baseLine);
//        }
//        _noiseLevel /= PeakList.size();
//        _noiseLevel = (float) Math.sqrt(_noiseLevel);
    }
    float _baseLine = -1f;

    public float GetBaseLine() {
        if (_baseLine == -1) {
            CalculateBaseLine();
            if (_baseLine == 0) {
                _baseLine = 1f;
            }
        }
        return _baseLine;
    }
    float _noiseLevel = -1f;

    public float GetNoiseLevel() {
        if (_noiseLevel == -1) {
            CalculateBaseLine();
        }
        return _noiseLevel;
    }

    public float EndRT() {
        if (endrt == -1f) {
            endrt = PeakList.get(PeakList.size() - 2).getX();
        }
        return endrt;
    }

    public float LastScanRT() {
        return PeakList.get(PeakList.size() - 1).getX();
    }

    public XYPointCollection GetPeakCollection() {
        XYPointCollection PtCollection = new XYPointCollection();

        for (int i = 0; i < SmoothData.Data.size(); i++) {
            PtCollection.AddPoint(SmoothData.Data.get(i).getX(), SmoothData.Data.get(i).getY());
        }
        return PtCollection;
    }

    public XYPointCollection GetSmoothPeakCollection(float startRT, float endRT) {
        XYPointCollection PtCollection = new XYPointCollection();

        for (int i = 0; i < SmoothData.PointCount(); i++) {
            XYData pt = SmoothData.Data.get(i);
            if (pt.getX() > endRT) {
                break;
            } else if (pt.getX() >= startRT && pt.getX() <= endRT) {
                PtCollection.AddPoint(pt.getX(), pt.getY());
            }
        }
        return PtCollection;
    }

    public float DetermineIntByRTRange(float StartRT, float EndRT) {
        float Intensity = 0f;
        for (int j = 0; j < GetSmoothedList().PointCount(); j++) {
            XYData pt = GetSmoothedList().Data.get(j);
            if (pt.getX() >= StartRT && pt.getX() <= EndRT) {
                if (pt.getY() > Intensity) {
                    Intensity = pt.getY();
                }
            }
        }
        return Intensity;
    }

//    public float GetGaussianFitness() {
//        if (fitness == -1) {
//            fitness = (float) GetGaussianFitter().getFitGoodness();
//        }
//        return fitness;
//    }
//    public float GetRsquired() {
//        return (float) GetGaussianFitter().getRSquared();
//    }
//    public ArrayList<XYPoint> GetGaussianFittingData() {
//        if (FitData == null) {
//            FitData = new ArrayList<>();
//            FitData.add(new XYPoint(SmoothData.Data.get(0).X - 0.1f, (float) GetGaussianFitter().f(GetGaussianFitter().getParams(), SmoothData.Data.get(0).X - 0.1f)));
//            for (float rt = SmoothData.Data.get(0).X; rt < SmoothData.Data.get(SmoothData.Data.size() - 1).X; rt += 0.05f) {
//                FitData.add(new XYPoint(rt, (float) GetGaussianFitter().f(GetGaussianFitter().getParams(), rt)));
//            }
//            FitData.add(new XYPoint(SmoothData.Data.get(SmoothData.Data.size() - 1).X + 0.1f, (float) GetGaussianFitter().f(GetGaussianFitter().getParams(), SmoothData.Data.get(SmoothData.Data.size() - 1).X + 0.1f)));
//        }
//        return FitData;
//    }
    public float RTWidth() {

        float Width = 0f;
        if (PeakList.size() > 0) {
            Width = PeakList.get(PeakList.size() - 1).getX() - PeakList.get(0).getX();
        } else if (SmoothData.PointCount() > 0) {
            Width = SmoothData.Data.get(SmoothData.PointCount() - 1).getX() - SmoothData.Data.get(0).getX();
        }
        return Width;
    }
//    public CurveFitter GetGaussianFitter() {
//        if (fitter == null) {
//            double[] xData = new double[SmoothData.Data.size() + 2];
//            double[] yData = new double[SmoothData.Data.size() + 2];
//
//            xData[0] = SmoothData.Data.get(0).X - 0.1f;
//            yData[0] = 0f;
//            for (int i = 0; i < SmoothData.Data.size(); i++) {
//                xData[i] = SmoothData.Data.get(i).X;
//                yData[i] = SmoothData.Data.get(i).Y;
//            }
//            xData[SmoothData.Data.size()] = SmoothData.Data.get(SmoothData.Data.size() - 1).X + 0.1;
//            yData[SmoothData.Data.size()] = 0f;
//            
//            fitter = new CurveFitter(xData, yData);
//            
//            fitter.doFit(CurveFitter.GAUSSIAN);
//        }
//        return fitter;
//    }

    public ArrayList<XYZData> GetPeakList() {
//        if(PeakList==null)
//        {
//            //ReadPeakResult();
//            ReadPeakResultMySQL(connection, Filename);
//        }
        return PeakList;
    }

    public XYPointCollection GetSmoothedList() {
//        if(PeakList==null)
//        {
//            //ReadPeakResult();
//            ReadPeakResultMySQL(connection, Filename);
//        }        
        return SmoothData;
    }

    public ArrayList<XYZData> GetPeakRegionList() {
        return PeakRegionList;
    }

    public void ReleasePeakData() {

        this.PeakList = null;
        this.SmoothData.dispose();
        this.SmoothData = null;
        this.PeakRegionList = null;
        this.PeakRidgeList = null;
        this.waveletMassDetector.DataPoint.clear();
//        for (int i = 0; i < this.waveletMassDetector.waveletCWT.length; i++) {
//            if (this.waveletMassDetector.waveletCWT[i] != null) {
//                this.waveletMassDetector.waveletCWT[i].clear();
//                this.waveletMassDetector.waveletCWT[i] = null;
//            }
//        }
        this.waveletMassDetector = null;
    }

    public void ReleaseRawPeak() {
        this.PeakList = null;
        this.PeakRegionList = null;
        this.PeakRidgeList = null;
        this.waveletMassDetector = null;
    }

    public XYSeries GetChartXYDatasetRAW() {
        XYSeries series1 = new XYSeries("Curve:" + Index + "_RAW");
        for (XYZData xYZPoint : PeakList) {
            series1.add(xYZPoint.getX(), xYZPoint.getZ());
        }
        return series1;
    }

    public XYSeries GetChartXYDatasetSmooth(String titleString) {
        XYSeries series1 = new XYSeries("Curve(" + TargetMz + "):" + Index);
        if (titleString != null) {
            series1.setKey(titleString);
        }
        for (int i = 0; i < SmoothData.PointCount(); i++) {
            XYData xYPoint = SmoothData.Data.get(i);
            series1.add(xYPoint.getX(), xYPoint.getY());
        }
        return series1;
    }
//    public XYSeries GetChartXYDatasetFitting()
//    {
//         XYSeries series1 = new XYSeries("Curve("+TargetMz+"):" + Index + "_Fitting");
//        for (XYPoint xYPoint : FitData) {
//            series1.add(xYPoint.X, xYPoint.Y);
//        }
//        return series1;
//    }

    public void ReadPeakResultMySQL(ResultSet rs) throws SQLException {
        //String PeakString = rs.getString("RAW_Peak");
        String SmoothedPeak = rs.getString("Smoothed_Peak");
        //String PeakRegion = rs.getString("PeakRegion");
//        String[] Peak = PeakString.split("#");
//        for (int i = 0; i < Peak.length; i++) {
//            if (Peak[i] != "") {
//                AddPeak(new XYZData(Float.parseFloat(Peak[i].split("_")[0]), Float.parseFloat(Peak[i].split("_")[1]), Float.parseFloat(Peak[i].split("_")[2])));
//            }
//        }
        SmoothData = new XYPointCollection();
        String[] Peak = SmoothedPeak.split("#");
        for (int i = 0; i < Peak.length; i++) {
            if (Peak[i] != "") {
                SmoothData.AddPoint(new XYData(Float.parseFloat(Peak[i].split("_")[0]), Float.parseFloat(Peak[i].split("_")[1])));
            }
        }
        SmoothData.Data.Finalize();
//        PeakRegionList = new ArrayList<>();
//        if (PeakRegion != null) {
//            Peak = PeakRegion.split("#");
//            for (int i = 0; i < Peak.length; i++) {
//                if (Peak[i] != "") {
//                    PeakRegionList.add(new XYZPoint(Float.parseFloat(Peak[i].split("_")[0]), Float.parseFloat(Peak[i].split("_")[1]), Float.parseFloat(Peak[i].split("_")[2])));
//                }
//            }
//        }
    }

    public void ReadPeakResultMySQL(Connection connection, String Filename) throws SQLException {
        Statement state = connection.createStatement();
        ResultSet rs = state.executeQuery("SELECT * FROM " + FilenameUtils.getBaseName(Filename) + "_PeakCurve WHERE Curve_index='" + Index + "'");

        rs.next();
        ReadPeakResultMySQL(rs);
        state.close();
        state = null;
    }

    public void ExportPeakResultMySQL(Connection connection, String Filename) throws SQLException {

        Statement state = connection.createStatement();
        String PeakString = "";
        for (int j = 0; j < PeakList.size(); j++) {
            PeakString += PeakList.get(j).getX() + "_" + PeakList.get(j).getY() + "_" + PeakList.get(j).getZ() + "#";
        }
        String SmoothedPeak = "";
        for (int j = 0; j < SmoothData.PointCount(); j++) {
            SmoothedPeak += SmoothData.Data.get(j).getX() + "_" + SmoothData.Data.get(j).getY() + "#";
        }
        String PeakRegion = "";
        for (int j = 0; j < PeakRegionList.size(); j++) {
            PeakRegion += PeakRegionList.get(j).getX() + "_" + PeakRegionList.get(j).getY() + "_" + PeakRegionList.get(j).getZ() + "#";
        }

//        String GaussianPeak = "";
//        ArrayList<XYPoint> FitData = GetGaussianFittingData();
//        for (int j = 0; j < FitData.size(); j++) {
//            GaussianPeak += FitData.get(j).X + "-" + FitData.get(j).Y + "_";
//        }                
        state.execute("INSERT INTO " + FilenameUtils.getBaseName(Filename) + "_PeakCurve (Curve_index, StartRT, EndRT, mz, RAW_Peak, Smoothed_Peak, PeakRegion) VALUES ('" + Index + "','" + StartRT() + "','" + EndRT() + "','" + TargetMz + "','" + PeakString + "','" + SmoothedPeak + "','" + PeakRegion + "')");
        state.close();
        //ReleasePeakData();
    }

    public void ExportPeakResultCSV(FileWriter writer) throws IOException {

//        String PeakString = "";
//        for (int j = 0; j < PeakList.size(); j++) {
//            PeakString += PeakList.get(j).getX() + "_" + PeakList.get(j).getY() + "_" + PeakList.get(j).getZ() + "#";
//        }
        String SmoothedPeak = "";
        for (int j = 0; j < SmoothData.PointCount(); j++) {
            SmoothedPeak += SmoothData.Data.get(j).getX() + "_" + SmoothData.Data.get(j).getY() + "#";
        }
        String Ridges = "";
        for (Float ridge : RegionRidge) {
            Ridges += ridge + "#";
        }
//        String PeakRegion = "";
//        for (int j = 0; j < PeakRegionList.size(); j++) {
//            PeakRegion += PeakRegionList.get(j).getX() + "_" + PeakRegionList.get(j).getY() + "_" + PeakRegionList.get(j).getZ() + "#";
//        }
//        writer.write(Index + "," + StartRT() + "," + EndRT() + "," + StartScan + "," + EndScan + "," + TargetMz + "," + PeakString + "," + SmoothedPeak + "," + PeakRegion + "\n");
        writer.write(Index + "," + StartRT() + "," + EndRT() + "," + TargetMz + "," + ApexRT + "," + ApexInt + "," + Ridges + "," + SmoothedPeak + "\n");
    }

    public void ExportPeakResultCSV_V2(FileWriter writer) throws IOException {

        String SmoothedPeak = "";
        for (int j = 0; j < SmoothData.PointCount(); j++) {
            SmoothedPeak += SmoothData.Data.get(j).getX() + "_" + SmoothData.Data.get(j).getY() + "#";
        }
        String Ridges = "";
        if (RegionRidge != null) {
            for (Float ridge : RegionRidge) {
                Ridges += ridge + "#";
            }
        }
        writer.write(Index + "," + StartRT() + "," + EndRT() + "," + TargetMz + "," + MzVar + "," + ApexRT + "," + ApexInt + "," + Ridges/*+ "," + PeakString*/ + "," + SmoothedPeak + "\n");
    }

    public void AddPeak(XYZData xYZPoint) {

        PeakList.add(xYZPoint);
        TotalIntMzF += xYZPoint.getY() * xYZPoint.getZ() * xYZPoint.getZ();
        TotalIntF += xYZPoint.getZ() * xYZPoint.getZ();
        if (xYZPoint.getZ() > ApexInt) {
            ApexInt = xYZPoint.getZ();
            ApexRT = xYZPoint.getX();
        }
        if (xYZPoint.getZ() < minIntF) {
            minIntF = xYZPoint.getZ();
        }
        TargetMz = TotalIntMzF / TotalIntF;
    }

    public void CalculateMzVar() {
        MzVar = 0f;
        for (int j = 0; j < PeakList.size(); j++) {
            MzVar += (PeakList.get(j).getX() - TargetMz) * (PeakList.get(j).getX() - TargetMz);
        }
        MzVar /= PeakList.size();
    }

}
