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
package MSUmpire.SpectraLib;

import MSUmpire.SpectralProcessingModule.ScoreFunction;
import MSUmpire.BaseDataStructure.InstrumentParameter;
import MSUmpire.BaseDataStructure.ScanData;
import MSUmpire.BaseDataStructure.XYData;
import MSUmpire.BaseDataStructure.XYPointCollection;
import MSUmpire.SpectralProcessingModule.Binning;
import java.util.ArrayList;

/**
 *
 * @author Chih-Chiang Tsou <chihchiang.tsou@gmail.com>
 */
public class SpectraST {

    public float MatchScore;
    public ScanData ScanA;
    public ScanData ScanB;
    public ArrayList<XYData> MatchPeakMzInt;//X: mz, Y: intensity
    

    public void GeneratePeakMatching(float intsitythreshold, float PPM) {

        MatchPeakMzInt = new ArrayList<>();

        for (XYData peakB : ScanB.Data) {
            float lowmz = InstrumentParameter.GetMzByPPM(peakB.getX(), 1, PPM);
            int startidx = ScanA.GetLowerIndexOfX(lowmz);
            XYData closetPeak = null;

            for (int idx = startidx; idx < ScanA.Data.size(); idx++) {
                XYData peakA = ScanA.Data.get(idx);
                if (InstrumentParameter.CalcPPM(peakB.getX(), peakA.getX()) <= PPM) {
                    if (closetPeak == null || peakA.getY() > closetPeak.getY()) {
                        closetPeak = peakA;
                    }
                } else if (peakA.getX() > peakB.getX()) {
                    break;
                }
            }
            if (closetPeak != null) {
                MatchPeakMzInt.add(new XYData(peakB.getX(), peakB.getY()));
            }
        }
    }

    public void SetDataByScanData(ScanData ScanA, ScanData ScanB, float intensitythreshold) {
        this.ScanA = ScanA;
        this.ScanB = ScanB;
        MatchScore = CalcDotProductByScan(ScanA, ScanB, intensitythreshold);
    }

    public float CalcDotProductByScan(XYPointCollection scanA, XYPointCollection scanB, float IntThreshold) {
        if (scanA.PointCount() == 0 || scanB.PointCount() == 0) {
            return 0f;
        }
        Binning binning=new Binning();
                
        XYPointCollection normalizedA = ScoreFunction.SpectralNormalizationForScan(binning.Binning(scanA, IntThreshold, null));

        //if(normalizedA.PointCount()==0)
        //normalizedA = SpectralSTNormalization(Binning(scanA, IntThreshold, null));
        XYPointCollection normalizedB = ScoreFunction.SpectralNormalizationForScan(binning.Binning(scanB, IntThreshold, null));
        //if(normalizedB.PointCount()==0)
        //normalizedB = SpectralSTNormalization(Binning(scanB, IntThreshold, null));
        return ScoreFunction.CalcDotProductForScan(normalizedA, normalizedB);
    }

    public XYPointCollection GenerateConsensusSpec(ArrayList<XYPointCollection> Specs, float minratio) {
        int min = Integer.MAX_VALUE;
        int max = 0;
        for (XYPointCollection spec : Specs) {
            if (spec.Data.get(0).getX() < min) {
                min = (int) spec.Data.get(0).getX();
            }
            if (spec.Data.get(spec.PointCount() - 1).getX() > max) {
                max = (int) spec.Data.get(spec.PointCount() - 1).getX();
            }
        }

        int size = (max - min + 1);
        float[][] mzarray = new float[size][3];
        //dim0=mz;
        //dim1=intensity;
        //dim2=count;

        for (int i = 0; i < size; i++) {
            mzarray[i][0] = min + i;
        }
        for (XYPointCollection spec : Specs) {
            for (int i = 0; i < spec.PointCount(); i++) {
                mzarray[(int) spec.Data.get(i).getX() - min][1] += spec.Data.get(i).getY();
                mzarray[(int) spec.Data.get(i).getX() - min][2] += 1;
            }
        }

        XYPointCollection ConsensusSpec = new XYPointCollection();

        for (int i = 0; i < size; i++) {
            if (mzarray[i][2] >= Specs.size() * minratio) {
                ConsensusSpec.Data.add(new XYData(mzarray[i][0], mzarray[i][1] / mzarray[i][2]));
            }
        }
        return ConsensusSpec;
    }
}
