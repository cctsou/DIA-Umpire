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
package MSUmpire.DIA;

import java.io.Serializable;
import weka.core.Attribute;
import weka.core.FastVector;
import weka.core.Instance;

/**
 *
 * @author Chih-Chiang Tsou <chihchiang.tsou@gmail.com>
 */
public class MatchSubscore implements Serializable{
    private static final long serialVersionUID = -6347236539892953039L;    
    public subscore[] Subscores;    
    public int NoEnableSubscores;
    public double[] SubSCoeff;
    public String[] SubSName;
    
    FastVector Features=null;
    
    public FastVector GetFeatureFastVectorSVR(){
        if (Features == null) {
            NoEnableSubscores = 0;
            for (int i = 0; i < Subscores.length; i++) {
                if (Subscores[i].enable) {
                    NoEnableSubscores++;
                }
            }
            Features = new FastVector(NoEnableSubscores);

            for (int i = 0; i < Subscores.length; i++) {
                if (Subscores[i].enable) {
                    Attribute attributeF = new Attribute(Subscores[i].name);
                    Features.addElement(attributeF);
                }
            }
            Attribute attributeF = new Attribute("Type");
            Features.addElement(attributeF);
        }
        return Features;
    }
    
    public FastVector GetFeatureFastVectorSVM(){
        if (Features == null) {
            NoEnableSubscores = 0;
            for (int i = 0; i < Subscores.length; i++) {
                if (Subscores[i].enable) {
                    NoEnableSubscores++;
                }
            }
            Features = new FastVector(NoEnableSubscores);

            for (int i = 0; i < Subscores.length; i++) {
                if (Subscores[i].enable) {
                    Attribute attributeF = new Attribute(Subscores[i].name);
                    Features.addElement(attributeF);
                }
            }
            FastVector fvClassVal = new FastVector(2);
            fvClassVal.addElement("ID");
            fvClassVal.addElement("Decoy");
            Attribute attributeF = new Attribute("Type", fvClassVal);
            Features.addElement(attributeF);
        }
        return Features;
    }
    
    public class subscore implements Serializable{
        private static final long serialVersionUID = 2801998695766081693L;
        public String name;
        public boolean enable;
        public double initialweight;

        public subscore(String name, boolean enable,double iniweight) {
            this.name = name;
            this.enable = enable;
            this.initialweight=iniweight;
        }
    }
    public MatchSubscore() {
        Subscores = new subscore[18];
        Subscores[0] = new subscore("SpecDotProduct", true, 25d);
        Subscores[1] = new subscore("SpecCorrelation", false, 25d);
        Subscores[2] = new subscore("SpecContrastAngle", false, 25d);
        Subscores[3] = new subscore("CorrScore", true, -2d);
        Subscores[4] = new subscore("IntScore", false, -0.2d);
        Subscores[5] = new subscore("PPMScore", true, 4d);
        Subscores[6] = new subscore("ApexDeltaScore", true, -2d);
        Subscores[7] = new subscore("RTOverlapScore", true, 1d);
        Subscores[8] = new subscore("NoMatchB", false, 1d);
        Subscores[9] = new subscore("NoMatchY", false, 1d);
        Subscores[10] = new subscore("PrecursorCorr", true, 1d);
        Subscores[11] = new subscore("RTDiff", true, 0.1d);
        Subscores[12] = new subscore("PrecursorPPM", true, 1d);
        Subscores[13] = new subscore("WeightTotalMatch", false, 1d);
        Subscores[14] = new subscore("NorTotalMatch", false, 1d);
        Subscores[15] = new subscore("PrecursorIsotope", true, 1d);
        Subscores[16] = new subscore("SumCorr", false, 1d);
        Subscores[17] = new subscore("SumPPMCorr", true, 5d);
        //Subscores[18] = new subscore("DeltaPPMCorr", false, 1d);
    }
    
    public void InitializeLDACoeff() {
        NoEnableSubscores = 0;
        for (int i = 0; i < Subscores.length; i++) {
            if (Subscores[i].enable) {
                NoEnableSubscores++;
            }
        }
        SubSCoeff=new double[NoEnableSubscores];
        SubSName=new String[NoEnableSubscores];
        int idx=0;
        for (int i = 0; i < Subscores.length; i++) {
            if (Subscores[i].enable) {
                SubSCoeff[idx] = Subscores[i].initialweight;
                SubSName[idx]=Subscores[i].name;
                idx++;
            }
        }
    }
  
    
    public float [] GetSubScoreArray(PeakGroupScore peakgroup){
        float[] Subs = new float[18];
        int idx=0;
        Subs[idx++] = peakgroup.SpecDotProduct;
        Subs[idx++] = peakgroup.SpecCorrelation;
        Subs[idx++] = peakgroup.ContrastAngle;
        Subs[idx++] = peakgroup.AveCorrScore;
        Subs[idx++] = peakgroup.FragIntAvgScore;
        Subs[idx++] = peakgroup.PPMScore;
        Subs[idx++] = peakgroup.ApexDeltaScore;
        Subs[idx++] = peakgroup.RTOverlapScore;
        Subs[idx++] = peakgroup.NoMatchB;
        Subs[idx++] = peakgroup.NoMatchY;
        Subs[idx++] = peakgroup.PrecursorCorr;
        Subs[idx++] = peakgroup.RTDiff;
        Subs[idx++] = peakgroup.PrecursorPPM;
        Subs[idx++] = peakgroup.MaxMatchCorr*(peakgroup.NoMatchY+peakgroup.NoMatchB);
        Subs[idx++] = (peakgroup.MaxMatchCorr*(peakgroup.NoMatchY+peakgroup.NoMatchB))/peakgroup.Peplength;
        Subs[idx++] = peakgroup.PrecursorIsoPattern;
        Subs[idx++] = peakgroup.SumCorrScore;
        Subs[idx++] = peakgroup.SumCorrPPMScore;
        //Subs[idx++] = peakgroup.PrecursorCentralRank;
        return Subs;
    }
    
    public Instance GetFeatureInstance(PeakGroupScore peakgroup){
        Instance feature=new Instance(NoEnableSubscores+1);
        float[] Subs=GetSubScoreArray(peakgroup);
        int idx=0;
        for (int i = 0; i < Subscores.length; i++) {
            if (Subscores[i].enable) {
                feature.setValue((Attribute)GetFeatureFastVectorSVR().elementAt(idx), Subs[i]);
                idx++;
            }
        }
        return feature;
    }
    
    public double[] GetEnableSubScoreArray(PeakGroupScore peakgroup){
        float[] Subs=GetSubScoreArray(peakgroup);
        double[] enablesubsore=new double[NoEnableSubscores];
        int idx=0;
        for (int i = 0; i < Subscores.length; i++) {
            if (Subscores[i].enable) {
                enablesubsore[idx++] = Subs[i];
            }
        }    
        return enablesubsore;
    }
}
