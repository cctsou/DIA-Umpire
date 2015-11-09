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
import MSUmpire.MathPackage.PearsonCorr;
import MSUmpire.PeakDataStructure.PeakCurve;

/**
 *
 * @author Chih-Chiang Tsou <chihchiang.tsou@gmail.com>
 */
public class CurveCorrUnit implements Runnable {

    public PeakCurve peakA;
    public PeakCurve peakB;
    public int i;
    public int j;
    public float corr = 0f;
    public InstrumentParameter parameter;

    public CurveCorrUnit(PeakCurve peakA, PeakCurve peakB, int i, int j, InstrumentParameter parameter) {
        this.peakA = peakA;
        this.peakB = peakB;
        this.i = i;
        this.j = j;
        this.parameter = parameter;
    }

    private float CalPeakCorr(PeakCurve peakA, PeakCurve peakB) {
        //System.out.print("Doing "+i+"_"+j+"\n");
        PearsonCorr corr = new PearsonCorr();
        float corre = corr.CalcCorr(peakA.GetPeakCollection(), peakB.GetPeakCollection(), parameter.NoPeakPerMin);
        corr = null;
        return corre;
    }

    @Override
    public void run() {
        boolean overlap = false;
        if (peakA.StartRT() >= peakB.StartRT() && peakA.StartRT() <= peakB.EndRT()) {
            overlap = true;
        } else if (peakA.EndRT() >= peakB.StartRT() && peakA.EndRT() <= peakB.EndRT()) {
            overlap = true;
        } else if (peakB.StartRT() >= peakA.StartRT() && peakB.StartRT() <= peakA.EndRT()) {
            overlap = true;
        } else if (peakB.EndRT() >= peakA.StartRT() && peakB.EndRT() <= peakA.EndRT()) {
            overlap = true;
        }
        if (overlap) {
            corr = CalPeakCorr(peakA, peakB);
        }
    }
}
