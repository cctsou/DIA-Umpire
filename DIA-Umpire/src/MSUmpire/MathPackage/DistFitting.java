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
package MSUmpire.MathPackage;

import MSUmpire.BaseDataStructure.XYData;
import ij.measure.CurveFitter;
import java.util.ArrayList;

/**
 *
 * @author Chih-Chiang Tsou <chihchiang.tsou@gmail.com>
 */


public class DistFitting {
       public CurveFitter GetGaussianFitter(double[] xData, double[] yData) {
        CurveFitter fitter = new CurveFitter(xData, yData);
        fitter.doFit(CurveFitter.GAUSSIAN);
        return fitter;
    }
       public ArrayList<XYData> GetGaussianFittingData() {
        
            ArrayList<XYData>  FitData = new ArrayList<>();
//            FitData.add(new XYData(SmoothData.Data.get(0).X - 0.1f, (float) GetGaussianFitter().f(GetGaussianFitter().getParams(), SmoothData.Data.get(0).X - 0.1f)));
//            for (float rt = SmoothData.Data.get(0).X; rt < SmoothData.Data.get(SmoothData.Data.size() - 1).X; rt += 0.05f) {
//                FitData.add(new XYData(rt, (float) GetGaussianFitter().f(GetGaussianFitter().getParams(), rt)));
//            }
//            FitData.add(new XYData(SmoothData.Data.get(SmoothData.Data.size() - 1).X + 0.1f, (float) GetGaussianFitter().f(GetGaussianFitter().getParams(), SmoothData.Data.get(SmoothData.Data.size() - 1).X + 0.1f)));
//        
        return FitData;
    }
}
