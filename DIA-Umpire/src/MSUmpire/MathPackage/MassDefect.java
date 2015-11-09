/* 
 * Author: Chih-Chiang Tsou <chihchiang.tsou@gmail.com>
 *             Nesvizhskii Lab, Department of Computational Medicine and Bioinformatics, 
 *             University of Michigan, Ann Arbor
 *
 * Copyright 2015 University of Michigan, Ann Arbor, MI
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

/**
 *
 * @author Chih-Chiang Tsou
 */
public class MassDefect {

    public boolean InMassDefectRange(float mass){
        //upper = 0.00052738*x + 0.066015 +0.1 
        //lower = 0.00042565*x + 0.00038210 -0.1

        double u = GetMassDefect(0.00052738d*mass + 0.066015d +0.1d);
        double l = GetMassDefect(0.00042565d*mass + 0.00038210d -0.1d);
        
        double defect=GetMassDefect(mass);
        if (u > l) {
            return (defect>=l && defect<=u);
        }
        return (defect>=l || defect<=u);
    }
    
    public double GetMassDefect(double mass){
        return mass-Math.floor(mass);
    }
    
}
