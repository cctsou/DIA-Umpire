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
package MSUmpire.SpectralProcessingModule;

import MSUmpire.BaseDataStructure.ScanData;
import java.util.ArrayList;
import java.util.Collections;

/**
 *
 * @author Chih-Chiang Tsou
 */
public class BackgroundDetector {

    public ScanData Scan;
    private float Ratio = 2f;

    public BackgroundDetector(ScanData Scan) {
        this.Scan = Scan;
    }

    public void DetermineBackgroundDynamic(float boundary) {
        if (Scan.Data.isEmpty()) {
            return;
        }
        ArrayList<Float> IntList = new ArrayList<>();
        for (int i = 0; i < Scan.Data.size(); i++) {
            IntList.add(Scan.Data.get(i).getY());
        }
        Collections.sort(IntList);
        float lower = IntList.get((int) (IntList.size() * 0.6f));
        int idx = 1;
        int count = 1;
        int remainingNo = IntList.size();
        while (idx < IntList.size()) {
            if (IntList.get(idx).equals(IntList.get(idx - 1))) {
                count++;
            } else {
                if (count > remainingNo / 5) {
                    Scan.background = IntList.get(idx - 1);
                    remainingNo = Math.max(200, IntList.size() - idx);
                    if (IntList.get(idx) > boundary) {
                        break;
                    }
                } else {
                    if (IntList.get(idx - 1) <= lower) {
                        Scan.background = IntList.get(idx - 1);
                    } else {
                        break;
                    }
                }
                count = 1;
            }
            idx++;
        }
        IntList.clear();
        IntList = null;
    }

    public void DetermineConstantBackground() {
        if (Scan.Data.isEmpty()) {
            return;
        }
        ArrayList<Float> IntList = new ArrayList<>();
        IntList.add(1f);
        for (int i = 0; i < Scan.Data.size(); i++) {
            IntList.add(Scan.Data.get(i).getY());
        }
        Collections.sort(IntList);
        int idx = 1;
        while (idx < IntList.size()) {
            if (IntList.get(idx).equals(IntList.get(idx - 1))) {
            } else {
                Scan.background = IntList.get(idx - 1);
            }
            idx++;
        }

        IntList.clear();
        IntList = null;
    }

    public void AdjacentPeakHistogram() {
        if (Scan.PointCount() < 10) {
            return;
        }
        ArrayList<Float> IntList = new ArrayList<>();
        for (int i = 0; i < Scan.Data.size(); i++) {
            IntList.add(Scan.Data.get(i).getY());
        }        
        Collections.sort(IntList);
        float upper = IntList.get((int) (IntList.size() * 0.7f));
        float lower = IntList.get(0);

        //FileWriter writer = new FileWriter("C:\\Umich\\Box Sync\\Default Sync Folder\\Background\\Test\\" + Scan.Num + "_hist.xls");
        //FileWriter writer2 = new FileWriter("C:\\Umich\\Box Sync\\Default Sync Folder\\Background\\Test\\" + Scan.Num + "_count.xls");
        //writer2.write(0+"\t"); 
        //writer.write(0 + "\t");
        if(upper<=lower+0.001){
            return;
        }
        int count1 = 0;
        int count2 = 0;
        int count3 = 0;
        int count4 = 0;
        int noise = 0;

        float bk = 0f;
        float interval = (upper - lower) / 20f;
        
        for (bk = lower; bk < upper; bk += interval) {
            count1 = 0;
            count2 = 0;
            count3 = 0;
            count4 = 0;
            noise = 0;
            //writer.write(bk + "\t");
            //writer2.write(bk+"\t");
            int preidx = -1;
            for (int i = 1; i < Scan.Data.size(); i++) {
                if (Scan.Data.get(i).getY() > bk) {
                    if (preidx != -1) {
                        float dist = Scan.Data.get(i).getX() - Scan.Data.get(preidx).getX();
                        //writer.write(dist + "\t");
                        if (dist > 0.95 && dist < 1.05 && Scan.Data.get(preidx).getY() > Scan.Data.get(i).getY()) {
                            count1++;
                        } else if (dist > 0.45 && dist < 0.55 && Scan.Data.get(preidx).getY() > Scan.Data.get(i).getY()) {
                            count2++;
                        } else if (dist > 0.3 && dist < 0.36 && Scan.Data.get(preidx).getY() > Scan.Data.get(i).getY()) {
                            count3++;
                        } else if (dist > 0.24 && dist < 0.26 && Scan.Data.get(preidx).getY() > Scan.Data.get(i).getY()) {
                            count4++;
                        } else if (dist < 0.23f) {
                            noise++;
                        }                        
                    }
                    preidx=i;
                }                
            }
            if (noise < (count1 + count2 + count3 + count4) * Ratio) {
                break;
            }
            //writer2.write(count1+"\t"+count2+"\t"+count3+"\t"+noise+"\n");
            //writer.write("\n");
        }
        if (bk > 0f) {
            Scan.background = bk;
            Scan.RemoveSignalBelowBG();
        }
        //writer.close();
        //writer2.close();
    }
}
