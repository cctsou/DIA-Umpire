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
package Utility;

import java.text.DecimalFormat;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 *
 * @author Chih-Chiang Tsou <chihchiang.tsou@gmail.com>
 */
public class UpdateProcess implements Runnable {

    public int processed = 0;
    int TotalJobs;

    public void Update() {
        synchronized (this) {
            processed++;
        }
    }

    public void SetTotal(int TotalJobs) {
        this.TotalJobs = TotalJobs;
    }
    String msg = "";
    char back = '\b';

    private void OutputProgressPercentage() throws InterruptedException {
        DecimalFormat df = new DecimalFormat("#.#");
        while (processed < TotalJobs) {
            Thread.sleep(500);
            for (int i = 0; i < msg.length(); i++) {
                System.out.print(back);
            }
            double progressPercentage = (float) processed / (float) TotalJobs * 100;
            msg = df.format(progressPercentage) + " %";
            System.out.print(msg);
        }
    }

    public void ClearMSG() {
        for (int i = 0; i < msg.length(); i++) {
            System.out.print(back);
        }
    }

    private void OutputProgressBar() throws InterruptedException {
        double gap = 0.01d;
        while (processed <= TotalJobs) {
            Thread.sleep(500);
            double progressPercentage = (float) processed / (float) TotalJobs;
            if (progressPercentage > gap) {
                updateProgress(progressPercentage);
                gap += 0.01d;
            }
        }
        System.out.print("\r");
    }

    static void updateProgress(double progressPercentage) {
        final int width = 50; // progress bar width in chars
        System.out.print("\r[");
        int i = 0;
        for (; i <= (int) (progressPercentage * width); i++) {
            System.out.print(".");
        }
        for (; i < width; i++) {
            System.out.print(" ");
        }
        System.out.print("]" + Math.round(progressPercentage * 100) + " %");
    }

    @Override
    public void run() {        
        try {
            //OutputProgressBar();
            OutputProgressPercentage();        
        } catch (InterruptedException ex) {
            Logger.getLogger(UpdateProcess.class.getName()).log(Level.SEVERE, null, ex);
        }
    }
}
