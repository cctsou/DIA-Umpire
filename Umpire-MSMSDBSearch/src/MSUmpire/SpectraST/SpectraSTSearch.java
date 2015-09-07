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
package MSUmpire.SpectraST;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.InputStreamReader;
import org.apache.commons.io.FilenameUtils;

/**
 *
 * @author Chih-Chiang Tsou <chihchiang.tsou@gmail.com>
 */
public class SpectraSTSearch {

    public void RunSpectraST(String mzxml) throws IOException, InterruptedException {

        String[] cmd = {"spectrast", "-sF" + FilenameUtils.getFullPath(mzxml) + "/Spectrast.params", "-sD" + FilenameUtils.getFullPath(mzxml) + "/UPS_1_2_Ecoli_PlusRevTag2.fa", "-sL" + FilenameUtils.getFullPath(mzxml) + "/UPS_Ecoli_Consensus_decoy.splib", "-sO" + FilenameUtils.getFullPath(mzxml) + "/", mzxml};
        Process p = Runtime.getRuntime().exec(cmd);
        System.out.print("Searching using spectrast...." + mzxml + "\n");
        p.waitFor();
        PrintOutput(p);

        boolean finished = false;
        while (!finished) {
            finished = true;
            String[] xinteractcmd = {"xinteract", "-OpdEA", "-PPM", "-drev", "-p0", FilenameUtils.getFullPath(mzxml) + "/" + FilenameUtils.getBaseName(mzxml) + ".pep.xml", "-N" + FilenameUtils.getFullPath(mzxml) + "/interact-" + FilenameUtils.getBaseName(mzxml) + ".pep.xml"};

            p = Runtime.getRuntime().exec(xinteractcmd);
            System.out.print("processing xinteract....\n");
            //p.waitFor();
            PrintOutput(p);
            Thread.sleep(10000);

            int count = 0;
            while (!(new File(FilenameUtils.getFullPath(mzxml) + "/interact-" + FilenameUtils.getBaseName(mzxml) + ".pep.xml")).exists() || !(new File(FilenameUtils.getFullPath(mzxml) + "/interact-" + FilenameUtils.getBaseName(mzxml) + ".prot.xml")).exists()) {
                if (!(new File(FilenameUtils.getFullPath(mzxml) + "/interact-" + FilenameUtils.getBaseName(mzxml) + ".pep.xml")).exists()) {
                    System.out.print(FilenameUtils.getFullPath(mzxml) + "/interact-" + FilenameUtils.getBaseName(mzxml) + ".pep.xml" + " cannot be found. keep waiting....\n");
                    //PrintOutput(p);
                }
                if (!(new File(FilenameUtils.getFullPath(mzxml) + "/interact-" + FilenameUtils.getBaseName(mzxml) + ".prot.xml")).exists()) {
                    System.out.print(FilenameUtils.getFullPath(mzxml) + "/interact-" + FilenameUtils.getBaseName(mzxml) + ".prot.xml" + " cannot be found keep waiting....\n");
                    //PrintOutput(p);
                }
                Thread.sleep(10000);
                count++;
                if (count > 50) {
                    System.out.print("");
                    break;
                }
            }
        }

    }

    private void PrintOutput(Process p) throws IOException {
        BufferedReader is;  // reader for output of process
        String line;
        is = new BufferedReader(new InputStreamReader(p.getInputStream()));
        while ((line = is.readLine()) != null) {
            if (!"".equals(line.trim())) {
                System.out.println(line);
            }
        }
    }
}
