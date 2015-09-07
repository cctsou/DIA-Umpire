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
package MSUmpire.spectrumparser;

import MSUmpire.BaseDataStructure.ScanData;
import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import org.apache.commons.io.FilenameUtils;

/**
 *
 * @author Chih-Chiang Tsou <chihchiang.tsou@gmail.com>
 */
public class PKLScanParser {

    public String filename;
    public ScanData scan;

    public PKLScanParser(String filename) throws IOException {
        this.filename = filename;
        Parse();
    }

    private void Parse() throws FileNotFoundException, IOException {
        //806.080993652344,8429.974609375,1
        //832.287536621094,7226.927734375,1
        //854.039978027344,6682.37646484375,1
        //861.061340332031,8370.4716796875,1
        BufferedReader reader = new BufferedReader(new FileReader(filename));
        String line = "";
        String[] Values = null;
        scan = new ScanData();
        scan.MGFTitle = FilenameUtils.getBaseName(filename);
        while ((line = reader.readLine()) != null) {
            if ((Values = line.split(",")).length == 3) {
                scan.AddPoint(Float.parseFloat(Values[0]), Float.parseFloat(Values[1]));
            }
        }
        reader.close();
    }
}
