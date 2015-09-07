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
package MSUmpire.mzXMLWrapper;

import MSUmpire.BaseDataStructure.XYData;
import MSUmpire.BaseDataStructure.ScanCollection;
import MSUmpire.BaseDataStructure.ScanData;
import java.io.FileWriter;
import java.io.IOException;

/**
 *
 * @author Chih-Chiang Tsou <chihchiang.tsou@gmail.com>
 */
public class MGFWrapper {

    ScanCollection scanCollection;

    public MGFWrapper(ScanCollection scanCollection) {
        this.scanCollection = scanCollection;
    }

    public void Export(String Filename) throws IOException {
        FileWriter mgfWriter = new FileWriter(Filename);
        for (ScanData scan : scanCollection.ScanHashMap.values()) {
            if (scan.MsLevel == 2) {
                StringBuilder mgfString = new StringBuilder();
                mgfString.append("BEGIN IONS\n");
                mgfString.append("PEPMASS=" + scan.PrecursorMz + "\n");
                mgfString.append("CHARGE=" + scan.PrecursorCharge + "+\n");
                mgfString.append("RTINSECONDS=" + scan.RetentionTime * 60f + "\n");
                mgfString.append("TITLE=" + scan.MGFTitle + ")\n");

                if (scan.PointCount() > 0) {
                    mgfWriter.write(mgfString.toString());
                }
                for (int i = 0; i < scan.Data.size(); i++) {
                    XYData point = scan.Data.get(i);
                    mgfWriter.write(point.getX() + " " + point.getY() + "\n");
                }
                if (scan.PointCount() > 0) {
                    mgfWriter.write("END IONS\n\n");
                } else {
                    System.out.print("");
                }
            }
        }
        mgfWriter.close();
    }
}
