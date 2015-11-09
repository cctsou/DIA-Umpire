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
package MSUmpire.SpectrumParser;

import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.RandomAccessFile;
import org.apache.log4j.Logger;

/**
 *
 * @author Chih-Chiang Tsou
 */
public class mzXMLIndexCheck {

    public static boolean CheckIndex(String filename) throws FileNotFoundException, IOException {

        Logger.getRootLogger().info("Checking index :"+filename);
        //System.out.println("Checking index :"+filename);
        try (RandomAccessFile fileHandler = new RandomAccessFile(filename, "r")) {
            StringBuilder sb = new StringBuilder();
            long preidx = 0;

            String CurrentLine = "";
            long currentLastPt = fileHandler.length() - 1;
            boolean indexexist = false;
            int linecount = 0;
            while (!(CurrentLine.trim().startsWith("<index name=") | CurrentLine.trim().startsWith("</msRun>"))) {
                //Read backward
                for (long filePointer = currentLastPt; filePointer != -1; filePointer--) {
                    fileHandler.seek(filePointer);
                    int readByte = fileHandler.readByte();
                    if (readByte == 0xA) {
                        if (filePointer == currentLastPt) {
                            continue;
                        } else {
                            currentLastPt = filePointer;
                            break;
                        }
                    } else if (readByte == 0xD) {
                        if (filePointer == currentLastPt - 1) {
                            continue;
                        } else {
                            currentLastPt = filePointer;
                            break;
                        }
                    }
                    sb.append((char) readByte);
                }
                linecount++;
                CurrentLine = sb.reverse().toString();
                sb = new StringBuilder();

                if (CurrentLine.trim().startsWith("</index>")) {
                    indexexist = true;
                }

                if (CurrentLine.trim().startsWith("</scan>")) {
                    Logger.getRootLogger().debug("File : " + filename + " index exists but the spectra section is not complete");                        
                    return false;
                }
                
                if (!indexexist && linecount > 10) {
                    fileHandler.close();
                    return false;
                }

                if (CurrentLine.trim().startsWith("<offset id")) {
                    int scanNo = Integer.parseInt(CurrentLine.substring(CurrentLine.indexOf("<offset id=\"") + 12).split("\"")[0]);
                    long index = (long) Long.parseLong(CurrentLine.substring(CurrentLine.indexOf(">") + 1, CurrentLine.indexOf("</offset>")));
                    if (index < 0) {
                        index = index + 2147483647l + 2147483648l;
                    }
                    if (preidx == index) {
                        Logger.getRootLogger().debug("File : " + filename + " index is not correct, ScanNo:" + scanNo + " and " + scanNo + 1 + " have same index");                        
                        sb = null;
                        fileHandler.close();
                        return false;
                    }
                    preidx=index;
                }
            }
            sb = null;
            fileHandler.close();
            return true;
        }
    }
}
