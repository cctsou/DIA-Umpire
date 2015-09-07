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

package Test;

import MSUmpire.PSMDataStructure.LCMSID;
import MSUmpire.SearchResultParser.PepXMLParser;
import Utility.PrintThread;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.logging.Level;
import javax.xml.parsers.ParserConfigurationException;
import org.apache.commons.io.FilenameUtils;
import org.apache.log4j.Logger;
import org.xml.sax.SAXException;
import org.xmlpull.v1.XmlPullParserException;

/**
 *
 * @author Chih-Chiang Tsou
 */
public class SpectraSTBuild implements Runnable{

    String filename;
    String Decoytag="reverse";
    public String ResultName;
    public SpectraSTBuild(String filename, String Decoytag){
        this.filename=filename;
        this.Decoytag=Decoytag;
    }
    
    @Override
    public void run() {
        Build();
    }
    
    public void Build(){
        try {
            File fileEntry=new File(filename);
            ResultName= fileEntry.getParent() + "/" + FilenameUtils.getBaseName(fileEntry.getName()) + "_filtered" + ".splib";
            
            if(new File(ResultName).exists()){
                return;
            }
            
            LCMSID IDsummary = new LCMSID(fileEntry.getName(),Decoytag,"");
            PepXMLParser pepxml = new PepXMLParser(IDsummary, filename, 0f);
            IDsummary.FDR = 0.01f;
            IDsummary.FindPepProbThresholdByFDR();
            float threshold = IDsummary.PepProbThreshold;
                        
            Process p = null;
            ArrayList<String> cmdlist = new ArrayList<>();
            cmdlist.add("spectrast");
            cmdlist.add("-c_BIN!");
            cmdlist.add("-c_NPK1");
            cmdlist.add("-cf'Protein!~" + IDsummary.DecoyTag + "'");
            cmdlist.add("-cN\"" + fileEntry.getParent() + "/" + FilenameUtils.getBaseName(fileEntry.getName()) + "_filtered\"");
            cmdlist.add("-cP" + threshold);
            cmdlist.add(fileEntry.getAbsolutePath());
            String[] cmd = new String[cmdlist.size()];
            cmd = cmdlist.toArray(cmd);
            p = Runtime.getRuntime().exec(cmd);
            Logger.getRootLogger().debug("Command: " + Arrays.toString(cmd));
            PrintThread printThread = new PrintThread(p);
            printThread.start();
            p.waitFor();
            if (p.exitValue() != 0) {
                Logger.getRootLogger().info("spectrast : " + filename + " failed");
                System.exit(1);
            }
            IDsummary=null;
        } catch (ParserConfigurationException | SAXException | IOException | XmlPullParserException | InterruptedException ex) {
            java.util.logging.Logger.getLogger(SpectraSTBuild.class.getName()).log(Level.SEVERE, null, ex);
        }
    }

    
   
    
}
