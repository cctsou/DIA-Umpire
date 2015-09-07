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

import MSUmpire.FragmentLib.FragmentLibManager;
import Utility.ConsoleLogger;
import org.apache.commons.io.FilenameUtils;
import org.apache.log4j.Level;
import org.xmlpull.v1.XmlPullParserException;

/**
 *
 * @author Chih-Chiang Tsou
 */
public class GenerateLibFSBySPtxt {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) throws XmlPullParserException, Exception {
        ConsoleLogger logger=new ConsoleLogger();
        logger.SetConsoleLogger(Level.DEBUG);
        String WorkFolder="F:\\Data\\Pan-Human\\";
        String TraML="phl004_consensus.sptxt";        
        FragmentLibManager ExlibManager = new FragmentLibManager(FilenameUtils.getBaseName(TraML), null);
        ExlibManager.ImportFragLibBySPTXT(WorkFolder+TraML);
        ExlibManager.WriteFragmentLibSerialization(WorkFolder);
    }
}
