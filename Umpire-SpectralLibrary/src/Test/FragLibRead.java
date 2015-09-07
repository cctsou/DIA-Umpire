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
import jaligner.matrix.MatrixLoaderException;
import java.io.IOException;
import java.sql.SQLException;
import java.util.Arrays;
import org.apache.commons.io.FilenameUtils;
import org.apache.log4j.ConsoleAppender;
import org.apache.log4j.Level;
import org.apache.log4j.Logger;
import org.apache.log4j.PatternLayout;
import org.xmlpull.v1.XmlPullParserException;

/**
 *
 * @author Chih-Chiang Tsou
 */
public class FragLibRead {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) throws SQLException, MatrixLoaderException, IOException, XmlPullParserException, Exception {
         Logger logger = Logger.getRootLogger();
        ConsoleAppender ca = new ConsoleAppender();
        ca.setThreshold(Level.INFO);
        ca.setName("ConsoleLogger_Info");
        ca.setLayout(new PatternLayout("%d %-5p [%c{1}] %m%n"));
        ca.activateOptions();

        Logger.getRootLogger().info("Command: " + Arrays.toString(args));
        logger.getLoggerRepository().resetConfiguration();
        logger.addAppender(ca);
        
        String  sptxt="F:\\splib\\NIST_human_QTOF_2012-04-20_7AA.sptxt";
        sptxt="D:\\Projects\\NIST_human_QTOF_2012-04-20_7AA.sptxt";
        String traml="F:\\Data\\ETH_SEC_HEK\\Spectral_library_shotgun\\ConvertTSVToTraML.TraML";
     FragmentLibManager fragmentLibManager=new FragmentLibManager("Test", null);
     //fragmentLibManager.ImportFragLibBySPTXT(sptxt);
     fragmentLibManager.ImportFragLibByTraML(traml,"DECOY");
     
     fragmentLibManager.WriteFragmentLibSerialization(FilenameUtils.getFullPath(sptxt));
      //HashMap<String, PepFragmentLib> decoy=fragmentLibManager.GetDecoyLib();
    }

}
