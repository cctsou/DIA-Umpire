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

import MSUmpire.PSMDataStructure.EnzymeManager;
import MSUmpire.PSMDataStructure.PTMManager;
import Utility.ConsoleLogger;
import Utility.PrintThread;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;
import javax.xml.parsers.ParserConfigurationException;
import org.apache.log4j.Level;
import org.apache.log4j.Logger;
import org.xml.sax.SAXException;
import org.xmlpull.v1.XmlPullParserException;

/**
 *
 * @author Chih-Chiang Tsou
 */
public class GenerateLibUsingSpectraST {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) throws ParserConfigurationException, SAXException, IOException, XmlPullParserException, InterruptedException {
        ConsoleLogger.SetConsoleLogger(Level.DEBUG);
        String Path = "F:\\Data\\Pedro_Benchmark\\DIA_SpecLib\\";
        String combinelib = "CombineLib";
        File folder = new File(Path);
        ArrayList<SpectraSTBuild> tasklist = new ArrayList<>();
        
        ExecutorService executorPool = null;

        executorPool = Executors.newFixedThreadPool(10);
        EnzymeManager.GetInstance();
        PTMManager.GetInstance();
        
        for (final File fileEntry : folder.listFiles()) {
            if (fileEntry.getName().endsWith(".pep.xml")) {
                SpectraSTBuild build=new SpectraSTBuild(fileEntry.getAbsolutePath(), "reverse");
                tasklist.add(build);
                executorPool.execute(build);
            }
        }
        
        executorPool.shutdown();
        try {
            executorPool.awaitTermination(Long.MAX_VALUE, TimeUnit.NANOSECONDS);
        } catch (InterruptedException e) {
            System.out.println("interrupted..");
        }

        Process p = null;
        ArrayList<String> cmdlist = new ArrayList<>();
        cmdlist.add("spectrast");
        cmdlist.add("-cAC");
        cmdlist.add("-cJU");
        cmdlist.add("-c_DIS!");
        cmdlist.add("-c_QUO0.1");
        cmdlist.add("-cN" + Path + combinelib);
        for (SpectraSTBuild lib : tasklist) {
            cmdlist.add(lib.ResultName);
        }
        String[] cmd = new String[cmdlist.size()];
        cmd = cmdlist.toArray(cmd);
        p = Runtime.getRuntime().exec(cmd);
        Logger.getRootLogger().debug("Command: " + Arrays.toString(cmd));
        PrintThread printThread = new PrintThread(p);
        printThread.start();
        p.waitFor();
        if (p.exitValue() != 0) {
            Logger.getRootLogger().info("spectrast : combine library failed");
        }
    }
}
