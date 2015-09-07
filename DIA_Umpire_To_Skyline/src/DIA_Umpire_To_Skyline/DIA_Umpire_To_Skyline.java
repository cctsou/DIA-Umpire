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
package DIA_Umpire_To_Skyline;

/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
import MSUmpire.BaseDataStructure.UmpireInfo;
import Utility.ConsoleLogger;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.HashMap;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;
import org.apache.commons.io.FilenameUtils;
import org.apache.log4j.Level;
import org.apache.log4j.Logger;

/**
 *
 * @author Chih-Chiang Tsou
 */
public class DIA_Umpire_To_Skyline {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) throws FileNotFoundException, IOException, Exception {
        System.out.println("=================================================================================================");
        System.out.println("DIA-Umpire_To_Skyline (version: " + UmpireInfo.GetInstance().Version + ")");
        if (args.length < 1) {
            System.out.println("command format error, it should be like: java -jar -Xmx20G DIA_Umpire_To_Skyline.jar Path NoThreads");
            System.out.println("command : java -jar -Xmx20G DIA_Umpire_To_Skyline.jar Path [Option]\n");
            System.out.println("\nOptions");
            System.out.println("\t-t\tNo. of threads, Ex: -t4 (using four threads, default value)");
            System.out.println("\t-cP\tPath of msconvert.exe for mzXML conversion, Ex: -cP (using four threads, default value)");
            return;
        }
        try {
            ConsoleLogger.SetConsoleLogger(Level.DEBUG);
            ConsoleLogger.SetFileLogger(Level.DEBUG, FilenameUtils.getFullPath(args[0]) + "diaumpire_to_skyline.log");
        } catch (Exception e) {
            System.out.println("Logger initialization failed");
        }

        Logger.getRootLogger().info("Path:" + args[0]);
        String msconvertpath="C:/inetpub/tpp-bin/msconvert";

        String WorkFolder = args[0];
        int NoCPUs = 4;
        
        for (int i = 1; i < args.length; i++) {
            if (args[i].startsWith("-")) {
                if (args[i].startsWith("-cP")) {
                    msconvertpath = args[i].substring(3);
                    Logger.getRootLogger().info("MSConvert path: " + msconvertpath);
                }
                if (args[i].startsWith("-t")) {
                    NoCPUs = Integer.parseInt(args[i].substring(2));
                    Logger.getRootLogger().info("No. of threads: " + NoCPUs);
                }
            }
        }
        

        HashMap<String, File> AssignFiles = new HashMap<>();

        try {
            File folder = new File(WorkFolder);
            if (!folder.exists()) {
                Logger.getRootLogger().info("Path: " + folder.getAbsolutePath() + " cannot be found.");
            }
            for (final File fileEntry : folder.listFiles()) {
                if (fileEntry.isFile() && fileEntry.getAbsolutePath().toLowerCase().endsWith(".mzxml")
                        && !fileEntry.getAbsolutePath().toLowerCase().endsWith("q1.mzxml")
                        && !fileEntry.getAbsolutePath().toLowerCase().endsWith("q2.mzxml")
                        && !fileEntry.getAbsolutePath().toLowerCase().endsWith("q3.mzxml")) {
                    AssignFiles.put(fileEntry.getAbsolutePath(), fileEntry);
                }
                if (fileEntry.isDirectory()) {
                    for (final File fileEntry2 : fileEntry.listFiles()) {
                        if (fileEntry2.isFile() && fileEntry2.getAbsolutePath().toLowerCase().endsWith(".mzxml")
                                && !fileEntry2.getAbsolutePath().toLowerCase().endsWith("q1.mzxml")
                                && !fileEntry2.getAbsolutePath().toLowerCase().endsWith("q2.mzxml")
                                && !fileEntry2.getAbsolutePath().toLowerCase().endsWith("q3.mzxml")) {
                            AssignFiles.put(fileEntry2.getAbsolutePath(), fileEntry2);
                        }
                    }
                }
            }

            Logger.getRootLogger().info("No. of files assigned :" + AssignFiles.size());
            for (File fileEntry : AssignFiles.values()) {
                Logger.getRootLogger().info(fileEntry.getAbsolutePath());
            }

            ExecutorService executorPool = null;
            executorPool = Executors.newFixedThreadPool(3);

            for (File fileEntry : AssignFiles.values()) {
                String mzXMLFile = fileEntry.getAbsolutePath();
                FileThread thread = new FileThread(mzXMLFile, NoCPUs,msconvertpath);
                executorPool.execute(thread);
            }
            executorPool.shutdown();
            try {
                executorPool.awaitTermination(Long.MAX_VALUE, TimeUnit.NANOSECONDS);
            } catch (InterruptedException e) {
                Logger.getRootLogger().info("interrupted..");
            }
        } catch (Exception e) {
            Logger.getRootLogger().error(e.getMessage());
            throw e;
        }
        Logger.getRootLogger().info("Job done");
        Logger.getRootLogger().info("=================================================================================================");

    }
}
