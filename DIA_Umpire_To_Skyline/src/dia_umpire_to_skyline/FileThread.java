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

import MSUmpire.DIA.DIAPack;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.concurrent.TimeUnit;
import org.apache.commons.io.FilenameUtils;
import org.apache.commons.lang.exception.ExceptionUtils;
import org.apache.log4j.Logger;

/**
 *
 * @author Chih-Chiang Tsou
 */
public class FileThread implements Runnable{

    String mzXMLFile;
    String msconvertpath;
    int NoCPUs;
    public FileThread(String mzXMLFile,int NoCPUs,String msconvertpath){
        this.mzXMLFile=mzXMLFile;
        this.NoCPUs=NoCPUs;
        this.msconvertpath=msconvertpath;
    }
    @Override
    public void run() {
        GenerateSkylineFiles();
    }
    
    public void GenerateSkylineFiles(){
        try {
            long time = System.currentTimeMillis();
            
            DIAPack DiaFile = new DIAPack(mzXMLFile, NoCPUs);
            if (!new File(FilenameUtils.getFullPath(DiaFile.Filename) + DiaFile.GetQ1Name() + ".mzXML").exists()
                    | !new File(FilenameUtils.getFullPath(DiaFile.Filename) + DiaFile.GetQ2Name() + ".mzXML").exists()
                    | !new File(FilenameUtils.getFullPath(DiaFile.Filename) + DiaFile.GetQ3Name() + ".mzXML").exists()) {
                return;
            }
            Logger.getRootLogger().info("=================================================================================================");
            Logger.getRootLogger().info("Processing " + mzXMLFile);
            
            if (!DiaFile.RawMGFExist()) {
                if (!DiaFile.LoadDIASetting()) {
                    Logger.getRootLogger().info("Loading DIA setting failed, job is incomplete");
                    System.exit(1);
                }
                if (!DiaFile.LoadParams()) {
                    Logger.getRootLogger().info("Loading parameters failed, job is incomplete");
                    System.exit(1);
                }
                DiaFile.BuildStructure();
                if (!DiaFile.ms1lcms.ReadPeakCluster()) {
                    Logger.getRootLogger().info("Loading peak and structure failed, job is incomplete");
                    System.exit(1);
                }
                DiaFile.CreateSkylingImportFolder();
                DiaFile.GenerateRawMGF();
                DiaFile.ClearStructure();
            }
            
            DiaFile.ConvertRawMGF(msconvertpath);
            ChangeScanTitlePepXML();
            DiaFile = null;
            System.gc();
            time = System.currentTimeMillis() - time;
            Logger.getRootLogger().info(mzXMLFile + " processed time:" + String.format("%d hour, %d min, %d sec", TimeUnit.MILLISECONDS.toHours(time), TimeUnit.MILLISECONDS.toMinutes(time) - TimeUnit.HOURS.toMinutes(TimeUnit.MILLISECONDS.toHours(time)), TimeUnit.MILLISECONDS.toSeconds(time) - TimeUnit.MINUTES.toSeconds(TimeUnit.MILLISECONDS.toMinutes(time))));
        } catch (Exception ex) {
            Logger.getRootLogger().error(ExceptionUtils.getStackTrace(ex));
        }
    }
    
    private void ChangeScanTitlePepXML() throws FileNotFoundException, IOException{        
        File fileEntry =new File(FilenameUtils.getFullPath(mzXMLFile));        
        String basename=FilenameUtils.getBaseName(mzXMLFile);
        for(File file : fileEntry.listFiles()){
            if(file.isFile() && file.getAbsoluteFile().toString().toLowerCase().endsWith("pep.xml")){
                String pepxmlbase=file.getName().split("\\.")[0];
                if(pepxmlbase.equals(basename+"_Q1") || pepxmlbase.equals(basename+"_Q2") ||pepxmlbase.equals(basename+"_Q3")){                    
                    BufferedReader reader=new BufferedReader(new FileReader(file));
                    String outputname=file.getName().replace("_Q", ".ForLibQ");
                    Logger.getRootLogger().info("Writing new pepXML files and correct the scan titles: "+outputname);
                    FileWriter writer=new FileWriter(FilenameUtils.getFullPath(mzXMLFile) + FilenameUtils.getBaseName(mzXMLFile) + "_Skyline/"+outputname);
                    String line="";
                    while ((line=reader.readLine())!=null) {
                        writer.write(line.replaceAll(basename+"_Q",basename+".ForLibQ")+"\n");
                    }
                    writer.close();
                }                
            }
        }
    }
}
