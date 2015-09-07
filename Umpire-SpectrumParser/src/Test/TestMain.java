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

import MSUmpire.BaseDataStructure.InstrumentParameter;
import MSUmpire.BaseDataStructure.ScanData;
import MSUmpire.BaseDataStructure.SpectralDataType;
import MSUmpire.spectrumparser.mzXMLParser;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Arrays;
import java.util.concurrent.ExecutionException;
import java.util.zip.DataFormatException;
import javax.xml.parsers.ParserConfigurationException;
import org.apache.log4j.ConsoleAppender;
import org.apache.log4j.Level;
import org.apache.log4j.Logger;
import org.apache.log4j.PatternLayout;
import org.xml.sax.SAXException;

/**
 *
 * @author Chih-Chiang Tsou <chihchiang.tsou@gmail.com>
 */
public class TestMain {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) throws IOException, DataFormatException, ParserConfigurationException, InterruptedException, FileNotFoundException, ExecutionException, SAXException, Exception {

         Logger logger = Logger.getRootLogger();
        ConsoleAppender ca = new ConsoleAppender();
        ca.setThreshold(Level.DEBUG);
        ca.setName("ConsoleLogger_Info");
        ca.setLayout(new PatternLayout("%d %-5p [%c{1}] %m%n"));
        ca.activateOptions();

        Logger.getRootLogger().info("Command: " + Arrays.toString(args));
        logger.getLoggerRepository().resetConfiguration();
        logger.addAppender(ca);
        String Filename = "F:\\b1906_293T_proteinID_01A_QE3_122212.mzXML";
        //String Filename = "F:\\Data\\ETH_SEC_HEK\\F57_heuselm_J130809_013_SW\\F57_heuselm_J130809_013_SW_Q1.mzXML";
        //String Filename=args[0];
        //String Filename="/Test/tiny1.mzXML3.0.mzXML";
        //InstrumentParameter parameter= new InstrumentParameter(InstrumentParameter.InstrumentType.Q_TOF);
        InstrumentParameter parameter = new InstrumentParameter(InstrumentParameter.InstrumentType.TOF5600_ABcentroid);
        int NoCPUs = 10;
        mzXMLParser mzxml = new mzXMLParser(Filename, parameter, SpectralDataType.DataType.DDA, null, NoCPUs);
        long start=System.nanoTime();
        //mzXMLParser mzxml = new mzXMLParser(Filename, parameter, SpectralDataType.DataType.DIA_F_Window, new DIA_Setting(), NoCPUs);
        //parameter.DetermineBGByID=false;
        //parameter.Denoise=false;        
        //parameter.EstimateBG=false;
        //mzxml.GetAllScanCollectionByMSLabel(true, true, true, true);
        mzxml.GetAllScanCollectionByMSLabel(true, true, true, false);
        mzxml.GetAllScanCollectionMS2Only(true, true);
        FileWriter writer=new FileWriter(Filename.replace(".mzXML", "PeakNo.txt"));
        for(ScanData scan : mzxml.scanCollection.ScanHashMap.values()){
            if(scan.MsLevel==2){
            writer.write(scan.Num+"\t"+scan.PointCount()+"\n");
            }
        }
        writer.close();
        //ScanData scan = mzxml.GetSingleRawScan(66918);
        //scan.Centroiding(parameter.Resolution, 0f);
//        for (ScanData scan : mzxml.scanCollection.ScanHashMap.values()) {
//            BackgroundDetector detector = new BackgroundDetector(scan);
//            detector.AdjacentPeakHistogram();
//        }
    }
}
 