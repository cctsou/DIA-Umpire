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
import MSUmpire.spectrumparser.DIA_Setting;
import MSUmpire.spectrumparser.mzXMLParser;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.IOException;
import java.util.concurrent.ExecutionException;
import java.util.zip.DataFormatException;
import javax.xml.parsers.ParserConfigurationException;
import org.xml.sax.SAXException;

/**
 *
 * @author Chih-Chiang Tsou <chihchiang.tsou@gmail.com>
 */
public class BKTest {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) throws IOException, DataFormatException, ParserConfigurationException, InterruptedException, FileNotFoundException, ExecutionException, SAXException, Exception {

        String Filename = "F:\\Data\\ETH_SEC_HEK_DDA\\F57_heuselm_J130809_014\\F57_heuselm_J130809_014.mzXML";
        //String Filename = "F:\\Data\\ETH_SEC_HEK\\F57_heuselm_J130809_013_SW\\F57_heuselm_J130809_013_SW_Q1.mzXML";
        //String Filename=args[0];
        //String Filename="/Test/tiny1.mzXML3.0.mzXML";
        //InstrumentParameter parameter= new InstrumentParameter(InstrumentParameter.InstrumentType.Q_TOF);
        InstrumentParameter parameter = new InstrumentParameter(InstrumentParameter.InstrumentType.TOF5600_ABcentroid);
        int NoCPUs = 10;
        //mzXMLParser mzxml = new mzXMLParser(Filename, parameter, SpectralDataType.DataType.IDA, NoCPUs);
        long start=System.nanoTime();
        mzXMLParser mzxml = new mzXMLParser(Filename, parameter, SpectralDataType.DataType.DDA, new DIA_Setting(), NoCPUs);
        //parameter.DetermineBGByID=false;
        parameter.Denoise=false;
        parameter.EstimateBG=false;
       
        //ScanData scan = mzxml.GetSingleRawScan(66918);
        //scan.Centroiding(parameter.Resolution, 0f);
        for (ScanData scan : mzxml.scanCollection.ScanHashMap.values()) {
            FileWriter writer = new FileWriter(Filename.replace(".mzXML", "_" + scan.Num + ".txt"));
        }
    }
}
 