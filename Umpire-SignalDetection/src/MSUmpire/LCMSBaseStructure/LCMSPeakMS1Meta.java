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
package MSUmpire.LCMSBaseStructure;

import MSUmpire.BaseDataStructure.InstrumentParameter;
import MSUmpire.BaseDataStructure.ScanCollection;
import MSUmpire.MathPackage.ChiSquareGOF;
import MSUmpire.PeptidePeakClusterDetection.PDHandlerMS1Meta;
import MSUmpire.spectrumparser.mzXMLParser;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.sql.SQLException;
import java.util.concurrent.ExecutionException;
import java.util.logging.Level;
import java.util.logging.Logger;
import java.util.zip.DataFormatException;
import javax.xml.parsers.ParserConfigurationException;
import org.xml.sax.SAXException;
import org.xmlpull.v1.XmlPullParserException;

/**
 *
 * @author Chih-Chiang Tsou <chihchiang.tsou@gmail.com>
 */
public class LCMSPeakMS1Meta extends LCMSPeakBase {

    public LCMSPeakMS1Meta(String Filename, InstrumentParameter parameter) {
        this.ScanCollectionName = Filename;
        this.ParentmzXMLName = Filename;
        //this.connectionManager = connectionManager;
        this.parameter = parameter;
        this.MaxNoPeakCluster = parameter.MaxNoPeakCluster;
        ChiSquareGOF.GetInstance(this.MaxNoPeakCluster);
        this.MinNoPeakCluster = parameter.MinNoPeakCluster;
        this.StartCharge = parameter.StartCharge;
        this.EndCharge = parameter.EndCharge;
        this.MiniIntensity = parameter.MinMSIntensity;
        this.SNR = parameter.SNThreshold;
    }

    public mzXMLParser GetmzXML() {
        if (mzxml == null) {
            try {
                mzxml = new mzXMLParser(ScanCollectionName, parameter, datattype, null, NoCPUs);
            } catch (Exception ex) {
                Logger.getLogger(LCMSPeakMS1Meta.class.getName()).log(Level.SEVERE, null, ex);
            } 
        }
        return mzxml;
    }

    public void PeakClusterDetection() throws FileNotFoundException, IOException, InterruptedException, ExecutionException, ParserConfigurationException, SAXException, DataFormatException, SQLException, XmlPullParserException {
        GetmzXML().GetAllScanCollectionByMSLabel(true, true, true, false);
        ScanCollection scanCollection = GetmzXML().scanCollection;
        parameter.NoPeakPerMin = (int) (5f / GetmzXML().GetMS1CycleTime());
        PDHandlerMS1Meta detection = new PDHandlerMS1Meta(this, NoCPUs, parameter.MS1PPM);
        detection.DetectPeakCurves(scanCollection);
        ExportPeakCurveResult();
        ExportPeakCluster();
        scanCollection.dispose();
    }

}
