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
package MSUmpire.spectrumparser;

import MSUmpire.BaseDataStructure.InstrumentParameter;
import MSUmpire.BaseDataStructure.ScanCollection;
import MSUmpire.BaseDataStructure.ScanData;
import MSUmpire.BaseDataStructure.SpectralDataType;
import MSUmpire.SpectralProcessingModule.BackgroundDetector;
import MSUmpire.SpectralProcessingModule.Deisotoping;
import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;

/**
 *
 * @author Chih-Chiang Tsou <chihchiang.tsou@gmail.com>
 */
public class MGFParser {

    public ScanCollection scanCollection = null;
    public String filename;
    public InstrumentParameter parameter;
    public int TotalScan;
    public int NoCPUs = 4;
    public SpectralDataType.DataType datatype;

    public MGFParser(String filename, InstrumentParameter parameter, SpectralDataType.DataType datatype, int NoCPUs) throws IOException {
        this.filename = filename;
        scanCollection = new ScanCollection(parameter.Resolution);
        scanCollection.Filename = filename;
        this.parameter = parameter;
        this.datatype = datatype;
        this.NoCPUs = NoCPUs;
        //Parse();
    }

    public MGFParser(String filename) throws IOException {
        parameter = new InstrumentParameter(InstrumentParameter.InstrumentType.Q_TOF);
        this.filename = filename;
        scanCollection = new ScanCollection(parameter.Resolution);
        scanCollection.Filename = filename;
        this.datatype = SpectralDataType.DataType.DDA;
        this.NoCPUs = 12;
        //Parse();
    }

    public void GetAllScanCollectionDDA() throws FileNotFoundException, IOException {
        BufferedReader reader = new BufferedReader(new FileReader(filename));
        String line = "";
        String[] Values = null;
        int ScanNum = 1;

        while ((line = reader.readLine()) != null) {
            if (line.trim().startsWith("BEGIN IONS")) {
                ScanData scan = new ScanData();
                while (!(line = reader.readLine()).trim().startsWith("END IONS")) {
                    if (line.trim().startsWith("PEPMASS=")) {
                        scan.PrecursorMz = Float.parseFloat(line.trim().subSequence(8, line.trim().length()).toString());
                    }
                    if (line.trim().startsWith("CHARGE=")) {
                        scan.PrecursorCharge = Integer.parseInt(line.trim().subSequence(7, line.trim().length() - 1).toString());
                    }
                    if (line.trim().startsWith("RTINSECONDS=")) {
                        scan.RetentionTime = Float.parseFloat(line.trim().subSequence(12, line.trim().length()).toString()) / 60f;
                    }                    
                    if (line.trim().startsWith("TITLE=")) {
                        scan.MGFTitle = line.trim().subSequence(6, line.trim().length()).toString().replace(",", "_");
                        if (scan.MGFTitle.contains(" RT:")) {
                            scan.RetentionTime = Float.parseFloat(scan.MGFTitle.substring(scan.MGFTitle.indexOf(" RT:") + 4).split(" ")[0]);
                        }
                        if (scan.MGFTitle.contains("=Scan:")) {
                            scan.Num = Integer.parseInt(scan.MGFTitle.substring(scan.MGFTitle.indexOf("=Scan:") + 6).split(" ")[0]);
                            ScanNum = scan.Num;
                        }
                    }
                    if ((Values = line.split(" ")).length == 2) {
                        scan.AddPoint(Float.parseFloat(Values[0]), Float.parseFloat(Values[1]));
                    } else if ((Values = line.split("\t")).length == 2) {
                        scan.AddPoint(Float.parseFloat(Values[0]), Float.parseFloat(Values[1]));
                    }
                }
                if (scan.Num == 0) {
                    scan.Num = ScanNum++;
                }
                if (scan.PrecursorCharge == 0) {
                    scan.MsLevel = 1;
                } else {
                    scan.MsLevel = 2;
                }
                scan.centroided = false;
                scan.Data.Finalize();
                scan.background = 0f;
                if (parameter.EstimateBG) {
                    BackgroundDetector detector=new BackgroundDetector(scan);
                    detector.DetermineConstantBackground();
                } else {
                    if (scan.MsLevel == 1) {
                        scan.background = parameter.MinMSIntensity;
                    }
                    if (scan.MsLevel == 2) {
                        scan.background = parameter.MinMSMSIntensity;
                    }
                }

                if (parameter.Denoise) {
                    scan.RemoveSignalBelowBG();
                }

                if (!scan.centroided) {
                    scan.Centroiding(parameter.Resolution, scan.background);
                }
                if (parameter.Deisotoping && scan.MsLevel == 1) {
                    new Deisotoping(scan, parameter);
                }
                scanCollection.AddScan(scan);
            }
        }
        reader.close();
    }

    private void Parse() throws FileNotFoundException, IOException {
        //        BEGIN IONS
//PEPMASS=820.998855732003
//CHARGE=1+
//RTINSECONDS=200
//TITLE=Elution from: 0.14 to 0.14   period: 0   experiment: 2 cycles:  1
//200.9942 2.3857
//354.9856 2.3857
//370.9314 5.1571
//388.9714 9.6857
//390.9608 2.7429
//END IONS
        BufferedReader reader = new BufferedReader(new FileReader(filename));
        String line = "";
        String[] Values = null;
        int ScanNum = 1;

        while ((line = reader.readLine()) != null) {
            if (line.trim().startsWith("BEGIN IONS")) {
                ScanData scan = new ScanData();
                scan.Num = ScanNum++;
                scan.MsLevel = 2;
                while (!(line = reader.readLine()).trim().startsWith("END IONS")) {
                    if (line.trim().startsWith("PEPMASS=")) {
                        scan.PrecursorMz = Float.parseFloat(line.trim().subSequence(8, line.trim().length()).toString());
                    }
                    if (line.trim().startsWith("CHARGE=")) {
                        scan.PrecursorCharge = Integer.parseInt(line.trim().subSequence(7, line.trim().length() - 1).toString());
                    }
                    if (line.trim().startsWith("RTINSECONDS=")) {
                        scan.RetentionTime = Float.parseFloat(line.trim().subSequence(12, line.trim().length()).toString());
                    }
                    if (line.trim().startsWith("TITLE=")) {
                        scan.MGFTitle = line.trim().subSequence(6, line.trim().length()).toString().replace(",", "_");
                    }
                    if ((Values = line.split(" ")).length == 2) {
                        scan.AddPoint(Float.parseFloat(Values[0]), Float.parseFloat(Values[1]));
                    } else if ((Values = line.split("\t")).length == 2) {
                        scan.AddPoint(Float.parseFloat(Values[0]), Float.parseFloat(Values[1]));
                    }
                }
                scan.Data.Finalize();
                scanCollection.AddScan(scan);
            }
        }
        reader.close();
    }
}
