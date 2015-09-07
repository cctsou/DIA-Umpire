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
import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.TreeMap;

/**
 *
 * @author Chih-Chiang Tsou <chihchiang.tsou@gmail.com>
 */
public class MALDIDataParser {

    public ScanCollection scanCollection = null;
    public String filename;
    public InstrumentParameter parameter;
    public int TotalScan;
    public SpectralDataType.DataType datatype;
    public float cycletime=0.2f;
    public MALDIDataParser(String filename) throws IOException {
        parameter = new InstrumentParameter(InstrumentParameter.InstrumentType.Q_TOF);
        this.filename = filename;
        scanCollection = new ScanCollection(parameter.Resolution);
        scanCollection.Filename = filename;
        this.datatype = SpectralDataType.DataType.MALDI;                
    }

    public void Parse() throws FileNotFoundException, IOException {
        TreeMap<String, ScanData> SortedMap=new TreeMap<>();
        File folder =new File(filename);
                
        File[] listOfFiles = folder.listFiles();
        if (listOfFiles == null) {
            return;
        }
        int ScanNo=1;
        
        for (int i = 0; i < listOfFiles.length; i++) {
            if (listOfFiles[i].isFile()) {
                String files = listOfFiles[i].getName();
                if (files.toLowerCase().endsWith(".mgf")) {
                    MGFParser mgf = new MGFParser(folder + "/" + files);
                    for (Integer scanNo : mgf.scanCollection.GetMS2DescendingArray()) {
                        ScanData Scan = mgf.scanCollection.GetScan(scanNo);
                        Scan.MsLevel = 2;
                        Scan.MGFTitle = Scan.MGFTitle.split("_")[Scan.MGFTitle.split("_").length - 2];
                        Scan.Normalization();
                        Scan.Data.Finalize();
                        scanCollection.AddScan(Scan);
                    }
                } else if (files.toLowerCase().endsWith(".pkl")) {
                    PKLScanParser pkl = new PKLScanParser(folder + "/" + files);
                    pkl.scan.MsLevel = 1;
                    pkl.scan.MGFTitle = pkl.scan.MGFTitle.split("_")[pkl.scan.MGFTitle.split("_").length - 2];
                    pkl.scan.Num=ScanNo;
                    //pkl.scan.Normalization();
                    pkl.scan.RetentionTime=cycletime*ScanNo;
                    ScanNo++;
                    pkl.scan.Data.Finalize();
                    SortedMap.put(pkl.scan.MGFTitle, pkl.scan);
                    //scanCollection.AddScan(pkl.scan);
                }
            }
        }
        
        for(ScanData scan : SortedMap.values()){
            scanCollection.AddScan(scan);
        }        
    }
}
