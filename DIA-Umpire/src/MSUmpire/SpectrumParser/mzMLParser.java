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
package MSUmpire.SpectrumParser;

import MSUmpire.BaseDataStructure.InstrumentParameter;
import MSUmpire.BaseDataStructure.ScanCollection;
import MSUmpire.BaseDataStructure.ScanData;
import MSUmpire.BaseDataStructure.SpectralDataType;
import MSUmpire.BaseDataStructure.XYData;
import java.io.*;
import java.util.ArrayList;
import java.util.concurrent.*;
import java.util.zip.DataFormatException;
import javax.xml.parsers.ParserConfigurationException;
import org.apache.log4j.Logger;
import org.xml.sax.SAXException;
import uk.ac.ebi.jmzml.model.mzml.Spectrum;
import uk.ac.ebi.jmzml.xml.io.MzMLObjectIterator;
import uk.ac.ebi.jmzml.xml.io.MzMLUnmarshaller;


/*
 * Using multi-threading to parse mzXML files Load all scans at once in the
 * constructor To change this template, choose Tools | Templates and open the
 * template in the editor.
 * modified from UCSD mzML spectrum parser
 */
/**
 *
 * @author Chih-Chiang Tsou
 */
public final class mzMLParser extends SpectrumParserBase{

    private final MzMLUnmarshaller unmarshaller;
    public ScanCollection scanCollection = null;
    
    public mzMLParser(String filename, InstrumentParameter parameter, SpectralDataType.DataType datatype, DIA_Setting dIA_Setting, int NoCPUs) throws FileNotFoundException, IOException, InterruptedException, ExecutionException, ParserConfigurationException, SAXException, DataFormatException {
        super(filename, parameter, datatype, dIA_Setting, NoCPUs);        
        scanCollection = InitializeScanCollection();
        unmarshaller = new MzMLUnmarshaller(new File(filename));
        ParseAllScans();
    }
    
    private void ParseAllScans() {
        MzMLObjectIterator<Spectrum> itr = unmarshaller.unmarshalCollectionFromXpath("/run/spectrumList/spectrum", Spectrum.class);        
        ArrayList<mzMLSpecConverter> ScanList=new ArrayList<>();        
        ExecutorService executorPool = null;
        executorPool = Executors.newFixedThreadPool(NoCPUs);
        Logger.getRootLogger().info("Starting to convert from jmzspec");
        while (itr.hasNext()) {
            Spectrum jmzSpec = itr.next();
            mzMLSpecConverter converter=new mzMLSpecConverter(jmzSpec,parameter);
            executorPool.execute(converter);    
            ScanList.add(converter);
        }
        executorPool.shutdown();
        try {
            executorPool.awaitTermination(Long.MAX_VALUE, TimeUnit.NANOSECONDS);
        } catch (InterruptedException e) {
            Logger.getRootLogger().info("interrupted..");
        }
        
        for (mzMLSpecConverter converter : ScanList) {
            ScanData spec = converter.spec;
            scanCollection.AddScan(spec);
            ElutionTimeToScanNoMap.put(spec.RetentionTime, spec.ScanNum);
            ScanToElutionTime.put(spec.ScanNum, spec.RetentionTime);
            MsLevelList.put(spec.ScanNum, spec.MsLevel);
            if (datatype != SpectralDataType.DataType.DDA) {
                if (spec.isolationWindowLoffset > 0f && spec.isolationWindowRoffset > 0f) {
                    if (!dIA_Setting.DIAWindows.containsKey(new XYData(spec.isolationWindowTargetMz - spec.isolationWindowLoffset, spec.isolationWindowTargetMz + spec.isolationWindowRoffset))) {
                        ArrayList<Integer> scanList2 = new ArrayList<>();
                        dIA_Setting.DIAWindows.put(new XYData(spec.isolationWindowTargetMz - spec.isolationWindowLoffset, spec.isolationWindowTargetMz + spec.isolationWindowRoffset), scanList2);
                    }
                    dIA_Setting.DIAWindows.get(new XYData(spec.isolationWindowTargetMz - spec.isolationWindowLoffset, spec.isolationWindowTargetMz + spec.isolationWindowRoffset)).add(spec.ScanNum);
                } else {
                    if (datatype == SpectralDataType.DataType.DIA_F_Window) {
                        spec.isolationWindowLoffset = (dIA_Setting.F_DIA_WindowSize + 1) * 0.2f;
                        spec.isolationWindowRoffset = (dIA_Setting.F_DIA_WindowSize + 1) * 0.8f;
                    } else if (datatype == SpectralDataType.DataType.DIA_V_Window) {
                        for (XYData window : dIA_Setting.DIAWindows.keySet()) {
                            if (window.getX() <= spec.isolationWindowTargetMz && window.getY() >= spec.isolationWindowTargetMz) {
                                dIA_Setting.DIAWindows.get(window).add(spec.ScanNum);
                                break;
                            }
                        }
                    }
                }
            }
        }        
    }
    
    
    public ScanCollection GetScanCollectionDIAMS2(XYData DIAWindow, boolean IncludePeak,float startRT, float endRT) throws InterruptedException, ExecutionException, IOException {
        if (dIA_Setting == null) {
            Logger.getRootLogger().error("This is not DIA data" + filename);
            return null;
        }
        return GetScanDIAMS2(DIAWindow, IncludePeak, startRT, endRT);
    }

    public ScanCollection GetScanDIAMS2(XYData DIAWindow, boolean IncludePeak, float startTime, float endTime) throws InterruptedException, ExecutionException, IOException {
        if (dIA_Setting == null) {
            Logger.getRootLogger().error(filename + " is not DIA data");
            return null;
        }
        ScanCollection swathScanCollection = new ScanCollection(parameter.Resolution);

        int StartScanNo = 0;
        int EndScanNo = 0;

        StartScanNo = GetStartScan(startTime);
        EndScanNo = GetEndScan(endTime);
        for (Integer scannum : dIA_Setting.DIAWindows.get(DIAWindow)) {
            if (scannum >= StartScanNo && scannum <= EndScanNo) {
                ScanData scan = scanCollection.ScanHashMap.get(scannum);
                swathScanCollection.AddScan(scan);
                swathScanCollection.ElutionTimeToScanNoMap.put(scan.RetentionTime, scan.ScanNum);
            }
        }
        return swathScanCollection;
    }

    public ScanCollection GetAllScanCollectionByMSLabel(boolean MS1Included, boolean MS2Included, boolean MS1Peak, boolean MS2Peak) throws InterruptedException, ExecutionException, IOException {
        return GetAllScanCollectionByMSLabel(MS1Included, MS2Included, MS1Peak, MS2Peak, 0f, 999999f);
    }

    public ScanCollection GetAllScanCollectionByMSLabel(boolean MS1Included, boolean MS2Included, boolean MS1Peak, boolean MS2Peak, float startTime, float endTime) throws InterruptedException, ExecutionException, IOException {
        ScanCollection newscanCollection = InitializeScanCollection();
        Logger.getRootLogger().debug("Memory usage before loading scans:" + Math.round((Runtime.getRuntime().totalMemory() - Runtime.getRuntime().freeMemory()) / 1048576) + "MB (" + NoCPUs + " threads)");

        ArrayList<Integer> IncludedScans = new ArrayList<>();
        
        for(Integer ScanNum : MsLevelList.keySet()){
            if(MsLevelList.get(ScanNum)==1 && MS1Included){
                IncludedScans.add(ScanNum);
            }
            if(MsLevelList.get(ScanNum)==2 && MS2Included){
                IncludedScans.add(ScanNum);
            }
        }
        
        int StartScanNo = 0;
        int EndScanNo = 0;

        StartScanNo = GetStartScan(startTime);        
        EndScanNo = GetEndScan(endTime);
        
        for(Integer scannum : IncludedScans){
            if(scannum >= StartScanNo && scannum <= EndScanNo){
                ScanData scan=scanCollection.ScanHashMap.get(scannum);                
                newscanCollection.AddScan(scan);
                newscanCollection.ElutionTimeToScanNoMap.put(scan.RetentionTime, scan.ScanNum);
            }
        }        
        System.gc();
        //System.out.print(".....done\n");
        Logger.getRootLogger().debug("Memory usage after loading scans:" + Math.round((Runtime.getRuntime().totalMemory() - Runtime.getRuntime().freeMemory()) / 1048576) + "MB");
        return newscanCollection;
    }
    
}
