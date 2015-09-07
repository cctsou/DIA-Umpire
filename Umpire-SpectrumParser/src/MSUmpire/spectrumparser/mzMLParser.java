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
import MSUmpire.BaseDataStructure.XYData;
import Utility.UpdateProcess;
import java.io.*;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import java.util.Map.Entry;
import java.util.TreeMap;
import java.util.concurrent.*;
import java.util.zip.DataFormatException;
import javax.xml.parsers.ParserConfigurationException;
import org.apache.log4j.Logger;
import org.xml.sax.SAXException;
import uk.ac.ebi.jmzml.xml.io.MzMLUnmarshaller;


/*
 * Using multi-threading to parse mzXML files Load all scans at once in the
 * constructor To change this template, choose Tools | Templates and open the
 * template in the editor.
 */
/**
 *
 * @author Chih-Chiang Tsou
 */
public final class mzMLParser {

    public ScanCollection scanCollection = null;
    public String filename;
    public InstrumentParameter parameter;
    public int TotalScan;
    public TreeMap<Integer, Long> ScanIndex;
    public int NoCPUs = 4;
    public SpectralDataType.DataType datatype;
    public TreeMap<XYData, ArrayList<Integer>> DIAWindows = null;
    public TreeMap<XYData, ScanCollection> DIAMS2Scans = null;
    public TreeMap<Integer, Integer> MsLevelList;
    public TreeMap<Float, Integer> ElutionTimeToScanNoMap;
    public int NoMS1Scans = 0;

    public mzMLParser(String filename, InstrumentParameter parameter, SpectralDataType.DataType datatype, int NoCPUs) throws FileNotFoundException, IOException, InterruptedException, ExecutionException, ParserConfigurationException, SAXException, DataFormatException {
        this.filename = filename;
        ElutionTimeToScanNoMap = new TreeMap<>();
        MsLevelList = new TreeMap<>();
        scanCollection = new ScanCollection(parameter.Resolution);
        scanCollection.Filename = filename;
        this.parameter = parameter;
        this.datatype = datatype;
        this.NoCPUs = NoCPUs;        
    }

    public float GetMS1CycleTime() {
        return (ElutionTimeToScanNoMap.lastKey() - ElutionTimeToScanNoMap.firstKey()) / NoMS1Scans;
    }

    public ScanCollection GetScanCollectionDIAMS2(XYData DIAWindow, boolean IncludePeak) throws InterruptedException, ExecutionException, IOException {
        if (!DIAMS2Scans.containsKey(DIAWindow)) {
            DIAMS2Scans.put(DIAWindow, GetScanCollectionDIAMS2(DIAWindow, IncludePeak, 0f, 999999f));
        }
        return DIAMS2Scans.get(DIAWindow);
    }

    public ScanCollection GetScanCollectionDIAMS2(XYData DIAWindow, boolean IncludePeak, float startTime, float endTime) throws InterruptedException, ExecutionException, IOException {
        ScanCollection swathScanCollection = new ScanCollection(parameter.Resolution);
        //System.out.print("Multithreading: "+NoCPUs +" processors (Memory usage:"+ Math.round((Runtime.getRuntime().totalMemory() -Runtime.getRuntime().freeMemory())/1048576)+"MB)\n");
        //System.out.print("...Reading all scans of SWATH window:" + swathwin.X + " - " + swathwin.Y + "....");        
        ExecutorService executorPool = null;
        List<MzXMLthreadUnit> ScanList = new ArrayList<>();

        //UpdateProcess progress = new UpdateProcess();
        UpdateProcess progress = null;

        executorPool = Executors.newFixedThreadPool(NoCPUs);

        int StartScanNo = 0;
        int EndScanNo = 0;

        if (startTime <= ElutionTimeToScanNoMap.firstKey()) {
            StartScanNo = ElutionTimeToScanNoMap.firstEntry().getValue();
        } else {
            if (startTime >= ElutionTimeToScanNoMap.lastKey()) {
                StartScanNo = ElutionTimeToScanNoMap.lastEntry().getValue();
            } else {
                StartScanNo = ElutionTimeToScanNoMap.lowerEntry(startTime).getValue();
            }
        }
        if (endTime <= ElutionTimeToScanNoMap.firstKey()) {
            EndScanNo = ElutionTimeToScanNoMap.firstEntry().getValue();
        } else {
            if (endTime >= ElutionTimeToScanNoMap.lastKey()) {
                EndScanNo = ElutionTimeToScanNoMap.lastEntry().getValue();
            } else {
                EndScanNo = ElutionTimeToScanNoMap.higherEntry(endTime).getValue();
            }
        }
        Iterator<Entry<Integer, Long>> iter = ScanIndex.entrySet().iterator();

        Entry<Integer, Long> ent = iter.next();
        Long currentIdx = ent.getValue();
        int nextScanNo = ent.getKey();

        while (iter.hasNext()) {
            ent = iter.next();
            long startposition = currentIdx;
            long nexposition = ent.getValue();
            int currentScanNo = nextScanNo;
            nextScanNo = ent.getKey();
            currentIdx = nexposition;

            if (currentScanNo >= StartScanNo && currentScanNo <= EndScanNo && DIAWindows.get(DIAWindow).contains(currentScanNo)) {
                byte[] buffer = new byte[(int) (nexposition - startposition)];
                RandomAccessFile fileHandler = new RandomAccessFile(filename, "r");
                fileHandler.seek(startposition);
                fileHandler.read(buffer, 0, (int) (nexposition - startposition));
                fileHandler.close();
                String xmltext = new String(buffer);
                if (ent.getKey() == Integer.MAX_VALUE) {
                    xmltext = xmltext.replaceAll("</msRun>", "");
                    buffer = null;
                }
                boolean ReadPeak = true;
                MzXMLthreadUnit unit = new MzXMLthreadUnit(xmltext, parameter, datatype, progress, ReadPeak);
                ScanList.add(unit);
                buffer = null;
                xmltext = null;
                fileHandler = null;
            }
        }

        //progress.SetTotal(ScanList.size());
        //Thread thread = new Thread(progress);
        //thread.start();
        for (MzXMLthreadUnit unit : ScanList) {
            executorPool.execute(unit);
        }
        executorPool.shutdown();
//        while (!executorPool.isTerminated()) {
//        }
        try {
            executorPool.awaitTermination(Long.MAX_VALUE, TimeUnit.NANOSECONDS);
        } catch (InterruptedException e) {
            Logger.getRootLogger().info("interrupted..");
        }
        //thread.stop();
        //thread = null;
        //progress.ClearMSG();
        //progress = null;
        for (MzXMLthreadUnit result : ScanList) {
            result.scan.MsLevel = 1;
            swathScanCollection.AddScan(result.scan);
            swathScanCollection.ElutionTimeToScanNoMap.put(result.scan.RetentionTime, result.scan.Num);
        }
        executorPool = null;
        ScanList.clear();
        ScanList = null;
        ent = null;
        iter = null;
        //System.gc();
        //System.out.print(".....done\n");
        //System.out.print("Finished multithreading (Memory usage:"+ Math.round((Runtime.getRuntime().totalMemory() -Runtime.getRuntime().freeMemory())/1048576)+"MB)\n");

        return swathScanCollection;
    }

    public void GetAllScanCollectionByMSLabel(boolean MS1Included, boolean MS2Included, boolean MS1Peak, boolean MS2Peak) throws InterruptedException, ExecutionException, IOException {
        GetAllScanCollectionByMSLabel(MS1Included, MS2Included, MS1Peak, MS2Peak, 0f, 999999f);
    }
    
    public void Test(){
        File xmlFile = new File("path/to/your/mzml/file");
        MzMLUnmarshaller unmarshaller = new MzMLUnmarshaller(xmlFile);    
        
    }

    public void GetAllScanCollectionByMSLabel(boolean MS1Included, boolean MS2Included, boolean MS1Peak, boolean MS2Peak, float startTime, float endTime) throws InterruptedException, ExecutionException, IOException {
        //System.out.print("Multithreading: "+NoCPUs +" processors (Memory usage:"+ Math.round((Runtime.getRuntime().totalMemory() -Runtime.getRuntime().freeMemory())/1048576)+"MB)\n");
        //System.out.print("...Reading all scans....");
        System.out.print("Memory usage before loading scans:" + Math.round((Runtime.getRuntime().totalMemory() - Runtime.getRuntime().freeMemory()) / 1048576) + "MB (" + NoCPUs + " threads)\n");

        ArrayList<Integer> IncludedMSlevel = new ArrayList<>();
        if (MS1Included) {
            IncludedMSlevel.add(1);
        }
        if (MS2Included) {
            IncludedMSlevel.add(2);
        }
        ExecutorService executorPool = null;
        List<MzXMLthreadUnit> ScanList = new ArrayList<>();

        //UpdateProcess progress = new UpdateProcess();
        UpdateProcess progress = null;

        executorPool = Executors.newFixedThreadPool(NoCPUs);

        int StartScanNo = 0;
        int EndScanNo = 0;

        if (startTime <= ElutionTimeToScanNoMap.firstKey()) {
            StartScanNo = ElutionTimeToScanNoMap.firstEntry().getValue();
        } else {
            if (startTime >= ElutionTimeToScanNoMap.lastKey()) {
                StartScanNo = ElutionTimeToScanNoMap.lastEntry().getValue();
            } else {
                StartScanNo = ElutionTimeToScanNoMap.lowerEntry(startTime).getValue();
            }
        }
        if (endTime <= ElutionTimeToScanNoMap.firstKey()) {
            EndScanNo = ElutionTimeToScanNoMap.firstEntry().getValue();
        } else {
            if (endTime >= ElutionTimeToScanNoMap.lastKey()) {
                EndScanNo = ElutionTimeToScanNoMap.lastEntry().getValue();
            } else {
                EndScanNo = ElutionTimeToScanNoMap.higherEntry(endTime).getValue();
            }
        }
        Iterator<Entry<Integer, Long>> iter = ScanIndex.entrySet().iterator();

        Entry<Integer, Long> ent = iter.next();
        Long currentIdx = ent.getValue();
        int nextScanNo = ent.getKey();

        while (iter.hasNext()) {
            ent = iter.next();
            long startposition = currentIdx;
            long nexposition = ent.getValue();
            int currentScanNo = nextScanNo;
            nextScanNo = ent.getKey();
            currentIdx = nexposition;

            if (currentScanNo >= StartScanNo && currentScanNo <= EndScanNo && IncludedMSlevel.contains(MsLevelList.get(currentScanNo))) {
                byte[] buffer = new byte[(int) (nexposition - startposition)];
                RandomAccessFile fileHandler = new RandomAccessFile(filename, "r");
                fileHandler.seek(startposition);
                fileHandler.read(buffer, 0, (int) (nexposition - startposition));
                fileHandler.close();
                String xmltext = new String(buffer);
                if (ent.getKey() == Integer.MAX_VALUE) {
                    xmltext = xmltext.replaceAll("</msRun>", "");
                    buffer = null;
                }
                boolean ReadPeak = true;
                if (MsLevelList.get(currentScanNo) == 1) {
                    ReadPeak = MS1Peak;
                }
                if (MsLevelList.get(currentScanNo) == 2) {
                    ReadPeak = MS2Peak;
                }
                MzXMLthreadUnit unit = new MzXMLthreadUnit(xmltext, parameter, datatype, progress, ReadPeak);
                ScanList.add(unit);
                buffer = null;
                xmltext = null;
                fileHandler = null;
            }
        }
        //progress.SetTotal(ScanList.size());
        //Thread thread = new Thread(progress);
        //thread.start();
        for (MzXMLthreadUnit unit : ScanList) {
            executorPool.execute(unit);
        }
        executorPool.shutdown();
////        while (!executorPool.isTerminated()) {
////        }
        try {
            executorPool.awaitTermination(Long.MAX_VALUE, TimeUnit.NANOSECONDS);
        } catch (InterruptedException e) {
            Logger.getRootLogger().info("interrupted..");
        }
        //thread.stop();
        //thread = null;
        //progress.ClearMSG();
        //progress = null;
        for (MzXMLthreadUnit result : ScanList) {
            scanCollection.AddScan(result.scan);
        }
        executorPool = null;
        ScanList.clear();
        ScanList = null;
        ent = null;
        iter = null;

        System.gc();
        //System.out.print(".....done\n");
        System.out.print("Memory usage after loading scans:" + Math.round((Runtime.getRuntime().totalMemory() - Runtime.getRuntime().freeMemory()) / 1048576) + "MB\n");
    }

    public void GetAllScanCollectionMS2Only(boolean MS2Included, boolean MS2Peak) throws InterruptedException, ExecutionException, IOException {

        //System.out.print("Multithreading: "+NoCPUs +" processors (Memory usage:"+ Math.round((Runtime.getRuntime().totalMemory() -Runtime.getRuntime().freeMemory())/1048576)+"MB)\n");
        //System.out.print("...Reading all scans....");
        System.out.print("Memory usage before loading scans:" + Math.round((Runtime.getRuntime().totalMemory() - Runtime.getRuntime().freeMemory()) / 1048576) + "MB\n");

        ArrayList<Integer> IncludedMSlevel = new ArrayList<>();

        if (MS2Included) {
            IncludedMSlevel.add(2);
        }

        ExecutorService executorPool = null;
        List<MzXMLthreadUnit> ScanList = new ArrayList<>();

        //UpdateProcess progress = new UpdateProcess();
        UpdateProcess progress = null;

        executorPool = Executors.newFixedThreadPool(NoCPUs);

        Iterator<Entry<Integer, Long>> iter = ScanIndex.entrySet().iterator();

        Entry<Integer, Long> ent = iter.next();
        Long currentIdx = ent.getValue();
        int nextScanNo = ent.getKey();

        while (iter.hasNext()) {
            ent = iter.next();
            long startposition = currentIdx;
            long nexposition = ent.getValue();
            int currentScanNo = nextScanNo;
            nextScanNo = ent.getKey();
            currentIdx = nexposition;

            if (IncludedMSlevel.contains(MsLevelList.get(currentScanNo))) {
                byte[] buffer = new byte[(int) (nexposition - startposition)];
                RandomAccessFile fileHandler = new RandomAccessFile(filename, "r");
                fileHandler.seek(startposition);
                fileHandler.read(buffer, 0, (int) (nexposition - startposition));
                fileHandler.close();
                String xmltext = new String(buffer);
                if (ent.getKey() == Integer.MAX_VALUE) {
                    xmltext = xmltext.replaceAll("</msRun>", "");
                    buffer = null;
                }
                boolean ReadPeak = true;
                if (MsLevelList.get(currentScanNo) == 2) {
                    ReadPeak = MS2Peak;
                }
                MzXMLthreadUnit unit = new MzXMLthreadUnit(xmltext, parameter, datatype, progress, ReadPeak);
                ScanList.add(unit);
                buffer = null;
                xmltext = null;
                fileHandler = null;
            }
        }
        //progress.SetTotal(ScanList.size());
        //Thread thread = new Thread(progress);
        //thread.start();
        for (MzXMLthreadUnit unit : ScanList) {
            executorPool.execute(unit);
        }
        executorPool.shutdown();
//        while (!executorPool.isTerminated()) {
//        }
        try {
            executorPool.awaitTermination(Long.MAX_VALUE, TimeUnit.NANOSECONDS);
        } catch (InterruptedException e) {
            Logger.getRootLogger().info("interrupted..");
        }
        //thread.stop();
        //thread = null;
        //progress.ClearMSG();
        //progress = null;
        for (MzXMLthreadUnit result : ScanList) {
            scanCollection.AddScan(result.scan);
        }
        executorPool = null;
        ScanList.clear();
        ScanList = null;
        ent = null;
        iter = null;
        System.gc();
        //System.out.print(".....done\n");
        System.out.print("Memory usage after loading scans:" + Math.round((Runtime.getRuntime().totalMemory() - Runtime.getRuntime().freeMemory()) / 1048576) + "MB\n");
    }

    public ScanData GetSingleRawScan(int scanNO) throws FileNotFoundException, IOException, Exception {
        long startposition = ScanIndex.get(scanNO);
        long nexposition = ScanIndex.ceilingEntry(scanNO + 1).getValue();
        byte[] buffer = new byte[(int) (nexposition - startposition)];
        RandomAccessFile fileHandler = new RandomAccessFile(filename, "r");
        fileHandler.seek(startposition);
        fileHandler.read(buffer, 0, (int) (nexposition - startposition));
        String xmltext = new String(buffer);
        xmltext = xmltext.replaceAll("</msRun>", "");

        mzXMLReadUnit read = new mzXMLReadUnit(xmltext);
        ScanData scan = read.Parse();
        xmltext = null;
        read = null;
        fileHandler.close();
        return scan;
    }

    public ScanData GetSingleScanByScanNumberAndRelease(int scanNO) throws FileNotFoundException, IOException, Exception {
        if (!scanCollection.ScanAdded(scanNO)) {
            long startposition = ScanIndex.get(scanNO);
            long nexposition = ScanIndex.ceilingEntry(scanNO + 1).getValue();
            byte[] buffer = new byte[(int) (nexposition - startposition)];
            RandomAccessFile fileHandler = new RandomAccessFile(filename, "r");
            fileHandler.seek(startposition);
            fileHandler.read(buffer, 0, (int) (nexposition - startposition));
            String xmltext = new String(buffer);
            xmltext = xmltext.replaceAll("</msRun>", "");

            MzXMLthreadUnit unit = new MzXMLthreadUnit(xmltext, parameter, datatype);
            unit.run();
            fileHandler.close();
            return unit.scan;
        }
        return scanCollection.GetScan(scanNO);
    }

    public void GetScanCollectionRawByScanNos(ArrayList<Integer> ScanNos) throws InterruptedException, ExecutionException, IOException {
        //System.out.print("Multithreading: "+NoCPUs +" processors (Memory usage:"+ Math.round((Runtime.getRuntime().totalMemory() -Runtime.getRuntime().freeMemory())/1048576)+"MB)\n");
        //System.out.print("...Reading all scans....");
        System.out.print("Memory usage before loading scans:" + Math.round((Runtime.getRuntime().totalMemory() - Runtime.getRuntime().freeMemory()) / 1048576) + "MB\n");

        ExecutorService executorPool = null;
        List<MzXMLthreadUnit> ScanList = new ArrayList<>();

        //UpdateProcess progress = new UpdateProcess();

        executorPool = Executors.newFixedThreadPool(NoCPUs);
        for (Integer scanNO : ScanNos) {
            long startposition = ScanIndex.get(scanNO);
            long nexposition = ScanIndex.ceilingEntry(scanNO + 1).getValue();
            byte[] buffer = new byte[(int) (nexposition - startposition)];
            RandomAccessFile fileHandler = new RandomAccessFile(filename, "r");
            fileHandler.seek(startposition);
            fileHandler.read(buffer, 0, (int) (nexposition - startposition));
            String xmltext = new String(buffer);
            if (nexposition == Integer.MAX_VALUE) {
                xmltext = xmltext.replaceAll("</msRun>", "");
            }
            MzXMLthreadUnit unit = new MzXMLthreadUnit(xmltext, parameter, datatype);
            ScanList.add(unit);
            buffer = null;
            xmltext = null;
            fileHandler.close();
        }
        //progress.SetTotal(ScanList.size());
        //Thread thread = new Thread(progress);
        //thread.start();
        for (MzXMLthreadUnit unit : ScanList) {
            executorPool.execute(unit);
        }
        executorPool.shutdown();
//        while (!executorPool.isTerminated()) {
//        }
        try {
            executorPool.awaitTermination(Long.MAX_VALUE, TimeUnit.NANOSECONDS);
        } catch (InterruptedException e) {
            Logger.getRootLogger().info("interrupted..");
        }
        //thread.stop();
        //thread = null;
        //progress.ClearMSG();
        //progress = null;
        for (MzXMLthreadUnit result : ScanList) {
            scanCollection.AddScan(result.scan);
        }
        executorPool = null;
        ScanList.clear();
        ScanList = null;

        System.gc();
        //System.out.print(".....done\n");
        System.out.print("Memory usage after loading scans:" + Math.round((Runtime.getRuntime().totalMemory() - Runtime.getRuntime().freeMemory()) / 1048576) + "MB\n");
    }

    public ScanData GetSingleScanByScanNumber(int scanNO) throws FileNotFoundException, IOException, Exception {
        if (!scanCollection.ScanAdded(scanNO)) {
            long startposition = ScanIndex.get(scanNO);
            long nexposition = ScanIndex.ceilingEntry(scanNO + 1).getValue();
            byte[] buffer = new byte[(int) (nexposition - startposition)];
            RandomAccessFile fileHandler = new RandomAccessFile(filename, "r");
            fileHandler.seek(startposition);
            fileHandler.read(buffer, 0, (int) (nexposition - startposition));
            String xmltext = new String(buffer);
            if (nexposition == Integer.MAX_VALUE) {
                xmltext = xmltext.replaceAll("</msRun>", "");
            }
            MzXMLthreadUnit unit = new MzXMLthreadUnit(xmltext, parameter, datatype);
            unit.run();
            scanCollection.AddScan(unit.scan);
            buffer = null;
            fileHandler.close();
        }
        return scanCollection.GetScan(scanNO);
    }
}
