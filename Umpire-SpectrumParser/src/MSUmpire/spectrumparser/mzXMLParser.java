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
import java.util.Comparator;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map.Entry;
import java.util.TreeMap;
import java.util.concurrent.*;
import java.util.zip.DataFormatException;
import javax.xml.parsers.ParserConfigurationException;
import org.apache.commons.io.FilenameUtils;
import org.apache.commons.lang.exception.ExceptionUtils;
import org.apache.log4j.Logger;
import org.nustaq.serialization.FSTObjectInput;
import org.nustaq.serialization.FSTObjectOutput;
import org.xml.sax.SAXException;


/* * 
 * constructor To change this template, choose Tools | Templates and open the
 * template in the editor.
 */
/**
 *
 * @author Chih-Chiang Tsou
 */
public final class mzXMLParser {

    public ScanCollection scanCollection = null;
    public String filename;
    public InstrumentParameter parameter;
    public int TotalScan;
    public TreeMap<Integer, Long> ScanIndex;
    public int NoCPUs = 4;
    public SpectralDataType.DataType datatype;
    public DIA_Setting dIA_Setting;
    public TreeMap<XYData, ScanCollection> DIAMS2Scans = new TreeMap<>(new Comparator<XYData>() {
        @Override
        public int compare(XYData o1, XYData o2) {
            return Float.compare(o1.getX(), o2.getX());
        }
    });
     public TreeMap<XYData, ScanCollection> MS1WindowScans = new TreeMap<>(new Comparator<XYData>() {
        @Override
        public int compare(XYData o1, XYData o2) {
            return Float.compare(o1.getX(), o2.getX());
        }
    });
    public TreeMap<Integer, Integer> MsLevelList;
    private TreeMap<Float, Integer> ElutionTimeToScanNoMap;
    private HashMap<Integer, Float> ScanToElutionTime;
    public int NoMS1Scans = 0;

    public mzXMLParser(String filename, InstrumentParameter parameter, SpectralDataType.DataType datatype, DIA_Setting dIA_Setting, int NoCPUs) throws Exception {
        this.filename = filename;
        this.ElutionTimeToScanNoMap = new TreeMap<>();
        this.dIA_Setting = dIA_Setting;
        this.MsLevelList = new TreeMap<>();
        this.scanCollection = new ScanCollection(parameter.Resolution);
        this.scanCollection.Filename = filename;
        this.parameter = parameter;
        this.datatype = datatype;
        this.NoCPUs = NoCPUs;
        ReadElutionAndScanIndex();
    }

    public float GetMS1CycleTime() {
        return (ElutionTimeToScanNoMap.lastKey() - ElutionTimeToScanNoMap.firstKey()) / NoMS1Scans;
    }

    public HashMap<Integer, Float> GetScanElutionTimeMap() {
        return ScanToElutionTime;
    }

    public void ReadElutionAndScanIndex() throws Exception {
        //long startRead =System.nanoTime();
        if (!FSScanPosRead()) {
            //long start =System.nanoTime();
            ParseIndex();
            //System.out.printf("ParseIndex() took: %.0f ms\n", (System.nanoTime() - start)/1e6);
            //start =System.nanoTime();
            WriteIndexSerialization();
            //System.out.printf("WriteIndexSerialization() took: %.0f ms\n", (System.nanoTime() - start)/1e6);
        }
        //System.out.printf("FSScanPosRead() took: %.0f ms\n", (System.nanoTime() - startRead)/1e6);
        //startRead =System.nanoTime();
        if (!FSElutionIndexRead()) {
            //long start =System.nanoTime();
            ParseElutionIndex();
            //System.out.printf("ParseElutionIndex() took: %.0f ms\n", (System.nanoTime() - start)/1e6);
            //start =System.nanoTime();
            FSElutionIndexWrite();
            //System.out.printf("FSElutionIndexWrite() took: %.0f ms\n", (System.nanoTime() - start)/1e6);
        }
        //System.out.printf("FSElutionIndexRead() took: %.0f ms\n", (System.nanoTime() - startRead)/1e6);
    }

    private void WriteIndexSerialization() {
        FSScanIdxWrite();
    }

    private void FSScanIdxWrite() {
        try {
            Logger.getRootLogger().debug("Writing ScanPos to file:" + FilenameUtils.removeExtension(filename) + ".ScanPosFS..");
            FileOutputStream fout = new FileOutputStream(FilenameUtils.removeExtension(filename) + ".ScanPosFS", false);
            FSTObjectOutput oos = new FSTObjectOutput(fout);
            oos.writeObject(ScanIndex);
            oos.close();
            fout.close();
        } catch (Exception ex) {
            Logger.getRootLogger().error(ExceptionUtils.getStackTrace(ex));
        }
    }
    
    private boolean FSScanPosRead() {
        if (!new File(FilenameUtils.removeExtension(filename) + ".ScanPosFS").exists()) {
            return false;
        }
        try {
            Logger.getRootLogger().debug("Reading ScanPos:" + FilenameUtils.removeExtension(filename) + ".ScanPosFS...");
            FileInputStream fileIn = new FileInputStream(FilenameUtils.removeExtension(filename) + ".ScanPosFS");
            FSTObjectInput in = new FSTObjectInput(fileIn);
            ScanIndex = (TreeMap<Integer, Long>) in.readObject();
            TotalScan = ScanIndex.size();
            in.close();
            fileIn.close();

        } catch (Exception ex) {
            Logger.getRootLogger().debug("ScanIndex serialization file failed");
            //Logger.getRootLogger().error(ExceptionUtils.getStackTrace(ex));
            return false;
        }
        return true;
    }
    
    public void ParseIndex() throws FileNotFoundException, IOException {
        TotalScan = 0;
        ScanIndex = new TreeMap<>();
        try (RandomAccessFile fileHandler = new RandomAccessFile(filename, "r")) {
            StringBuilder sb = new StringBuilder();

            String CurrentLine = "";
            long currentLastPt = fileHandler.length() - 1;
            boolean indexexist=false;
            int linecount=0;
            while (!(CurrentLine.trim().startsWith("<index name=") | CurrentLine.trim().startsWith("</msRun>"))) {
                //Read backward
                for (long filePointer = currentLastPt; filePointer != -1; filePointer--) {
                    fileHandler.seek(filePointer);
                    int readByte = fileHandler.readByte();
                    if (readByte == 0xA) {
                        if (filePointer == currentLastPt) {
                            continue;
                        } else {
                            currentLastPt = filePointer;
                            break;
                        }
                    } else if (readByte == 0xD) {
                        if (filePointer == currentLastPt - 1) {
                            continue;
                        } else {
                            currentLastPt = filePointer;
                            break;
                        }
                    }
                    sb.append((char) readByte);
                }
                linecount++;
                CurrentLine = sb.reverse().toString();
                sb = new StringBuilder();

                if (CurrentLine.trim().startsWith("</index>")) {
                    indexexist = true;
                }

                if (!indexexist && linecount > 10) {
                    fileHandler.close();
                    Logger.getRootLogger().debug("File : " + filename + " doesn't have index. the processing will stop.");
                    System.exit(1);
                }

                if (CurrentLine.trim().startsWith("<offset id")) {
                    int scanNo = Integer.parseInt(CurrentLine.substring(CurrentLine.indexOf("<offset id=\"") + 12).split("\"")[0]);
                    long index = (long) Long.parseLong(CurrentLine.substring(CurrentLine.indexOf(">") + 1, CurrentLine.indexOf("</offset>")));
                    if (index < 0) {
                        index = index + 2147483647l + 2147483648l;
                    }
                    if (ScanIndex.containsKey(scanNo + 1) && ScanIndex.get(scanNo + 1) == index) {
                        Logger.getRootLogger().debug("File : " + filename + " index is not correct, ScanNo:" + scanNo + " and " + scanNo + 1 + " have same index");
                        Logger.getRootLogger().debug("Please use indexmzXML from  TPP package to fix incorrect index of the mzXML file.");
                        Logger.getRootLogger().debug("command: indexmzXML filename.mzXML");
                        System.exit(1);
                    }
                    ScanIndex.put(scanNo, index);
                } else if (CurrentLine.trim().startsWith("<indexOffset>")) {
                    long IndexEnd = (long) Long.parseLong(CurrentLine.substring(CurrentLine.indexOf("<indexOffset>") + 13, CurrentLine.indexOf("</indexOffset>")));
                    if (IndexEnd < 0) {
                        IndexEnd = IndexEnd + 2147483647l + 2147483648l;
                    }
                    ScanIndex.put(Integer.MAX_VALUE, IndexEnd);
                }
            }
            TotalScan = ScanIndex.size();
            sb = null;
            fileHandler.close();
        }
    }

    public void ParseElutionIndex() throws Exception {
        ElutionTimeToScanNoMap = new TreeMap<>();
        scanCollection.ElutionTimeToScanNoMap = new TreeMap<>();
        ScanToElutionTime = new HashMap<>();
        MsLevelList = new TreeMap<>();
        if(ScanIndex==null | ScanIndex.isEmpty()){
            return;
        }

        boolean ReadDIAWindow = false;
        if (datatype != SpectralDataType.DataType.DDA) {
            if (datatype != SpectralDataType.DataType.DIA_V_Window) {
                dIA_Setting.DIAWindows = new TreeMap<>();
            }
            ReadDIAWindow = true;
        }
        try (RandomAccessFile fileHandler = new RandomAccessFile(filename, "r")) {
            Iterator<Entry<Integer, Long>> iter = ScanIndex.entrySet().iterator();
            Long currentIdx = iter.next().getValue();
            while (iter.hasNext()) {
                long startposition = currentIdx;
                long nexposition = iter.next().getValue();
                currentIdx = nexposition;
                fileHandler.seek(startposition);
                
//                if(startposition>=nexposition){
//                    System.out.println("");
//                }
                
                byte[] bufr = new byte[(int) (nexposition - startposition)];
                fileHandler.readFully(bufr, 0, (int) (nexposition - startposition));

                String temp = new String(bufr);

                float nowtime = 0f;
                int scanno = 0;
                int mslevel = 0;
                //float precursorF=0f;
                if (!temp.contains("<scan")) {
                    fileHandler.close();
                    return;
                }

                if (temp.contains("<scan num=") && (temp.contains("retentionTime=\"PT"))) {
                    String substr = temp.substring(temp.indexOf("<scan num=") + 11);
                    scanno = Integer.parseInt(substr.substring(0, substr.indexOf("\"")));
                    
                    
                    nowtime = Float.parseFloat(temp.substring(temp.indexOf("retentionTime=\"PT") + 17).split("S\"")[0]);
                    
                    mslevel = Integer.parseInt(temp.substring(temp.indexOf("msLevel=") + 9, temp.indexOf("msLevel=") + 10));
                    if (temp.contains("scanType=\"calibration\"")) {
                        mslevel = -1;
                    }
                    if (mslevel == 1) {
                        NoMS1Scans++;                        
                        if (temp.contains("scanType=\"SIM\"") && datatype == SpectralDataType.DataType.WiSIM) {
                            int startidx = temp.indexOf("lowMz=\"") + 7;
                            int stopidx = startidx + 1;
                            for (int i = startidx + 1; i < temp.length(); i++) {
                                if (temp.charAt(i) == '\"') {
                                    stopidx = i;
                                    break;
                                }
                            }
                            float lowmz = Float.parseFloat(temp.substring(startidx, stopidx));
                            startidx = temp.indexOf("highMz=\"") + 8;
                            stopidx = startidx + 1;
                            for (int i = startidx + 1; i < temp.length(); i++) {
                                if (temp.charAt(i) == '\"') {
                                    stopidx = i;
                                    break;
                                }
                            }
                            float highmz = Float.parseFloat(temp.substring(startidx, stopidx));
                            for (XYData MS1win : dIA_Setting.MS1Windows.keySet()) {
                                if (MS1win.getX() <= lowmz && MS1win.getY() >= highmz) {
                                    dIA_Setting.MS1Windows.get(MS1win).add(scanno);
                                }
                            }
                        }                        
                    }
                    //precursorF= Float.parseFloat(temp.substring(temp.indexOf("precursorIntensity=") + 20, temp.indexOf("\" precursorCharge")));
                    
                    if (datatype != SpectralDataType.DataType.DDA) {
                        if (ReadDIAWindow && mslevel == 2) {
                            if (datatype == SpectralDataType.DataType.MSX) {
                                substr = temp;
                                while (substr.contains("</precursorMz>")) {
                                    int stopidx = substr.indexOf("</precursorMz>");
                                    int startidx = 0;
                                    for (int i = stopidx; i > 0; i--) {
                                        if (substr.charAt(i) == '>') {
                                            startidx = i + 1;
                                            break;
                                        }
                                    }
                                    float precursormz = Float.parseFloat(substr.substring(startidx, stopidx));

                                    startidx = substr.indexOf("windowWideness=\"") + 16;
                                    stopidx = startidx + 1;
                                    for (int i = startidx + 1; i < substr.length(); i++) {
                                        if (substr.charAt(i) == '\"') {
                                            stopidx = i;
                                            break;
                                        }
                                    }
                                    float windowwideness = Float.parseFloat(substr.substring(startidx, stopidx));
                                    float Loffset = windowwideness / 2f;
                                    float Roffset = windowwideness / 2f;

                                    if (!dIA_Setting.DIAWindows.containsKey(new XYData(precursormz - Loffset, precursormz + Roffset))) {
                                        ArrayList<Integer> scanList = new ArrayList<>();
                                        dIA_Setting.DIAWindows.put(new XYData(precursormz - Loffset, precursormz + Roffset), scanList);
                                    }
                                    dIA_Setting.DIAWindows.get(new XYData(precursormz - Loffset, precursormz + Roffset)).add(scanno);
                                    substr = substr.substring(substr.indexOf("</precursorMz>") + 14);
                                }
                            } else if (datatype == SpectralDataType.DataType.DIA_F_Window || datatype == SpectralDataType.DataType.pSMART || datatype == SpectralDataType.DataType.WiSIM) {
                                int stopidx = temp.indexOf("</precursorMz>");
                                if (stopidx == -1) {
                                    Logger.getRootLogger().error("Parsing </precursorMz> failed. scan number :" + scanno);                                    
                                    System.exit(3);
                                }
                                int startidx = 0;
                                for (int i = stopidx; i > 0; i--) {
                                    if (temp.charAt(i) == '>') {
                                        startidx = i + 1;
                                        break;
                                    }
                                }
                                float precursormz = Float.parseFloat(temp.substring(startidx, stopidx));
                                float Loffset = (dIA_Setting.F_DIA_WindowSize + 1) * 0.2f;
                                float Roffset = (dIA_Setting.F_DIA_WindowSize + 1) * 0.8f;
//                                float Loffset = 5f;
//                                float Roffset = 21f;

                                if (temp.contains("windowWideness=\"")) {
                                    startidx = temp.indexOf("windowWideness=\"") + 16;
                                    stopidx = startidx + 1;
                                    for (int i = startidx + 1; i < temp.length(); i++) {
                                        if (temp.charAt(i) == '\"') {
                                            stopidx = i;
                                            break;
                                        }
                                    }
                                    float windowwideness = Float.parseFloat(temp.substring(startidx, stopidx));
                                    Loffset = windowwideness / 2f;
                                    Roffset = windowwideness / 2f;
                                }

                                if (!dIA_Setting.DIAWindows.containsKey(new XYData(precursormz - Loffset, precursormz + Roffset))) {
                                    ArrayList<Integer> scanList = new ArrayList<>();
                                    dIA_Setting.DIAWindows.put(new XYData(precursormz - Loffset, precursormz + Roffset), scanList);
                                }
                                dIA_Setting.DIAWindows.get(new XYData(precursormz - Loffset, precursormz + Roffset)).add(scanno);
                            } else if (datatype == SpectralDataType.DataType.DIA_V_Window) {
                                int stopidx = temp.indexOf("</precursorMz>");
                                int startidx = 0;
                                for (int i = stopidx; i > 0; i--) {
                                    if (temp.charAt(i) == '>') {
                                        startidx = i + 1;
                                        break;
                                    }
                                }
                                float precursormz = Float.parseFloat(temp.substring(startidx, stopidx));
                                for (XYData window : dIA_Setting.DIAWindows.keySet()) {
                                    if (window.getX() <= precursormz && window.getY() >= precursormz) {
                                        dIA_Setting.DIAWindows.get(window).add(scanno);
                                        break;
                                    }
                                }
                            } else if (datatype == SpectralDataType.DataType.MSe) {
                                float mzlowF = 0f;
                                float mzhighF = 10000f;
                                if (!dIA_Setting.DIAWindows.containsKey(new XYData(mzlowF, mzhighF))) {
                                    ArrayList<Integer> scanList = new ArrayList<>();
                                    dIA_Setting.DIAWindows.put(new XYData(mzlowF, mzhighF), scanList);
                                }
                                dIA_Setting.DIAWindows.get(new XYData(mzlowF, mzhighF)).add(scanno);
                            }
                        }
                    }
                } else {
                    Logger.getRootLogger().error("index of mzXML error");
                    System.exit(1);
                }
                ElutionTimeToScanNoMap.put(nowtime / 60f, scanno);
                scanCollection.ElutionTimeToScanNoMap.put(nowtime / 60f, scanno);
                ScanToElutionTime.put(scanno, nowtime / 60f);
                MsLevelList.put(scanno, mslevel);
            }
            //WriteElutionIndex();            
            fileHandler.close();
        }
    }

    private void FSElutionIndexWrite() throws IOException {

        try {
            Logger.getRootLogger().debug("Writing RTidx to file:" + FilenameUtils.removeExtension(filename) + ".RTidxFS..");
            FileOutputStream fout = new FileOutputStream(FilenameUtils.removeExtension(filename) + ".RTidxFS", false);
            FSTObjectOutput oos = new FSTObjectOutput(fout);
            oos.writeObject(ElutionTimeToScanNoMap);
            oos.close();
            fout.close();
        } catch (Exception ex) {
            Logger.getRootLogger().error(ExceptionUtils.getStackTrace(ex));
        }

        try {
            Logger.getRootLogger().debug("Writing Scanidx to file:" + FilenameUtils.removeExtension(filename) + ".ScanidxFS..");
            FileOutputStream fout = new FileOutputStream(FilenameUtils.removeExtension(filename) + ".ScanidxFS", false);
            FSTObjectOutput oos = new FSTObjectOutput(fout);
            oos.writeObject(MsLevelList);
            oos.close();
            fout.close();
        } catch (Exception ex) {
            Logger.getRootLogger().error(ExceptionUtils.getStackTrace(ex));
        }

        try {
            Logger.getRootLogger().debug("Writing ScanRT to file:" + FilenameUtils.removeExtension(filename) + ".ScanRTFS..");
            FileOutputStream fout = new FileOutputStream(FilenameUtils.removeExtension(filename) + ".ScanRTFS", false);
            FSTObjectOutput oos = new FSTObjectOutput(fout);
            oos.writeObject(ScanToElutionTime);
            oos.close();
            fout.close();
        } catch (Exception ex) {
            Logger.getRootLogger().error(ExceptionUtils.getStackTrace(ex));
        }

        if (datatype !=  SpectralDataType.DataType.DDA) {
            try {
                Logger.getRootLogger().debug("Writing DIAWindows to file:" + FilenameUtils.removeExtension(filename) + ".DIAWindowsFS..");
                FileOutputStream fout = new FileOutputStream(FilenameUtils.removeExtension(filename) + ".DIAWindowsFS", false);
                FSTObjectOutput oos = new FSTObjectOutput(fout);
                oos.writeObject(dIA_Setting.DIAWindows);
                oos.close();
                fout.close();
            } catch (Exception ex) {
                Logger.getRootLogger().error(ExceptionUtils.getStackTrace(ex));
            }
        }
         if (datatype ==  SpectralDataType.DataType.WiSIM) {
            try {
                Logger.getRootLogger().debug("Writing MS1 windows to file:" + FilenameUtils.removeExtension(filename) + ".MS1WindowsFS..");
                FileOutputStream fout = new FileOutputStream(FilenameUtils.removeExtension(filename) + ".MS1WindowsFS", false);
                FSTObjectOutput oos = new FSTObjectOutput(fout);
                oos.writeObject(dIA_Setting.MS1Windows);
                oos.close();
                fout.close();
            } catch (Exception ex) {
                Logger.getRootLogger().error(ExceptionUtils.getStackTrace(ex));
            }
        }
    }

    
//    //<editor-fold defaultstate="collapsed" desc="Old code for reading and writing index file">
  //  private void WriteElutionIndex() throws IOException {
//        BufferedWriter writer = new BufferedWriter(new FileWriter(FilenameUtils.removeExtension(filename) + ".RTidx"));
//        BufferedWriter writer2 = new BufferedWriter(new FileWriter(FilenameUtils.removeExtension(filename) + ".Scanidx"));
//
//        for (Float entry : ElutionTimeToScanNoMap.keySet()) {
//            writer.write(entry + "," + ElutionTimeToScanNoMap.get(entry) + "\n");
//        }
//        for (Integer entry : MsLevelList.keySet()) {
//            writer2.write(entry + "," + MsLevelList.get(entry) + "\n");
//        }
//
//        writer.close();
//        writer2.close();
//        if (datatype == SpectralDataType.DataType.DIA_V_Window || datatype == SpectralDataType.DataType.DIA_F_Window || datatype == SpectralDataType.DataType.MSe || datatype == SpectralDataType.DataType.MSX) {
//            BufferedWriter writer3 = new BufferedWriter(new FileWriter(FilenameUtils.removeExtension(filename) + ".DIAWindows"));
//            for (XYData window : dIA_Setting.DIAWindows.keySet()) {
//                String ScanString = "";
//                ArrayList<Integer> scanlist = dIA_Setting.DIAWindows.get(window);
//                for (Integer scan : scanlist) {
//                    ScanString += scan + "_";
//                }
//                writer3.write(window.getX() + "," + window.getY() + "," + ScanString + "\n");
//            }
//            writer3.close();
//        }
//    }


//    private void ReadElutionIndexFile() throws NumberFormatException, FileNotFoundException, IOException {
//        BufferedReader reader = new BufferedReader(new FileReader(FilenameUtils.removeExtension(filename) + ".RTidx"));
//        BufferedReader reader2 = new BufferedReader(new FileReader(FilenameUtils.removeExtension(filename) + ".Scanidx"));
//        String line = "";
//        while ((line = reader.readLine()) != null) {
//            ElutionTimeToScanNoMap.put(Float.parseFloat(line.split(",")[0]), Integer.parseInt(line.split(",")[1]));
//            scanCollection.ElutionTimeToScanNoMap.put(Float.parseFloat(line.split(",")[0]), Integer.parseInt(line.split(",")[1]));
//        }
//        line = "";
//        while ((line = reader2.readLine()) != null) {
//            int value = Integer.parseInt(line.split(",")[1]);
//            MsLevelList.put(Integer.parseInt(line.split(",")[0]), value);
//            if (value == 1) {
//                NoMS1Scans++;
//            }
//        }
//        if (datatype != SpectralDataType.DataType.DDA) {
//            dIA_Setting.DIAWindows = new TreeMap<>(new Comparator<XYData>() {
//                @Override
//                public int compare(XYData o1, XYData o2) {
//                    return Float.compare(o1.getX(), o2.getX());
//                }
//            });
//            BufferedReader reader3 = new BufferedReader(new FileReader(FilenameUtils.removeExtension(filename) + ".DIAWindows"));
//            line = "";
//            while ((line = reader3.readLine()) != null) {
//                XYData DIAwin = new XYData(Float.parseFloat(line.split(",")[0]), Float.parseFloat(line.split(",")[1]));
//                ArrayList<Integer> scanList = new ArrayList<>();
//                for (int i = 0; i < line.split(",")[2].split("_").length - 1; i++) {
//                    scanList.add(Integer.parseInt(line.split(",")[2].split("_")[i]));
//                }
//                dIA_Setting.DIAWindows.put(DIAwin, scanList);
//            }
//        }
//        line = null;
//        return;
//    }
//</editor-fold>
    
    private boolean FSElutionIndexRead() throws NumberFormatException, FileNotFoundException, IOException {
        if (!new File(FilenameUtils.removeExtension(filename) + ".RTidxFS").exists()) {
            return false;
        }
        if (!new File(FilenameUtils.removeExtension(filename) + ".ScanRTFS").exists()) {
            return false;
        }
        if (!new File(FilenameUtils.removeExtension(filename) + ".ScanidxFS").exists()) {
            return false;
        }
        
                
        try {
            Logger.getRootLogger().debug("Reading RTidx:" + FilenameUtils.removeExtension(filename) + ".RTidxFS...");
            FileInputStream fileIn = new FileInputStream(FilenameUtils.removeExtension(filename) + ".RTidxFS");
            FSTObjectInput in = new FSTObjectInput(fileIn);
            ElutionTimeToScanNoMap = (TreeMap<Float, Integer>) in.readObject();
            scanCollection.ElutionTimeToScanNoMap = ElutionTimeToScanNoMap;
            in.close();
            fileIn.close();
        } catch (Exception ex) {
            //i.printStackTrace();
            Logger.getRootLogger().debug("RTidxFS serialization file failed");
            //Logger.getRootLogger().error(ExceptionUtils.getStackTrace(ex));
            return false;
        }

        //////////
        try {
            Logger.getRootLogger().debug("Reading ScanRT:" + FilenameUtils.removeExtension(filename) + ".ScanRTFS...");
            FileInputStream fileIn = new FileInputStream(FilenameUtils.removeExtension(filename) + ".ScanRTFS");
            FSTObjectInput in = new FSTObjectInput(fileIn);
            ScanToElutionTime = (HashMap<Integer, Float>) in.readObject();
            in.close();
            fileIn.close();
        } catch (Exception ex) {
            //i.printStackTrace();
            Logger.getRootLogger().debug("ScanRTFS serialization file failed");
            //Logger.getRootLogger().error(ExceptionUtils.getStackTrace(ex));
            return false;
        }
        ///////////

        try {
            Logger.getRootLogger().debug("Reading Scanidx:" + FilenameUtils.removeExtension(filename) + ".ScanidxFS...");
            FileInputStream fileIn = new FileInputStream(FilenameUtils.removeExtension(filename) + ".ScanidxFS");
            FSTObjectInput in = new FSTObjectInput(fileIn);
            MsLevelList = (TreeMap<Integer, Integer>) in.readObject();
            for (Integer value : MsLevelList.values()) {
                if (value == 1) {
                    NoMS1Scans++;
                }
            }
            in.close();
            fileIn.close();
        } catch (Exception ex) {
            //i.printStackTrace();
            Logger.getRootLogger().debug("ScanidxFS serialization file failed");
            //Logger.getRootLogger().error(ExceptionUtils.getStackTrace(ex));
            return false;
        }

        if (datatype != SpectralDataType.DataType.DDA) {
            if (!new File(FilenameUtils.removeExtension(filename) + ".DIAWindowsFS").exists()) {
                return false;
            }
            try {
                Logger.getRootLogger().debug("Reading DIAWindows:" + FilenameUtils.removeExtension(filename) + ".DIAWindowsFS...");
                FileInputStream fileIn = new FileInputStream(FilenameUtils.removeExtension(filename) + ".DIAWindowsFS");
                FSTObjectInput in = new FSTObjectInput(fileIn);
                dIA_Setting.DIAWindows = (TreeMap<XYData, ArrayList<Integer>>) in.readObject();
                in.close();
                fileIn.close();
            } catch (Exception ex) {
                //i.printStackTrace();
                Logger.getRootLogger().debug("DIAWindowsFS serialization file failed");
                //Logger.getRootLogger().error(ExceptionUtils.getStackTrace(ex));
                return false;
            }
        }
        
         if (datatype == SpectralDataType.DataType.WiSIM) {
            if (!new File(FilenameUtils.removeExtension(filename) + ".MS1WindowsFS").exists()) {
                return false;
            }
            try {
                Logger.getRootLogger().debug("Reading MS1 Windows:" + FilenameUtils.removeExtension(filename) + ".MS1WindowsFS...");
                FileInputStream fileIn = new FileInputStream(FilenameUtils.removeExtension(filename) + ".MS1WindowsFS");
                FSTObjectInput in = new FSTObjectInput(fileIn);
                dIA_Setting.MS1Windows = (TreeMap<XYData, ArrayList<Integer>>) in.readObject();
                in.close();
                fileIn.close();
            } catch (Exception ex) {
                //i.printStackTrace();
                Logger.getRootLogger().debug("MS1WindowsFS serialization file failed");
                //Logger.getRootLogger().error(ExceptionUtils.getStackTrace(ex));
                return false;
            }
        }
        return true;
    }

    private void GetScanCollectionSingleThread(int msLevel) throws IOException, ParserConfigurationException, SAXException, DataFormatException {
        Iterator<Entry<Integer, Long>> iter = ScanIndex.entrySet().iterator();

        Long currentIdx = iter.next().getValue();
        while (iter.hasNext()) {
            Entry<Integer, Long> ent = iter.next();
            long startposition = currentIdx;
            long nexposition = ent.getValue();
            currentIdx = nexposition;
            byte[] buffer = new byte[(int) (nexposition - startposition)];
            RandomAccessFile fileHandler = new RandomAccessFile(filename, "r");
            fileHandler.seek(startposition);
            fileHandler.read(buffer, 0, (int) (nexposition - startposition));
            String xmltext = new String(buffer);
            if (ent.getKey() == Integer.MAX_VALUE) {
                xmltext = xmltext.replaceAll("</msRun>", "");
            }

            mzXMLReadUnit read = new mzXMLReadUnit(xmltext);
            try {
                ScanData scanData = read.Parse();
                scanCollection.AddScan(scanData);
                Logger.getRootLogger().debug(scanData.Num + ":" + scanData.Data.size() + "\n");
            } catch (Exception ex) {
                Logger.getRootLogger().error("ScanNo: " + ent.getKey() + " error. (" + ex.getMessage() + ")");
                Logger.getRootLogger().error(ExceptionUtils.getStackTrace(ex));
            }
            buffer = null;
            xmltext = null;
            fileHandler.close();
        }
    }

    
    public ScanCollection GetScanCollectionMS1Window(XYData MS1Window, boolean IncludePeak) throws InterruptedException, ExecutionException, IOException {
        if (dIA_Setting == null) {
            Logger.getRootLogger().error("This is not DIA data" + filename);
            return null;
        }
        if (!MS1WindowScans.containsKey(MS1Window)) {
            MS1WindowScans.put(MS1Window, GetScanCollectionMS1Window(MS1Window, IncludePeak, 0f, 999999f));
        }
        return MS1WindowScans.get(MS1Window);
    }

    public int GetScanNoByRT(float RT) {
        int ScanNo = 0;
        if (RT <= ElutionTimeToScanNoMap.firstKey()) {
            ScanNo = ElutionTimeToScanNoMap.firstEntry().getValue();
        } else if (RT >= ElutionTimeToScanNoMap.lastKey()) {
            ScanNo = ElutionTimeToScanNoMap.lastEntry().getValue();
        } else {
            ScanNo = ElutionTimeToScanNoMap.lowerEntry(RT).getValue();
        }
        return ScanNo;
    }
    
    public ScanCollection GetScanCollectionMS1Window(XYData MS1Window, boolean IncludePeak, float startTime, float endTime) throws InterruptedException, ExecutionException, IOException {
        if (dIA_Setting == null) {
            Logger.getRootLogger().error(filename + " is not DIA data");
            return null;
        }
        ScanCollection MS1WindowScanCollection = new ScanCollection(parameter.Resolution);
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

            if (currentScanNo >= StartScanNo && currentScanNo <= EndScanNo && dIA_Setting.MS1Windows.get(MS1Window).contains(currentScanNo)) {
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
            MS1WindowScanCollection.AddScan(result.scan);
            MS1WindowScanCollection.ElutionTimeToScanNoMap.put(result.scan.RetentionTime, result.scan.Num);
        }
        executorPool = null;
        ScanList.clear();
        ScanList = null;
        ent = null;
        iter = null;
        //System.gc();
        //System.out.print(".....done\n");
        //System.out.print("Finished multithreading (Memory usage:"+ Math.round((Runtime.getRuntime().totalMemory() -Runtime.getRuntime().freeMemory())/1048576)+"MB)\n");

        return MS1WindowScanCollection;
    }

    
    public ScanCollection GetScanCollectionDIAMS2(XYData DIAWindow, boolean IncludePeak,float startRT, float endRT) throws InterruptedException, ExecutionException, IOException {
        if (dIA_Setting == null) {
            Logger.getRootLogger().error("This is not DIA data" + filename);
            return null;
        }
        if (!DIAMS2Scans.containsKey(DIAWindow)) {
            DIAMS2Scans.put(DIAWindow, GetScanDIAMS2(DIAWindow, IncludePeak, startRT, endRT));
        }
        return DIAMS2Scans.get(DIAWindow);
    }

    public ScanCollection GetScanDIAMS2(XYData DIAWindow, boolean IncludePeak, float startTime, float endTime) throws InterruptedException, ExecutionException, IOException {
        if (dIA_Setting == null) {
            Logger.getRootLogger().error(filename + " is not DIA data");
            return null;
        }
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

            if (currentScanNo >= StartScanNo && currentScanNo <= EndScanNo && dIA_Setting.DIAWindows.get(DIAWindow).contains(currentScanNo)) {
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

    public void GetAllScanCollectionByMSLabel(boolean MS1Included, boolean MS2Included, boolean MS1Peak, boolean MS2Peak, float startTime, float endTime) throws InterruptedException, ExecutionException, IOException {
        //System.out.print("Multithreading: "+NoCPUs +" processors (Memory usage:"+ Math.round((Runtime.getRuntime().totalMemory() -Runtime.getRuntime().freeMemory())/1048576)+"MB)\n");
        //System.out.print("...Reading all scans....");
        Logger.getRootLogger().debug("Memory usage before loading scans:" + Math.round((Runtime.getRuntime().totalMemory() - Runtime.getRuntime().freeMemory()) / 1048576) + "MB (" + NoCPUs + " threads)");

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
                executorPool.execute(unit);
                buffer = null;
                xmltext = null;
                fileHandler = null;
            }
        }
        //progress.SetTotal(ScanList.size());
        //Thread thread = new Thread(progress);
        //thread.start();
//        for (MzXMLthreadUnit unit : ScanList) {
//            executorPool.execute(unit);
//        }
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
        Logger.getRootLogger().debug("Memory usage after loading scans:" + Math.round((Runtime.getRuntime().totalMemory() - Runtime.getRuntime().freeMemory()) / 1048576) + "MB");
    }

    public void GetAllScanCollectionMS2Only(boolean MS2Included, boolean MS2Peak) throws InterruptedException, ExecutionException, IOException {

        //System.out.print("Multithreading: "+NoCPUs +" processors (Memory usage:"+ Math.round((Runtime.getRuntime().totalMemory() -Runtime.getRuntime().freeMemory())/1048576)+"MB)\n");
        //System.out.print("...Reading all scans....");
        Logger.getRootLogger().debug("Memory usage before loading scans:" + Math.round((Runtime.getRuntime().totalMemory() - Runtime.getRuntime().freeMemory()) / 1048576) + "MB");

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
        Logger.getRootLogger().debug("Memory usage after loading scans:" + Math.round((Runtime.getRuntime().totalMemory() - Runtime.getRuntime().freeMemory()) / 1048576) + "MB");
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
        Logger.getRootLogger().debug("Memory usage before loading scans:" + Math.round((Runtime.getRuntime().totalMemory() - Runtime.getRuntime().freeMemory()) / 1048576) + "MB");

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
        Logger.getRootLogger().debug("Memory usage after loading scans:" + Math.round((Runtime.getRuntime().totalMemory() - Runtime.getRuntime().freeMemory()) / 1048576) + "MB");
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
