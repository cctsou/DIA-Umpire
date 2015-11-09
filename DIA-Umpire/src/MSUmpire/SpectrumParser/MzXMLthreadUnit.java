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
import MSUmpire.BaseDataStructure.ScanData;
import MSUmpire.BaseDataStructure.SpectralDataType;
import java.awt.Color;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.concurrent.locks.ReentrantReadWriteLock;
import java.util.zip.DataFormatException;
import javax.xml.parsers.ParserConfigurationException;
import org.apache.commons.lang.exception.ExceptionUtils;
import org.apache.log4j.Logger;
import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartUtilities;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.chart.plot.XYPlot;
import org.jfree.data.statistics.HistogramDataset;
import org.jfree.data.statistics.HistogramType;
import org.xml.sax.SAXException;

/**
 *
 * @author Chih-Chiang Tsou <chihchiang.tsou@gmail.com>
 */
public class MzXMLthreadUnit implements Runnable {

    ReentrantReadWriteLock Lock = new ReentrantReadWriteLock();
    public ScanData scan;
    private String XMLtext;
    private InstrumentParameter parameter;
    boolean ReadPeak = true;
    SpectralDataType.DataType dataType = SpectralDataType.DataType.DDA;

    public MzXMLthreadUnit(String XMLtext, InstrumentParameter parameter, SpectralDataType.DataType dataType,boolean ReadPeak) {
        this.XMLtext = XMLtext;
        this.parameter = parameter;
        this.ReadPeak = ReadPeak;
        this.dataType = dataType;
    }

    public MzXMLthreadUnit(String XMLtext, InstrumentParameter parameter, SpectralDataType.DataType dataType) {
        this.XMLtext = XMLtext;
        this.parameter = parameter;
        this.dataType = dataType;
    }

    private void DrawIntDis() {
        double[] PeakInt = new double[scan.Data.size()];
        for (int i = 0; i < scan.Data.size(); i++) {
            PeakInt[i] = scan.Data.get(i).getY();
        }
        HistogramDataset histogramDataset = new HistogramDataset();
        histogramDataset.setType(HistogramType.RELATIVE_FREQUENCY);
        histogramDataset.addSeries("Intensity", PeakInt, 100);

        JFreeChart chart = ChartFactory.createHistogram("Peak intensity distribution", "Intensity", "Frequency", histogramDataset, PlotOrientation.VERTICAL, true, false, false);
        XYPlot plot = chart.getXYPlot();
        plot.setBackgroundPaint(Color.white);
        plot.setDomainGridlinePaint(Color.white);
        plot.setRangeGridlinePaint(Color.white);
        plot.setForegroundAlpha(0.8f);
        chart.setBackgroundPaint(Color.white);
        try {
            ChartUtilities.saveChartAsPNG(new File(scan.MsLevel + "_ " + scan.ScanNum + ".PNG"), chart, 1000, 600);
        } catch (IOException e) {
        }
    }

    private void Read() throws FileNotFoundException, IOException, ParserConfigurationException, SAXException, DataFormatException {
        mzXMLReadUnit read = new mzXMLReadUnit(this.XMLtext);
        this.scan = read.Parse();
        this.XMLtext = null;
        read = null;
    }

    @Override
    public void run() {
        try {
            Read();
        } catch (Exception ex) {
            Logger.getRootLogger().error(ExceptionUtils.getStackTrace(ex));
        }
        
       scan.Preprocessing(parameter);
        
    }
}
