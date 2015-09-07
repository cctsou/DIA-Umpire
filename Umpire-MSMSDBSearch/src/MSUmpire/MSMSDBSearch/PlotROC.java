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
package MSUmpire.MSMSDBSearch;

import MSUmpire.PSMDataStructure.LCMSID;
import MSUmpire.PSMDataStructure.PSM;
import MSUmpire.PSMDataStructure.SortedPSMListProb;
import MSUmpire.SearchResultParser.PepXMLParser;
import java.awt.Color;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;

import javax.xml.parsers.ParserConfigurationException;
import org.apache.commons.io.FilenameUtils;
import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartUtilities;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.chart.plot.XYPlot;
import org.jfree.data.xy.XYDataItem;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;
import org.xml.sax.SAXException;
import org.xmlpull.v1.XmlPullParserException;

/**
 *
 * @author Chih-Chiang Tsou <chihchiang.tsou@gmail.com>
 */
public class PlotROC {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) throws IOException, InterruptedException, SAXException, XmlPullParserException, ParserConfigurationException {
        if (args == null || args.length < 2) {
            args = new String[7];
            args[0] = "IDA_UPS2_5pm_Ecoli_1ug_rep3.mzXML";
            args[1] = "10";
            args[2] = "101";
            args[3] = "10";
            args[4] = "101";
            args[5] = "10";
            args[6] = "IDA_UPS2_5pm_Ecoli_1ug_rep3_ROC.csv";
        }
        if (args.length < 7) {
            System.out.print("java -jar -d64 -Xms15G PlotROC.jar mzXML StartPrecursorPPM EndPrecursorPPM StartFragPPM EndFragPPM PPMInterval ExportCSV\n");
            return;
        }

        String Filepath = new java.io.File(".").getCanonicalPath() + "/";
        //Filepath = "C:/Umich/Java_projects/JQuant/UPS_Ecoli_SWATH/";
        Filepath = "C:/Umich/Data/SWATH_IDA_Project94/New_Rep/FindParameter/IDA_UPS2_5pm_Ecoli_1ug_rep3/";
        String mzxML = args[0];

        int StartPrecursorPPM = Integer.parseInt(args[1]);
        int EndPrecursorPPM = Integer.parseInt(args[2]);
        int StartFragPPM = Integer.parseInt(args[3]);
        int EndFragPPM = Integer.parseInt(args[4]);
        int PPMInterval = Integer.parseInt(args[5]);
        String ExportCSV = args[6];

        int maxpos = 0;
        FileWriter writer = new FileWriter(Filepath + ExportCSV);
        float bestprecursorPPM = 0f;
        float bestfragPPM = 0f;
        final XYSeriesCollection dataset = new XYSeriesCollection();
        LCMSID pepID = new LCMSID(mzxML,"rev","");
        for (int precursorppm = StartPrecursorPPM; precursorppm < EndPrecursorPPM; precursorppm += PPMInterval) {
            for (int fragppm = StartFragPPM; fragppm < EndFragPPM; fragppm += PPMInterval) {
                PepXMLParser pepXMLParser = new PepXMLParser(pepID, Filepath + "interact-" + FilenameUtils.getBaseName(mzxML) + "_" + precursorppm + "_" + fragppm + ".pep.xml", 0f);

                SortedPSMListProb sortedPSMs = new SortedPSMListProb();

                for (PSM psm : pepID.PSMList.values()) {
                    sortedPSMs.add(psm);
                }
                int pos = 0;
                int total = 0;
                int decoy = 0;
                final XYSeries series = new XYSeries(precursorppm + "_" + fragppm);

                for (int i = 0; i < sortedPSMs.size(); i++) {
                    PSM topscore = sortedPSMs.get(i);
                    total++;
                    if (topscore.IsDecoy("rev")) {
                        decoy++;
                        int fdr = decoy;
                        int sensitivity = (total - decoy);
                        series.add(fdr, sensitivity);
                        if (fdr > 700) {
                            break;
                        }
                    }
                }

                writer.write(precursorppm + "_" + fragppm + "_X\t");
                for (Object xy : series.getItems()) {
                    writer.write(((XYDataItem) xy).getXValue() + "\t");
                }
                writer.write("\n");
                writer.write(precursorppm + "_" + fragppm + "_Y\t");
                for (Object xy : series.getItems()) {
                    writer.write(((XYDataItem) xy).getYValue() + "\t");
                }
                writer.write("\n");
                dataset.addSeries(series);
            }
        }
        writer.close();
        JFreeChart chart = ChartFactory.createXYLineChart("", "decoy", "positive", dataset, PlotOrientation.VERTICAL, true, true, false);
        XYPlot plot = chart.getXYPlot();
        //plot.getRangeAxis().setRange(0, 100);
        plot.setBackgroundPaint(Color.white);
//        //plot.setForegroundAlpha(0.8f);
        plot.setDomainGridlinePaint(Color.white);
        plot.setRangeGridlinePaint(Color.white);
        chart.setBackgroundPaint(Color.white);
//        final XYLineAndShapeRenderer renderer = new XYLineAndShapeRenderer();
//        renderer.setSeriesLinesVisible(0, false);
//        renderer.setSeriesShapesVisible(1, false);
        //plot.setRenderer(renderer);
        ChartUtilities.saveChartAsPNG(new File(Filepath + FilenameUtils.getBaseName(mzxML) + "_ROC.PNG"), chart, 1200, 800);
        System.out.print(args[0] + " Best precursor PPM:" + bestprecursorPPM + " Best fragment PPM:" + bestfragPPM + " Positive matches:" + maxpos + "\n");
    }
}
