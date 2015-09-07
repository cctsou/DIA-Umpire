/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package MSUmpire.GUI;

import MSUmpire.PeakDataStructure.PeakCluster;
import MSUmpire.PeakDataStructure.PeakCurve;
import MSUmpire.PeakDataStructure.PrecursorFragmentPairEdge;
import java.awt.BasicStroke;
import java.awt.Color;
import org.jfree.chart.ChartFactory;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.axis.NumberAxis;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.chart.plot.XYPlot;
import org.jfree.chart.renderer.xy.XYLineAndShapeRenderer;
import org.jfree.data.xy.XYDataset;
import org.jfree.data.xy.XYSeriesCollection;

/**
 *
 * @author Chih-Chiang Tsou <chihchiang.tsou@gmail.com>
 */
public class ChartDataConverter {

    public static XYDataset createClusterDataset(PeakCluster cluster) {
        final XYSeriesCollection dataset = new XYSeriesCollection();
        if (cluster.IsoPeaksCurves[0] != null) {
            dataset.addSeries(cluster.IsoPeaksCurves[0].GetChartXYDatasetSmooth("1st peak(m/z:" + cluster.IsoPeaksCurves[0].TargetMz + ")"));
        }
        if (cluster.IsoPeaksCurves[1] != null) {
            dataset.addSeries(cluster.IsoPeaksCurves[1].GetChartXYDatasetSmooth("2nd peak(m/z:" + cluster.IsoPeaksCurves[1].TargetMz + ")"));
        }
        if (cluster.IsoPeaksCurves[2] != null) {
            dataset.addSeries(cluster.IsoPeaksCurves[2].GetChartXYDatasetSmooth("3rd peak(m/z:" + cluster.IsoPeaksCurves[2].TargetMz + ")"));
        }
        if (cluster.IsoPeaksCurves.length > 3 && cluster.IsoPeaksCurves[3] != null) {
            dataset.addSeries(cluster.IsoPeaksCurves[3].GetChartXYDatasetSmooth("4th peak(m/z:" + cluster.IsoPeaksCurves[3].TargetMz + ")"));
        }
        if (cluster.GroupedFragmentPeaks != null) {
            createSwathClusterDataset(dataset, cluster);
        }
        return dataset;
    }

    public static XYDataset createCurveDataset(PeakCurve peakCurve) {
        final XYSeriesCollection dataset = new XYSeriesCollection();
        if(peakCurve.GetSmoothedList().Data.isEmpty()){
            dataset.addSeries(peakCurve.GetChartXYDatasetRAW());
        }
        else {
            dataset.addSeries(peakCurve.GetChartXYDatasetSmooth(null));
        }
        return dataset;
    }

    public static void createSwathClusterDataset(XYSeriesCollection dataset, PeakCluster cluster) {
        for (PrecursorFragmentPairEdge fragment : cluster.GroupedFragmentPeaks) {
            //dataset.addSeries(fragment.peakCurve.GetChartXYDatasetSmooth(String.valueOf(" mz: " + fragment.fragmentinfo.FragmentMz + "(Index:" + fragment.fragmentinfo.PeakCurveIndexB + "; Corr:" + fragment.fragmentinfo.Correlation + ")")));
        }
    }

    public static JFreeChart createXICChart(PeakCurve peak) {
        XYDataset dataset = createCurveDataset(peak);
        JFreeChart chart = createXICChart(dataset, true);
        return chart;
    }
    
    public static JFreeChart createXICChart(PeakCluster cluster) {
        XYDataset dataset = createClusterDataset(cluster);
        JFreeChart chart = createXICChart(dataset, true);
        return chart;
    }

    public static JFreeChart createXICChart(final XYDataset dataset, boolean IncludeLengend) {

        // create the chart...
        final JFreeChart chart = ChartFactory.createXYLineChart(
                "", // chart title
                "RT", // x axis label
                "Intensity", // y axis label
                dataset, // data
                PlotOrientation.VERTICAL,
                IncludeLengend,//!(SWATHButton.isSelected()||SWATHButtonALL.isSelected()), // include legend
                true, // tooltips
                true // urls
        );

        // NOW DO SOME OPTIONAL CUSTOMISATION OF THE CHART...
        chart.setBackgroundPaint(Color.white);

//        final StandardLegend legend = (StandardLegend) chart.getLegend();
        //      legend.setDisplaySeriesShapes(true);
        // get a reference to the plot for further customisation...
        final XYPlot plot = chart.getXYPlot();
        plot.setBackgroundPaint(Color.white);
        //    plot.setAxisOffset(new Spacer(Spacer.ABSOLUTE, 5.0, 5.0, 5.0, 5.0));
        plot.setDomainGridlinePaint(Color.white);
        plot.setRangeGridlinePaint(Color.white);

        final XYLineAndShapeRenderer renderer = new XYLineAndShapeRenderer();

        for (int i = 0; i < plot.getSeriesCount(); i++) {
            renderer.setSeriesLinesVisible(i, true);
            renderer.setSeriesShapesVisible(i, false);
            //renderer.setSeriesItemLabelsVisible(i, false);
            renderer.setSeriesVisibleInLegend(i, false);
            renderer.setSeriesPaint(i, Color.GRAY);
        }

        renderer.setSeriesVisibleInLegend(0, true);
        renderer.setSeriesVisibleInLegend(1, true);
        renderer.setSeriesVisibleInLegend(2, true);
        renderer.setSeriesStroke(0, new BasicStroke(4.0f));
        renderer.setSeriesStroke(1, new BasicStroke(2.0f));
        renderer.setSeriesStroke(2, new BasicStroke(1.0f));
        renderer.setSeriesStroke(3, new BasicStroke(1.0f));
        renderer.setSeriesPaint(0, Color.BLACK);
        renderer.setSeriesPaint(1, Color.DARK_GRAY);
        renderer.setSeriesPaint(2, Color.GRAY);
        renderer.setSeriesPaint(3, Color.GRAY);
        plot.setRenderer(renderer);

        // change the auto tick unit selection to integer units only...
        final NumberAxis rangeAxis = (NumberAxis) plot.getRangeAxis();
        rangeAxis.setStandardTickUnits(NumberAxis.createIntegerTickUnits());
        // OPTIONAL CUSTOMISATION COMPLETED.
        return chart;

    }
}
