/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package CXL_PeakPairFinder;

import MSUmpire.BaseDataStructure.InstrumentParameter;
import MSUmpire.PeakDataStructure.PeakCluster;
import MSUmpire.PeptidePeakClusterDetection.PeakCurveCorrCalc;
import crosslinker.DC4;
import java.io.IOException;
import java.util.HashMap;
import net.sf.javaml.core.kdtree.KDTree;
import net.sf.javaml.core.kdtree.KeySizeException;
import org.apache.avalon.framework.ExceptionUtil;
import org.apache.log4j.Logger;

/**
 * This class finds precursor peak of crosslinked peptides (given two peak pairs)
 * @author Chih-Chiang Tsou
 */ 
public class PrecursorCrossPepFinder implements Runnable{
    public PairGroup LowMassPeakGroup;
    public PairGroup HighMassPeakGroup;
    public HashMap<Integer,CoElutePeak> PrecursorCrossPepPeaks;
    public float MaxLowMassPeakCorr=0f;
    public float MaxHighMassPeakCorr=0f;
    public float LowHighPeakCorr=0f;
    private float PrecCorssPeakPPM=1000f;
    private InstrumentParameter parameter;
    float lowrt = 0f;
    float highrt = 0f;
    private final KDTree PeakClusterSearchTree;
    private CoElutePeak BestPrecursorPeak=null;
    
    public PrecursorCrossPepFinder(PairGroup LowMassPeak, PairGroup HighMassPeak, KDTree PeakClusterSearchTree, InstrumentParameter parameter){
        this.LowMassPeakGroup=LowMassPeak;
        this.HighMassPeakGroup=HighMassPeak;
        this.PeakClusterSearchTree=PeakClusterSearchTree;
        this.parameter=parameter;
        try {
            LowHighPeakCorr = PeakCurveCorrCalc.CalPeakCorr(LowMassPeak.lowMassPeak.MonoIsotopePeak, HighMassPeak.lowMassPeak.MonoIsotopePeak, parameter.NoPeakPerMin);
            if (Float.isNaN(LowHighPeakCorr)) {
                LowHighPeakCorr = 0f;
            }
        } catch (IOException ex) {
           Logger.getRootLogger().error(ExceptionUtil.printStackTrace(ex));
        }
        
        lowrt = Math.max(LowMassPeak.lowMassPeak.PeakHeightRT[0],HighMassPeak.lowMassPeak.PeakHeightRT[0]) - parameter.RTtol;
        highrt = Math.min(LowMassPeak.lowMassPeak.PeakHeightRT[0],HighMassPeak.lowMassPeak.PeakHeightRT[0]) + parameter.RTtol;
    }
    public void FindPrecursorCrossPeak(){
        float IntactMW = LowMassPeakGroup.lowMassPeak.NeutralMass()+HighMassPeakGroup.lowMassPeak.NeutralMass()+DC4.DABCO;
        
        float lowMW = InstrumentParameter.GetMzByPPM(IntactMW, 1, parameter.MS1PPM);
        float highMW = InstrumentParameter.GetMzByPPM(IntactMW, 1, -parameter.MS1PPM);

        Object[] found = null;
        try {
            found = PeakClusterSearchTree.range(new double[]{lowrt, lowMW}, new double[]{highrt, highMW});
        } catch (KeySizeException ex) {
        }
        if (found == null || found.length == 0) {
            return;
        }

        PrecursorCrossPepPeaks = new HashMap<>();
        for (Object foundpeak : found) {
            PeakCluster peakB = (PeakCluster) foundpeak;
            float ppm = InstrumentParameter.CalcPPM(IntactMW, peakB.NeutralMass());
            if (ppm < parameter.MS1PPM) {
                float corrlow = 0.2f;
                float corrhigh = 0.2f;
                try {
                    corrlow = PeakCurveCorrCalc.CalPeakCorr(LowMassPeakGroup.lowMassPeak.MonoIsotopePeak, peakB.MonoIsotopePeak, parameter.NoPeakPerMin);
                } catch (IOException ex) {
                    Logger.getRootLogger().error(ex.getMessage());
                }
                if (Float.isNaN(corrlow)) {
                    corrlow = 0f;
                }
                try {
                    corrhigh = PeakCurveCorrCalc.CalPeakCorr(HighMassPeakGroup.lowMassPeak.MonoIsotopePeak, peakB.MonoIsotopePeak, parameter.NoPeakPerMin);
                } catch (IOException ex) {
                    Logger.getRootLogger().error(ex.getMessage());
                }
                if (Float.isNaN(corrhigh)) {
                    corrhigh = 0f;
                }

                //if (corrhigh>0.5f && corrlow>0.5f && corrhigh+corrlow > MaxLowMassPeakCorr+MaxHighMassPeakCorr) {
                if ((corrhigh>0.5f && corrlow>0.5f) && (PrecursorCrossPepPeaks.get(peakB.Charge)==null || ppm < PrecursorCrossPepPeaks.get(peakB.Charge).PPM)) {
                    CoElutePeak peak =new CoElutePeak(peakB, corrlow+corrhigh, ppm);
                    PrecursorCrossPepPeaks.put(peakB.Charge, peak);
                    MaxLowMassPeakCorr = corrlow;
                    MaxHighMassPeakCorr=corrhigh;
                    if (ppm < PrecCorssPeakPPM) {
                        PrecCorssPeakPPM = ppm;
                    }
                }
            }
        }
    }

     public CoElutePeak GetBestPrecursor() {
        if (BestPrecursorPeak == null && PrecursorCrossPepPeaks!=null) {
            for (CoElutePeak coElutePeak : PrecursorCrossPepPeaks.values()) {
                if (BestPrecursorPeak==null || (coElutePeak.Correlation > BestPrecursorPeak.Correlation || (coElutePeak.Correlation == BestPrecursorPeak.Correlation && coElutePeak.PPM < BestPrecursorPeak.PPM))) {
                    BestPrecursorPeak = coElutePeak;
                }
            }
        }
        return BestPrecursorPeak;
     }
    @Override
    public void run() {
        FindPrecursorCrossPeak();
    }    
}
