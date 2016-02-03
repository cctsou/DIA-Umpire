/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package CXL_PeakPairFinder;

import MSUmpire.BaseDataStructure.InstrumentParameter;
import MSUmpire.BaseDataStructure.XYData;
import MSUmpire.LCMSPeakStructure.LCMSPeakBase;
import MSUmpire.PeakDataStructure.PeakCluster;
import MSUmpire.SpectralProcessingModule.ScanPeakGroup;
import java.util.ArrayList;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;
import org.apache.log4j.Logger;

/**
 *
 * @author Chih-Chiang Tsou
 */
public class CrosslinkerPepFinder {

    public LCMSPeakBase LCMSPeakBase;
    public XYData mzRange = new XYData(Float.MIN_VALUE, Float.MAX_VALUE);
    public SortedFragFinder PairList;
    public ArrayList<PrecursorCrossPepFinder> IntactPepList;
    public ArrayList<MS2PeakPairFinder> IntactPepListMS2;
    public float medianIntensity=-1f;
    public CrosslinkerPepFinder(LCMSPeakBase LCMSPeakBase,float medianIntensity) {
        this.LCMSPeakBase = LCMSPeakBase;
        this.medianIntensity=medianIntensity;
    }

     public void FindPairPeakMS2(ArrayList<ScanPeakGroup> MS2PeakGroups, InstrumentParameter parameter, int NoCPUs) {
        Logger.getRootLogger().info("Finding MS2 peak pairs from MS/MS scans........");
        
        ExecutorService executorPool = null;
        executorPool = Executors.newFixedThreadPool(NoCPUs);
       
         IntactPepListMS2 = new ArrayList<>();
         //For each scan
         for (ScanPeakGroup scan : MS2PeakGroups) {
             MS2PeakPairFinder finder = new MS2PeakPairFinder(scan, parameter);
             executorPool.execute(finder);
             //finder.run();
             IntactPepListMS2.add(finder);
         }
        
        executorPool.shutdown();
        try {
            executorPool.awaitTermination(Long.MAX_VALUE, TimeUnit.NANOSECONDS);
        } catch (InterruptedException e) {
            Logger.getRootLogger().info("interrupted..");
        }

        //Logger.getRootLogger().info("No of ion clusters:" + LCMSPeakBase.PeakClusters.size() + " (Memory usage:" + Math.round((Runtime.getRuntime().totalMemory() - Runtime.getRuntime().freeMemory()) / 1048576) + "MB)");

    }


    public void FindAllPairPeaks(int NoCPUs) {
        Logger.getRootLogger().info("Finding MS1 peak pairs........");
                
        ExecutorService executorPool = null;
        executorPool = Executors.newFixedThreadPool(NoCPUs);
        PairList = new SortedFragFinder();
        ArrayList<PeakPairFinder> templist=new ArrayList<>();

        //For each MS1 peak cluster, find the peak pair
        for (int i = 0; i < LCMSPeakBase.PeakClusters.size(); i++) {
            PeakCluster peakCluster = LCMSPeakBase.PeakClusters.get(i);
            if (peakCluster.NeutralMass() >= mzRange.getX() && peakCluster.NeutralMass() <= mzRange.getY()) {
                PeakPairFinder unit = new PeakPairFinder(peakCluster, LCMSPeakBase.GetPeakClusterMassSearchTree(), LCMSPeakBase.parameter);
                templist.add(unit);
                executorPool.execute(unit);
            }
        }
        
        executorPool.shutdown();
        try {
            executorPool.awaitTermination(Long.MAX_VALUE, TimeUnit.NANOSECONDS);
        } catch (InterruptedException e) {
            Logger.getRootLogger().info("interrupted..");
        }

        for(PeakPairFinder finder : templist){
            if (finder.pairgroup.GetBestPeakPair() != null && (finder.pairgroup.lowMassPeak.PeakHeight[0]>medianIntensity || finder.pairgroup.GetBestPeakPair().peakpair.PeakHeight[0]>medianIntensity)) {
                PairList.add(finder);
            }
        }
        
        IntactPepList=new ArrayList<>();
        Logger.getRootLogger().info("Finding crosslinking groups");
        executorPool = Executors.newFixedThreadPool(NoCPUs);
        for (int i = 0; i < PairList.size(); i++) {
            if (PairList.get(i).pairgroup.GetBestPeakPair() != null ) {
                PairGroup lowmass = PairList.get(i).pairgroup;                
                for (int j = i; j < PairList.size(); j++) {                    
                    if (PairList.get(j).pairgroup.GetBestPeakPair() != null) {
                        PairGroup highmass = PairList.get(j).pairgroup;
                        //For any two peak pairs, check if they co-elute together (RT diff < RTtol)
                        if (lowmass.lowMassPeak.NeutralMass()<=highmass.lowMassPeak.NeutralMass() 
                                && Math.abs(lowmass.lowMassPeak.PeakHeightRT[0] - highmass.lowMassPeak.PeakHeightRT[0]) < LCMSPeakBase.parameter.RTtol) {
                            PrecursorCrossPepFinder intactLinkedPep = new PrecursorCrossPepFinder(lowmass, highmass, LCMSPeakBase.GetPeakClusterMassSearchTree(), LCMSPeakBase.parameter);
                            IntactPepList.add(intactLinkedPep);
                            executorPool.execute(intactLinkedPep);
                        }
                    }
                }
            }
        }

        executorPool.shutdown();
        try {
            executorPool.awaitTermination(Long.MAX_VALUE, TimeUnit.NANOSECONDS);
        } catch (InterruptedException e) {
            Logger.getRootLogger().info("interrupted..");
        }
        //Logger.getRootLogger().info("No of ion clusters:" + LCMSPeakBase.PeakClusters.size() + " (Memory usage:" + Math.round((Runtime.getRuntime().totalMemory() - Runtime.getRuntime().freeMemory()) / 1048576) + "MB)");

    }
}
