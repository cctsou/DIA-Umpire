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
package MSUmpire.PostProcessing;

import MSUmpire.BaseDataStructure.InstrumentParameter;
import MSUmpire.DIA.DIAPack;
import MSUmpire.PeakDataStructure.PeakCluster;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;
import org.apache.log4j.Logger;

/**
 *
 * @author Chih-Chiang Tsou <chihchiang.tsou@gmail.com>
 */
public class PrecursorGroupAlignment {
    ArrayList<DIAPack> DIAfiles;
    InstrumentParameter parameter;
    
    public PrecursorGroupAlignment(ArrayList<DIAPack> DIAFiles, InstrumentParameter parameter){
        this.DIAfiles=DIAFiles;
        this.parameter=parameter;
    }
    
    public void CreateScoreMatrix(String path) throws IOException{
        System.out.println("Calculate matrix");
        for (int i = 0; i < DIAfiles.size(); i++) {
            for (PeakCluster cluster : DIAfiles.get(i).ms1lcms.PeakClusters) {
                cluster.CreateLock();
            }
        }
        
        for (int i = 0; i < DIAfiles.size(); i++) {
            for (int j = i + 1; j < DIAfiles.size(); j++) {
                FileWriter writer = new FileWriter(path + i + "_" + j + "_score.txt");
                ArrayList<String> IncludedIndex=new ArrayList<>();                
                ExecutorService executorPool = null;
                executorPool = Executors.newFixedThreadPool(20);
                ArrayList<PeakClusterMatch> ResultList = new ArrayList<>();                
//                for (PeakCluster cluster : DIAfiles.get(i).ms1lcms.MZSortedClusters) {
//                    int startindex = DIAfiles.get(j).ms1lcms.MZSortedClusters.BinarySearchLower(InstrumentParameter.GetMzByPPM(cluster.TargetMz(), cluster.Charge, parameter.MS1PPM));
//                    int endindex = DIAfiles.get(j).ms1lcms.MZSortedClusters.BinarySearchLower(InstrumentParameter.GetMzByPPM(cluster.TargetMz(), cluster.Charge, -parameter.MS1PPM));
//                    for(int idx=startindex;idx<=endindex;idx++){
//                        PeakCluster cluster2=DIAfiles.get(j).ms1lcms.MZSortedClusters.get(idx);
//                        if (Math.abs(cluster2.PeakHeightRT[0] - cluster.PeakHeightRT[0]) < 10f) {
//                            if (!IncludedIndex.contains(cluster.Index + "_" + cluster2.Index)) {
//                                PeakClusterMatch match=new PeakClusterMatch(cluster, cluster2, parameter);
//                                ResultList.add(match);
//                                executorPool.execute(match);                                
//                                IncludedIndex.add(cluster.Index + "_" + cluster2.Index);
//                            }
//                        }
//                    }
//                }
//                for (PeakCluster cluster : DIAfiles.get(j).ms1lcms.MZSortedClusters) {
//                    int startindex = DIAfiles.get(i).ms1lcms.MZSortedClusters.BinarySearchLower(InstrumentParameter.GetMzByPPM(cluster.TargetMz(), cluster.Charge, parameter.MS1PPM));
//                    int endindex = DIAfiles.get(i).ms1lcms.MZSortedClusters.BinarySearchLower(InstrumentParameter.GetMzByPPM(cluster.TargetMz(), cluster.Charge, -parameter.MS1PPM));
//                    for (int idx = startindex; idx <= endindex; idx++) {
//                        PeakCluster cluster2 = DIAfiles.get(i).ms1lcms.MZSortedClusters.get(idx);
//                        if (Math.abs(cluster2.PeakHeightRT[0] - cluster.PeakHeightRT[0]) < 10f) {
//                            if (!IncludedIndex.contains(cluster2.Index + "_" + cluster.Index)) {
//                                PeakClusterMatch match=new PeakClusterMatch(cluster2, cluster, parameter);
//                                ResultList.add(match);
//                                executorPool.execute(match);    
//                                IncludedIndex.add(cluster2.Index + "_" + cluster.Index);
//                            }
//                        }
//                    }
//                }
                executorPool.shutdown();
//                while (!executorPool.isTerminated()) {
//                }
                try {
                    executorPool.awaitTermination(Long.MAX_VALUE, TimeUnit.NANOSECONDS);
                } catch (InterruptedException e) {
                    Logger.getRootLogger().info("interrupted..");
                }
                
                for(PeakClusterMatch match : ResultList){
                    writer.write(match.peakClusterA.PeakHeightRT[0] + "\t" + match.peakClusterB.PeakHeightRT[0]+ "\t" +match.similiarity+"\n");
                }
                writer.close();
            }
        }
    }
}
