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
package MSUmpire.DIA;

import MSUmpire.PeakDataStructure.PeakCluster;
import java.io.Serializable;

/**
 *
 * @author Chih-Chiang Tsou <chihchiang.tsou@gmail.com>
 */
public class PeakGroupScore implements Serializable{
    private static final long serialVersionUID = 164894161684L;

    public transient PeakCluster cluster;
    public int MSlevel = 0;
    public float PrecursorRT;

    public float SpecDotProduct=0f;
    public float SpecCorrelation=0f;
    public float ContrastAngle=0f;
    public float AveCorrScore=0f;
    public float SumCorrScore=0f;
    public float SumCorrPPMScore=0f;
    public float PrecursorCentralRank=0f;
    public float MaxMatchCorr=0f;
    public float FragIntAvgScore=0f;
    public float PPMScore=0f;
    public float ApexDeltaScore=0f;
    public float RTOverlapScore=0f;
    public int NoMatchB = 0;
    public int NoMatchY = 0;
    public float PrecursorCorr=0f;
    public float RTDiff=0f;
    public float PrecursorPPM=0f;
    public int Peplength=0;
    public float PrecursorScore=0f;    
    public float PrecursorIsoPattern=0f;    
    public int NoFragmentLib = 0;
    public float MixtureModelProb;
    public float MixtureModelLocalProb;
    public float UmpireScore;
    public float PrecursorNeutralMass;
    
    public PeakGroupScore(PeakCluster cluster) {
        this.cluster = cluster;
        this.PrecursorRT=cluster.PeakHeightRT[0];
        this.PrecursorNeutralMass=cluster.NeutralMass();
    }

    void AddScore(float UmpireScore) {
        cluster.AddScore(UmpireScore);
    }

    float GetScoreRank(float UmpireScore) {
        return cluster.GetScoreRank(UmpireScore);
    }
    
}
