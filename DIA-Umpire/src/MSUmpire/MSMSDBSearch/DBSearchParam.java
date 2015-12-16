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

/**
 *
 * @author Chih-Chiang Tsou <chihchiang.tsou@gmail.com>
 */
public abstract class DBSearchParam implements Cloneable{

    public float FragPPM;
    public float PrecursorPPM;
    public int MinNoPeaksScoring;
    public int MinNoPeaks;
    public int TotalPeaks;
    public String SpectrumPath;
    public String RawSearchResult;
    public String InteractPepXMLPath;
    public String PepXMLPath;
    public String ProtXMLPath;
    public String CombinedPepXML;
    public String CombinedProt;
    public String FastaPath;
    public String DecoyFasta;
    public String OutputSeqPath;
    public int MissCleavage = 1;
    public boolean SemiCleavage = false;
    public boolean NonSpecificCleavage=false;
    public boolean IsotopeError = false;    
    public String parameterPath;
    public String templateParamFile;
    public SearchInstrumentType defaultType;
    public int NoCPUs = 2;
    public float PepFDR = 0.01f;
    public float ProtFDR = 0.01f;
    public boolean Overwrite = false;
    public String DecoyPrefix="rev_";
    
    public enum SearchInstrumentType {
        Orbitrap,
        TOF5600,
        QExactive,
        };

    @Override
    public Object clone() throws CloneNotSupportedException {
        return super.clone();
    }

    public abstract void SetResultFilePath(String mzXMLfile);
    public abstract void SetCombineFileName(String filename, String tag);
            
    protected void SetParameter(SearchInstrumentType type) {
        switch (type) {
            case QExactive: {
                PrecursorPPM = 10;
                FragPPM = 20;
                MinNoPeaksScoring = 3;
                MinNoPeaks = 15;
                TotalPeaks = 140;
                MissCleavage = 1;
                break;
            }
            case Orbitrap: {
                PrecursorPPM = 10;
                FragPPM = 500;
                MinNoPeaksScoring = 3;
                MinNoPeaks = 15;
                TotalPeaks = 100;
                break;
            }

            case TOF5600: {
                PrecursorPPM = 30;
                FragPPM = 40;
                MinNoPeaksScoring = 3;
                MinNoPeaks = 3;
                TotalPeaks = 140;
                MissCleavage = 1;
                SemiCleavage = false;
                break;
            }
        }
    }
}
