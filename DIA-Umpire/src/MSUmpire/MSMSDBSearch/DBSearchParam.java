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
    
    //<editor-fold defaultstate="collapsed" desc="[-inst InstrumentID] (0: Low-res LCQ/LTQ (Default), 1: High-res LTQ, 2: TOF, 3: Q-Exactive)">
    // [-m FragmentMethodID] (0: As written in the spectrum or CID if no info (Default), 1: CID, 2: ETD, 3: HCD)
    //[-inst MS2DetectorID] (0: Low-res LCQ/LTQ (Default), 1: Orbitrap/FTICR, 2: TOF, 3: Q-Exactive)
    //[-e EnzymeID] (0: unspecific cleavage, 1: Trypsin (Default), 2: Chymotrypsin, 3: Lys-C, 4: Lys-N, 5: glutamyl endopeptidase, 6: Arg-C, 7: Asp-N, 8: alphaLP, 9: no cleavage)
    public int MSGFInstrumentID=2;
    public int MSGFFragmentMethodID=0;
    public int MSGFEnzymeID=1;
    
//</editor-fold>
    
    public String parameterPath;
    public String templateParamFile;
    public SearchInstrumentType defaultType;
    public int NoCPUs = 2;
    public float PepFDR = 0.01f;
    public float ProtFDR = 0.01f;
    public boolean Overwrite = false;
    public String xinteractpath = "C:/inetpub/tpp-bin/xinteract";
    public String msconvertpath = "C:/inetpub/tpp-bin/msconvert";
    public String DecoyPrefix="rev_";
    public String xinteractpara = "-OpdEAP -PPM -drev -p0.1";
   
    
    public enum SearchInstrumentType {
        Orbitrap,
        TOF5600,
        YJ_TOF5600,
        QExactive,
        dfermin_phospho,
        Orbit_Velos,
        Orbit_elite,
        Orbit_Velos_High_Field,
        MSe,
    };

    @Override
    public Object clone() throws CloneNotSupportedException {
        return super.clone();
    }

    public abstract void SetResultFilePath(String mzXMLfile);
    public abstract void SetCombineFileName(String filename, String tag);
    public abstract void GenerateParamFile();
        
    protected void SetParameter(SearchInstrumentType type) {
        switch (type) {
            case Orbit_Velos: {
                PrecursorPPM = 20;
                FragPPM = 100;
                MinNoPeaksScoring = 3;
                MinNoPeaks = 15;
                TotalPeaks = 100;
                break;
            }
            case Orbit_elite: {
                PrecursorPPM = 20;
                FragPPM = 50;
                MinNoPeaksScoring = 3;
                MinNoPeaks = 15;
                TotalPeaks = 100;
                break;
            }
            case QExactive: {
                PrecursorPPM = 10;
                FragPPM = 20;
                MinNoPeaksScoring = 3;
                MinNoPeaks = 15;
                TotalPeaks = 140;
                MissCleavage = 1;
                break;
            }
            case MSe: {
                PrecursorPPM = 25;
                FragPPM = 50;
                MinNoPeaksScoring = 3;
                MinNoPeaks = 15;
                TotalPeaks = 140;
                MissCleavage = 1;
                break;
            }

            case dfermin_phospho: {
                PrecursorPPM = 200;
                FragPPM = 200;
                MinNoPeaksScoring = 3;
                MinNoPeaks = 15;
                TotalPeaks = 100;
                break;
            }
            case Orbit_Velos_High_Field: {
                PrecursorPPM = 20;
                FragPPM = 200;
                MinNoPeaksScoring = 3;
                MinNoPeaks = 15;
                TotalPeaks = 100;
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
            case YJ_TOF5600: {
                PrecursorPPM = 30;
                FragPPM = 70;
                MinNoPeaksScoring = 3;
                MinNoPeaks = 15;
                TotalPeaks = 140;
                MissCleavage = 1;
                SemiCleavage = false;
                break;

            }
        }
    }
}
