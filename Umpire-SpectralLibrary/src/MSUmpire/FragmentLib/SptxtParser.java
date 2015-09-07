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
package MSUmpire.FragmentLib;

import MSUmpire.PSMDataStructure.ModStringConvert;
import MSUmpire.PSMDataStructure.PepFragmentLib;
import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.HashMap;
import org.xmlpull.v1.XmlPullParserException;

/**
 *
 * @author Chih-Chiang Tsou <chihchiang.tsou@gmail.com>
 */


//<editor-fold defaultstate="collapsed" desc="Example 1">
/*Name: n[43]AAAAAAGAGPEM[147]VR/2
LibID: 2
MW: 1301.6387
PrecursorMZ: 650.8193
Status: Normal
FullName: -.n[43]AAAAAAGAGPEM[147]VR.G/2
Comment: AvePrecursorMz=651.2202 BinaryFileOffset=10348 DUScorr=3.3/1.9/2.9 Dot_cons=0.887/0.060 Dotbest=0.89 Dotfull=0.814/0.049 Dottheory=0.89 Flags=1,0,0 FracUnassigned=0.12,1/5;0.20,7/20;0.07,37/152 Fullname=-.AAAAAAGAGPEM(O)VR.G/2 Inst=qtof Max2med_orig=100.0/0.0 Missing=0.1425/0.0856 Mods=2/0,A,Acetyl/11,M,Oxidation Mz_av=651.234 Mz_diff=-0.017 Mz_exact=650.8199 NAA=14 NISTProtein="sp|P28482|MK01_HUMAN Mitogen-activated protein kinase 1 [Homo sapiens]" NMC=0 NTT=2 Naa=14 Nreps=13/23 Organism=human Parent=650.820 Parent_med=650.8110/0.01 Pep=Tryptic Pfin=7.1e+011 Pfract=0 Prob=1.0000 Probcorr=1 Protein=1/sp|P28482|MK01_HUMAN Pseq=36/2 Sample=6/ca_hela2_bentscheff_cam,1,1/ca_k562_2_bantscheff_cam,2,3/gu_platelets_cofradic_none,0,3/mpi_a459_cam,4,6/mpi_hct116_cam,1,1/mpi_jurkat_cam,5,9 Se=3^X12:ex=1.25e-005/7.977e-005,td=0/4.033e+005,sd=0/0,hs=47.35/3.367,bs=2.9e-007,b2=6.8e-007,bd=1.9e+006^O13:ex=7.42e-011/1.894e-006,td=8150/9.962e+004,pr=5.17e-014/1.161e-009,bs=1.11e-012,b2=2.1e-012,bd=49500^P11:sc=33.5/2.091,dc=22.1/2.327,ps=3.3/0.4809,bs=26.5,b2=26.5,bd=16.8 Spec=Consensus Tfratio=5.4e+005 Unassign_all=0.11 Unassigned=0.023
NumPeaks: 152
86.1000	1215.0	a1/0.04	7/13 15.1
115.1000	148.0	a3^2/0.53,y2-44^2/-0.50,y2-46^2/0.51	8/12 0.5
143.1000	221.0	?	11/13 0.4*/
//</editor-fold>

//<editor-fold defaultstate="collapsed" desc="Example 2 (with retention time)">
/*Name: FVFSLVDAM[147]NGK/2
LibID: 2157
MW: 1344.6737
PrecursorMZ: 672.3368
Status: Normal
FullName: R.FVFSLVDAM[147]NGK.E/2 (CID)
Comment: AvePrecursorMz=672.7876 BestRawSpectrum=18299_REP2_500ng_HumanLysate_IDA_1.23998.23998 BinaryFileOffset=9236618 ConsFracAssignedPeaks=0.562 DotConsensus=0.77,0.13;1/2 FracUnassigned=0.40,2/5;0.56,12/20;0.36,35/84 Inst="1/mass analyzer type,2,2" MassDiff=0.0020 MassDiffCounts=1/0:2 MaxRepSN=2.9 Mods=1/8,M,Oxidation NAA=12 NMC=0 NTT=2 Nreps=2/2 OrigMaxIntensity=0.48 Parent=672.337 Pep=Tryptic PrecursorIntensity=0 Prob=0.9999 ProbRange=0.999999,0.9999,0.99265,0.9854 Protein=3/trG3XAL0|G3XAL0_HUMAN/trE9PDB2|E9PDB2_HUMAN/sp_P40926|MDHM_HUMAN RepFracAssignedPeaks=0.594 RepNumPeaks=63.5/20.5 RetentionTime=3696.3,3967.6,3425.1 SN=200.0 Sample=1/f__Data_SWATH_Umpire_Paper_OpenSWATH_Human_speclib_interact-18299_REP2_500ng_HumanLysate_IDA_1_FIXEDPROB_,2,2 Se=1^X2:ex=0.0000e+000/0.0000e+000,fv=4.2188/1.3055,hs=0.0000e+000/0.0000e+000,ns=0.0000e+000/0.0000e+000,pb=0.9927/0.0073 Spec=Consensus TotalIonCurrent=1.6e+004
NumPeaks: 112
110.0782	3054.2	a2^2/0.00	1/1 0.0000|0.00
120.0813	8219.9	IFA/0.00,IM_147_A/0.03,a1/0.00	2/1 0.0000|0.30
129.0000	5995.0	IKD/0.00,y1-18/-0.10	1/1 0.0000|0.00*/
//</editor-fold>


public class SptxtParser {

    public void Parse(String sptxt, HashMap<String, PepFragmentLib> PeptideFragmentLib) throws FileNotFoundException, IOException, XmlPullParserException {
        
        BufferedReader reader=new BufferedReader(new FileReader(sptxt));
        String line="";        
        while ((line = reader.readLine()) != null) {
            if (line.startsWith("Name: ")) {
                PepFragmentLib fraglib = new PepFragmentLib();
                String modseq = line.split(":")[1].split("/")[0].trim();
                ModStringConvert.ConvertTPPModString(modseq, null);
                fraglib.Charge=Integer.parseInt(line.split(":")[1].split("/")[1].trim());

                do {
                    line = reader.readLine();
                    if (line.startsWith("PrecursorMZ: ")) {
                        fraglib.PrecursorMz = Float.parseFloat(line.split(":")[1].trim());                        
                    }
                    else if (line.startsWith("PrecursorMZ: ")) {
                        fraglib.PrecursorMz = Float.parseFloat(line.split(":")[1].trim());                        
                    }
                } while ("".equals(line));

                if (!PeptideFragmentLib.containsKey(fraglib.GetKey())) {
                }
            }
              
        
//        if (!PeptideFragmentLib.containsKey(pepIonID.GetKey())) {
//            PepFragmentLib fraglib = new PepFragmentLib();
//            fraglib.Sequence = pepIonID.Sequence;
//            fraglib.ModificationString = pepIonID.GetModificationString();
//            fraglib.Charge = pepIonID.Charge;
//            fraglib.ModSequence = pepIonID.ModSequence;
//            fraglib.PrecursorMz = pepIonID.NeutralPrecursorMz();
//            fraglib.MS1Score = pepIonID.PeakClusterScore;
//            if (pepIonID.MaxProbability > fraglib.MaxProbability) {
//                fraglib.MaxProbability = pepIonID.MaxProbability;
//            }
//            if (pepIonID.PeakClusterScore > fraglib.MS1Score) {
//                fraglib.MS1Score = pepIonID.PeakClusterScore;
//            }
//            PeptideFragmentLib.put(pepIonID.GetKey(), fraglib);
//        }
//        if (pepIonID.FragmentPeaks != null && !pepIonID.FragmentPeaks.isEmpty()) {
//            PeptideFragmentLib.get(pepIonID.GetKey()).AddFragments(pepIonID.FragmentPeaks);
        }
    }    
}
