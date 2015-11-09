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

import java.io.*;
import java.util.ArrayList;
import org.apache.avalon.framework.ExceptionUtil;
import org.apache.commons.io.FilenameUtils;
import org.apache.log4j.Logger;

/**
 *
 * @author Chih-Chiang Tsou <chihchiang.tsou@gmail.com>
 */
public class TandemParam extends DBSearchParam {

    public String PotentialModification = "15.994915@M,57.021464@C";
    public String PotentialModMotif = "";
    public String FixModification = "";
    public ArrayList<String> LabelingModification = new ArrayList<>();
    public String tandempath = "C:/inetpub/tpp-bin/tandem.exe";
    public String tandem2XML = "C:/inetpub/tpp-bin/Tandem2XML";
    public String Scoring = "Native";
    public boolean SpectrumConditioning = true;

    public TandemParam(SearchInstrumentType type) {
        defaultType = type;
        SetParameter(type);
    }

    private String AddPotentialModification(String output, String Modinfo) {
        if (!"".equals(Modinfo)) {
            output = output.replace("<note type=\"input\" label=\"residue, potential modification mass\"></note>", "<note type=\"input\" label=\"residue, potential modification mass\">" + Modinfo + "</note>");
        }
        return output;
    }

    private String AddFixModification(String output, String Modinfo) {
        if (!"".equals(Modinfo)) {
            output = output.replace("<note type=\"input\" label=\"residue, modification mass\"></note>", "<note type=\"input\" label=\"residue, modification mass\">" + Modinfo + "</note>");
        }
        return output;
    }

    private String AddPotentialModMotif(String output, String ModMotif) {
        if (!"".equals(ModMotif)) {
            output = output.replace("<note type=\"input\" label=\"residue, potential modification motif\"></note>", "<note type=\"input\" label=\"residue, potential modification motif\">" + ModMotif + "</note>");
        }
        return output;
    }

    private String AddlableingModification(String output, int Index, String Modinfo) {
        output = output.replace("<note> Add labeling modification by replacing this line</note>", "<note type=\"input\" label=\"residue, modification mass " + Index + "\">" + Modinfo + "</note>\n<note> Add labeling modification by replacing this line</note>");
        return output;
    }

    public void SetCombineFileName(String filename, String tag) {
        CombinedPepXML = FilenameUtils.separatorsToUnix(FilenameUtils.getFullPath(filename) + "interact-" + FilenameUtils.getBaseName(filename) + tag + ".tandem.combine.pep.xml");
        CombinedProt = FilenameUtils.getFullPath(filename) + FilenameUtils.getBaseName(filename) + tag + ".tandem.Qcombine.prot.xml";
    }

    @Override
    public void SetResultFilePath(String mzXMLfile) {
        SpectrumPath = FilenameUtils.separatorsToUnix(FilenameUtils.getFullPath(mzXMLfile) + FilenameUtils.getName(mzXMLfile));
        PepXMLPath = FilenameUtils.separatorsToUnix(FilenameUtils.getFullPath(mzXMLfile) + FilenameUtils.getBaseName(mzXMLfile) + ".tandem.pep.xml");
        InteractPepXMLPath = FilenameUtils.separatorsToUnix(FilenameUtils.getFullPath(mzXMLfile) + "interact-" + FilenameUtils.getBaseName(mzXMLfile) + ".tandem.pep.xml");
        ProtXMLPath = InteractPepXMLPath.replace(".pep.xml", ".prot.xml");
        parameterPath = FilenameUtils.separatorsToUnix(FilenameUtils.getFullPath(mzXMLfile) + FilenameUtils.getBaseName(mzXMLfile) + ".tandem.param");
        RawSearchResult = FilenameUtils.separatorsToUnix(FilenameUtils.getFullPath(mzXMLfile) + FilenameUtils.getBaseName(mzXMLfile) + ".tandem");
    }

    @Override
    public void GenerateParamFile() {
        FileWriter writer = null;
        try {
            InputStream is = DBSearchParam.class.getClassLoader().getResourceAsStream("resource/tandem.xml");
            if (templateParamFile != null && new File(templateParamFile).exists()) {
                Logger.getRootLogger().info("Using X! Tandem parameter template: " + templateParamFile);
                is = new FileInputStream(templateParamFile);
            }
            BufferedReader reader = new BufferedReader(new InputStreamReader(is));
            String line = "";
            StringBuffer sb = new StringBuffer();
            while ((line = reader.readLine()) != null) {
                sb.append(line + "\n");
            }
            String output = sb.toString();
            if (Scoring.equals("Kscore")) {
                output = output.replace("<note type=\"input\" label=\"list path, default parameters\">isb_default_input_native.xml</note>", "<note type=\"input\" label=\"list path, default parameters\">" + FilenameUtils.getFullPath(parameterPath) + "/isb_default_input_kscore.xml</note>");
            } else if (Scoring.equals("Native")) {
                output = output.replace("<note type=\"input\" label=\"list path, default parameters\">isb_default_input_native.xml</note>", "<note type=\"input\" label=\"list path, default parameters\">" + FilenameUtils.getFullPath(parameterPath) + "/isb_default_input_native.xml</note>");
            }
            output = output.replace("<note type=\"input\" label=\"spectrum, fragment monoisotopic mass error\">15</note>", "<note type=\"input\" label=\"spectrum, fragment monoisotopic mass error\">" + FragPPM + "</note>");
            output = output.replace("<note type=\"input\" label=\"spectrum, parent monoisotopic mass error minus\">20</note>", "<note type=\"input\" label=\"spectrum, parent monoisotopic mass error minus\">" + PrecursorPPM + "</note>");
            output = output.replace("<note type=\"input\" label=\"spectrum, parent monoisotopic mass error plus\">20</note>", "<note type=\"input\" label=\"spectrum, parent monoisotopic mass error plus\">" + PrecursorPPM + "</note>");
            output = output.replace("<note type=\"input\" label=\"spectrum, path\"></note>", "<note type=\"input\" label=\"spectrum, path\">" + SpectrumPath + "</note>");
            if (NonSpecificCleavage) {
                output = output.replace("<note type=\"input\" label=\"protein, cleavage site\">[RK]|{P}</note>", "<note type=\"input\" label=\"protein, cleavage site\">[X]|[X]</note>");
            }
            if (!SpectrumConditioning) {
                output = output.replace("<note type=\"input\" label=\"spectrum, use conditioning\">yes</note>", "<note type=\"input\" label=\"spectrum, use conditioning\">no</note>");
            }
            output = output.replace("<note type=\"input\" label=\"output, path\"></note>", "<note type=\"input\" label=\"output, path\">" + RawSearchResult + "</note>");
            //output=output.replace("<note type=\"input\" label=\"output, sequence path\"></note>", "<note type=\"input\" label=\"output, sequence path\">" + OutputSeqPath + "</note>");
            output = output.replace("<note type=\"input\" label=\"scoring, minimum ion count\">3</note>", "<note type=\"input\" label=\"scoring, minimum ion count\">" + MinNoPeaksScoring + "</note>");
            output = output.replace("<note type=\"input\" label=\"scoring, maximum missed cleavage sites\">0</note>", "<note type=\"input\" label=\"scoring, maximum missed cleavage sites\">" + MissCleavage + "</note>");
            output = output.replace("<note type=\"input\" label=\"list path, taxonomy information\">taxnomy.xml</note>", "<note type=\"input\" label=\"list path, taxonomy information\">" + FilenameUtils.getFullPath(parameterPath) + "taxonomy.xml" + "</note>");
            output = output.replace("<note type=\"input\" label=\"spectrum, minimum peaks\">3</note>", "<note type=\"input\" label=\"spectrum, minimum peaks\">" + MinNoPeaks + "</note>");
            output = output.replace("<note type=\"input\" label=\"spectrum, total peaks\">100</note>", "<note type=\"input\" label=\"spectrum, total peaks\">" + TotalPeaks + "</note>");
            for (int i = 0; i < LabelingModification.size(); i++) {
                output = AddlableingModification(output, i + 1, LabelingModification.get(i));
            }
            output = AddPotentialModification(output, PotentialModification);
            output = AddFixModification(output, FixModification);
            output = AddPotentialModMotif(output, PotentialModMotif);
            if (SemiCleavage) {
                output = output.replace("<note type=\"input\" label=\"protein, cleavage semi\">no</note>", "<note type=\"input\" label=\"protein, cleavage semi\">yes</note>");
            }
            if (IsotopeError) {
                output = output.replace("<note type=\"input\" label=\"spectrum, parent monoisotopic mass isotope error\">no</note>", "<note type=\"input\" label=\"spectrum, parent monoisotopic mass isotope error\">yes</note>");
            }
            output = output.replace("<note type=\"input\" label=\"spectrum, threads\">6</note>", "<note type=\"input\" label=\"spectrum, threads\">" + NoCPUs + "</note>");
            writer = new FileWriter(parameterPath);
            writer.write(output);
            writer.close();
            is = DBSearchParam.class.getClassLoader().getResourceAsStream("resource/taxonomy.xml");
            reader = new BufferedReader(new InputStreamReader(is));
            line = "";
            sb = new StringBuffer();
            while ((line = reader.readLine()) != null) {
                sb.append(line);
            }
            output = sb.toString();
            output = output.replace("<file format=\"peptide\" URL=\"fasta\"/>", "<file format=\"peptide\" URL=\"" + FastaPath + "\"/>");
            output = output.replace("<file format=\"peptide\" URL=\"fastadecoy\"/>", "<file format=\"peptide\" URL=\"" + DecoyFasta + "\"/>");
            writer = new FileWriter(FilenameUtils.getFullPath(parameterPath) + "taxonomy.xml");
            writer.write(output);
            writer.close();
            is = DBSearchParam.class.getClassLoader().getResourceAsStream("resource/isb_default_input_kscore.xml");
            reader = new BufferedReader(new InputStreamReader(is));
            line = "";
            sb = new StringBuffer();
            while ((line = reader.readLine()) != null) {
                sb.append(line);
            }
            output = sb.toString();
            if (!(new File(FilenameUtils.getFullPath(parameterPath) + "isb_default_input_kscore.xml")).exists()) {
                writer = new FileWriter(FilenameUtils.getFullPath(parameterPath) + "isb_default_input_kscore.xml");
                writer.write(output);
                writer.close();
            }
            is = DBSearchParam.class.getClassLoader().getResourceAsStream("resource/isb_default_input_native.xml");
            reader = new BufferedReader(new InputStreamReader(is));
            line = "";
            sb = new StringBuffer();
            while ((line = reader.readLine()) != null) {
                sb.append(line);
            }
            output = sb.toString();
            if (!(new File(FilenameUtils.getFullPath(parameterPath) + "isb_default_input_native.xml")).exists()) {
                writer = new FileWriter(FilenameUtils.getFullPath(parameterPath) + "isb_default_input_native.xml");
                writer.write(output);
                writer.close();
            }
            writer.close();
        } catch (Exception ex) {
            Logger.getRootLogger().error(ExceptionUtil.printStackTrace(ex));
        }
    }
}
