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
import org.apache.avalon.framework.ExceptionUtil;
import org.apache.commons.io.FilenameUtils;
import org.apache.log4j.Logger;

/**
 *
 * @author Chih-Chiang Tsou <chihchiang.tsou@gmail.com>
 */
public class CometParam extends DBSearchParam {

    public String cometpath = "";    

    public CometParam(SearchInstrumentType type) {
        defaultType = type;
        SetParameter(type);
    }

    @Override
    public void GenerateParamFile() {
        try {
            InputStream is = DBSearchParam.class.getClassLoader().getResourceAsStream("resource/comet.params");
            
            if(templateParamFile!=null && new File(templateParamFile).exists()){
                Logger.getRootLogger().info("Using Comet parameter template: "+templateParamFile);
                is = new FileInputStream(templateParamFile);
            }  
            
            BufferedReader reader = new BufferedReader(new InputStreamReader(is));
            String line = "";
            StringBuffer sb = new StringBuffer();
            while ((line = reader.readLine()) != null) {
                sb.append(line + "\n");
            }
            String output = sb.toString();
            output = output.replace("database_name = /some/path/db.fasta", "database_name =" + FastaPath);
            output = output.replace("peptide_mass_tolerance = 30.0", "peptide_mass_tolerance =" + PrecursorPPM);
            output = output.replace("allowed_missed_cleavage = 1", "allowed_missed_cleavage = "+MissCleavage);
            output = output.replace("num_threads = 0", "num_threads = "+NoCPUs);
            if (IsotopeError) {
                output = output.replace("isotope_error = 0", "isotope_error = 1");
            }
            if (NonSpecificCleavage) {
                output = output.replace("search_enzyme_number = 1", "search_enzyme_number = 0");
            }
            if (FragPPM > 200) {
                output = output.replace("fragment_bin_tol = 0.02", "fragment_bin_tol = 1.0005");
                output = output.replace("theoretical_fragment_ions = 0", "theoretical_fragment_ions = 1");
                output = output.replace("fragment_bin_offset = 0.0", "fragment_bin_offset = 0.4");
            }
            
            FileWriter writer = new FileWriter(parameterPath);
            writer.write(output);
            writer.close();
        } catch (IOException ex) {
            Logger.getRootLogger().error(ExceptionUtil.printStackTrace(ex));
        }
    }
    
    public void SetCombineFileName(String filename, String tag){
        CombinedPepXML=FilenameUtils.separatorsToUnix(FilenameUtils.getFullPath(filename) + "interact-" + FilenameUtils.getBaseName(filename) + tag+".comet.combine.pep.xml");        
        CombinedProt = FilenameUtils.getFullPath(filename) + FilenameUtils.getBaseName(filename) + tag+ ".comet.Qcombine.prot.xml";
    }
    
    @Override
    public void SetResultFilePath(String mzXMLfile) {
        SpectrumPath = FilenameUtils.separatorsToUnix(FilenameUtils.getFullPath(mzXMLfile) + FilenameUtils.getName(mzXMLfile));
        PepXMLPath = FilenameUtils.separatorsToUnix(FilenameUtils.getFullPath(mzXMLfile)  + FilenameUtils.getBaseName(mzXMLfile) + ".comet.pep.xml");
        InteractPepXMLPath = FilenameUtils.separatorsToUnix(FilenameUtils.getFullPath(mzXMLfile) + "interact-" + FilenameUtils.getBaseName(mzXMLfile) + ".comet.pep.xml");
        ProtXMLPath = InteractPepXMLPath.replace(".pep.xml", ".prot.xml");        
        parameterPath = FilenameUtils.separatorsToUnix(FilenameUtils.getFullPath(mzXMLfile) + FilenameUtils.getBaseName(mzXMLfile) + ".comet.param");        
    }
}
