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
public class MSGFParam extends DBSearchParam {
    public String idconvert = "";    
    public String msgfpath="";
    public String JavaPath="java";
    public String msgfcommand="";
    
    public MSGFParam(SearchInstrumentType type) {
        defaultType = type;
        SetParameter(type);
    }

    @Override
    public void GenerateParamFile() {
        try {
            InputStream is = DBSearchParam.class.getClassLoader().getResourceAsStream("resource/MSGFMods.txt");
            
            if(templateParamFile!=null && new File(templateParamFile).exists()){
                Logger.getRootLogger().info("Using MSGF+ modofication template file: "+templateParamFile);
                is = new FileInputStream(templateParamFile);
            }
            BufferedReader reader = new BufferedReader(new InputStreamReader(is));
            String line = "";
            StringBuffer sb = new StringBuffer();
            while ((line = reader.readLine()) != null) {
                sb.append(line + "\n");
            }
            String output = sb.toString();
            FileWriter writer = new FileWriter(parameterPath);
            writer.write(output);
            writer.close();
        } catch (IOException ex) {
            Logger.getRootLogger().error(ExceptionUtil.printStackTrace(ex));
        }
    }
    public void SetCombineFileName(String filename, String tag){
        CombinedPepXML=FilenameUtils.separatorsToUnix(FilenameUtils.getFullPath(filename) + "interact-" + FilenameUtils.getBaseName(filename) + tag+".msgf.combine.pep.xml");        
        CombinedProt = FilenameUtils.getFullPath(filename) + FilenameUtils.getBaseName(filename) + tag+ ".msgf.Qcombine.prot.xml";
    }
    
    @Override
    public void SetResultFilePath(String mzXMLfile) {
        SpectrumPath = FilenameUtils.separatorsToUnix(FilenameUtils.getFullPath(mzXMLfile) + FilenameUtils.getName(mzXMLfile));
        RawSearchResult=FilenameUtils.separatorsToUnix(FilenameUtils.getFullPath(mzXMLfile) + FilenameUtils.getBaseName(mzXMLfile) + ".msgf.mzid");
        PepXMLPath = FilenameUtils.separatorsToUnix(FilenameUtils.getFullPath(mzXMLfile)  + FilenameUtils.getBaseName(mzXMLfile) + ".msgf.pepXML");
        InteractPepXMLPath = FilenameUtils.separatorsToUnix(FilenameUtils.getFullPath(mzXMLfile) + "interact-" + FilenameUtils.getBaseName(mzXMLfile) + ".msgf.pep.xml");
        ProtXMLPath = InteractPepXMLPath.replace(".pep.xml", ".prot.xml");
        parameterPath = FilenameUtils.separatorsToUnix(FilenameUtils.getFullPath(mzXMLfile) + FilenameUtils.getBaseName(mzXMLfile) + ".msgf.param");          
    }
}
