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
package MSUmpire.BaseDataStructure;

import org.apache.commons.io.FilenameUtils;

/**
 *
 * @author Chih-Chiang Tsou <chihchiang.tsou@gmail.com>
 */
public class TandemParam extends DBSearchParam {
    
    public TandemParam(SearchInstrumentType type) {
        defaultType = type;
        SetParameter(type);
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
}
