/*
 * Copyright 2014 Chih-Chiang Tsou.
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
package diaumpire;

import MSUmpire.DIA.DIAPack;
import MSUmpire.MSMSDBSearch.MSMSDBSearch;
import java.util.ArrayList;
import org.apache.commons.lang.exception.ExceptionUtils;
import org.apache.log4j.Logger;

/**
 *
 * @author Chih-Chiang Tsou
 */
public class TPPThread implements Runnable {

    DIAPack DiaFile;
    ArrayList<MSMSDBSearch> searches=new ArrayList<>();
    String MGFtag="";

    public TPPThread(DIAPack DIAfile, ArrayList<MSMSDBSearch> searches, String MGFtag) throws CloneNotSupportedException {
        this.DiaFile = DIAfile;
        this.MGFtag=MGFtag;
        for (MSMSDBSearch search : searches) {
            this.searches.add(search.Clone());
        }
    }

    @Override
    public void run() {
        try {
            for (MSMSDBSearch dbsearch : searches) {
                DiaFile.QcombineProteinProphet(dbsearch,MGFtag);                
            }
            //DiaFile.iProphet(searches);
            DiaFile.Qsplit_iProphet(searches);
        } catch (Exception ex) {
            Logger.getRootLogger().error(ExceptionUtils.getStackTrace(ex));
        }
    }
}
