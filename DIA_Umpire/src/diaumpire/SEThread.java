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
import java.util.concurrent.TimeUnit;
import org.apache.commons.lang.exception.ExceptionUtils;
import org.apache.log4j.Logger;

/**
 *
 * @author Chih-Chiang Tsou
 */
public class SEThread implements Runnable {

    DIAPack DiaFile;

    public SEThread(DIAPack DIAfile) {
        this.DiaFile = DIAfile;
    }

    @Override
    public void run() {
        try {
            long time = System.currentTimeMillis();
            Logger.getRootLogger().info("Module A: Signal extraction");
            DiaFile.process();
            DiaFile.ms1lcms.ClearAllPeaks();
            DiaFile.ClearStructure();
            time = System.currentTimeMillis() - time;
            Logger.getRootLogger().info(DiaFile.Filename + " processed time:" + String.format("%d hour, %d min, %d sec", TimeUnit.MILLISECONDS.toHours(time), TimeUnit.MILLISECONDS.toMinutes(time) - TimeUnit.HOURS.toMinutes(TimeUnit.MILLISECONDS.toHours(time)), TimeUnit.MILLISECONDS.toSeconds(time) - TimeUnit.MINUTES.toSeconds(TimeUnit.MILLISECONDS.toMinutes(time))));
        } catch (Exception ex) {
            Logger.getRootLogger().error(ExceptionUtils.getStackTrace(ex));
            System.exit(2);
        }
    }
}
