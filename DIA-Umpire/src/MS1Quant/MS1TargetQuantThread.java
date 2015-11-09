/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package MS1Quant;

import MSUmpire.BaseDataStructure.InstrumentParameter;
import MSUmpire.LCMSBaseStructure.LCMSPeakMS1;
import MSUmpire.PSMDataStructure.LCMSID;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.sql.SQLException;
import java.util.HashMap;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.TimeUnit;
import java.util.zip.DataFormatException;
import javax.xml.parsers.ParserConfigurationException;
import org.apache.commons.io.FilenameUtils;
import org.apache.commons.lang.exception.ExceptionUtils;
import org.apache.log4j.FileAppender;
import org.apache.log4j.Level;
import org.apache.log4j.Logger;
import org.apache.log4j.PatternLayout;
import org.apache.log4j.Priority;
import org.xml.sax.SAXException;
import org.xmlpull.v1.XmlPullParserException;

/**
 *
 * @author Chih-Chiang Tsou
 */
public class MS1TargetQuantThread implements Runnable{

    InstrumentParameter param=null;
    File mzxmlfile=null;
    int NoCPUs=5;
    LCMSID id=null;
    String outputfolder="";
    public MS1TargetQuantThread(File mzxmlfile, LCMSID id, int NoCPUs, String outputfolder, InstrumentParameter param){
        this.param=param;
        this.mzxmlfile=mzxmlfile;
        this.NoCPUs=NoCPUs;
        this.id=id;
        this.outputfolder=outputfolder;
    }
    
    @Override
    public void run() {
        try {
//            Logger logger=Logger.getLogger(mzxmlfile.getAbsolutePath());
//            FileAppender fa = new FileAppender();
//            fa.setName("FileLogger_Debug");
//            fa.setFile(mzxmlfile.getAbsolutePath()+"_quant.log");
//            fa.setLayout(new PatternLayout("%d %-5p [%c{1}] %m%n"));
//            fa.setThreshold(Level.DEBUG);
//            fa.setAppend(false);
//            fa.activateOptions();
//            logger.addAppender(fa);
//            
//            logger.info("Processing file " + mzxmlfile.getAbsolutePath() + "....");
            Logger.getRootLogger().info("Processing file " + mzxmlfile.getAbsolutePath() + "....");
            
            long time = System.currentTimeMillis();
            LCMSPeakMS1 LCMS1 = new LCMSPeakMS1(mzxmlfile.getAbsolutePath(), NoCPUs);
            LCMS1.SetParameter(param);
            
            LCMS1.Resume = false;
            if (!param.TargetIDOnly) {
                LCMS1.CreatePeakFolder();
            }
            LCMS1.ExportPeakClusterTable = false;
            
            
            if(id.PSMList.isEmpty()){
                Logger.getRootLogger().warn("There is no PSM mapped to the file:"+mzxmlfile.getName()+", skipping the file.");
                return;
            }
            LCMS1.IDsummary=id;
            LCMS1.IDsummary.mzXMLFileName=mzxmlfile.getAbsolutePath();
            
            if(param.TargetIDOnly){
                LCMS1.SaveSerializationFile=false;
            }
            
            if (param.TargetIDOnly || !LCMS1.ReadPeakCluster()) {
                LCMS1.PeakClusterDetection();
            }
            
            LCMS1.AssignQuant(false);
            LCMS1.IDsummary.ExportPepID(outputfolder);
            
            time = System.currentTimeMillis() - time;
            //logger.info(LCMS1.ParentmzXMLName + " processed time:" + String.format("%d hour, %d min, %d sec", TimeUnit.MILLISECONDS.toHours(time), TimeUnit.MILLISECONDS.toMinutes(time) - TimeUnit.HOURS.toMinutes(TimeUnit.MILLISECONDS.toHours(time)), TimeUnit.MILLISECONDS.toSeconds(time) - TimeUnit.MINUTES.toSeconds(TimeUnit.MILLISECONDS.toMinutes(time))));
            Logger.getRootLogger().info(LCMS1.ParentmzXMLName + " processed time:" + String.format("%d hour, %d min, %d sec", TimeUnit.MILLISECONDS.toHours(time), TimeUnit.MILLISECONDS.toMinutes(time) - TimeUnit.HOURS.toMinutes(TimeUnit.MILLISECONDS.toHours(time)), TimeUnit.MILLISECONDS.toSeconds(time) - TimeUnit.MINUTES.toSeconds(TimeUnit.MILLISECONDS.toMinutes(time))));
            LCMS1.BaseClearAllPeaks();
            LCMS1.SetmzXML(null);
            LCMS1.IDsummary = null;
            LCMS1 = null;
            id.ReleaseIDs();
            id=null;
            System.gc();        
        } catch (Exception ex) {
            Logger.getRootLogger().error(ExceptionUtils.getStackTrace(ex));            
        }
    }

}
