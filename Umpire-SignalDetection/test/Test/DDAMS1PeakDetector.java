/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package Test;

import MSUmpire.BaseDataStructure.InstrumentParameter;
import MSUmpire.LCMSBaseStructure.LCMSPeakMS1;
import MSUmpire.MySQLTool.ConnectionManager;
import MSUmpire.PSMDataStructure.LCMSID;
import MSUmpire.PSMDataStructure.PSM;
import MSUmpire.SearchResultParser.PepXMLParser;
import java.io.*;
import java.sql.SQLException;
import java.util.ArrayList;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.TimeUnit;
import java.util.zip.DataFormatException;
import javax.xml.parsers.ParserConfigurationException;
import org.apache.commons.io.FilenameUtils;
import org.xml.sax.SAXException;
import org.xmlpull.v1.XmlPullParserException;

/**
 *
 * @author Chih-Chiang Tsou <chihchiang.tsou@gmail.com>
 */
public class DDAMS1PeakDetector {

    /**
     * @param args the command line arguments
     *
     * Path mzxmlname pepxmlname protxmlname startRT endRT
     */
    public static void main(String[] args) throws InterruptedException, FileNotFoundException, ExecutionException, IOException, ParserConfigurationException, DataFormatException, SAXException, SQLException, XmlPullParserException {
        System.out.print("=================================================================================================\n");
        System.out.print("Parameter file:" + args[0] + "\n");
        BufferedReader reader = new BufferedReader(new FileReader(args[0]));
        String line = "";
        String Dir = "";
        String PepXML;
        InstrumentParameter para = null;
        ArrayList<String[]> FileArrayList = new ArrayList();
        int NoCPUs = 2;
        ConnectionManager connectionManager = null;
        float pepthreshold = 0.8f;
        while ((line = reader.readLine()) != null) {
            if (line.startsWith("####Folder:")) {
                Dir = reader.readLine().trim();
                System.out.print("Working folder: " + Dir + "\n");
            } else if (line.startsWith("####Files:")) {
                while (!(line = reader.readLine()).startsWith("####End Files")) {
                    String[] Pair = new String[2];
                    Pair[0] = line.split(" ")[0];
                    if (line.split(" ").length > 1) {
                        Pair[1] = line.split(" ")[1];
                    }
                    FileArrayList.add(Pair);
                }
            } else if (line.startsWith("####Instrument:")) {
                String parastring = reader.readLine().trim();
                switch (parastring) {
                    case "TOF5600": {
                        para = new InstrumentParameter(InstrumentParameter.InstrumentType.TOF5600_msconvert);
                        break;
                    }
                    case "TOF5600_ABConverted": {
                        para = new InstrumentParameter(InstrumentParameter.InstrumentType.TOF5600_ABcentroid);
                        break;
                    }
                    case "Orbitrap": {
                        para = new InstrumentParameter(InstrumentParameter.InstrumentType.Orbitrap);
                        break;
                    }
                }
                System.out.print("Instrument type: " + para.InsType.toString() + "\n");
            } else if (line.startsWith("####No of CPUs")) {
                NoCPUs = Integer.parseInt(reader.readLine().trim());
                System.out.print("No of CPUs: " + NoCPUs + "\n");
            } else if (line.startsWith("####peptide prob threshold")) {
                pepthreshold = Float.parseFloat(reader.readLine().trim());
                System.out.print("PepXML prob threshold: " + pepthreshold + "\n");
            } else if (line.startsWith("####MYSQL DB INFO")) {
                String address = reader.readLine().split(":")[1];
                int port = Integer.parseInt(reader.readLine().split(":")[1]);
                String db = reader.readLine().split(":")[1];
                String user = reader.readLine().split(":")[1];
                String pwd = reader.readLine().split(":")[1];
                connectionManager = new ConnectionManager(address, port, db, user, pwd);
            }
        }

        new File(Dir).mkdirs();
        System.out.print("IDA / DDA peak cluster detection\n");
        for (String[] Filepair : FileArrayList) {
            System.out.print("=================================================================================================\n");
            long time = System.currentTimeMillis();
            LCMSID summary = new LCMSID(FilenameUtils.getBaseName(Filepair[0]),"rev_","");

            if (Filepair[1] != "") {
                PepXMLParser parser = new PepXMLParser(summary, Dir + Filepair[1], pepthreshold);
            }
            if (summary.PSMList.size() > 0) {
                int decoy = 0;
                for (PSM psm : summary.PSMList.values()) {
                    boolean isdecoy = true;
                    for (String ProID : psm.ParentProtIDs) {
                        if (!ProID.startsWith("rev_")) {
                            isdecoy = false;
                            break;
                        }
                    }
                    if (isdecoy) {
                        decoy++;
                    }
                }
                System.out.print("No. of identified spectra:" + summary.PSMList.size() + "(" + summary.GetPepIonList().size() + " peptide ions) decoy:" + decoy + "\n");
            } else {
                summary = null;
            }

            LCMSPeakMS1 QuantSummary = new LCMSPeakMS1(Dir + Filepair[0], NoCPUs);
            
        QuantSummary.SetParameter(para);
        //QuantSummary.SetMySQLConnection(connectionManager);
            QuantSummary.AssignIDResult(summary);
            QuantSummary.PeakClusterDetection();
            time = System.currentTimeMillis() - time;

            System.out.print("Processed time:" + String.format("%d hour, %d min, %d sec\n", TimeUnit.MILLISECONDS.toHours(time), TimeUnit.MILLISECONDS.toMinutes(time) - TimeUnit.HOURS.toMinutes(TimeUnit.MILLISECONDS.toHours(time)), TimeUnit.MILLISECONDS.toSeconds(time) - TimeUnit.MINUTES.toSeconds(TimeUnit.MILLISECONDS.toMinutes(time))));
            QuantSummary = null;
            summary = null;
            System.gc();
        }
    }
}
