/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package MSUmpire.DIA;

import MSUmpire.MSMSDBSearch.DBSearchParam;
import MSUmpire.PSMDataStructure.LCMSID;
import java.io.IOException;
import java.util.logging.Level;
import java.util.logging.Logger;
import javax.xml.parsers.ParserConfigurationException;
import org.xml.sax.SAXException;
import org.xmlpull.v1.XmlPullParserException;

/**
 *
 * @author Chih-Chiang Tsou
 */
public class DIA_ParsePeptXMLThread implements Runnable{

    DIAPack Diafile;
    DBSearchParam param;
    
    public DIA_ParsePeptXMLThread(DIAPack Diafile,DBSearchParam param){
        this.Diafile=Diafile;
        this.param=param;
    }

    @Override
    public void run() {
        try {
            Diafile.ParseiProphetPepXML(param);
            //Diafile.IDsummary.ClearPSMs();            
        } catch (ParserConfigurationException | SAXException | IOException | XmlPullParserException ex) {
            Logger.getLogger(DIA_ParsePeptXMLThread.class.getName()).log(Level.SEVERE, null, ex);
        } catch (InterruptedException ex) {
            Logger.getLogger(DIA_ParsePeptXMLThread.class.getName()).log(Level.SEVERE, null, ex);
        } catch (Exception ex) {
            Logger.getLogger(DIA_ParsePeptXMLThread.class.getName()).log(Level.SEVERE, null, ex);
        }
    }
}
