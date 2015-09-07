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
public class DIA_ParseiProphet implements Runnable{

    DIAPack Diafile;
    DBSearchParam param;
    LCMSID combineID;
    
    public DIA_ParseiProphet(DIAPack Diafile,DBSearchParam param,LCMSID combineID){
        this.Diafile=Diafile;
        this.param=param;
        this.combineID=combineID;
    }

    @Override
    public void run() {
        try {
            Diafile.ParseSearchEngineResultiProphet(param,combineID);
            //Diafile.IDsummary.ClearPSMs();            
        } catch (ParserConfigurationException | SAXException | IOException | XmlPullParserException ex) {
            Logger.getLogger(DIA_ParseiProphet.class.getName()).log(Level.SEVERE, null, ex);
        } catch (InterruptedException ex) {
            Logger.getLogger(DIA_ParseiProphet.class.getName()).log(Level.SEVERE, null, ex);
        } catch (Exception ex) {
            Logger.getLogger(DIA_ParseiProphet.class.getName()).log(Level.SEVERE, null, ex);
        }
    }
}
