/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package Test;

import java.io.IOException;
import MSUmpire.MSMSDBSearch.DBSearchParam;
import MSUmpire.MSMSDBSearch.TandemParam;
import MSUmpire.MSMSDBSearch.TandemSearch;
import org.apache.commons.io.FilenameUtils;

/**
 *
 * @author Chih-Chiang Tsou <chihchiang.tsou@gmail.com>
 */
public class TandemDBSearchTest {

    /**
     * @param args the command line arguments
     * @throws java.io.IOException
     * @throws java.lang.InterruptedException
     */
    public static void main(String[] args) throws IOException, InterruptedException {

        String SpectraFile = "ec01.raw.mzXML";
//        String FastaFile = "UPS_1_2_Ecoli.plusREV.fa";
        String Folder = "C:/Umich/Data/wessels_071713/";
        //String SpectraFile = "14343_UPS1_400fm_Ecolilysate_IDA_5600_ABConvert.mzXML";
        String FastaFile = "UPS_1_2_Ecoli_PlusRevTag2.fa";
        //String DecoyFasta="UPS_1_2_Ecoli_ShufTag.fa";
        //String Folder = "C:/inetpub/wwwroot/ISB/data/UPS1_Ecoli/";
        TandemParam para = new TandemParam(DBSearchParam.SearchInstrumentType.Orbit_Velos);
        para.FastaPath = FilenameUtils.separatorsToUnix(Folder + FastaFile);
        //para.DecoyFasta=FilenameUtils.separatorsToUnix(Folder+DecoyFasta);
        para.NoCPUs = 3;
        para.parameterPath = FilenameUtils.separatorsToUnix(Folder + FilenameUtils.getBaseName(SpectraFile) + "_tandem.para");
        //para.decoyparamterPath= FilenameUtils.separatorsToUnix(Folder + FilenameUtils.getBaseName(SpectraFile) + "_decoytandem.para");
        para.RawSearchResult = FilenameUtils.separatorsToUnix(Folder + FilenameUtils.getBaseName(SpectraFile) + ".tandem");
        //para.OutputTandemDecoyPath= FilenameUtils.separatorsToUnix(Folder + FilenameUtils.getBaseName(SpectraFile) + "_decoy.tandem");
        para.InteractPepXMLPath = FilenameUtils.separatorsToUnix(Folder + "interact-" + FilenameUtils.getBaseName(SpectraFile) + ".pep.xml");
para.ProtXMLPath=para.InteractPepXMLPath.replace(".pep.xml", ".prot.xml");        
//para.OutputSeqPath=Filepath+FilenameUtils.getBaseName(mzxML)+"_"+precursorppm+"_"+fragppm+".seq";
        para.SpectrumPath = FilenameUtils.separatorsToUnix(Folder + SpectraFile);
        TandemSearch search = new TandemSearch(para);
        search.RunTandem();
    }

}
