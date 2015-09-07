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
package MSUmpire.mzXMLWrapper;

import MSUmpire.BaseDataStructure.ScanCollection;
import MSUmpire.BaseDataStructure.ScanData;
import MSUmpire.BaseDataStructure.XYData;
import java.io.FileWriter;
import java.io.IOException;
import java.nio.ByteBuffer;
import org.apache.commons.codec.binary.Base64;
import org.apache.log4j.Logger;

/**
 *
 * @author Chih-Chiang Tsou
 *
 */
public class mzXMLWrapper {

    ScanCollection scanCollection;
    String Filename;

    public mzXMLWrapper(ScanCollection scanCollection, String Filename) throws IOException {
        this.scanCollection = scanCollection;
        this.Filename = Filename;
        Wrap();
    }

    private void Wrap() throws IOException {
        FileWriter writer = new FileWriter(Filename);
        writer.write("<?xml version=\"1.0\" encoding=\"ISO-8859-1\"?>\n"
                + "<mzXML xmlns=\"http://sashimi.sourceforge.net/schema_revision/mzXML_3.2\nxmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\"\n"
                + "xsi:schemaLocation=\"http://sashimi.sourceforge.net/schema_revision/mzXML_3.2 http://sashimi.sourceforge.net/schema_revision/mzXML_3.2/mzXML_idx_3.2.xsd\">");
        int lastMS1ScanNo = 1;
        for (Integer ScanNo : scanCollection.ScanHashMap.keySet()) {
            ScanData scanData = scanCollection.ScanHashMap.get(ScanNo);
            writer.write("<scan num=" + ScanNo + "\"\n");
            writer.write("scanType=\"Full\"\n");
            String line = "0";
            if (scanData.centroided) {
                line = "1";
            }
            writer.write("centroided=\"" + line + "\"\n");
            writer.write("msLevel=\"" + scanData.MsLevel + "\"\n");
            writer.write("peaksCount=\"" + scanData.PointCount() + "\"\n");
            writer.write("polarity=\"+\"\n");
            writer.write("retentionTime=\"PT" + scanData.RetentionTime * 60f + "S\"\n");
            writer.write("lowMz=\"" + scanData.StartMz + "\"\n");
            writer.write("highMz=\"" + scanData.EndMz + "\"\n");
            writer.write("basePeakMz=\"" + scanData.BasePeakMz + "\"\n");
            writer.write("basePeakIntensity=\"" + scanData.BasePeakIntensity + "\"\n");
            writer.write("totIonCurrent=\"" + scanData.TotIonCurrent() + "\"\n");

            if (scanData.MsLevel == 2) {
                writer.write("<precursorMz precursorScanNum=\"" + lastMS1ScanNo + "\" precursorIntensity=\"" + scanData.PrecursorIntensity + "\" precursorCharge=\"" + scanData.PrecursorCharge + "\" activationMethod=\"" + scanData.ActivationMethod + "\">" + scanData.PrecursorMz + "</precursorMz>\n");
            } else {
                lastMS1ScanNo = ScanNo;
            }
            writer.write("<peaks compressionType=\"none\"\n");
            writer.write("compressedLen=\"0\"\n");
            writer.write("precision=\"32\"\n");
            writer.write("byteOrder=\"network\"\n");
            writer.write("contentType=\"m/z-int\">");
            writer.write(ScanDataToBase64PeakString(scanData));
            writer.write("</peaks>\n");
            writer.write("</scan>\n");
        }
        //IndexBuilder(Filename);
    }

//    private void IndexBuilder(String newmzXMLPath)
//       {           
//           StringBuilder sb = new StringBuilder();
//           
//           FileStream fs = new FileStream(newmzXMLPath, FileMode.Open);
//           byte[] arr = new byte[1];
//           byte[] arr2 = new byte[10];
//           double byteindex = 0;
//           int scancount = 0;
//           double indexOffsetindex = 0;
//
//           while (byteindex != -1)
//           {
//               fs.Read(arr, 0, 1);
//               string arrstring = ByteArrToString(arr);
//
//               if (arrstring.Contains("<"))
//               {
//                   fs.Read(arr2, 0, 6);
//                   string arr2string = ByteArrToString(arr2);
//                   if (arr2string.Contains("scan n"))
//                   {
//                       scancount++;
//                       sb.append("    <offset id=\"" + scancount.ToString() + "\">" + (byteindex-2).ToString() + "</offset>\n");
//                   }
//                   if (arr2string.Contains("index "))
//                   {
//                       indexOffsetindex = byteindex;
//
//                       break;
//                   }
//                   byteindex = byteindex + 6;
//               }
//
//               byteindex++;
//           }
//
//           fs.Close();
//           sb.Append("  </index>\n");
//           sb.Append("  <indexOffset>" + indexOffsetindex.ToString() + "</indexOffset>\n");
//           sb.Append("  <sha1>");
//
//           StreamWriter sw = new StreamWriter(newmzXMLPath, true);
//           sw.Write(sb.ToString());
//           sw.Close();
//
//           fs = new FileStream(newmzXMLPath, FileMode.Open);
//           byte[] resultarr = sha1.ComputeHash(fs);
//           fs.Close();
//
//           string resultarrstring = "";
//           for (int i = 0; i < resultarr.Length;i++ )
//           {
//               resultarrstring = resultarrstring + resultarr[i].ToString("X");
//           }
//
//           sw = new StreamWriter(newmzXMLPath, true);
//           sw.Write(resultarrstring.ToLower() + "</sha1>\n");
//           sw.WriteLine("</mzXML>");
//           sw.Close();
//       }
    private String ScanDataToBase64PeakString(ScanData scan) {
        int offset;
        int peaksCount = scan.PointCount();

        byte[] decoded = new byte[peaksCount * 8];

        ByteBuffer byteBuffer = ByteBuffer.wrap(decoded);

        for (int i = 0; i < peaksCount; i++) {
            XYData xy = scan.Data.get(i);
            byteBuffer.putFloat(xy.getX());
            byteBuffer.putFloat(xy.getY());
        }
        Base64 encoder2 = new Base64();
        //sun.misc.BASE64Encoder encoder=new BASE64Encoder();
        byteBuffer = null;
        Logger.getRootLogger().info(encoder2.encode(decoded).toString() + "\n");
        //System.out.print(encoder.encode(decoded)+"\n");
        //return encoder.encode(decoded);
        return encoder2.encode(decoded).toString();
    }
}
