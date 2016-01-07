/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package crosslinker;

import java.util.ArrayList;

/**
 *
 * @author Chih-Chiang Tsou <chihchiang.tsou@gmail.com>
 */
public class MALDI_Spot {

    public String SpotTag;
    public CrossLinkerScanResult MS1result;
    public ArrayList<CrossLinkerScanResult> MS2result = new ArrayList<>();

}
