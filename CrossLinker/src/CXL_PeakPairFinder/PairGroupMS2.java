/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package CXL_PeakPairFinder;

import MSUmpire.SpectralProcessingModule.IsotopePeakGroup;

/**
 *
 * @author Chih-Chiang Tsou
 */
public class PairGroupMS2 {
    public IsotopePeakGroup LowMassPeak;
    public IsotopePeakGroup HighMassPeak;
    public IsotopePeakGroup DeadEndpairs;
    
    public PairGroupMS2(IsotopePeakGroup lowmasspeak){
        this.LowMassPeak=lowmasspeak;
    }
}
