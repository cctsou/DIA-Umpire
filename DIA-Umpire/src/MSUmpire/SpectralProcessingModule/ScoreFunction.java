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
package MSUmpire.SpectralProcessingModule;

import MSUmpire.BaseDataStructure.XYData;
import MSUmpire.BaseDataStructure.XYPointCollection;

/**
 * Score function for comparing two spectra
 * @author Chih-Chiang Tsou
 */
public class ScoreFunction {

    
    public static XYPointCollection SpectralNormalizationForScan(XYPointCollection PointCollection) {
        float total = 0f;
        for (int i = 0; i < PointCollection.PointCount(); i++) {
            total += PointCollection.Data.get(i).getY() * PointCollection.Data.get(i).getY();
        }
        total = (float) Math.sqrt(total);

        XYPointCollection result=new XYPointCollection();
        for (int i = 0; i < PointCollection.PointCount(); i++) {
            result.AddPoint(PointCollection.Data.get(i).getX(),PointCollection.Data.get(i).getY() / total);
        }
        return result;
    }
    
    public static float CalcDotProductForScan(XYPointCollection pointset, XYPointCollection pointset2) {
        int index_1 = 0;
        int index_2 = 0;
        float InnerProduct = 0;
        if (pointset == null || pointset2 == null) {
            return 0f;
        }
        while (index_1 < pointset.PointCount() && index_2 < pointset2.PointCount()) {
            float diff = pointset.Data.get(index_1).getX() - pointset2.Data.get(index_2).getX();
            if (diff == 0) {
                InnerProduct += pointset.Data.get(index_1).getY() * pointset2.Data.get(index_2).getY();
                index_2++;
                index_1++;
            } else if (diff > 0)//index_1.X > index_2.X
            {
                index_2++;
            } else if (diff < 0)//index_2.X > index_1.X
            {
                index_1++;
            }
        }
        return InnerProduct;
    }
    
    public static XYPointCollection SpectralNormalizationForPairCollection(XYPointCollection PointCollection) {
        float totalX = 0f;
        float totalY = 0f;
        for (int i = 0; i < PointCollection.PointCount(); i++) {
            totalY += PointCollection.Data.get(i).getY() * PointCollection.Data.get(i).getY();
            totalX += PointCollection.Data.get(i).getX() * PointCollection.Data.get(i).getX();
        }
        totalX = (float) Math.sqrt(totalX);
        totalY = (float) Math.sqrt(totalY);

        XYPointCollection result=new XYPointCollection();
        for (int i = 0; i < PointCollection.PointCount(); i++) {
            result.AddPoint(PointCollection.Data.get(i).getX() / totalX, PointCollection.Data.get(i).getY() / totalY);
        }
        return result;
    }
    
    public static float CalcDotProductForPairPointCollection(XYPointCollection pointset) {        
        if (pointset == null || pointset.PointCount()<=1) {
            return 0f;
        }
        float InnerProduct = 0;
        for(XYData point : pointset.Data){
            InnerProduct += point.getY() * point.getX();
        }
        return InnerProduct;
    }
    
    public static float CalcSpecCorrForPairPointCollection(XYPointCollection pointset) {        
        if (pointset == null || pointset.PointCount()<=1) {
            return 0f;
        }
        float Score = 0f;
        float sumX = pointset.GetSumX();
        float sumY = pointset.GetSumY();
        float IProduct = 0f;
        float NormX = 0f;
        float NormY = 0f;
        float xdiff = 0f;
        float ydiff = 0f;
        for (XYData point : pointset.Data) {
            xdiff = point.getX() - sumX;
            ydiff = point.getY() - sumY;
            IProduct += xdiff * ydiff;
            NormX += xdiff * xdiff;
            NormY += ydiff * ydiff;
        }
        NormX = (float) Math.sqrt(NormX);
        NormY = (float) Math.sqrt(NormY);
        Score = IProduct / (NormX * NormY);
        Score=0.5f*(1+Score);           
        return Score;
    }
    
    public static float CalcSpecContrastAngleForPairPointCollection(XYPointCollection pointset) {
        if (pointset == null || pointset.PointCount()<=1) {
            return 0f;
        }
        float Score = 0f;
        float NormX = 0f;
        float NormY = 0f;
        for (XYData point : pointset.Data) {
            NormX += point.getX() *point.getX();
            NormY += point.getY()*point.getY();            
        }
        double cosine=CalcDotProductForPairPointCollection(pointset)/(Math.sqrt(NormX)*Math.sqrt(NormY));
        if(cosine>1d){
            cosine=1d;
        }
        float angle=(float)Math.acos(cosine);        
        Score = (float)(1-(2*angle/Math.PI));            
        return Score;
    }
        
}
