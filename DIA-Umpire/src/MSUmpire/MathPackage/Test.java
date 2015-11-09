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
package MSUmpire.MathPackage;

import ExtPackages.jMEF.MixtureModel;
import ExtPackages.jMEF.PVector;
import ExtPackages.jMEF.ExpectationMaximization1D;
import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Vector;

/**
 *
 * @author Chih-Chiang Tsou <chihchiang.tsou@gmail.com>
 */
public class Test {

    public static void main(String[] args) throws FileNotFoundException, IOException {

        BufferedReader reader = new BufferedReader(new FileReader("C:/Umich/Box Sync/Default Sync Folder/MS1_Score/a_study6_6qc1_peakcluster.csv"));
        String line = "";
        reader.readLine();

        // Variables
        int n = 2;

        // Initial mixture model
        ArrayList<Double> IDList = new ArrayList<>();
        ArrayList<Double> UnIDList = new ArrayList<>();
        double meanID = 0d;
        double meanUnID = 0d;
        double SDID = 0d;
        double SDUnID = 0d;

        while ((line = reader.readLine()) != null) {
            String identified = line.split(",")[3];
            if (identified.equals("1")) {
                IDList.add(Double.parseDouble(line.split(",")[38]));
                meanID += Double.parseDouble(line.split(",")[38]);
            } else {
                UnIDList.add(Double.parseDouble(line.split(",")[38]));
                meanUnID += Double.parseDouble(line.split(",")[38]);
            }
        }
        meanID /= IDList.size();
        meanUnID /= UnIDList.size();
        double[] obs=new double[UnIDList.size()];
        PVector[] points = new PVector[UnIDList.size()];
        for (double value : IDList) {
            SDID += (value - meanID) * (value - meanID);
        }
        SDID /= IDList.size();
        int i = 0;
        for (double value : UnIDList) {
            SDUnID += (value - meanUnID) * (value - meanUnID);
            points[i] = new PVector(1);
            points[i].array[0] = value;
            obs[i]=value;
            i++;
        }
        SDUnID /= UnIDList.size();       
                
        FileWriter writer = new FileWriter("C:\\Users\\Tsou\\Downloads\\kerneltest.xls");
        FileWriter writer2 = new FileWriter("C:\\Users\\Tsou\\Downloads\\kerneltest2.xls");
        for (double value : obs) {
            writer.write(value+"\n");
        }
        
        writer.close();
        writer2.close();
        Vector<PVector>[] clusters = KMeans.run(points, n);
        MixtureModel mmc = ExpectationMaximization1D.initialize(clusters);
        mmc = ExpectationMaximization1D.run(points, mmc);

//        MixtureModel mm = new MixtureModel(n);
//        mm.EF = new UnivariateGaussian();
//        PVector param = new PVector(2);
//        param.array[0] = meanID;
//        param.array[1] = SDID;
//        mm.param[0] = param;
//        mm.weight[0] = 1;
//        PVector param2 = new PVector(2);
//        param2.array[0] = meanUnID;
//        param2.array[1] = SDUnID;
//        mm.param[1] = param2;
//        mm.weight[1] = 2;
//        mm.normalizeWeights();
//        System.out.println("Initial mixure model \n" + mm + "\n");
        // Draw points from initial mixture model and compute the n clusters
//        PVector[] points = mm.drawRandomPoints(m);
//        Vector<PVector>[] clusters = KMeans.run(points, n);
//
//        // Classical EM
//        MixtureModel mmc;
//        mmc = ExpectationMaximization1D.initialize(clusters);
        //MixtureModel mmc;
        //mmc = ExpectationMaximization1D.run(points, mm);
        System.out.println("Mixure model estimated using classical EM \n" + mmc + "\n");

        PVector point = new PVector(2);
        System.out.print("Mean:" + ((PVector) mmc.param[0]).array[0] + "\nSD:" + ((PVector) mmc.param[0]).array[1] + "\n");
        point.array[0] = ((PVector) mmc.param[0]).array[0] - 2 * ((PVector) mmc.param[0]).array[1];
        point.array[1] = 0f;
        point.array[1] = mmc.EF.density(point, mmc.param[0]);
        System.out.print("Mean:" + point.array[0] + "\nSD:" + point.array[1]);

    }

}
