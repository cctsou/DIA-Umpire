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
package MSUmpire.PeakDataStructure;

import MSUmpire.BaseDataStructure.XYData;
import MSUmpire.BaseDataStructure.XYPointCollection;

/**
 * B-spline smoothing
 * @author Chih-Chiang Tsou <chihchiang.tsou@gmail.com>
 */
public class Bspline {

    private float[] bspline_T_ = null;

    public XYPointCollection Run(XYPointCollection data, int PtNum, int smoothDegree) {
        XYPointCollection bsplineCollection = new XYPointCollection();
        int p = smoothDegree;
        int n = data.Data.size() - 1;
        int m = data.Data.size() + p;
        bspline_T_ = new float[m + p];

        if (data.Data.size() <= p) {
            return data;
        }

        for (int i = 0; i <= n; i++) {
            bspline_T_[i] = 0;
            bspline_T_[m - i] = 1;
        }
        float intv = 1.0f / (m - 2 * p);
        for (int i = 1; i <= (m - 1); i++) {
            bspline_T_[p + i] = bspline_T_[p + i - 1] + intv;
        }

        float t;
        for (int i = 0; i <= PtNum; i++) {
            t = ((float) i / PtNum);
            XYData pt = getbspline(data, t, n, p);
            bsplineCollection.AddPoint(pt);
        }
        if (bsplineCollection.Data.get(bsplineCollection.PointCount() - 1).getX() < data.Data.get(data.PointCount() - 1).getX()) {
            bsplineCollection.AddPoint(data.Data.get(data.PointCount() - 1));
        }
        if (bsplineCollection.Data.get(0).getX() > data.Data.get(0).getX()) {
            bsplineCollection.AddPoint(data.Data.get(0));
        }
        bsplineCollection.Data.Finalize();
        return bsplineCollection;
    }

    XYData getbspline(XYPointCollection data, float t, int n, int p) {
        XYData pt = new XYData(0, 0);

        int itp = 0;
        for (int i = 0; i <= n; i++) {
            pt.setX(pt.getX() + data.Data.get(itp).getX() * bspline_base(i, p, t));
            pt.setY(pt.getY() + data.Data.get(itp).getY() * bspline_base(i, p, t));
            itp++;
        }
        return pt;
    }

    float bspline_base(int i, int p, float t) {
        float n, c1, c2;
        float tn1 = 0;
        float tn2 = 0;
        if (p == 0) {
            if (bspline_T_[i] <= t && t < bspline_T_[i + 1] && bspline_T_[i] < bspline_T_[i + 1]) {
                n = 1;
            } else {
                n = 0;
            }
        } else {
            if ((bspline_T_[i + p] - bspline_T_[i]) == 0) {
                c1 = 0;
            } else {
                tn1 = bspline_base(i, (p - 1), t);
                c1 = (t - bspline_T_[i]) / (bspline_T_[i + p] - bspline_T_[i]);
            }
            if ((bspline_T_[i + p + 1] - bspline_T_[i + 1]) == 0) {
                c2 = 0;
            } else {
                tn2 = bspline_base((i + 1), (p - 1), t);
                c2 = (bspline_T_[i + p + 1] - t) / (bspline_T_[i + p + 1] - bspline_T_[i + 1]);
            }
            n = (c1 * tn1) + (c2 * tn2);
        }
        return n;
    }
}
