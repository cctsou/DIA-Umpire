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
package MSUmpire.BaseDataStructure;

import java.io.Serializable;

/*
 * Two dimensional data
 */
/**
 *
 * @author Chih-Chiang Tsou <chihchiang.tsou@gmail.com>
 */
public class XYData implements Comparable<XYData>, Serializable {

    private static final long serialVersionUID = 973492749274921L;

    private float X;
    private float Y;

//    private synchronized void writeObject(java.io.ObjectOutputStream stream) throws java.io.IOException {        
//        stream.defaultWriteObject();
//        stream.writeFloat(xydata[0]);        
//        stream.writeFloat(xydata[1]);        
//    }
//    private void readObject(java.io.ObjectInputStream in) throws ClassNotFoundException, IOException {
//        in.defaultReadObject();
//        xydata[0]=in.readFloat();
//        xydata[1]=in.readFloat();
//    }
    public XYData(float x, float y) {
        setX(x);
        setY(y);
    }

    @Override
    public int compareTo(XYData o) {
        return Float.compare(o.getX(), getX());
    }

    /**
     * @return the X
     */
    public float getX() {
        return X;
    }

    /**
     * @param X the X to set
     */
    public void setX(float X) {
        this.X = X;
    }

    /**
     * @return the Y
     */
    public float getY() {
        return Y;
    }

    /**
     * @param Y the Y to set
     */
    public void setY(float Y) {
        this.Y = Y;
    }

    public XYData cloneXYData() {
        return new XYData(getX(), getY());
    }
}
