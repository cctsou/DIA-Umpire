/**
 * Copyright (c) 2010 Institute for Systems Biology
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 *
 * Contact info:
 * Mi-Youn K. Brusniak
 * Insitute for Systems Biology
 * 1441 North 34th St.
 * Seattle, WA  98103  USA
 * mbrusniak@systemsbiology.org
 *
 */
package ExternalPackages.org.systemsbiology.apps.tramlcreator;


import ExternalPackages.org.hupo.psi.ms.traml.TraMLType;
import ExternalPackages.org.systemsbiology.constants.JTRAML_URL;

import javax.xml.bind.JAXBContext;
import javax.xml.bind.JAXBElement;
import javax.xml.bind.JAXBException;
import javax.xml.bind.Marshaller;
import javax.xml.namespace.QName;
import java.io.StringWriter;

public class TraMLCreator {

    // Define instance variables.
    private TraMLType traML;


    // Constructor
    public TraMLCreator() {

    }

    /* File parse method.
    *
    */
    public void setTraML(TraMLType newTraML) {
        this.traML = newTraML;
    }

    public String asString() throws JAXBException {

        java.io.StringWriter sw = new StringWriter();
        JAXBContext ctx = JAXBContext.newInstance("org.hupo.psi.ms.traml");

        Marshaller m = ctx.createMarshaller();
        m.setProperty(Marshaller.JAXB_SCHEMA_LOCATION, JTRAML_URL.TRAML_XSD_LOCATION);
        m.setProperty(Marshaller.JAXB_FORMATTED_OUTPUT, true);
        m.setProperty(Marshaller.JAXB_ENCODING, "UTF-8");        
        
        JAXBElement<TraMLType> tramlWrap =
                new JAXBElement<TraMLType>(new QName(JTRAML_URL.TRAML_URI, "TraML"), TraMLType.class, traML);

        m.marshal(tramlWrap, sw);

        return sw.toString();
    }

} // End class

