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
/**
 * Copyright 2008 The European Bioinformatics Institute, and others.
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
 */
package ExternalPackages.org.systemsbiology.apps.tramlvalidator;

import javax.xml.bind.JAXBContext;
import javax.xml.bind.JAXBElement;
import javax.xml.bind.JAXBException;
import javax.xml.bind.UnmarshalException;
import javax.xml.bind.Unmarshaller;

import psidev.psi.tools.ontology_manager.impl.local.OntologyLoaderException;
import psidev.psi.tools.validator.MessageLevel;
import psidev.psi.tools.validator.ValidatorException;
import psidev.psi.tools.validator.ValidatorMessage;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Collection;

import ExternalPackages.org.hupo.psi.ms.traml.TraMLType;

public class TraMLValidatorRunner {

    public static void main(String[] args) throws Exception {

        int iarg = 0;
        String option;
        String tramlFileName = null;
        String ataqsBin = null;
        String errFile = null;
        String temp_str[];

        while (iarg < args.length) {
            option = args[iarg];
            if (option.indexOf("inputfile") != -1) {
                temp_str = option.split("=");
                tramlFileName = temp_str[1];
            } else if (option.indexOf("ataqs_bin") != -1) {
                temp_str = option.split("=");
                ataqsBin = temp_str[1];
            } else if (option.indexOf("errfile") != -1) {
                temp_str = option.split("=");
                errFile = temp_str[1];
            } else {
                printOptions();
            }
            iarg++;
        }

        // check for all required options
        if ((tramlFileName == null) || (ataqsBin == null) || (errFile == null))
            printOptions();

        String errorMessage = null;
        try {

            File ontologyFile = new File(ataqsBin + File.separator + "ontologies.xml");
            File mappingRules = new File(ataqsBin + File.separator + "TraML-mapping.xml");
            File objectRules = new File(ataqsBin + File.separator + "object_rules.xml");

            TraMLValidator validator = new TraMLValidator(new FileInputStream(ontologyFile),
                    new FileInputStream(mappingRules),
                    new FileInputStream(objectRules));

            System.out.println("Create Context\n");

            JAXBContext ctx = JAXBContext.newInstance("org.hupo.psi.ms.traml");
            Unmarshaller um = ctx.createUnmarshaller();
            JAXBElement<TraMLType> jaxb_traml = (JAXBElement<TraMLType>) um.unmarshal(new File(tramlFileName));
            TraMLType traml = jaxb_traml.getValue();

            final Collection<ValidatorMessage> messages = validator.validate(traml);

            System.out.println("Validation run collected " + messages.size() + " message(s):");
            for (ValidatorMessage message : messages) {
                if (message.getLevel().isHigher(MessageLevel.INFO)) {
                    errorMessage += message.getMessage() + "\n";
                }
            }

        } catch (FileNotFoundException e) {

            errorMessage += "TraML Validation: FileNotFoundException " + e.toString() + "\n";

        } catch (ValidatorException e) {

            errorMessage += "TraML Validation: ValidatorException " + e.toString() + "\n";

        } catch (OntologyLoaderException e) {

            errorMessage += "TraML Validation: Ontology loading failed " + e.toString() + "\n";

        } catch (UnmarshalException e) {

            errorMessage += "TraML Validation: Unmarshall failed " + e.toString() + "\n";

        } catch (JAXBException e) {

            errorMessage += "TraML Validation failed " + e.toString() + "\n";
        }

        // only write out error messages if we have some and errFile has been passed as a parameter
        if ((errorMessage != null) && (errFile != null)) {

            try {
                BufferedWriter out = new BufferedWriter(new FileWriter(errFile));
                out.write(errorMessage);
                out.close();
            } catch (IOException e) {

                System.out.println("Exception writing out messages for TraMLValidator");
            }
        }
    }

    private static void printOptions() {
        System.err.println("Usage: " + TraMLValidatorRunner.class.getName());
        System.err.println("\t inputfile");
        System.err.println("\t\t Example: inputfile=Example.TraML");
        System.err.println("\t ataqs_bin");
        System.err.println("\t\t Example: ataqs_bin=/serum/ATAQS/bin");
        System.err.println("\t logfile");
        System.err.println("\t\t Example: logfile=TraMLValidator.log");
        System.exit(0);
    }

}





