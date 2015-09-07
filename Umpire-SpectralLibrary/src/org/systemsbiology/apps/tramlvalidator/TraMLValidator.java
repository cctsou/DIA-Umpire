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
package org.systemsbiology.apps.tramlvalidator;

import psidev.psi.tools.ontology_manager.impl.local.OntologyLoaderException;
import psidev.psi.tools.validator.Validator;
import psidev.psi.tools.validator.ValidatorException;
import psidev.psi.tools.validator.ValidatorMessage;

import java.io.InputStream;
import java.util.ArrayList;
import java.util.Collection;

public class TraMLValidator extends Validator {

    public TraMLValidator(InputStream ontoConfig,
                          InputStream cvRuleConfig,
                          InputStream objectRuleConfig) throws ValidatorException, OntologyLoaderException {
        super(ontoConfig, cvRuleConfig, objectRuleConfig);
    }

    public TraMLValidator(InputStream ontoConfig,
                          InputStream cvRuleConfig) throws ValidatorException, OntologyLoaderException {
        super(ontoConfig, cvRuleConfig);
    }

    public TraMLValidator(InputStream ontoConfig) throws OntologyLoaderException {
        super(ontoConfig);
    }


    public Collection<ValidatorMessage> validate(org.hupo.psi.ms.traml.TraMLType traml) throws ValidatorException {

        final Collection<ValidatorMessage> messages = new ArrayList<ValidatorMessage>();

        // Validate CV Mapping RUles
        messages.addAll(checkCvMappingRules());

        // Run CV Mapping rules
        messages.addAll(super.checkCvMapping(traml, "/TraMLType"));

        // Run Object Rules
        messages.addAll(super.validate(traml));

        return messages;
    }
}

