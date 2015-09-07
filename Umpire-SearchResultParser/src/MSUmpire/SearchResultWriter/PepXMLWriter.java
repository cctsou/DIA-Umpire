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
package MSUmpire.SearchResultWriter;

import MSUmpire.BaseDataStructure.UmpireInfo;
import MSUmpire.FastaParser.FastaParser_V2;
import MSUmpire.PSMDataStructure.PepIonID;
import Utility.DateTimeTag;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import org.apache.commons.io.FilenameUtils;
import org.apache.log4j.Logger;
import org.xmlpull.v1.XmlPullParserException;

/**
 *
 * @author Chih-Chiang Tsou
 */
public class PepXMLWriter {

    ArrayList<PepIonID> PepList;
    StringBuilder sb = new StringBuilder();
    String Filename;
    String Fasta;

    public class ParentProtein {

        String Acc;
        char PreAA = '-';
        char NextAA = '-';
        String Des;
    }

    public PepXMLWriter(String Filename, String Fasta) throws IOException, XmlPullParserException {
        this.Filename = Filename;
        this.Fasta = Fasta;        
    }
    
    public void SetPepList(HashMap<String, PepIonID> Import){
        this.PepList=new ArrayList<>();
        for(PepIonID pepIonID : Import.values()){
            PepList.add(pepIonID);
        }                
    }

    public void SetPepList(ArrayList<PepIonID> Import){
        this.PepList=Import;
    }
        
    public void Write() throws IOException, XmlPullParserException {
        int minlength=Integer.MAX_VALUE;
        int maxlength=0;
        int maxmiss=0;
        for(PepIonID pepIonID :PepList){
            if(pepIonID.Sequence.length()>maxlength){
                maxlength=pepIonID.Sequence.length();
            }
            if(pepIonID.Sequence.length()<minlength){
                minlength=pepIonID.Sequence.length();                
            }
            if(pepIonID.getNMissedCleavages()>maxmiss){
                maxmiss=pepIonID.getNMissedCleavages();
            }
        }
        
        FastaParser_V2 fastaparser = FastaParser_V2.FasterSerialzationRead(Fasta);
        if (fastaparser == null) {
            fastaparser = new FastaParser_V2(Fasta);
            fastaparser.digestion(maxmiss, minlength, maxlength,"DECOY");
        }
        Header();
        sb.append("<msms_run_summary base_name=\"" + FilenameUtils.getBaseName(Filename) + "\" search_engine=\"Umpire\" >\n");
        SearchSummary();
        int index = 0;
        for (PepIonID pepIonID : PepList) {
            ArrayList<ParentProtein> parentprots = new ArrayList<>();
            if (!fastaparser.PeptideList.containsKey(pepIonID.Sequence)) {
                Logger.getRootLogger().warn(pepIonID.Sequence + " is not found as tryptic peptides in " + Fasta + ". Searching protein sequences..");
                for (String acc : fastaparser.ProteinList.keySet()) {
                    String ProtSeq = fastaparser.ProteinList.get(acc).Seq;
                    if (ProtSeq.contains(pepIonID.Sequence)) {
                        ParentProtein prot = new ParentProtein();
                        prot.Acc = acc;
                        prot.Des = fastaparser.ProteinList.get(acc).Des;
                        int aaindex = ProtSeq.indexOf(pepIonID.Sequence) - 1;
                        if (aaindex >= 0) {
                            prot.PreAA = ProtSeq.charAt(aaindex);
                        }
                        aaindex = ProtSeq.indexOf(pepIonID.Sequence) + pepIonID.Sequence.length();
                        if (aaindex < ProtSeq.length()) {
                            prot.NextAA = ProtSeq.charAt(aaindex);
                        }
                        parentprots.add(prot);
                    }
                }
            } else {
                for (String acc : fastaparser.PeptideList.get(pepIonID.Sequence).Proteins) {
                    String ProtSeq = fastaparser.ProteinList.get(acc).Seq;
                    if (ProtSeq.contains(pepIonID.Sequence)) {
                        ParentProtein prot = new ParentProtein();
                        prot.Acc = acc;
                        prot.Des = fastaparser.ProteinList.get(acc).Des;
                        int aaindex = ProtSeq.indexOf(pepIonID.Sequence) - 1;
                        if (aaindex >= 0) {
                            prot.PreAA = ProtSeq.charAt(aaindex);
                        }
                        aaindex = ProtSeq.indexOf(pepIonID.Sequence) + pepIonID.Sequence.length();
                        if (aaindex < ProtSeq.length()) {
                            prot.NextAA = ProtSeq.charAt(aaindex);
                        }
                        parentprots.add(prot);
                    }
                }
            }

            if (parentprots.isEmpty()) {
                Logger.getRootLogger().warn(pepIonID.Sequence + " is not found in " + Fasta );
            }
            else{
                sb.append("<spectrum_query spectrum=\"" + pepIonID.GetKey() + "\" precursor_neutral_mass=\"" + pepIonID.CalcNeutralPepMass() + "\" assumed_charge=\"" + pepIonID.Charge + "\" index=\"" + (index++) + "\" retention_time_sec=\"" + pepIonID.PeakRT * 60f + "\">\n"
                        + "<search_result>\n"
                        + "<search_hit hit_rank=\"1\" peptide=\"" + pepIonID.Sequence + "\" peptide_prev_aa=\"" + parentprots.get(0).PreAA + "\" peptide_next_aa=\"" + parentprots.get(0).NextAA + "\" protein=\"" + parentprots.get(0).Acc + "\" protein_descr=\"" + parentprots.get(0).Des + "\" num_tot_proteins=\"" + parentprots.size() + "\" num_matched_ions=\"" + pepIonID.GetFragCount()+ "\" tot_num_ions=\"" + 2 * (pepIonID.Sequence.length() - 1) + "\" calc_neutral_pep_mass=\"" + pepIonID.CalcNeutralPepMass() + "\" massdiff=\"" + (pepIonID.CalcNeutralPepMass() - pepIonID.ObservedMass()) + "\" num_missed_cleavages=\"" + pepIonID.getNMissedCleavages() + "\" is_rejected=\"0\">\n");

                for (int i = 1; i < parentprots.size(); i++) {
                    sb.append("<alternative_protein protein=\"" + parentprots.get(i).Acc + "\" protein_descr=\"" + parentprots.get(i).Des + "\" num_tol_term=\"2\" peptide_prev_aa=\"" + parentprots.get(i).PreAA + "\" peptide_next_aa=\"" + parentprots.get(i).NextAA + "\"/>");
                }
                        //                        + "<modification_info modified_peptide=\"KITIADCGQLE\">\n"
                //                        + "<mod_aminoacid_mass position=\"7\" mass=\"160.0306\"/>\n"
                //                        + "</modification_info>\n"                        
                //                        + "<analysis_result analysis=\"peptideprophet\">\n"
                //                        + "<peptideprophet_result probability=\"0.8297\" all_ntt_prob=\"(0.0059,0.3657,0.8297)\">\n"
                //                        + "<search_score_summary>\n"
                //                        + "<parameter name=\"fval\" value=\"1.6872\"/>\n"
                //                        + "<parameter name=\"ntt\" value=\"2\"/>\n"
                //                        + "<parameter name=\"nmc\" value=\"0\"/>\n"
                //                        + "<parameter name=\"massd\" value=\"4.011\"/>\n"
                //                        + "<parameter name=\"isomassd\" value=\"0\"/>\n"
                //                        + "</search_score_summary>\n"
                //                        + "</peptideprophet_result>\n"
                //                        + "</analysis_result>\n"
                sb.append("<analysis_result analysis=\"interprophet\">\n"
                        + "<interprophet_result probability=\"" + pepIonID.TargetedProbability() + "\">\n"
                        //                        + "<search_score_summary>\n"
                        //                        + "<parameter name=\"nrs\" value=\"0\"/>\n"
                        //                        + "<parameter name=\"nsi\" value=\"0\"/>\n"
                        //                        + "<parameter name=\"nsm\" value=\"0\"/>\n"
                        //                        + "</search_score_summary>\n"
                        + "</interprophet_result>\n"
                        + "</analysis_result>\n"
                        + "</search_hit>\n"
                        + "</search_result>\n"
                        + "</spectrum_query>\n");
            }

        }
        sb.append("</msms_run_summary>\n");
        sb.append("</msms_pipeline_analysis>\n");
        FileWriter writer = new FileWriter(Filename);
        writer.write(sb.toString());
        writer.close();
    }

    private void SearchSummary() {
        sb.append("<sample_enzyme name=\"trypsin\">\n"
                + "<specificity cut=\"KR\" no_cut=\"P\" sense=\"C\"/>\n"
                + "</sample_enzyme>\n"
                + "<search_summary base_name=\"" + FilenameUtils.getBaseName(Filename) + "\" search_engine=\"Umpire\" precursor_mass_type=\"monoisotopic\" fragment_mass_type=\"monoisotopic\" search_id=\"1\">\n"
                + "<search_database local_path=\"" + Fasta + "\" type=\"AA\"/>\n"
                + "<enzymatic_search_constraint enzyme=\"trypsin\" max_num_internal_cleavages=\"1\" min_number_termini=\"2\"/>\n"
                + "<aminoacid_modification aminoacid=\"C\" massdiff=\"57.0215\" mass=\"160.0306\" variable=\"N\"/>\n"
                + "<aminoacid_modification aminoacid=\"C\" massdiff=\"-17.0265\" mass=\"143.0041\" variable=\"Y\" symbol=\"^\"/>\n"
                + "<!--X! Tandem n-terminal AA variable modification-->\n"
                + "<aminoacid_modification aminoacid=\"E\" massdiff=\"-18.0106\" mass=\"111.0320\" variable=\"Y\" symbol=\"^\"/>\n"
                + "<!--X! Tandem n-terminal AA variable modification-->\n"
                + "<aminoacid_modification aminoacid=\"M\" massdiff=\"15.9949\" mass=\"147.0354\" variable=\"Y\"/>\n"
                + "<aminoacid_modification aminoacid=\"Q\" massdiff=\"-17.0265\" mass=\"111.0321\" variable=\"Y\" symbol=\"^\"/>\n"
                + "<!--X! Tandem n-terminal AA variable modification-->\n"
                + "<terminal_modification terminus=\"n\" massdiff=\"42.0106\" mass=\"43.0184\" protein_terminus=\"N\" variable=\"Y\" symbol=\"^\"/>"
                + "    </search_summary>\n");
    }

    private void Header() {
        sb.append("<?xml version=\"1.0\" encoding=\"ISO-8859-1\"?>\n");
        sb.append("<msms_pipeline_analysis date=\"" + DateTimeTag.GetTag() + "\" summary_xml=\"" + Filename + "\" xmlns=\"http://regis-web.systemsbiology.net/pepXML\" xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\" xsi:schemaLocation=\"http://sashimi.sourceforge.net/schema_revision/pepXML/pepXML_v18.xsd\">\n");
        sb.append("<analysis_summary analysis=\"Umpire\" version=\"" + UmpireInfo.GetInstance().Version + "\" time=\"" + DateTimeTag.GetTag() + "\"/>\n");
    }
}
