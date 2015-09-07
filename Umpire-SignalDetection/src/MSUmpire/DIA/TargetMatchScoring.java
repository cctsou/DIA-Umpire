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
package MSUmpire.DIA;

import MSUmpire.BaseDataStructure.XYPointCollection;
import MSUmpire.MathPackage.KMeans;
import MSUmpire.MathPackage.KernelDensityEstimator;
import MSUmpire.MathPackage.Regression;
import MSUmpire.PSMDataStructure.PepIonID;
import jMEF.ExpectationMaximization1D;
import jMEF.MixtureModel;
import jMEF.PVector;
import java.awt.Color;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;
import java.io.Serializable;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Vector;
import javastat.multivariate.DiscriminantAnalysis;
import org.apache.commons.io.FilenameUtils;
import org.apache.commons.lang.exception.ExceptionUtils;
import org.apache.commons.math3.distribution.NormalDistribution;
import org.apache.log4j.Logger;
import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartUtilities;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.axis.NumberAxis;
import org.jfree.chart.plot.DatasetRenderingOrder;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.chart.plot.XYPlot;
import org.jfree.chart.renderer.xy.XYLineAndShapeRenderer;
import org.jfree.data.statistics.HistogramDataset;
import org.jfree.data.statistics.HistogramType;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;
import org.nustaq.serialization.FSTObjectInput;
import org.nustaq.serialization.FSTObjectOutput;
import umontreal.iro.lecuyer.probdist.WeibullDist;
import weka.classifiers.Classifier;
import weka.classifiers.Evaluation;
import weka.classifiers.functions.LibSVM;
import weka.core.Attribute;
import weka.core.Debug.Random;
import weka.core.FastVector;
import weka.core.Instance;
import weka.core.Instances;
import weka.core.SelectedTag;
import weka.core.Utils;

/**
 *
 * @author Chih-Chiang Tsou <chihchiang.tsou@gmail.com>
 */
public class TargetMatchScoring implements Serializable {

    private static final long serialVersionUID = 3845269874265L;
    public ArrayList<UmpireSpecLibMatch> libTargetMatches;
    public ArrayList<UmpireSpecLibMatch> libIDMatches;
    private boolean Terminate = false;
    String Filename;
    String LibID;
    public transient int NoBinPoints = 1000;
    MatchSubscore matchSubscore = new MatchSubscore();
    Float[][] MixtureModelProb;
    ArrayList<PeakGroupScore> decoyModelingList = new ArrayList<>();
    public transient TrainingModel model = TrainingModel.LDA;
    transient Instances TrainingData = null;
    transient LibSVM svm = null;

    public enum TrainingModel {

        SVR,
        SVM,
        LDA,
    }

    public TargetMatchScoring(String Filename, String LibID) {
        libTargetMatches = new ArrayList<>();
        libIDMatches = new ArrayList<>();
        this.Filename = Filename;
        this.LibID = LibID;
    }

    public void RemoveLowScoreEntry() {
        if (Terminate) {
            return;
        }
        Logger.getRootLogger().info("Removing low probability (<0.2) hits: total target entries: " + libTargetMatches.size());
        ArrayList<UmpireSpecLibMatch> newlist = new ArrayList<>();
        for (UmpireSpecLibMatch match : libTargetMatches) {
            if ((match.BestDecoyHit != null && match.BestDecoyHit.MixtureModelLocalProb > 0.2) || (match.BestHit != null && match.BestHit.MixtureModelLocalProb > 0.2)) {
                newlist.add(match);
            }
        }
        libTargetMatches = newlist;
        Logger.getRootLogger().info("Remaining entries: " + libTargetMatches.size());
        decoyModelingList.clear();
    }

    public void Process() throws IOException {
        //ClearGroupFragments();
        //LibraryMatchWrite();

        if (model == TrainingModel.SVM) {
            SVMTraining();
            AssignProbSVM();
        } else {
            if (model == TrainingModel.LDA) {
                EM_LDATraining();
            } else if (model == TrainingModel.SVR) {
                SVRTraining();
            }
            //MixtureModeling();
            MixtureModelingSemiSupervised();
            //GeneratePrecursorCentrolRanking();
            AssignMixtureModelProb();
        }
        ExportAlignmentResult();
        RemoveLowScoreEntry();
        ClearPepFrag();
        LibraryMatchWrite();
    }

    public void ClearPepFrag() {
        if (Terminate) {
            return;
        }
        for (UmpireSpecLibMatch match : libTargetMatches) {
            match.pepIonID.ClearPepFragFactory();
        }
    }

    public void GeneratePrecursorCentrolRanking() {
        for (UmpireSpecLibMatch match : libTargetMatches) {
            for (PeakGroupScore peakGroupScore : match.TargetHits) {
                peakGroupScore.AddScore(peakGroupScore.UmpireScore);
            }
            for (PeakGroupScore peakGroupScore : match.DecoyHits) {
                peakGroupScore.AddScore(peakGroupScore.UmpireScore);
            }
        }
        for (UmpireSpecLibMatch match : libIDMatches) {
            for (PeakGroupScore peakGroupScore : match.TargetHits) {
                peakGroupScore.AddScore(peakGroupScore.UmpireScore);
            }
            for (PeakGroupScore peakGroupScore : match.DecoyHits) {
                peakGroupScore.AddScore(peakGroupScore.UmpireScore);
            }
        }
        /////////////////////////////////calculate delta
        for (UmpireSpecLibMatch match : libTargetMatches) {
            for (PeakGroupScore peakGroupScore : match.TargetHits) {
                peakGroupScore.PrecursorCentralRank = peakGroupScore.GetScoreRank(peakGroupScore.UmpireScore);
            }
            for (PeakGroupScore peakGroupScore : match.DecoyHits) {
                peakGroupScore.PrecursorCentralRank = peakGroupScore.GetScoreRank(peakGroupScore.UmpireScore);
            }
        }
        for (UmpireSpecLibMatch match : libIDMatches) {
            for (PeakGroupScore peakGroupScore : match.TargetHits) {
                peakGroupScore.PrecursorCentralRank = peakGroupScore.GetScoreRank(peakGroupScore.UmpireScore);
            }
            for (PeakGroupScore peakGroupScore : match.DecoyHits) {
                peakGroupScore.PrecursorCentralRank = peakGroupScore.GetScoreRank(peakGroupScore.UmpireScore);
            }
        }
    }

    public void AssignModelingDecoy() {
        for (UmpireSpecLibMatch match : libTargetMatches) {
            for (PeakGroupScore peakGroupScore : match.DecoyHits) {
                decoyModelingList.add(peakGroupScore);
            }
        }
        for (UmpireSpecLibMatch match : libIDMatches) {
            for (PeakGroupScore peakGroupScore : match.DecoyHits) {
                decoyModelingList.add(peakGroupScore);
            }
        }
    }

    private void ExportAlignmentResult() throws IOException {
        try {
            FileWriter writer = new FileWriter(FilenameUtils.getFullPath(Filename) + "/" + FilenameUtils.getBaseName(Filename) + "_" + LibID + "_TargetReExtraction.xls");
            writer.write("Modseq\tPredictRT\tPrecursorMz\tType\tClusterRT\tMS level\tApexDeltaScore\tCorrScore\tSpecDotProduct\tSpecCorrelation\tSpecConstrastAngle\tPPMScore\tRTOverlap\tIntScore\tMS1Corr\tNoB\tNoY\tRTdiff\tPrcursorPPM\tSumCorrScore\tSumCorrPPMScore\tPrecursorIsoPattern\tPrecursorCentralRank\tRank\tU-Score\tProbability\n");
            for (UmpireSpecLibMatch match : libTargetMatches) {
                int Rank = 1;
                for (PeakGroupScore target : match.TargetHits) {
                    writer.write(match.pepIonID.ModSequence + "\t" + match.pepIonID.PredictRTString() + "\t" + match.pepIonID.NeutralPrecursorMz() + "\tTarget\t" + target.PrecursorRT + "\t" + target.MSlevel + "\t" + target.ApexDeltaScore + "\t" + target.AveCorrScore + "\t" + target.SpecDotProduct + "\t" + target.SpecCorrelation + "\t" + target.ContrastAngle + "\t" + target.PPMScore + "\t" + target.RTOverlapScore + "\t" + target.FragIntAvgScore + "\t" + target.PrecursorCorr + "\t" + target.NoMatchB + "\t" + target.NoMatchY + "\t" + target.RTDiff + "\t" + target.PrecursorPPM + "\t" + target.SumCorrScore + "\t" + target.SumCorrPPMScore + "\t" + target.PrecursorIsoPattern + "\t" + target.PrecursorCentralRank + "\t" + (Rank++) + "\t" + target.UmpireScore + "\t" + target.MixtureModelLocalProb + "\n");
                }
                Rank = 1;
                for (PeakGroupScore target : match.DecoyHits) {
                    writer.write(match.pepIonID.ModSequence + "\t" + match.pepIonID.PredictRTString() + "\t" + match.pepIonID.NeutralPrecursorMz() + "\tTarget_Decoy\t" + target.PrecursorRT + "\t" + target.MSlevel + "\t" + target.ApexDeltaScore + "\t" + target.AveCorrScore + "\t" + target.SpecDotProduct + "\t" + target.SpecCorrelation + "\t" + target.ContrastAngle + "\t" + target.PPMScore + "\t" + target.RTOverlapScore + "\t" + target.FragIntAvgScore + "\t" + target.PrecursorCorr + "\t" + target.NoMatchB + "\t" + target.NoMatchY + "\t" + target.RTDiff + "\t" + target.PrecursorPPM + "\t" + target.SumCorrScore + "\t" + target.SumCorrPPMScore + "\t" + target.PrecursorIsoPattern + "\t" + target.PrecursorCentralRank + "\t" + (Rank++) + "\t" + target.UmpireScore + "\t" + target.MixtureModelLocalProb + "\n");
                }
            }
            for (UmpireSpecLibMatch match : libIDMatches) {
                int Rank = 1;
                for (PeakGroupScore target : match.TargetHits) {
                    writer.write(match.pepIonID.ModSequence + "\t" + match.pepIonID.PredictRTString() + "\t" + match.pepIonID.NeutralPrecursorMz() + "\tID\t" + target.PrecursorRT + "\t" + target.MSlevel + "\t" + target.ApexDeltaScore + "\t" + target.AveCorrScore + "\t" + target.SpecDotProduct + "\t" + target.SpecCorrelation + "\t" + target.ContrastAngle + "\t" + target.PPMScore + "\t" + target.RTOverlapScore + "\t" + target.FragIntAvgScore + "\t" + target.PrecursorCorr + "\t" + target.NoMatchB + "\t" + target.NoMatchY + "\t" + target.RTDiff + "\t" + target.PrecursorPPM + "\t" + target.SumCorrScore + "\t" + target.SumCorrPPMScore + "\t" + target.PrecursorIsoPattern + "\t" + target.PrecursorCentralRank + "\t" + (Rank++) + "\t" + target.UmpireScore + "\t" + target.MixtureModelLocalProb + "\n");
                }
                for (PeakGroupScore target : match.DecoyHits) {
                    writer.write(match.pepIonID.ModSequence + "\t" + match.pepIonID.PredictRTString() + "\t" + match.pepIonID.NeutralPrecursorMz() + "\tID_Decoy\t" + target.PrecursorRT + "\t" + target.MSlevel + "\t" + target.ApexDeltaScore + "\t" + target.AveCorrScore + "\t" + target.SpecDotProduct + "\t" + target.SpecCorrelation + "\t" + target.ContrastAngle + "\t" + target.PPMScore + "\t" + target.RTOverlapScore + "\t" + target.FragIntAvgScore + "\t" + target.PrecursorCorr + "\t" + target.NoMatchB + "\t" + target.NoMatchY + "\t" + target.RTDiff + "\t" + target.PrecursorPPM + "\t" + target.SumCorrScore + "\t" + target.SumCorrPPMScore + "\t" + target.PrecursorIsoPattern + "\t" + target.PrecursorCentralRank + "\t" + (Rank++) + "\t" + target.UmpireScore + "\t" + target.MixtureModelLocalProb + "\n");
                }
            }
            writer.close();
        } catch (Exception e) {
            Logger.getRootLogger().error(ExceptionUtils.getStackTrace(e));
        }
    }

    private void CalcUmpireScore() {
        Logger.getRootLogger().info("Calculating Umpire Scores");
        for (UmpireSpecLibMatch match : libTargetMatches) {
            for (PeakGroupScore target : match.TargetHits) {
                AssignUmpireScore(target);
            }
            for (PeakGroupScore target : match.DecoyHits) {
                AssignUmpireScore(target);
            }
            match.Ranking();
        }
        for (UmpireSpecLibMatch match : libIDMatches) {
            for (PeakGroupScore target : match.TargetHits) {
                AssignUmpireScore(target);
            }
            for (PeakGroupScore target : match.DecoyHits) {
                AssignUmpireScore(target);
            }
            match.Ranking();
        }
    }

    public void AssignUmpireScore(PeakGroupScore peakgroup) {
        float[] Subs = matchSubscore.GetSubScoreArray(peakgroup);
        peakgroup.UmpireScore = 0f;
        int idx = 0;
        if (model == TrainingModel.LDA) {
            for (int i = 0; i < Subs.length; i++) {
                if (matchSubscore.Subscores[i].enable) {
                    peakgroup.UmpireScore += (float) (matchSubscore.SubSCoeff[idx++] * Subs[i]);
                }
            }
        } else if (model == TrainingModel.SVR) {
            Instance ins = matchSubscore.GetFeatureInstance(peakgroup);
            ins.setDataset(TrainingData);
            try {
                peakgroup.UmpireScore = (float) svm.classifyInstance(ins);
            } catch (Exception ex) {
                Logger.getRootLogger().error(ExceptionUtils.getStackTrace(ex));
            }
        } else if (model == TrainingModel.SVM) {
            Instance ins = matchSubscore.GetFeatureInstance(peakgroup);
            ins.setDataset(TrainingData);
            try {
                double[] temp = svm.distributionForInstance(ins);
                peakgroup.MixtureModelLocalProb = (float) temp[0];
                peakgroup.UmpireScore = (float) temp[0];
            } catch (Exception ex) {
                Logger.getRootLogger().error(ExceptionUtils.getStackTrace(ex));
            }
        }
    }

    private void EM_LDATraining() throws IOException {

        Regression regression = new Regression();
        XYPointCollection points = new XYPointCollection();
        matchSubscore.InitializeLDACoeff();
        int NoLDAComp = matchSubscore.NoEnableSubscores;

        CalcUmpireScore();

        float LDASimialrity = 0f;
        float StopThreshold = 0.95f;
        int MaxIterations = 50;
        int iteration = 0;
        float samplingratio = 0.5f;

        ArrayList<PeakGroupScore> IDList = new ArrayList<>();
        ArrayList<PeakGroupScore> decoyList = new ArrayList<>();

        for (UmpireSpecLibMatch match : libTargetMatches) {
//            for (PeakGroupScore peakGroupScore : match.DecoyHits) {
//                decoyList.add(peakGroupScore);
//            }
            if (match.BestDecoyHit != null) {
                decoyList.add(match.BestDecoyHit);
            }
            if (match.BestMS2DecoyHit != null) {
                decoyList.add(match.BestMS2DecoyHit);
            }
        }
        for (UmpireSpecLibMatch match : libIDMatches) {
            if (match.BestHit != null) {
                IDList.add(match.BestHit);
            }
//            for (PeakGroupScore peakGroupScore : match.DecoyHits) {
//                //decoyList.add(peakGroupScore);
//            }
        }

//        FileWriter writer=new FileWriter(FilenameUtils.getFullPath(Filename) + "/" + FilenameUtils.getBaseName(Filename) + "_LDATraining.xls");
//        writer.write("File:"+Filename+"\n");
//        writer.write("No. of identified groups:"+idcount+"\n");
//        writer.write("No. of target groups:"+targetcount+"\n");
//        writer.write("No. of decoy groups:"+decoyList.size()+"\n");
        Collections.shuffle(decoyList);
        ArrayList<PeakGroupScore> decoyTList = new ArrayList<>();
        for (int i = 0; i < decoyList.size() / 2; i++) {
            decoyModelingList.add(decoyList.get(i));
        }
        for (int i = decoyList.size() / 2; i < decoyList.size(); i++) {
            decoyTList.add(decoyList.get(i));
        }
        decoyList = decoyTList;
//        for (PeakGroupScore score : decoyModelingList) {
//            decoyList.remove(score);
//        }

        int targetNo = (int) (IDList.size() * samplingratio);
        int TrainNo = Math.min(targetNo, decoyList.size());

        if (TrainNo < 5) {
            Terminate = true;
            return;
        }
        //writer.write("Training set size:"+TrainNo+"\n\n");            
        //writer.write("Iteration\t");
//        for (int i = 0; i < NoFeatures; i++) {
//            writer.write(LDACompName[i] + "\t");
//        }
//        writer.write(iteration + "\t");
//        for (int i = 0; i < NoFeatures; i++) {
//            writer.write(SubSCoeff[i] + "\t");
//        }
//        writer.write(LDASimialrity + "\n");
        while (LDASimialrity < StopThreshold && iteration < MaxIterations) {
            Collections.shuffle(decoyList);
            Collections.shuffle(IDList);

            double[][] testdata = new double[NoLDAComp][2 * TrainNo];
            double[] testgroup = new double[2 * TrainNo];

            int idx = 0;
            for (int i = 0; i < TrainNo; i++) {
                PeakGroupScore peakgroup = IDList.get(i);
                double[] enablesubscore = matchSubscore.GetEnableSubScoreArray(peakgroup);
                for (int j = 0; j < enablesubscore.length; j++) {
                    testdata[j][idx] = enablesubscore[j];
                }
                testgroup[idx] = 1;
                idx++;
            }

            for (int i = 0; i < TrainNo; i++) {
                PeakGroupScore peakgroup = decoyList.get(i);
                double[] enablesubscore = matchSubscore.GetEnableSubScoreArray(peakgroup);
                for (int j = 0; j < enablesubscore.length; j++) {
                    testdata[j][idx] = enablesubscore[j];
                }
                testgroup[idx] = 0;
                idx++;
            }
            DiscriminantAnalysis LDA = new DiscriminantAnalysis();
            int[] group = LDA.predictedGroup(testgroup, testdata, testdata);
            boolean modelvalid = true;
            for (int i = 0; i < NoLDAComp; i++) {
                if (Float.isNaN((float) LDA.linearDiscriminants[0][i])) {
                    modelvalid = false;
                    break;
                }
                points.AddPoint((float) matchSubscore.SubSCoeff[i], (float) -LDA.linearDiscriminants[0][i]);
                matchSubscore.SubSCoeff[i] = -LDA.linearDiscriminants[0][i];
            }

            if (!modelvalid) {
                Logger.getRootLogger().debug("LDA failed at iteration:" + iteration);
                break;
            }
            regression.SetData(points);

            LDASimialrity = regression.GetR2();
            if (Float.isNaN(LDASimialrity)) {
                LDASimialrity = -1f;
            }

            DecimalFormat df = new DecimalFormat("#.####");
            Logger.getRootLogger().debug("----------------------------------------------------------------------------------------");
            Logger.getRootLogger().debug("Iteration:" + (iteration++));
            Logger.getRootLogger().debug("No of target hits:" + targetNo);
            Logger.getRootLogger().debug("No of decoy hits:" + decoyList.size());
            Logger.getRootLogger().debug("Training set size:" + TrainNo * 2);
            Logger.getRootLogger().debug("Trained weights:");
            for (int i = 0; i < NoLDAComp; i++) {
                Logger.getRootLogger().debug(matchSubscore.SubSName[i] + ":" + df.format(matchSubscore.SubSCoeff[i]));
            }
            Logger.getRootLogger().debug("LDA weight similarity to previous iteration:" + df.format(LDASimialrity));

            CalcUmpireScore();
//            writer.write(iteration + "\t");++) {
//                writer.write(SubSCoeff[i] + "\t");
//            }
//            writer.write(LDASimialrity + "\n");
//            for (int i = 0; i < NoFeatures; i
        }
        //writer.close();
    }

    private void SVRTraining() throws IOException {
        try {
            float samplingratio = 0.5f;
            model = TrainingModel.LDA;
            matchSubscore.InitializeLDACoeff();
            CalcUmpireScore();
            model = TrainingModel.SVR;
            ArrayList<PeakGroupScore> IDList = new ArrayList<>();
            ArrayList<PeakGroupScore> decoyList = new ArrayList<>();

            for (UmpireSpecLibMatch match : libTargetMatches) {
                if (match.BestDecoyHit != null) {
                    decoyList.add(match.BestDecoyHit);
                }
                if (match.BestMS2DecoyHit != null) {
                    decoyList.add(match.BestMS2DecoyHit);
                }
            }
            for (UmpireSpecLibMatch match : libIDMatches) {
                if (match.BestHit != null) {
                    IDList.add(match.BestHit);
                }
            }

            svm = new LibSVM();
            //svm.setShrinking(false);
            svm.setSVMType(new SelectedTag(LibSVM.SVMTYPE_NU_SVR, LibSVM.TAGS_SVMTYPE));
            svm.setKernelType(new SelectedTag(LibSVM.KERNELTYPE_RBF, LibSVM.TAGS_KERNELTYPE));

            FastVector Features = matchSubscore.GetFeatureFastVectorSVR();

            Collections.shuffle(decoyList);
            ArrayList<PeakGroupScore> decoyTList = new ArrayList<>();
            for (int i = 0; i < decoyList.size() / 2; i++) {
                decoyModelingList.add(decoyList.get(i));
            }
            for (int i = decoyList.size() / 2; i < decoyList.size(); i++) {
                decoyTList.add(decoyList.get(i));
            }
            decoyList = decoyTList;

            int targetNo = (int) (IDList.size() * samplingratio);
            int TrainNo = Math.min(targetNo, decoyList.size());

            if (TrainNo < 5) {
                Terminate = true;
                return;
            }

            Collections.shuffle(decoyList);
            Collections.shuffle(IDList);
            TrainingData = new Instances("Training", Features, TrainNo * 2);
            TrainingData.setClassIndex(TrainingData.numAttributes() - 1);

            for (int i = 0; i < TrainNo; i++) {
                PeakGroupScore peakgroup = IDList.get(i);
                Instance feature = matchSubscore.GetFeatureInstance(peakgroup);
                feature.setValue((Attribute) matchSubscore.GetFeatureFastVectorSVR().elementAt(TrainingData.numAttributes() - 1), 100);
                TrainingData.add(feature);
            }

            for (int i = 0; i < TrainNo; i++) {
                PeakGroupScore peakgroup = decoyList.get(i);
                Instance feature = matchSubscore.GetFeatureInstance(peakgroup);
                feature.setValue((Attribute) matchSubscore.GetFeatureFastVectorSVR().elementAt(TrainingData.numAttributes() - 1), 0);
                TrainingData.add(feature);
            }
            svm.buildClassifier(TrainingData);

            Evaluation eTest = new Evaluation(TrainingData);
            eTest.evaluateModel(svm, TrainingData);

            DecimalFormat df = new DecimalFormat("#.####");
            Logger.getRootLogger().debug("----------------------------------------------------------------------------------------");
            Logger.getRootLogger().debug("SVR model:");
            Logger.getRootLogger().debug(eTest.toSummaryString());

            CalcUmpireScore();
        } catch (Exception ex) {
            Logger.getRootLogger().error(ExceptionUtils.getStackTrace(ex));
        }
    }

    private void SVMTraining() throws IOException {
        try {
            ArrayList<PeakGroupScore> IDList = new ArrayList<>();
            ArrayList<PeakGroupScore> decoyList = new ArrayList<>();

            for (UmpireSpecLibMatch match : libTargetMatches) {
                for (PeakGroupScore hit : match.DecoyHits) {
                    decoyList.add(hit);
                }
            }
            for (UmpireSpecLibMatch match : libIDMatches) {
                for (PeakGroupScore hit : match.TargetHits) {
                    IDList.add(hit);
                }
            }

            svm = new LibSVM();
            svm.setShrinking(false);
            svm.setSVMType(new SelectedTag(LibSVM.SVMTYPE_NU_SVC, LibSVM.TAGS_SVMTYPE));
            svm.setKernelType(new SelectedTag(LibSVM.KERNELTYPE_RBF, LibSVM.TAGS_KERNELTYPE));
            svm.setProbabilityEstimates(true);

            FastVector Features = matchSubscore.GetFeatureFastVectorSVM();

            int TrainNo = Math.min(IDList.size(), decoyList.size());

            Collections.shuffle(decoyList);
            Collections.shuffle(IDList);
            
            TrainingData = new Instances("Training", Features, TrainNo);
            TrainingData.setClassIndex(TrainingData.numAttributes() - 1);

            for (int i = 0; i < TrainNo; i++) {
                PeakGroupScore peakgroup = IDList.get(i);
                Instance feature = matchSubscore.GetFeatureInstance(peakgroup);
                feature.setValue((Attribute) matchSubscore.GetFeatureFastVectorSVM().elementAt(TrainingData.numAttributes() - 1), "ID");
                TrainingData.add(feature);
            }

            for (int i = 0; i < TrainNo; i++) {
                PeakGroupScore peakgroup = decoyList.get(i);
                Instance feature = matchSubscore.GetFeatureInstance(peakgroup);
                feature.setValue((Attribute) matchSubscore.GetFeatureFastVectorSVM().elementAt(TrainingData.numAttributes() - 1), "Decoy");
                TrainingData.add(feature);
            }

            int seed = 319;
            int folds = 5;
            Random rand = new Random(seed);   // create seeded number generator
            Instances randData = new Instances(TrainingData);   // create copy of original data
            randData.randomize(rand);         // randomize data with number generator
            if (randData.classAttribute().isNominal()) {
                randData.stratify(folds);
            }
           //////////////////////////
            // perform cross-validation
            
            Evaluation eval = new Evaluation(randData);
            for (int n = 0; n < folds; n++) {
                Instances train = randData.trainCV(folds, n);
                Instances test = randData.testCV(folds, n);
                Classifier clsCopy = Classifier.makeCopy(svm);
                clsCopy.buildClassifier(train);
                eval.evaluateModel(clsCopy, test);
            }

            // output evaluation
            System.out.println();
            System.out.println("=== Setup ===");
            System.out.println("Classifier: " + svm.getClass().getName() + " " + Utils.joinOptions(svm.getOptions()));
            System.out.println("Dataset: " + TrainingData.relationName());
            System.out.println("Folds: " + folds);
            System.out.println("Seed: " + seed);
            System.out.println();
            System.out.println(eval.toSummaryString("=== " + folds + "-fold Cross-validation ===", false));
            ////////////////////////////
            CalcUmpireScore();
        } catch (Exception ex) {
            Logger.getRootLogger().error(ExceptionUtils.getStackTrace(ex));
        }
    }

    public void MixtureModelingByNormalKDE() {

        if (libTargetMatches.isEmpty()) {
            return;
        }
        int targetNo = 0;
        int decoyNo = 0;
        for (UmpireSpecLibMatch match : libTargetMatches) {
            targetNo += match.TargetHits.size();
            decoyNo += match.DecoyHits.size();
        }

        String pngfile = FilenameUtils.getFullPath(Filename) + "/" + FilenameUtils.getBaseName(Filename) + "_LibMatchModel.png";

        double Decoymean = 0d;
        double Decoysigma = 0d;

        for (UmpireSpecLibMatch match : libTargetMatches) {
            for (PeakGroupScore peakscore : match.DecoyHits) {
                Decoymean += peakscore.UmpireScore;
            }
        }

        Decoymean /= decoyNo;

        for (UmpireSpecLibMatch match : libTargetMatches) {
            for (PeakGroupScore peakscore : match.DecoyHits) {
                Decoysigma += (peakscore.UmpireScore - Decoymean) * (peakscore.UmpireScore - Decoymean);
            }
        }
        Decoysigma /= decoyNo;
        Decoysigma = Math.sqrt(Decoysigma);

        int idx = 0;
        int didx = 0;
        double[] TargetObs = new double[targetNo];
        double[] DecoyObs = new double[decoyNo];
        for (UmpireSpecLibMatch match : libTargetMatches) {
            for (PeakGroupScore peakscore : match.TargetHits) {
                TargetObs[idx++] = peakscore.UmpireScore;
            }
            for (PeakGroupScore peakscore : match.DecoyHits) {
                DecoyObs[didx++] = peakscore.UmpireScore;
            }
        }
        //Arrays.sort(IDObs);

        KernelDensityEstimator kde = new KernelDensityEstimator();
        List<double[]> density = null;

        NormalDistribution nor = new NormalDistribution(Decoymean, Decoysigma);

        double DecoyDensityRatio = 0f;
        for (idx = 0; idx < density.get(0).length - 1; idx++) {
            if (density.get(0)[idx] < Decoymean && density.get(0)[idx + 1] > Decoymean) {
                DecoyDensityRatio = density.get(1)[idx] * targetNo / nor.probability(Decoymean);
                break;
            }
        }

        XYSeries model1 = new XYSeries("Decoy model");
        XYSeries model2 = new XYSeries("Kernel density");
        XYSeriesCollection dataset = new XYSeriesCollection();

        HistogramDataset histogramDataset = new HistogramDataset();
        histogramDataset.setType(HistogramType.RELATIVE_FREQUENCY);
        histogramDataset.addSeries("Target matches", TargetObs, 100);
        histogramDataset.addSeries("Decoy matches", DecoyObs, 100);

        for (idx = 0; idx < density.get(0).length - 1; idx++) {
            model2.add(density.get(0)[idx], density.get(1)[idx] * targetNo);
            model1.add(density.get(0)[idx], nor.probability(density.get(0)[idx]) * DecoyDensityRatio);
        }

        dataset.addSeries(model1);
        dataset.addSeries(model2);
        JFreeChart chart = ChartFactory.createXYLineChart(FilenameUtils.getBaseName(pngfile), "Score", "Frequency", dataset, PlotOrientation.VERTICAL, true, true, false);
        XYPlot plot = chart.getXYPlot();
        plot.setBackgroundPaint(Color.white);
        plot.setDomainGridlinePaint(Color.white);
        plot.setRangeGridlinePaint(Color.white);
        plot.setForegroundAlpha(0.8f);
        chart.setBackgroundPaint(Color.white);
        try {
            ChartUtilities.saveChartAsPNG(new File(pngfile), chart, 1000, 600);
        } catch (IOException e) {
        }
    }

    public void MixtureModelingSemiSupervised() throws IOException {

        if (libTargetMatches.isEmpty() || Terminate) {
            return;
        }
        Logger.getRootLogger().info("Semi-parametric mixture modeling");
        int IDNo = 0;
        int decoyNo = 0;
        int modelNo = 0;

        for (UmpireSpecLibMatch match : libIDMatches) {
            if (match.BestHit != null) {
                IDNo++;
            }
        }
        decoyNo = decoyModelingList.size();

        for (UmpireSpecLibMatch match : libTargetMatches) {
            //modelNo+= match.TargetHits.size();
            if (match.BestMS1Hit != null) {
                modelNo++;
            }
            if (match.BestMS2Hit != null) {
                modelNo++;
            }
        }

        double[] IDObs = new double[IDNo];
        double[] DecoyObs = new double[decoyNo];
        double[] ModelObs = new double[modelNo];
        int idx = 0;
        int didx = 0;
        int midx = 0;
        for (UmpireSpecLibMatch match : libIDMatches) {
            if (match.BestHit != null) {
                IDObs[idx++] = match.BestHit.UmpireScore;
            }
        }
        for (PeakGroupScore peakGroupScore : decoyModelingList) {
            DecoyObs[didx++] = peakGroupScore.UmpireScore;
        }

        for (UmpireSpecLibMatch match : libTargetMatches) {
            if (match.BestMS1Hit != null) {
                ModelObs[midx++] = match.BestMS1Hit.UmpireScore;
            }
            if (match.BestMS2Hit != null) {
                ModelObs[midx++] = match.BestMS2Hit.UmpireScore;
            }
        }

        MixtureModelKDESemiParametric mixkde = new MixtureModelKDESemiParametric();
        mixkde.NoBinPoints = NoBinPoints;
        mixkde.SetData(ModelObs, DecoyObs, IDObs);
        mixkde.Modeling();
        mixkde.GeneratePlot(FilenameUtils.getFullPath(Filename) + "/" + FilenameUtils.getBaseName(Filename) + "_" + LibID + "_LibMatchModel.png");
        MixtureModelProb = mixkde.MixtureModelProb;
    }

    public void MixtureModelingByGammaKDE() {

        if (libTargetMatches.isEmpty()) {
            return;
        }
        int IDNo = 0;
        int decoyNo = 0;
        int modelNo = 0;
        double IDmean = 0d;
        double Decoymean = 0d;

        for (UmpireSpecLibMatch match : libIDMatches) {
            if (match.BestHit != null) {
                IDNo++;
                IDmean += match.BestHit.UmpireScore;
            }
        }

        decoyNo = decoyModelingList.size();
        double[] DecoyObs = new double[decoyNo];

        int didx = 0;
        for (PeakGroupScore peakGroupScore : decoyModelingList) {
            Decoymean += peakGroupScore.UmpireScore;
            DecoyObs[didx++] = peakGroupScore.UmpireScore;
        }

        float max = Float.NEGATIVE_INFINITY;
        float min = Float.POSITIVE_INFINITY;
        for (UmpireSpecLibMatch match : libTargetMatches) {
            for (PeakGroupScore peakscore : match.TargetHits) {
                if (peakscore.UmpireScore > max) {
                    max = peakscore.UmpireScore;
                }
                if (peakscore.UmpireScore < min) {
                    min = peakscore.UmpireScore;
                }
            }
            for (PeakGroupScore peakscore : match.DecoyHits) {
                if (peakscore.UmpireScore > max) {
                    max = peakscore.UmpireScore;
                }
                if (peakscore.UmpireScore < min) {
                    min = peakscore.UmpireScore;
                }
            }

            if (match.BestMS1Hit != null) {
                modelNo++;
            }
            if (match.BestMS2Hit != null) {
                modelNo++;
            }
        }

        Decoymean /= decoyNo;
        IDmean /= IDNo;

        double[] IDObs = new double[IDNo];
        double[] ModelObs = new double[modelNo];

        int midx = 0;
        int idx = 0;
        for (UmpireSpecLibMatch match : libIDMatches) {
            if (match.BestHit != null) {
                IDObs[idx++] = match.BestHit.UmpireScore;
            }
        }
        for (UmpireSpecLibMatch match : libTargetMatches) {
            if (match.BestMS1Hit != null) {
                ModelObs[midx++] = match.BestMS1Hit.UmpireScore;
            }
            if (match.BestMS2Hit != null) {
                ModelObs[midx++] = match.BestMS2Hit.UmpireScore;
            }
        }

        String pngfile = FilenameUtils.getFullPath(Filename) + "/" + FilenameUtils.getBaseName(Filename) + "_LibMatchModel.png";
        XYSeries model1 = new XYSeries("Incorrect matches");
        XYSeries model2 = new XYSeries("Correct matches");
        XYSeries model3 = new XYSeries("All target hits");

        KernelDensityEstimator Modelkde = new KernelDensityEstimator();
        Modelkde.SetData(ModelObs);

        //CurveFitter fitter=new CurveFitter(IDObs, IDObs)
        WeibullDist DecoyDensity = WeibullDist.getInstanceFromMLE(DecoyObs, DecoyObs.length);

        double[] model_kde_x = new double[NoBinPoints];

        float intv = (max - min) / NoBinPoints;
        for (int i = 0; i < NoBinPoints; i++) {
            model_kde_x[i] = max - i * intv;
        }

        double[] model_kde_y = Modelkde.Density(model_kde_x);
        double[] decoy_kde_y = new double[model_kde_x.length];
        for (int i = 0; i < model_kde_x.length; i++) {
            decoy_kde_y[i] = DecoyDensity.density(model_kde_x[i]);
        }
        float maxdecoy = 0f;
        int maxindex = 1;
        for (int i = 0; i < NoBinPoints; i++) {
            if (decoy_kde_y[i] > maxdecoy) {
                maxindex = i;
                maxdecoy = (float) decoy_kde_y[i];
            }
        }
        float weight = (float) (model_kde_y[maxindex] / maxdecoy);
        float positiveaccu = 0f;
        float negativeaccu = 0f;
        MixtureModelProb = new Float[NoBinPoints + 1][3];

        for (int i = 0; i < NoBinPoints; i++) {
            float positiveNumber = (float) (model_kde_y[i] - decoy_kde_y[i] * weight);
            float negativeNumber = (float) (decoy_kde_y[i] * weight);
            model2.add(model_kde_x[i], positiveNumber);
            model1.add(model_kde_x[i], negativeNumber);
            model3.add(model_kde_x[i], model_kde_y[i]);
            positiveaccu += positiveNumber;
            negativeaccu += negativeNumber;
            MixtureModelProb[i][0] = (float) model_kde_x[i];
            MixtureModelProb[i][2] = positiveNumber / (negativeNumber + positiveNumber);
            MixtureModelProb[i][1] = positiveaccu / (negativeaccu + positiveaccu);
        }
        XYSeriesCollection dataset = new XYSeriesCollection();
        dataset.addSeries(model1);
        dataset.addSeries(model2);
        dataset.addSeries(model3);

        HistogramDataset histogramDataset = new HistogramDataset();
        histogramDataset.setType(HistogramType.SCALE_AREA_TO_1);
        histogramDataset.addSeries("ID hits", IDObs, 100);
        histogramDataset.addSeries("Decoy hits", DecoyObs, 100);
        //histogramDataset.addSeries("Model hits", ModelObs, 100);

        JFreeChart chart = ChartFactory.createHistogram(FilenameUtils.getBaseName(pngfile), "Score", "Hits", histogramDataset, PlotOrientation.VERTICAL, true, false, false);
        XYPlot plot = chart.getXYPlot();

        plot.setBackgroundPaint(Color.white);
        plot.setDomainGridlinePaint(Color.white);
        plot.setRangeGridlinePaint(Color.white);
        plot.setForegroundAlpha(0.8f);
        chart.setBackgroundPaint(Color.white);

        XYLineAndShapeRenderer render = new XYLineAndShapeRenderer();
        plot.setDataset(1, dataset);
        plot.setRenderer(1, render);
        plot.setDatasetRenderingOrder(DatasetRenderingOrder.FORWARD);
        try {
            ChartUtilities.saveChartAsPNG(new File(pngfile), chart, 1000, 600);
        } catch (IOException e) {
        }
    }

    public void MixtureModeling() {
        if (libTargetMatches.isEmpty()) {
            return;
        }
        int IDNo = 0;
        int decoyNo = 0;
        int modelNo = 0;
        double IDmean = 0d;
        double Decoymean = 0d;

        for (UmpireSpecLibMatch match : libIDMatches) {
            if (match.BestHit != null) {
                IDNo++;
                IDmean += match.BestHit.UmpireScore;
            }
        }

        decoyNo = decoyModelingList.size();
        for (PeakGroupScore peakGroupScore : decoyModelingList) {
            Decoymean += peakGroupScore.UmpireScore;
        }

        for (UmpireSpecLibMatch match : libTargetMatches) {
            //modelNo+= match.TargetHits.size();
            if (match.BestMS1Hit != null) {
                modelNo++;
            }
            if (match.BestMS2Hit != null) {
                modelNo++;
            }
        }

        Decoymean /= decoyNo;
        IDmean /= IDNo;

        PVector[] points = new PVector[modelNo];
        PVector[] centroids = new PVector[2];

        int idx = 0;
        for (UmpireSpecLibMatch match : libTargetMatches) {
            if (match.BestMS1Hit != null) {
                points[idx] = new PVector(1);
                points[idx].array[0] = match.BestMS1Hit.UmpireScore;
                idx++;
            }
            if (match.BestMS2Hit != null) {
                points[idx] = new PVector(1);
                points[idx].array[0] = match.BestMS2Hit.UmpireScore;
                idx++;
            }
//            for(PeakGroupScore peakGroupScore : match.TargetHits){
//                points[idx] = new PVector(1);
//                points[idx].array[0] = match.BestMS2Hit.UmpireScore;
//                idx++;
//            }
        }

        MixtureModel mmc;
        centroids[0] = new PVector(1);
        centroids[0].array[0] = Decoymean;
        centroids[1] = new PVector(1);
        centroids[1].array[0] = IDmean;
        Vector<PVector>[] clusters = KMeans.run(points, 2, centroids);
        MixtureModel mm = ExpectationMaximization1D.initialize(clusters);
        mmc = ExpectationMaximization1D.run(points, mm);
        DecimalFormat df = new DecimalFormat("#.####");
        Logger.getRootLogger().debug("----------------------------------------------------------------------------------------");
        Logger.getRootLogger().debug("No. of modeling points=" + modelNo);
        Logger.getRootLogger().debug("ID hits mean=" + df.format(IDmean));
        Logger.getRootLogger().debug("Decoy hits mean=" + df.format(Decoymean));
        //System.out.print("T-test: p-value=" + df.format(model.ttest.pValue).toString() + "\n");
        Logger.getRootLogger().debug("Incorrect hits model mean=" + df.format(((PVector) mmc.param[0]).array[0]) + " variance=" + df.format(((PVector) mmc.param[0]).array[1]) + " weight=" + df.format(mmc.weight[0]));
        Logger.getRootLogger().debug("Correct hits model mean=" + df.format(((PVector) mmc.param[1]).array[0]) + " variance=" + df.format(((PVector) mmc.param[1]).array[1]) + " weight=" + df.format(mmc.weight[1]));

        if (((PVector) mmc.param[0]).array[0] > ((PVector) mmc.param[1]).array[0]) {
            return;
        }

        float max = (float) (((PVector) mmc.param[1]).array[0] + 4 * Math.sqrt(((PVector) mmc.param[1]).array[1]));
        float min = (float) (((PVector) mmc.param[0]).array[0] - 4 * Math.sqrt(((PVector) mmc.param[0]).array[1]));

        IDNo = 0;
        decoyNo = 0;
        modelNo = 0;

        for (PeakGroupScore peakGroupScore : decoyModelingList) {
            if (peakGroupScore.UmpireScore > min && peakGroupScore.UmpireScore < max) {
                decoyNo++;
            }
        }

        for (UmpireSpecLibMatch match : libIDMatches) {
            if (match.BestHit != null && match.BestHit.UmpireScore > min && match.BestHit.UmpireScore < max) {
                IDNo++;
            }
        }

        for (UmpireSpecLibMatch match : libTargetMatches) {
            //targetNo += match.TargetHits.size();
            //decoyNo += match.DecoyHits.size();
            if (match.BestMS1Hit != null && match.BestMS1Hit.UmpireScore > min && match.BestMS1Hit.UmpireScore < max) {
                modelNo++;
            }
            if (match.BestMS2Hit != null && match.BestMS2Hit.UmpireScore > min && match.BestMS2Hit.UmpireScore < max) {
                modelNo++;
            }
            //modelNo += match.TargetHits.size();            
        }

        double[] IDObs = new double[IDNo];
        double[] DecoyObs = new double[decoyNo];
        double[] ModelObs = new double[modelNo];
        idx = 0;
        int didx = 0;
        int midx = 0;
        for (UmpireSpecLibMatch match : libIDMatches) {
            if (match.BestHit != null && match.BestHit.UmpireScore > min && match.BestHit.UmpireScore < max) {
                IDObs[idx++] = match.BestHit.UmpireScore;
            }
        }
        for (PeakGroupScore peakGroupScore : decoyModelingList) {
            if (peakGroupScore.UmpireScore > min && peakGroupScore.UmpireScore < max) {
                DecoyObs[didx++] = peakGroupScore.UmpireScore;
            }
        }

        for (UmpireSpecLibMatch match : libTargetMatches) {
//            for(PeakGroupScore peak : match.TargetHits){
//                ModelObs[midx++]=peak.UmpireScore;
//            }
            if (match.BestMS1Hit != null && match.BestMS1Hit.UmpireScore > min && match.BestMS1Hit.UmpireScore < max) {
                ModelObs[midx++] = match.BestMS1Hit.UmpireScore;
            }
            if (match.BestMS2Hit != null && match.BestMS2Hit.UmpireScore > min && match.BestMS2Hit.UmpireScore < max) {
                ModelObs[midx++] = match.BestMS2Hit.UmpireScore;
            }
        }

        String pngfile = FilenameUtils.getFullPath(Filename) + "/" + FilenameUtils.getBaseName(Filename) + "_" + LibID + "_LibMatchModel.png";
        XYSeries model1 = new XYSeries("Incorrect matches");
        XYSeries model2 = new XYSeries("Correct matches");
        XYSeries model3 = new XYSeries("All target hits");

        int NoPoints = 1000;
        double[] model_kde_x = new double[NoPoints];
        float intv = (max - min) / NoPoints;
        PVector point = new PVector(2);
        for (int i = 0; i < NoPoints; i++) {
            point.array[0] = max - i * intv;
            model_kde_x[i] = point.array[0];
            point.array[1] = mmc.EF.density(point, mmc.param[0]) * mmc.weight[0];
            model1.add(point.array[0], point.array[1]);
            point.array[1] = mmc.EF.density(point, mmc.param[1]) * mmc.weight[1];
            model2.add(point.array[0], point.array[1]);
        }

        KernelDensityEstimator kde = new KernelDensityEstimator();
        kde.SetData(ModelObs);
        double[] model_kde_y = kde.Density(model_kde_x);

        for (int i = 0; i < NoPoints; i++) {
            if (model_kde_x[i] > min && model_kde_x[i] < max) {
                model3.add(model_kde_x[i], model_kde_y[i]);
            }
        }
        MixtureModelProb = new Float[NoPoints + 1][3];
        float positiveaccu = 0f;
        float negativeaccu = 0f;

        MixtureModelProb[0][0] = (float) model2.getMaxX() + Float.MIN_VALUE;
        MixtureModelProb[0][1] = 1f;
        MixtureModelProb[0][2] = 1f;

        for (int i = 1; i < NoPoints + 1; i++) {
            float positiveNumber = model2.getY(NoPoints - i).floatValue();
            float negativeNumber = model1.getY(NoPoints - i).floatValue();
            MixtureModelProb[i][0] = model2.getX(NoPoints - i).floatValue();
            positiveaccu += positiveNumber;
            negativeaccu += negativeNumber;
            MixtureModelProb[i][2] = positiveNumber / (negativeNumber + positiveNumber);
            MixtureModelProb[i][1] = positiveaccu / (negativeaccu + positiveaccu);
        }

        XYSeriesCollection dataset = new XYSeriesCollection();
        dataset.addSeries(model1);
        dataset.addSeries(model2);
        dataset.addSeries(model3);

        HistogramDataset histogramDataset = new HistogramDataset();
        histogramDataset.setType(HistogramType.SCALE_AREA_TO_1);
        histogramDataset.addSeries("ID hits", IDObs, 100);
        histogramDataset.addSeries("Decoy hits", DecoyObs, 100);
        //histogramDataset.addSeries("Model hits", ModelObs, 100);

        JFreeChart chart = ChartFactory.createHistogram(FilenameUtils.getBaseName(pngfile), "Score", "Hits", histogramDataset, PlotOrientation.VERTICAL, true, false, false);
        XYPlot plot = chart.getXYPlot();

        NumberAxis domain = (NumberAxis) plot.getDomainAxis();
        domain.setRange(min, max);
        plot.setBackgroundPaint(Color.white);
        plot.setDomainGridlinePaint(Color.white);
        plot.setRangeGridlinePaint(Color.white);
        plot.setForegroundAlpha(0.8f);
        chart.setBackgroundPaint(Color.white);

        XYLineAndShapeRenderer render = new XYLineAndShapeRenderer();
//        render.setSeriesPaint(0, Color.DARK_GRAY);
//        render.setSeriesPaint(1, Color.DARK_GRAY); 
//        render.setSeriesPaint(2, Color.GREEN); 
//        render.setSeriesShape(0, new Ellipse2D.Double(0, 0, 2, 2));
//        render.setSeriesShape(1, new Ellipse2D.Double(0, 0, 2, 2));
//        render.setSeriesShape(2, new Ellipse2D.Double(0, 0, 2.5f, 2.5f));
//        render.setSeriesStroke(1, new BasicStroke(1.0f));
//        render.setSeriesStroke(0, new BasicStroke(1.0f));
//        render.setSeriesStroke(2, new BasicStroke(2.0f));
        plot.setDataset(1, dataset);
        plot.setRenderer(1, render);
        plot.setDatasetRenderingOrder(DatasetRenderingOrder.FORWARD);
        try {
            ChartUtilities.saveChartAsPNG(new File(pngfile), chart, 1000, 600);
        } catch (IOException e) {
        }
    }

    public void MixtureModelingMixDecoy() throws IOException {
        if (libTargetMatches.isEmpty()) {
            return;
        }
        int IDNo = 0;
        int decoyNo = 0;
        int modelNo = 0;
        double IDmean = 0d;
        double Decoymean = 0d;

        for (UmpireSpecLibMatch match : libIDMatches) {
            if (match.BestHit != null) {
                IDNo++;
                IDmean += match.BestHit.UmpireScore;
            }
        }

        decoyNo = decoyModelingList.size();
        for (PeakGroupScore peakGroupScore : decoyModelingList) {
            Decoymean += peakGroupScore.UmpireScore;
        }

        for (UmpireSpecLibMatch match : libTargetMatches) {
            modelNo += match.TargetHits.size();
        }
        modelNo += decoyModelingList.size();

        Decoymean /= decoyNo;
        IDmean /= IDNo;

        PVector[] points = new PVector[modelNo];
        PVector[] centroids = new PVector[2];

        int idx = 0;
        for (UmpireSpecLibMatch match : libTargetMatches) {
            for (PeakGroupScore peak : match.TargetHits) {
                points[idx] = new PVector(1);
                points[idx].array[0] = peak.UmpireScore;
                idx++;
            }
        }
        for (PeakGroupScore peak : decoyModelingList) {
            points[idx] = new PVector(1);
            points[idx].array[0] = peak.UmpireScore;
            idx++;
        }

        MixtureModel mmc;
        centroids[0] = new PVector(1);
        centroids[0].array[0] = Decoymean;
        centroids[1] = new PVector(1);
        centroids[1].array[0] = IDmean;
        Vector<PVector>[] clusters = KMeans.run(points, 2, centroids);
        MixtureModel mm = ExpectationMaximization1D.initialize(clusters);
        mmc = ExpectationMaximization1D.run(points, mm);
        DecimalFormat df = new DecimalFormat("#.####");
        Logger.getRootLogger().debug("----------------------------------------------------------------------------------------");
        Logger.getRootLogger().debug("No. of modeling points=" + modelNo);
        Logger.getRootLogger().debug("ID hits mean=" + df.format(IDmean));
        Logger.getRootLogger().debug("Decoy hits mean=" + df.format(Decoymean));
        //System.out.print("T-test: p-value=" + df.format(model.ttest.pValue).toString() + "\n");
        Logger.getRootLogger().debug("Incorrect hits model mean=" + df.format(((PVector) mmc.param[0]).array[0]) + " variance=" + df.format(((PVector) mmc.param[0]).array[1]) + " weight=" + df.format(mmc.weight[0]));
        Logger.getRootLogger().debug("Correct hits model mean=" + df.format(((PVector) mmc.param[1]).array[0]) + " variance=" + df.format(((PVector) mmc.param[1]).array[1]) + " weight=" + df.format(mmc.weight[1]));

        if (((PVector) mmc.param[0]).array[0] > ((PVector) mmc.param[1]).array[0]) {
            return;
        }

        float max = (float) (((PVector) mmc.param[1]).array[0] + 4 * Math.sqrt(((PVector) mmc.param[1]).array[1]));
        float min = (float) (((PVector) mmc.param[0]).array[0] - 4 * Math.sqrt(((PVector) mmc.param[0]).array[1]));

        ///////////////////////////////////////////////////////////////////////////////////////////////
        IDNo = 0;
        decoyNo = 0;
        modelNo = 0;

        for (PeakGroupScore peakGroupScore : decoyModelingList) {
            if (peakGroupScore.UmpireScore > min && peakGroupScore.UmpireScore < max) {
                decoyNo++;
            }
        }

        for (UmpireSpecLibMatch match : libIDMatches) {
            if (match.BestHit != null && match.BestHit.UmpireScore > min && match.BestHit.UmpireScore < max) {
                IDNo++;
            }
        }

        for (UmpireSpecLibMatch match : libTargetMatches) {
            for (PeakGroupScore peak : match.TargetHits) {
                if (peak.UmpireScore > min && peak.UmpireScore < max) {
                    modelNo++;
                }
            }
        }

        double[] IDObs = new double[IDNo];
        double[] DecoyObs = new double[decoyNo];
        double[] ModelObs = new double[modelNo];
        idx = 0;
        int didx = 0;
        int midx = 0;
        for (UmpireSpecLibMatch match : libIDMatches) {
            if (match.BestHit != null && match.BestHit.UmpireScore > min && match.BestHit.UmpireScore < max) {
                IDObs[idx++] = match.BestHit.UmpireScore;
            }
        }
        for (PeakGroupScore peakGroupScore : decoyModelingList) {
            if (peakGroupScore.UmpireScore > min && peakGroupScore.UmpireScore < max) {
                DecoyObs[didx++] = peakGroupScore.UmpireScore;
            }
        }

        for (UmpireSpecLibMatch match : libTargetMatches) {
            for (PeakGroupScore peak : match.TargetHits) {
                if (peak.UmpireScore > min && peak.UmpireScore < max) {
                    ModelObs[midx++] = peak.UmpireScore;
                }
            }
        }

        int NoPoints = 1000;
        double[] model_kde_x = new double[NoPoints];
        float intv = (max - min) / NoPoints;
        for (int i = 0; i < NoPoints; i++) {
            model_kde_x[i] = max - i * intv;
        }

        KernelDensityEstimator kde_decoy = new KernelDensityEstimator();
        kde_decoy.SetData(DecoyObs);
        double[] model_kde_decoy_y = kde_decoy.Density(model_kde_x);

        String pngfile = FilenameUtils.getFullPath(Filename) + "/" + FilenameUtils.getBaseName(Filename) + "_" + LibID + "_LibMatchModel.png";
        String modelfile = FilenameUtils.getFullPath(Filename) + "/" + FilenameUtils.getBaseName(Filename) + "_" + LibID + "_ModelPoints.txt";
        FileWriter writer = new FileWriter(modelfile);
        XYSeries model1 = new XYSeries("Incorrect matches");
        XYSeries model2 = new XYSeries("Correct matches");
        XYSeries model3 = new XYSeries("All target hits");
        double x = 0f;
        double y_c = 0f;
        double y_i = 0f;
        PVector point = new PVector(2);
        for (int i = 0; i < NoPoints; i++) {
            point.array[0] = model_kde_x[i];
            x = point.array[0];
            point.array[1] = (mmc.EF.density(point, mmc.param[0]) - model_kde_decoy_y[i]) * mmc.weight[0];
            model1.add(point.array[0], point.array[1]);
            y_i = point.array[1];
            point.array[1] = (mmc.EF.density(point, mmc.param[1]) - model_kde_decoy_y[i]) * mmc.weight[1];
            model2.add(point.array[0], point.array[1]);
            y_c = point.array[1];
            writer.write(x + "\t" + y_i + "\t" + y_c + "\n");
        }
        writer.close();

        KernelDensityEstimator kde = new KernelDensityEstimator();
        kde.SetData(ModelObs);
        double[] model_kde_y = kde.Density(model_kde_x);

        for (int i = 0; i < NoPoints; i++) {
            if (model_kde_x[i] > min && model_kde_x[i] < max) {
                model3.add(model_kde_x[i], model_kde_y[i]);
            }
        }
        MixtureModelProb = new Float[NoPoints + 1][3];
        float positiveaccu = 0f;
        float negativeaccu = 0f;

        MixtureModelProb[0][0] = (float) model2.getMaxX() + Float.MIN_VALUE;
        MixtureModelProb[0][1] = 1f;
        MixtureModelProb[0][2] = 1f;

        for (int i = 1; i < NoPoints + 1; i++) {
            float positiveNumber = model2.getY(NoPoints - i).floatValue();
            float negativeNumber = model1.getY(NoPoints - i).floatValue();
            MixtureModelProb[i][0] = model2.getX(NoPoints - i).floatValue();
            positiveaccu += positiveNumber;
            negativeaccu += negativeNumber;
            MixtureModelProb[i][2] = 0.999999f * positiveNumber / (negativeNumber + positiveNumber);
            MixtureModelProb[i][1] = 0.999999f * positiveaccu / (negativeaccu + positiveaccu);
        }

        XYSeriesCollection dataset = new XYSeriesCollection();
        dataset.addSeries(model1);
        dataset.addSeries(model2);
        dataset.addSeries(model3);

        HistogramDataset histogramDataset = new HistogramDataset();
        histogramDataset.setType(HistogramType.SCALE_AREA_TO_1);
        histogramDataset.addSeries("ID hits", IDObs, 100);
        histogramDataset.addSeries("Decoy hits", DecoyObs, 100);
        //histogramDataset.addSeries("Model hits", ModelObs, 100);

        JFreeChart chart = ChartFactory.createHistogram(FilenameUtils.getBaseName(pngfile), "Score", "Hits", histogramDataset, PlotOrientation.VERTICAL, true, false, false);
        XYPlot plot = chart.getXYPlot();

        NumberAxis domain = (NumberAxis) plot.getDomainAxis();
        domain.setRange(min, max);
        plot.setBackgroundPaint(Color.white);
        plot.setDomainGridlinePaint(Color.white);
        plot.setRangeGridlinePaint(Color.white);
        plot.setForegroundAlpha(0.8f);
        chart.setBackgroundPaint(Color.white);

        XYLineAndShapeRenderer render = new XYLineAndShapeRenderer();
//        render.setSeriesPaint(0, Color.DARK_GRAY);
//        render.setSeriesPaint(1, Color.DARK_GRAY); 
//        render.setSeriesPaint(2, Color.GREEN); 
//        render.setSeriesShape(0, new Ellipse2D.Double(0, 0, 2, 2));
//        render.setSeriesShape(1, new Ellipse2D.Double(0, 0, 2, 2));
//        render.setSeriesShape(2, new Ellipse2D.Double(0, 0, 2.5f, 2.5f));
//        render.setSeriesStroke(1, new BasicStroke(1.0f));
//        render.setSeriesStroke(0, new BasicStroke(1.0f));
//        render.setSeriesStroke(2, new BasicStroke(2.0f));
        plot.setDataset(1, dataset);
        plot.setRenderer(1, render);
        plot.setDatasetRenderingOrder(DatasetRenderingOrder.FORWARD);
        try {
            ChartUtilities.saveChartAsPNG(new File(pngfile), chart, 1000, 600);
        } catch (IOException e) {
        }
    }

    public void AssignMixtureModelProb() {
        if (Terminate) {
            return;
        }

        for (UmpireSpecLibMatch match : libTargetMatches) {
            for (PeakGroupScore peakscore : match.TargetHits) {
                for (int i = 0; i < NoBinPoints - 1; i++) {
                    if (peakscore.UmpireScore >= MixtureModelProb[i][0]) {
                        peakscore.MixtureModelProb = MixtureModelProb[i][1];
                        peakscore.MixtureModelLocalProb = MixtureModelProb[i][2];
                        break;
                    }
                }
            }
            for (PeakGroupScore peakscore : match.DecoyHits) {
                for (int i = 0; i < NoBinPoints - 1; i++) {
                    if (peakscore.UmpireScore >= MixtureModelProb[i][0]) {
                        peakscore.MixtureModelProb = MixtureModelProb[i][1];
                        peakscore.MixtureModelLocalProb = MixtureModelProb[i][2];
                        break;
                    }
                }
            }
            match.AssignProbToPepIon();
        }
    }

    public void AssignProbSVM() {
        if (Terminate) {
            return;
        }

        for (UmpireSpecLibMatch match : libTargetMatches) {
            PepIonID pepIonID = match.pepIonID;
            pepIonID.MS1AlignmentLocalProbability = 0f;
            pepIonID.MS1AlignmentProbability = 0f;
            pepIonID.MS2AlignmentLocalProbability = 0f;
            pepIonID.MS2AlignmentProbability = 0f;

            for (PeakGroupScore peakscore : match.TargetHits) {
                if (peakscore.MSlevel == 1 && peakscore.UmpireScore > pepIonID.MS1AlignmentLocalProbability) {
                    pepIonID.MS1AlignmentLocalProbability = peakscore.MixtureModelLocalProb;
                }
                if (peakscore.MSlevel == 2 && peakscore.UmpireScore > pepIonID.MS2AlignmentLocalProbability) {
                    pepIonID.MS2AlignmentLocalProbability = peakscore.MixtureModelLocalProb;
                }
            }
        }
    }

    public static TargetMatchScoring LibraryMatchRead(String Filename, String LibID) throws FileNotFoundException {

        if (!new File(FilenameUtils.getFullPath(Filename) + FilenameUtils.getBaseName(Filename) + "_" + LibID + "_LibMatch.serFS").exists()) {
            return null;
        }
        TargetMatchScoring match = null;
        try {
            Logger.getRootLogger().info("Loading Target library match results to file:" + FilenameUtils.getFullPath(Filename) + FilenameUtils.getBaseName(Filename) + "_" + LibID + "_LibMatch.serFS...");
            FileInputStream fileIn = new FileInputStream(FilenameUtils.getFullPath(Filename) + FilenameUtils.getBaseName(Filename) + "_" + LibID + "_LibMatch.serFS");
            FSTObjectInput in = new FSTObjectInput(fileIn);
            match = (TargetMatchScoring) in.readObject();
            in.close();
            fileIn.close();
        } catch (Exception ex) {
            Logger.getRootLogger().error(ExceptionUtils.getStackTrace(ex));
            return null;
        }
        return match;
    }

    public static TargetMatchScoring LibraryMatchReadJS(String Filename, String LibID) throws FileNotFoundException {

        if (!new File(FilenameUtils.getFullPath(Filename) + FilenameUtils.getBaseName(Filename) + "_" + LibID + "_LibMatch.serFS").exists()) {
            return null;
        }
        TargetMatchScoring match = null;
        try {
            Logger.getRootLogger().info("Loading Target library match results to file:" + FilenameUtils.getFullPath(Filename) + FilenameUtils.getBaseName(Filename) + "_" + LibID + "_LibMatch.ser...");
            FileInputStream fileIn = new FileInputStream(FilenameUtils.getFullPath(Filename) + FilenameUtils.getBaseName(Filename) + "_" + LibID + "_LibMatch.ser");
            ObjectInputStream in = new ObjectInputStream(fileIn);
            match = (TargetMatchScoring) in.readObject();
            in.close();
            fileIn.close();
        } catch (Exception ex) {
            Logger.getRootLogger().error(ExceptionUtils.getStackTrace(ex));
            return null;
        }
        return match;
    }

    public boolean LibraryMatchWriteJS() throws FileNotFoundException {
        try {
            Logger.getRootLogger().info("Writing Target match results to file:" + FilenameUtils.getFullPath(Filename) + FilenameUtils.getBaseName(Filename) + "_" + LibID + "_LibMatch.ser...");
            FileOutputStream fout2 = new FileOutputStream(FilenameUtils.getFullPath(Filename) + FilenameUtils.getBaseName(Filename) + "_" + LibID + "_LibMatch.ser", false);
            ObjectOutputStream oos = new ObjectOutputStream(fout2);
            oos.writeObject(this);
            oos.close();
            fout2.close();
        } catch (Exception ex) {
            Logger.getRootLogger().error(ExceptionUtils.getStackTrace(ex));
            return false;
        }
        return true;
    }

    public boolean LibraryMatchWrite() throws FileNotFoundException {
        try {
            Logger.getRootLogger().info("Writing Target match results to file:" + FilenameUtils.getFullPath(Filename) + FilenameUtils.getBaseName(Filename) + "_" + LibID + "_LibMatch.serFS...");
            FileOutputStream fout = new FileOutputStream(FilenameUtils.getFullPath(Filename) + FilenameUtils.getBaseName(Filename) + "_" + LibID + "_LibMatch.serFS", false);
            FSTObjectOutput out = new FSTObjectOutput(fout);
            out.writeObject(this);
            out.close();
            fout.close();
        } catch (Exception ex) {
            Logger.getRootLogger().error(ExceptionUtils.getStackTrace(ex));
            return false;
        }
        return true;
    }

}
