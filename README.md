**DIA-Umpire**: 

Computational analysis for mass spectrometry-based proteomics data

The followings are the brief descriptions for each NetBeans project

**DIA-Umpire_SE**

The first step of DIA-Umpire analysis. It is a signal extraction module to generate pseudo MS/MS spectra given a DIA file

**DIA-Umpire_Quant:** 

DIA-Umpire quantification module which performs the following steps for a list of DIA files:

1. Generate untargeted peptide IDs
2. Internal library search
3. External library search
4. Generate a master protein list given an FDR threshold and map peptides to the master protein list and do final protein-level quantification

Alternatively, one run the four modules separately using the following Netbeans projects

**DIA_Umpire_LCMSIDGen**

After untargeted identification analysis, DIA-Umpire generate untargeted peptide IDs using pepXML files as *LCMSID file, which is the input file format for internal and external library search.

**DIA_Umpire_IntLibSearch**

DIA-Umpire internal library search using *LCMSID file as input file

**DIA_Umpire_ExtLibSearch**

DIA-Umpire external library search using *LCMSID file and a spectral library as input files

**DIA_Umpire_ProtQuant**

DIA-Umpire protein quantification module, the output is csv format


**DIA-Umpire_To_Skyline**

Module to generate raw-intensity pseudo MS/MS spectra. (without intensity adjustments) 

**DIA-Umpire**

Main DIA-Umpire class libraries

* ExternalPackages: external packages, currently including JAligner, SortedListLib, JMEF, and a traML parser developed by ISB.
* FDREstimator: wrapper to generate FDR filtered protein and peptide lists
* MSUmpire: Umpire libraries
  * BaseDataStructure : basic data structure, including x-y pair value data, scan data, and parameter setting class etc.
  * DIA : DIA specific classes
  * FragmentLib : Fragment library manager, basically it's a spectral library
  * LCMSPeakStructure : Data structure classes for MS1 or MS2 peak for to a LCMS run
  * MathPackage : Math calculation classes
  * PSMDataStructure : Data structure classes related to Mass-spec based identifications, ranging from PSM, peptide ion, protein and identification data structure for a LC-MS run (LCMSID.java). In addition, this package includes processing manager to extract PTM and peptide fragment information generated from Compomics library. 
  * PeakDataStructure : Data structure classes related to peak data, from peak curve, peak isotope cluster, and peak smoothing algorithms.
  * PeptidePeakClusterDetection : Processing classes to detect peak features 
  * SearchResultParser : PepXML, ProtXML parsers
  * SeqUtility : Classes for sequence processing, including FastaParser and shuffled sequence generator
  * SpectralProcessingModule : Spectrum peak processing classes
  * SpectrumParser: mzXML, mzML parsers
  * Utility: Other classes
  
* resource: all resource files in text format


**MS1Quant**

DDA-based MS1 quantification tool based on the feature detection algorithm. 
