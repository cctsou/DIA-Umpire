**DIA-Umpire**: 

Computational analysis for mass spectrometry-based proteomics data

**DIA-Umpire_Quant:** 

DIA-Umpire quantification module which performs the following steps for a list of DIA files:

1. Generate a master protein list given an FDR threshold (using prot.xml file)
2. Generate untargeted peptide IDs
3. Internal library search
4. External library search
5. Map peptides to the master protein list and do final protein-evel quantification

**DIA-Umpire_SE**

DIA-Umpire signal extraction module to generate pseudo MS/MS spectra given a DIA file

**DIA-Umpire_To_Skyline**

Module to generate raw-intensity pseudo MS/MS spectra. (without intensity adjustments) 

**DIA-Umpire**

Main DIA-Umpire class libraries

* ExternalPackages: external packages, currently including JAligner, SortedListLib, JMEF, and a traML parser developed by ISB.
* FDREstimator: wrapper to generate FDR filtered protein and peptide list
* MSUmpire: Umpire libraries
* resource: all resource files in txet format

**MS1Quant**

DDA-based MS1 quantification tool based on the feature detection algorithm. 
