
# Change Log
All notable changes to this project will be documented in this file.
 
The format is based on [Keep a Changelog](http://keepachangelog.com/)
and this project adheres to [Semantic Versioning](http://semver.org/).

## [Unreleased] - 2023-06-16
 - Added multi-sample MultiQC wf to run after stats files are available.

## [Unreleased] - 2023-06-09
 - Working on the Qualimap task, docker is not working, it needs more tests.
 
### Added
- [Trello ticket](https://trello.com/c/akf50Qlf/35-qc-pipeline-for-wgs): Original trello ticket
- [Fastp](https://github.com/OpenGene/fastp): A tool designed to provide fast all-in-one preprocessing for FastQ files. This tool is developed in C++ with multithreading supported to afford high performance.
- [Picard's CollectMutipleMetrics](https://gatk.broadinstitute.org/hc/en-us/articles/360037594031-CollectMultipleMetrics-Picard-): Collect multiple classes of metrics.This 'meta-metrics' tool runs one or more of the metrics collection modules at the same time to cut down on the time spent reading in data from input files. Available modules include CollectAlignmentSummaryMetrics, CollectInsertSizeMetrics, QualityScoreDistribution, MeanQualityByCycle, CollectBaseDistributionByCycle, CollectGcBiasMetrics, RnaSeqMetrics, CollectSequencingArtifactMetrics and CollectQualityYieldMetrics. The tool produces outputs of '.pdf' and '.txt' files for each module, except for the CollectAlignmentSummaryMetrics module, which outputs only a '.txt' file. Output files are named by specifying a base name (without any file extensions).
- [Picard's CollectWgsMetrics](https://gatk.broadinstitute.org/hc/en-us/articles/360037269351-CollectWgsMetrics-Picard-): Collect metrics about coverage and performance of whole genome sequencing (WGS) experiments. This tool collects metrics about the fractions of reads that pass base- and mapping-quality filters as well as coverage (read-depth) levels for WGS analyses. Both minimum base- and mapping-quality values as well as the maximum read depths (coverage cap) are user-defined.
- [MultiQC](https://multiqc.info/): MultiQC searches a given directory for analysis logs and compiles a HTML report. It's a general-use tool, perfect for summarising the output from numerous bioinformatics tools.
 
### Changed
 
### Fixed
 
## [0.0.1] - aaaa-mm-dd
   
### Added
 
### Changed
  
### Fixed
