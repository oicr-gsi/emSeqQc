# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]
### Added
- BAM file merging doesn't happen if only a single BAM file is produced.

### Changed
- Replace bamqc with smallBamQc.

## [Unreleased] - 2024-06-25
### Added
-[GRD-797](https://jira.oicr.on.ca/browse/GRD-797) - add vidarr labels to outputs (changes to medata only).

## [1.0.0]
### Added
- EM-Seq pipeline that takes in multiple FastQ pairs and outputs bamqc, samtools stats, and MethylDackel files.
