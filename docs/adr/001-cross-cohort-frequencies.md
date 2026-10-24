# ADR 001: Calculating allele frequencies across multiple cohorts

## Status

*Proposed*

## Context

We want to expand existing cohort allele frequency calculation tooling in the VRS-AnVIL toolkit to produce accurate cross-cohort counts. A few important considerations

* Different experiments/cohorts profile different genomic regions; "whole genome" isn't always exactly the same. Absent something like a gVCF, we want to be careful about calling an empty row in a VCF a ref-agree.
* Sample deduplication: double-counts are really bad in the context of small numbers (ie rare disease). Patients might get sampled twice across studies, and we want to be able to spot those cases if they happen.
* Be thoughtful about abstractions; different studies will package data in different ways, so we need to be smart about defining interface structures so that it's easy for us to develop "plugins" for each study. We've done this work already for GREGoR but expanding to other cohorts may be instructive and require some changes.

## Decision [proposed]

Workflow looks like this

1. Get individual cohort artifacts

  * Coverage file (BED)
  * Variant index
    * VRS ID/allele data
    * Sufficient data to tabix back to source VCF: a file URI + VCF `chr` and `pos` values

2. Prepare cross-cohort artifacts

  * Unified coverage file, consisting of the intersection of all cohorts' coverage (BED)
  * Deduplicated sample list, produced via an all-vs-all comparison of samples within regions of shared coverage
    * bloom filter first, exhaustive comparison second
  * (Probably) unified variant index, combining indices from individual cohorts

3. For an allele of interest:

  * First, get the VRS ID for the allele, either via some search/translation function or by already knowing it
  * Then, get the dereferencing params (file URI, VCF CHR, VCF POS) for each cohort
  * For each cohort, use those params to slice back into the original VCFs with tabix, returning sample counts
  * Add up counts, structure as VA-Spec CAF, and return


## Consequences

**Synthetic coverage calculation**: BED files aren't always available. One conservative approach could be to produce an artificial BED file based on called variants (+/- some # of BPs?).

**Sample deduplication**: need to work out bloom filter things, equality thresholds, etc

**Define deduplicated sample list format**: Probably a simple line delimited text file with sample names is good enough

**Update library interfaces**: Need to update function signatures, CLI etc to handle multiple cohorts/multiple input sources
