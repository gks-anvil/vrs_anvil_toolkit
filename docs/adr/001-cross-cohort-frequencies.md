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

  * Coverage file (BED), with reference assembly and provenance for any coverage estimation
  * Variant index
    * VRS ID/allele data
    * Sufficient data to tabix back to source VCF: a file URI + original, unnormalized VCF `chr` and `pos` values

2. Prepare cross-cohort artifacts

  * Unified coverage file, consisting of the intersection of all cohorts' coverage (BED)
  * Deduplicated sample list, produced via an all-vs-all comparison of samples within regions of shared coverage
    * bloom filter first, exhaustive comparison second
  * (Probably) unified variant index, combining indices from individual cohorts
  * Reproducible cohort definition recording source datasets, index versions, participant/sample selections, filters, and normalization/liftover, coverage, and phenotype-matching policies

3. For an allele of interest:

  * First, get the VRS ID for the allele, either via some search/translation function or by already knowing it
  * Then, get the dereferencing params (file URI, VCF CHR, VCF POS) for each cohort
  * For each cohort, use those params to slice back into the original VCFs with tabix, returning sample counts
  * Apply the cohort/subcohort sample selection and coverage policy to calculate coverage-adjusted allele numbers; exclude samples without adequate coverage rather than treating missing records as reference-agree calls
  * Add up counts, structure as VA-Spec CAF, and return, retaining per-study results and provenance alongside combined results

## Consequences

**Synthetic coverage calculation**: BED files aren't always available. One conservative approach could be to produce an artificial BED file based on called variants (+/- some # of BPs?). This remains a heuristic, not proof of callability; its method and window parameters need to be retained. The AnVIL proposal also calls for estimating coverage from alignment data when supplied profiles are unavailable. Variants and coverage ranges need a common reference assembly (GRCh38 in the proposal), with any liftover policy recorded.

**Sample deduplication**: need to work out bloom filter things, equality thresholds, etc. Signatures should use sets of called VRS variants to assess similarity both within and across studies without collecting unnecessary personal health information. A similarity flag is not confirmed participant identity; retain methods, thresholds, and review/exclusion decisions.

**Define deduplicated sample list format**: Probably a simple line delimited text file with sample names is good enough. The accompanying mappings must distinguish study-specific participant records from samples and retain study/dataset provenance so multiple samples do not silently become multiple participants.

**Update library interfaces**: Need to update function signatures, CLI etc to handle multiple cohorts/multiple input sources

**Metadata and evidence interfaces**: Study plugins should map Terra tables, JSON, or delimited metadata into a shared internal representation. Phenotype-defined subcohorts need explicit matching policies, including exact or HPO hierarchy-aware matching and vocabulary mappings where needed. CAF construction and schema validation should also accept aggregate genotype/phenotype counts from other tools, independently of this indexing workflow.
