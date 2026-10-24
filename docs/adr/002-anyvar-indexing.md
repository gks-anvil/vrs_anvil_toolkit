# ADR 002: Use AnyVar as indexing backend

## Status

*Proposed*

## Context

In the previous GREGoR work, we developed a variant index system to enable cohort allele frequency calculation without downloading the entirety of original VCFs. The solution involved pulling out VRS alleles from a VCF, storing them into a sqlite file, and packaging that together with a URI for the source file(s). Then, to calculate frequency, we get a VRS ID for the allele of interest, look it up in the index, and then use the information contained there to slice (ie via tabix) into the source file, pull out sample columns etc, maybe do further filtering based on patient-associated phenotypes, and return the final result.

There are some deficiencies here

* The index producer and consumer libs are totally separate, there's no shared schema
* All of this tooling is fairly out-in-the-wind; it adds a new corner pile of things for us to maintain
* Sqlite is simple and fine but doesn't scale incredibly well; there's also probably not a need for an OLTP db given that we are basically only doing writes up-front, well before any reads need to happen. It is nice that it's a file, though.

So, we'd like to consider alternatives. For the index itself, we want to support at least the following two main query modes:

1) Individual allele lookup by VRS ID, joined w/ associated source URI(s)
2) Genomic range-based lookup (i.e. retrieve all alleles with coordinates between `(start, end)`)



## Decision [proposed]

Use [AnyVar](https://github.com/biocommons/anyvar/) to ingest variants. Use the already-supported duckdb backend for the index.

**Benefits**

* Enables workflow consolidation -- AnyVar can annotate the dataset VCFs *and* build an index simultaneously
* Maintenance effort -- employs a project that we're working on for other reasons, so some maintenance/bugfixes/etc come for free and we have other eyes on library health
* Opens up some other interoperability features; there are some planned data sharing/export options that help bolster the case for employing the index
* Duckdb is better-suited for our use-case than sqlite -- it's file-based, which plays nicely with AnVIL, and has some better optimization in terms of storage size and reads for large genomic datasets. I wouldn't claim it is the single most optimal storage solution, but it's pretty versatile, very established, and the duckdb backing for AnyVar is largely written already.

**Possible drawbacks**

* AnyVar is an in-progress application; we're at a 1.x release, but further schema alterations are possible, which would make it complicated to re-use older indices in the future
  * Some kind of migration scripting might work (could possibly even build this into AnyVar directly if necessary)
* A critical part of the workflow is packaging in some kind of URI + file positioning to point back to the storage location of the original VCF to get sample counts for the allele. AnyVar doesn't have a particular efficient way of annotating this; we could include URI as an extension on every single allele, but that would represent a considerable size increase for a very simple (and almost always uniform or near-uniform) annotation.
  * Not unthinkable that we could build a solution into AnyVar for this, if it makes sense for other corners of the AnyVar userbase too
  * We could also inject this into the anyvar table within the vrs-anvil toolkit; it would still enable VRS-AnVIL Toolkit-authored indices to be used in other AnyVar contexts but would limit utility in the other direction


## Consequences

* Need to do some building in the vrs-anvil toolkit itself, both to consolidate index construction and
* Update existing annotator/index workflows -- we can now do this in a single pass
