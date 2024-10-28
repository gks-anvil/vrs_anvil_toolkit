# Cohort Allele Frequency Generation

## Testing

To run integrations tests in `tests/integration/gregor`:

```bash
export GCS_OAUTH_TOKEN=$(gcloud auth application-default print-access-token)

# TOOD: download phenotype table as tsv from AnVIL data on platform
export PHENOTYPE_TABLE=<PATH_TO_DOWNLOADED_PHENOTYPE_TABLE>

# TODO: get VCF index from collaborator
export VRS_VCF_INDEX=<PATH_TO_VRS_VCF_INDEX>
```
