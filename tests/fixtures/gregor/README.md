# Cohort Allele Frequency Generation

## Testing

To run integrations tests in `tests/integration/gregor`:

1. **Configure Google Cloud token**: create token to access remote files (eg remote VCF)
```bash
export GCS_OAUTH_TOKEN=$(gcloud auth application-default print-access-token)
```
2. **Locate phenotype table**: Download a GREGoR formatted phenotype table as tsv from AnVIL platform
   1. Then, either copy that file to `tests/gregor/fixtures/` as  `phenotypes.tsv` or `phenotypes.csv` or...
   2. `export PHENOTYPE_TABLE=<ABSOLUTE_PATH_TO_PHENO_TABLE>`
3. **Locate fixture**: Get VCF index from collaborators and either
   1. Copy that file to `tests/gregor/fixtures/chr3_chrY_index.db`
   2. `export VRS_VCF_INDEX=<ABSOLUTE_PATH_TO_VRS_VCF_INDEX>`
