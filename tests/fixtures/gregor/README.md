# Cohort Allele Frequency Generation

## Testing

To run integrations tests in `tests/integration/gregor`:

1. **Configure Google Cloud token**: create token to access remote files (eg remote VCF)
```bash
export GCS_OAUTH_TOKEN=$(gcloud auth application-default print-access-token)
```
2. **Locate phenotype table**: Download a GREGoR formatted phenotype table as tsv from AnVIL platform, then either
   1. Copy that file to `tests/gregor/fixtures/` as  `phenotypes.tsv` or `phenotypes.csv` or...
   2. `export PHENOTYPE_TABLE=<ABSOLUTE_PATH_TO_PHENO_TABLE>`
3. **Locate fixture**: Request access to the [VCF index](https://ohsuitg-my.sharepoint.com/my?id=%2Fpersonal%2Fwongq%5Fohsu%5Fedu%2FDocuments%2Fgregor%5Fcaf%5Fintegration%5Ftest) from collaborators, then either
   1. Copy that file to `tests/gregor/fixtures/chr3_chrY_index.db` or...
   2. `export VRS_VCF_INDEX=<ABSOLUTE_PATH_TO_VRS_VCF_INDEX>`
