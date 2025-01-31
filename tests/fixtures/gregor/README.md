# Cohort Allele Frequency Generation

## Testing

To run integrations tests in `tests/integration/gregor`:

1. **Configure Google Cloud token**: create token to access remote files (eg remote VCF)
```bash
export GCS_OAUTH_TOKEN=$(gcloud auth application-default print-access-token)
```
2. **Locate phenotype table**: There are two ways to do this depending on your testing environment
   1. If working in a Terra workspace,
      1. Import the phenotype table from the [GREGoR U08 Workspace](https://app.terra.bio/#workspaces/gregor-dcc/GREGOR_COMBINED_CONSORTIUM_U08/data) in your workspace.
      2. Then, the tests will automatically look for the table titled "phenotype" and load it.
   2. If not in Terra,
      1. Download the phenotype table from the [GREGoR U08 Workspace](https://app.terra.bio/#workspaces/gregor-dcc/GREGOR_COMBINED_CONSORTIUM_U08/data) as a tsv from AnVIL platform
      2. Ensure the file name is `phenotype.tsv`
      3. Locate the file by either
         1. Moving that file to `tests/fixtures/gregor` or...
         2. `export PHENOTYPE_TABLE=<ABSOLUTE_PATH_TO_PHENO_TABLE>`
3. **Locate fixture**: Request access to the [VCF index](https://ohsuitg-my.sharepoint.com/my?id=%2Fpersonal%2Fwongq%5Fohsu%5Fedu%2FDocuments%2Fgregor%5Fcaf%5Fintegration%5Ftest) from collaborators, then either
   1. Copy that file to `tests/fixtures/gregor/chr3_chrY_index.db` or...
   2. `export VRS_VCF_INDEX=<ABSOLUTE_PATH_TO_VRS_VCF_INDEX>`
