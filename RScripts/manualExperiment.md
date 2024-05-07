# Manual 'Experiment'

```bash
Rscript RScripts/initialize.R YYYYMMDD
Rscript RScripts/new_transfer.R 1
Rscript RScripts/mkgwl_setup.R 1

Rscript RScripts/new_transfer.R 1
Rscript RScripts/mkgwl_transfer.R 1 1
Rscript RScripts/make_barcodes.R 1 1 1
Rscript RScripts/generate_bc_files_transfer.R 1 1 1
```

