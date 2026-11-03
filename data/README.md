# Prepared lesson inputs

- `targets.rds`: metadata for the six LNCaP/PREC samples in the introduction.
- `rgSet.rds`: their raw EPICv2 red/green intensities, before QC or normalisation.

Source: GEO GSE240469, samples GSM7698438, GSM7698446, GSM7698462,
GSM7698435, GSM7698443 and GSM7698459 (500 ng DNA input).

The introduction executes its preparation chunks during rendering: it calls
GEOquery to download the IDAT files, imports the data, and saves both RDS files
here. The next episode loads the saved objects independently.

Keep both RDS files available for builds that render only the next episode,
including in the GitHub checkout used for deployment. To regenerate manually,
run the introduction's sample-table, download-data, read-idats and
save-input-objects chunks in order with `episodes/` as the working directory.
Replace both RDS files together when changing samples.

The saved Basename/filenames metadata records the original import paths. Those
paths are not required for reading the saved intensity object on another machine.
