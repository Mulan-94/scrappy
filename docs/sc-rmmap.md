### Introduction

<!-- ------
### Example -->


------
### Output



------
### Help Menu
```
usage: sc-rmmap [-h] [-q QFITS] [-u UFITS] [-i IFITS] [-f FREQ] [-ncore NUMPROCESSOR] [-mask MASKFITS] [-o PREFIX] [-niters NITERS] [-md MAX_DEPTH] [--depth-step DEPTH_STEP] [-debug] [-snr SNR]
                [-nd]

Performs linear least squares fitting to Q and U image.

options:
  -h, --help            show this help message and exit
  -q QFITS, --qcube QFITS
                        Stokes Q cube (fits)
  -u UFITS, --ucube UFITS
                        Stokes U cube (fits)
  -i IFITS, --icube IFITS
                        Stokes I cube (fits)
  -f FREQ, --freq FREQ  Frequency file (text)
  -ncore NUMPROCESSOR, --ncore NUMPROCESSOR
                        number of cores to use. Default 60.
  -mask MASKFITS, --maskfits MASKFITS
                        A mask image (fits)
  -o PREFIX, --prefix PREFIX
                        This is a prefix for output files.
  -niters NITERS, --niters NITERS
                        Number of clean iterations. Default 1000
  -md MAX_DEPTH, --max-depth MAX_DEPTH
                        Maximum Farady depth to fit for. Default 200
  --depth-step DEPTH_STEP
                        Depth stepping. Default 1
  -debug, --debug       Enable debug mode, disable parallel processing
  -snr SNR, --snr-threshold SNR
                        Threshold to mask out data. Default is 10
  -nd, --no-derotate    Use this switch to NOT derotate the RMTF by the mean lambda squared.
```