

Autogenerate lines-of-sight and their corresponding data.


### Help Menu
```
usage: sc-los [-h] [-nri NOISE_REF] [-todo] [--noverwrite] [-o] [-t] [-j] [--debug] [-ro] [-mrn] [-rf] [-rs] [-m] [-lo] [-idir IMAGE_DIR | -cubes CUBES CUBES CUBES] [-freqs] [-ref-image WCS_REF]
              [-nrf NRFILE] [-mft] [-psnr] [-po] [--plot-grid] [-piqu] [-pfp] [-plp] [--ymax] [--ymin] [--xmax] [--xmin] [-p [...]] [-ps [...]]

________________________________________________________________________________

 Description

 Generate I, Q and U data for various LoS from image cubes.

 This script uses the total intesity MFS image to determine regions
 of the source with enough SNR and selects those as the lines of sight.
 For now, this is tuned specifically for Pictor A, but it may be
 extended and revised in the future.

 The following could be generated from this script:

 . 1. Region file containing the selected lines of sight
 . 2. Plot showing locations of the regions on the source
 . 3. Pickle files containing the following keys for the data:

 .   - I      - Q     - U
 .   - i_mfs_noise    - q_mfs_noise   - u_mfs_noise
 .   - I_err    - Q_err   - U_err
 .   - lpol     - lpol_err
 .   - fpol     - fpol_err
 .   - pangle   - pangle_err
 .   - mask     - freqs
 . i.e for a single line of sight. Each LoS contains data for all
 . the available channels
 .
 . 4. Optionaly, generate plots for fractional polarisation vs
 lambda squared for each LoS

________________________________________________________________________________

options:
  -h, --help            show this help message and exit
  -todo , --todo        A string containing to do items. Specify using: (r): generate regions, (l): generate LOS data, (p): generate plots. Default is 'rl'

Options:
  -nri NOISE_REF, --noise-ref-image NOISE_REF
                        The total intensity image used to get the global noise reference.

General arguments:
  --noverwrite          Do not ovewrite everything along the way. Default is overwrite
  -o , -odir , --output-dir 
                        where to dump output
  -t , --threshold      If SNR below which data will not be considered. Default is 3
  -j , --nworkers       How many workers to use for processing
  --debug               Disble parallel processing and enables sequential mode.

Region generation arguments:
  -ro, --regions-only   Only generate the region files. This requires the following options:   --ref-image --mask --region-size 
  -mrn , --minimum-region-noise 
                        Specific noise floor to generate the regions.
  -rf , --region-file   An input region file. Otherwise, one will be auto-generated. Genereated regions will be stored here
  -rs , --region-size   Create regions of this circle radius and perform analyses on them. If you want to set the data threshold, please use --threshold.
  -m , --mask           Mask containing the area where LoS should be restricted. This can be a FITS file or a region file (*.reg). It's REQUIRED for automatically making regions.

LoS Data generation arguments:
  -lo, --los-only       Only generate the line of sight data files. The following options should be specified:   --idir/--cubes --ref-image --mask/--region-file --noise-ref-image --noise-
                        ref-file 
  -idir IMAGE_DIR, --image_dir IMAGE_DIR
                        Where the channelised I, Q and U images are
  -cubes CUBES CUBES CUBES, --cubes CUBES CUBES CUBES
                        The I, Q, U image cubes (in this specific order) to be used. This will require specification of --freq-file
  -freqs , --freq-file 
                        Text file containing frequencies to be used. This is only active when the input FITS images are cubes and is particularly useful when the frequencies of images that form the
                        cube do not increase monotonically.
  -ref-image WCS_REF, --ref-image WCS_REF
                        The reference image that will be used to generate the default region file. Must be the stokes I MFS image. This image will also be used to get the reference WCS for region
                        file generation.
  -nrf NRFILE, --noise-region-file NRFILE
                        A region file containing region to be used as the noise reference.
  -mft , --minimum-flag-threshold 
                        Fraction of flags above which lines of sight should be ignored. Can be useful if you want to plot all the generated LOS. Otherwise, they will be filtered out. The simple
                        filter is that where values of fractional polarisation i are >1 or <0, this data is 'flagged'. behaviour: flag size > 0.7 of total data is flagged, ingore. Max is 1, min is
                        >0. Default 0.7
  -psnr, --use-polzd-snr
                        Use to elect use of polarised SNR to determine valid LoS. Default S/N used is total intensity: rms noise

Plotting Arguments:
  -po, --plots-only     Only do plots
  --plot-grid           Enable to make gridded plots
  -piqu, --plot-iqu     Plot Q and U values
  -pfp, --plot-frac-pol
                        Plot Fractional polarization
  -plp, --plot-linear-pol
                        Plot linear polarization power
  --ymax                Y axis max limit
  --ymin                Y axis min limit
  --xmax                Y axis max limit
  --xmin                Y axis min limit
  -p [ ...], --plot [ ...]
                        Make plots for these region sizes manually. These will be linearly scaled
  -ps [ ...], --plot-scales [ ...]
                        Scales for the plots. Can be a space separated list of different scales. Options are linear or log.
```