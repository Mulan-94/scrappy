## sc-ifit

Fit a simple taylor polynomial on Stokes I cube, generating a model I


### Help Menu
```
usage: sc-ifit [-h] [-mask MASK] [-ex EXCLUDE] [-u] [-deg DEG] [-o ODIR] cube freqs

Fit a simple taylor polynomial on Stokes I cube, generating a model I

positional arguments:
  cube                  Name of the input I cube to be modelled
  freqs                 Text file containing respective frequencies. The freqs should be in Hz

options:
  -h, --help            show this help message and exit
  -mask MASK, --mask-name MASK
                        Region within which to perform fits. This is a FITS MASK
  -ex EXCLUDE, --exclude-channels EXCLUDE
                        Channels that should be excluded while fitting the image cube. This should be a text file containing those channel numbers. Data corresponding to the specified channels will
                        be set to zero.
  -u, --unstack         Unstack the output model cube as single channelised images. If -ex was specified those channels will also be excluded.
  -deg DEG, --degree DEG
                        Degree of the taylor polynimial to be fitted.
  -o ODIR, --output-dir ODIR
                        Directory where to dump the output files.
```