## sc-houdini

Create FITS masks

### Help Menu
```
usage: sc-houdini [-h] [-o  [...]] [-above] [-below] [-rb  [...]] [-er] iname

positional arguments:
  iname                 Input image from which to generate mask

options:
  -h, --help            show this help message and exit
  -o  [ ...], --outname  [ ...]
                        Where to dump the output(s). We can iterate over multiple reg files. See '-rb'
  -above , --above      Include everything above this value. ie. > above. Aka the lower limit
  -below , --below      Include everything below this value. i.e < below. Aka the upper limit
  -rb  [ ...], --region-bound  [ ...]
                        DS9 region file(s) within which to make our mask
  -er , --exclude-region 
                        DS9 region file(s) containing the regions that should be excluded
```