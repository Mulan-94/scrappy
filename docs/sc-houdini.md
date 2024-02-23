### Introduction

This tool creates a simple FITS mask given an input image and thresholds containing the values that should not be masked out. `Sc-houdini` also accepts region files specifying a region which its efforts should be concentrated within (see [DS9 region files](https://ds9.si.edu/doc/ref/region.html)). 


<!-- ------
### Example -->


------
### Outputs
Sc-houdini outputs:

- A FITS image mask
- A corresponding `*.png` file showing how the mask appears for quick mask inspection (*this is sometimes buggy and is being investigated*).

With a typical run using [showrunner][output-directory-structure], the output is stored in the `masks` directory and will appear like:

<pre><code>
├── masks                 # Location of the generated FITS masks
|   ├── <b style="color:red">*.fits </b>
|   ├── <b style="color:red">*.png</b>
</code></pre>


------
### Help Menu
The help menu is available using the command:

```
sc-houdini -h
```

which results in the output below. However, as the tool is under active development, the help menu may change in the future and should be **checked before use using the help command** after installation.

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