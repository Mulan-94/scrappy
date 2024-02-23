Some sub-band images generated during multi-frequency imaging may become irrevocably unusable because of corruption  by RFI or the PSF. Consequently, such images are excluded from spectropolarimetric analysis. Since `WSClean` assigns a weight to each of the sub-band images -- where the better ones have more weight -- such images end up with lower weights. Sc-beam-plot excludes sub-band images whose weight is below a certain threshold. When it is impossible to reconstruct an image, its beam size may also be set to zero; such images are also excluded.

As a diagnostic, sc-beam-plot generates a plot showing beam dimensions read from the `BMAJ & BMIN` FITS header keywords before and after "valid" images are selected; these are stored in `00-beam_vs_freq.png` and `selected_beam_vs_freq.png`, respectively. This way, strange beam variations can easily be identified and further action taken. This tool also generates other `*.txt` outputs that are used within other tools in Scrappy. These are described in the next section.


### Outputs

| Output                            | Description                           |
|-----------------------------------|---------------------------------------|
| `00-beams.npz`                    | Contains beam dimensions in radians as a NumPy pickle file. |
| `00-beam_vs_freq.png`             | Figure with all the images' (selected + unselected) beam dimensions. |
| `beam-dims.txt`                   | Contains the LARGEST beam size found amongst the selected images. |
| `frequencies.txt`                 | Frequencies of the selected sub-band images. |
| `not-selected-channels.txt`       | Contains channel numbers (`WSClean` syntax) of the excluded images. |
| `orig-wsums.txt`                  | Contains weights of all the input images. |
| `selected_beam_vs_freq.png`       | Figure showing the excluded channels. |
| `selected-channels.txt`           | Contains channel numbers of the selected sub-band images. |
| `selected-freq-images.txt`        | Contains names of the selected images. |
| `wsums.txt`                       | Weights of the selected images. |



### Help Menu
```
usage: sc-beam-plot [-h] [-o] [-search SEARCH] [-t] [-sc SEL] [-as] [-pre PREFIX] [-ec EXCLUDE]

positional arguments:
                        Directory where to find the input images

options:
  -h, --help            show this help message and exit
  -o , --output         Directory where to dump outputs. Default is current directory.
  -search SEARCH, --search SEARCH
                        Regex of what string to search in the directory
  -t , --threshold      Channels with WSUM/WSUM.max() below this value will be excluded. Only used if -s is active
  -sc SEL, --select-channels SEL
                        Custom channel selections. This can be a file containing those channels, or a comma separated list e.g. '0,1,2,3'
  -as, --auto-select    Try and suggest valid channels for selection. Use in conjuction with '-t'.
  -pre PREFIX, --prefix PREFIX
                        Prefix to append on the output images
  -ec EXCLUDE, --exclude-channels EXCLUDE
                        Channels to be excluded during autoselection. This can be a file containing those channels, or a comma separated list e.g. '0,1,2,3'.
```