
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