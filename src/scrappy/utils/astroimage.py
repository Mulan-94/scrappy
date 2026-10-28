from astropy.io import fits
import typer
import numpy as np

class AstroImage:
    _image_id = 0

    def __init__(self, fname):
        self.fname = fname
        self._is_cube = None
        self._hdr = self._get_header()

        if not self._is_cube:
            AstroImage._image_id += 1

    def __str__(self):
        return f"A set of accesibilties for FITS file {self.fname}"

    def __repr__(self):
        return f"AstroImage(fname={self.fname!r})"

    @property
    def is_cube_file(self):
        return self._is_cube

    @property
    def is_single_file(self):
        return not self._is_cube

    @property
    def command(self):
        """Return wsclean command that generated image"""
        return "".join(str(self._hdr["HISTORY"]).split("\n"))

    @property
    def data(self):
        _data = fits.getdata(self.fname).squeeze()
        return _data

    @property
    def frequency(self):
        """The central frequency of the image"""
        return self._hdr["CRVAL3"]

    @property
    def header(self):
        return self._hdr

    def _get_header(self):
        hdr = fits.getheader(self.fname)
        if hdr["NAXIS3"] >1:
            self._is_cube = True
        return hdr

    @property
    def wsum(self):
        return self._hdr["WSCVWSUM"]


    def get_beam_dims(self) -> dict:
        """
        Returns
        -------
        dims: dict[int: dict[key: val]]
            Of the format {chan: {x|y|z: }}
        """
        hkeys = "BMAJ BMIN BPA".split()

        dims = dict(chan=[], bmaj=[], bmin=[], bpa=[])
        if self.is_single_file:
            dims["chan"] = AstroImage._image_id
            dims["bmaj"] = self._hdr["BMAJ"]
            dims["bmin"] = self._hdr["BMIN"]
            dims["bpa"] = self._hdr["BPA"]
        else:
            nchans = self._hdr['NAXIS3']
            typer.echo(f"Found {nchans} frequencies in file")
            chans = range(1, nchans+1)

            dims = dict()
            for chan_id in chans:
                dims["chan"].append(chan_id)
                dims["bmaj"].append(self._hdr[f"BMAJ{chan_id}"])
                dims["bmin"].append(self._hdr[f"BMIN{chan_id}"])
                dims["bpa"].append(self._hdr[f"BPA{chan_id}"])

        return dims
