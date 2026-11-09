from copy import deepcopy

import typer
import xarrayfits as xf
from astropy.io import fits
from astropy.wcs import WCS

from scrappy.logging.config import configure_logger

snitch = configure_logger()


class AstroImage:
    _image_id = 0

    def __init__(self, fname, chunked=False, chunk_size=1024):
        snitch.info(f"Loading image: {fname.resolve()}")
        self.fname = fname
        self.chunked = chunked
        self.chunk_size = chunk_size
        self._is_cube = None
        self._hdr = self._get_header()
        self._data = None
        self._freq = None

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
        if self._data is None:
            if self.chunked:
                (self._data,) = xf.xds_from_fits(
                    str(self.fname), chunks={2: self.chunk_size, 3: self.chunk_size}
                )
                self._data = self._data.hdu
            else:
                self._data = fits.getdata(self.fname).squeeze()
        return self._data

    @data.setter
    def data(self, new_data=None):
        self._data = new_data

    @property
    def frequency(self):
        """The central frequency of the image"""
        if self._freq is None:
            self._freq = self._hdr["CRVAL3"]
        return self._freq

    @frequency.setter
    def frequency(self, value):
        self._freq = value

    @property
    def header(self):
        return self._hdr

    def _get_header(self):
        hdr = fits.getheader(self.fname)
        if "NAXIS3" in hdr.keys() and hdr["NAXIS3"] > 1:
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

        dims = {"chan": [], "bmaj": [], "bmin": [], "bpa": []}
        if self.is_single_file:
            dims["chan"] = 1
            dims["bmaj"] = self._hdr["BMAJ"]
            dims["bmin"] = self._hdr["BMIN"]
            dims["bpa"] = self._hdr["BPA"]
        else:
            nchans = self._hdr["NAXIS3"]
            typer.echo(f"Found {nchans} frequencies in file")
            chans = range(1, nchans + 1)

            dims = {}
            for chan_id in chans:
                dims["chan"].append(chan_id)
                dims["bmaj"].append(self._hdr[f"BMAJ{chan_id}"])
                dims["bmin"].append(self._hdr[f"BMIN{chan_id}"])
                dims["bpa"].append(self._hdr[f"BPA{chan_id}"])

        return dims

    def get_clean_header(self):
        """
        A header with only the important stuff
        """
        nhdr = deepcopy(self._hdr)
        nhdr.remove("HISTORY", remove_all=True)
        deletes = []
        for word in "wsc bma bmi bpa".upper().split():
            for key in nhdr.keys():
                if word in key:
                    deletes.append(key)

        for _ in deletes:
            snitch.debug(f"Deleting {_} key from header")
            nhdr.remove(_, remove_all=True)
        return nhdr

    def get_clean_2d_header(self):
        """
        Only the important stuff but for 2-D fits images
        """
        nhdr = self.get_clean_header()

        main_keys = [
            "SIMPLE",
            "BITPIX",
            "NAXIS",
            "NAXIS1",
            "NAXIS2",
            # 'EXTEND', 'BUNIT','BMAJ','BMIN','BPA',
            "EQUINOX",
            "BTYPE",
            "TELESCOP",
            "OBSERVER",
            "OBJECT",
            "SPECSYS",
            "DATE-OBS",
            # 'ORIGIN', "WSCDATAC", "WSCVDATE","WSCVERSI"
        ]
        extra = []

        naxes = 3
        for i in range(1, 4):
            for key in ["CTYPE", "CRPIX", "CRVAL", "CDELT", "CUNIT"]:
                if f"{key}{i}" in "CRVAL3 CDELT3".split():
                    continue
                extra.append(f"{key}{i}")
        main_keys.extend(extra)

        main_keys = {key: nhdr.get(key) for key in main_keys}
        main_keys.update({"NAXIS": naxes, "CUNIT3": None, "CTYPE3": None})

        nhdr = WCS(nhdr, preserve_units=True)
        if nhdr.naxis > naxes:
            for _ in range(nhdr.naxis - naxes):
                nhdr.dropaxis(-1)
        nhdr = nhdr.to_header()

        temp = nhdr.fromkeys(main_keys.keys())
        temp.update(main_keys)
        temp.add_history("Made by scrap.py")
        temp.add_history("https://www.lexyandati.com/scrappy/sc-rmmap/")
        return temp
