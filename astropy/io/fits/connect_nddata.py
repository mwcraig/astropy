# Licensed under a 3-clause BSD style license - see LICENSE.rst

from __future__ import print_function

import os
import warnings

from .. import registry as io_registry
from ...extern import six
from ...nddata import NDData
from ...utils import OrderedDict
from ...utils.exceptions import AstropyUserWarning

from . import HDUList, PrimaryHDU, ImageHDU, Header
from .hdu.hdulist import fitsopen as fits_open
from .util import first

from .connect_common import _put_meta_in_header, _put_header_in_meta, is_fits


def read_nddata_fits(input, hdu=None):
    """
    Read FITS file into object that implements NDDataBase interface.

    Parameters
    ----------

    input : str or file-like object or compatible `astropy.io.fits` HDU object
        If a string, the filename to read the table from. If a file object, or
        a compatible HDU object, the object to extract the table from. The
        following `astropy.io.fits` HDU objects can be used as input:
        - :class:`~astropy.io.fits.hdu.table.PrimaryHDU`
        - :class:`~astropy.io.fits.hdu.table.ImageHDU`
        - :class:`~astropy.io.fits.hdu.hdulist.HDUList`
    hdu : int or str, optional
        The HDU to read the table from.
    """

    if isinstance(input, HDUList):

        # Parse all nddata objects
        nddata_hdus = OrderedDict()
        for idx, an_hdu in enumerate(input):
            if isinstance(an_hdu, (ImageHDU, PrimaryHDU)) and an_hdu.header['naxis'] > 0:
                nddata_hdus[idx] = an_hdu

        if len(nddata_hdus) > 1:
            if hdu is None:
                warnings.warn("hdu= was not specified but multiple tables"
                              " are present, reading in first available"
                              " table (hdu={0})".format(first(nddata_hdus)),
                              AstropyUserWarning)
                hdu = first(nddata_hdus)

                # Grabbed the stuff below directly from table fits reader
                # hdu might not be an integer, so we first need to convert it
                # to the correct HDU index
                hdu = input.index_of(hdu)

                if hdu in nddata_hdus:
                    nddata = nddata_hdus[hdu]
                else:
                    raise ValueError("No table found in hdu={0}".format(hdu))

        elif len(nddata_hdus) == 1:
            nddata = nddata_hdus[first(nddata_hdus)]
        else:
            raise ValueError("No NDData-like extension found")

    elif isinstance(input, (ImageHDU, PrimaryHDU)) and input.header['naxis'] > 0:

        nddata = input

    else:

        # input is either a string or a file-like object, hopefully.
        hdulist = fits_open(input)

        try:
            return read_nddata_fits(hdulist, hdu=hdu)
        finally:
            hdulist.close()

    nd = NDData(nddata.data)

    _put_header_in_meta(nddata.header, nd.meta)

    return nd


def write_nddata_fits(input, output, overwrite=False):
    """
    Write an NDData-like object to a FITS file.

    Parameters
    ----------
    input : object that implements NDData interface
        The object whose data is to be written.
    output : str
        The filename to write the table to.
    overwrite : bool
        Whether to overwrite any existing file without warning.
    """

    # Thank you, io.fits.connect.write, for this snippet.
    # Check if output file already exists
    if isinstance(output, six.string_types) and os.path.exists(output):
        if overwrite:
            os.remove(output)
        else:
            raise IOError("File exists: {0}".format(output))

    hdu = PrimaryHDU(data=input.data, header=Header(input.meta))

    _put_meta_in_header(input.meta, hdu.header)

    hdu.writeto(output)

io_registry.register_reader('fits', NDData, read_nddata_fits)
io_registry.register_writer('fits', NDData, write_nddata_fits)
io_registry.register_identifier('fits', NDData, is_fits)
