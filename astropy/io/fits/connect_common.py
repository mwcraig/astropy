# Licensed under a 3-clause BSD style license - see LICENSE.rst

from __future__ import print_function

import re
import warnings

from .. import registry as io_registry
from ...utils.exceptions import AstropyUserWarning

from . import HDUList, TableHDU, BinTableHDU, GroupsHDU, PrimaryHDU, ImageHDU


# FITS file signature as per RFC 4047
FITS_SIGNATURE = (b"\x53\x49\x4d\x50\x4c\x45\x20\x20\x3d\x20\x20\x20\x20\x20"
                  b"\x20\x20\x20\x20\x20\x20\x20\x20\x20\x20\x20\x20\x20\x20"
                  b"\x20\x54")

# Keywords to remove for all tables that are read in
REMOVE_KEYWORDS = ['XTENSION', 'BITPIX', 'NAXIS', 'NAXIS1', 'NAXIS2',
                   'PCOUNT', 'GCOUNT', 'TFIELDS']

# Column-specific keywords
COLUMN_KEYWORDS = ['TFORM[0-9]+',
                   'TBCOL[0-9]+',
                   'TSCAL[0-9]+',
                   'TZERO[0-9]+',
                   'TNULL[0-9]+',
                   'TTYPE[0-9]+',
                   'TUNIT[0-9]+',
                   'TDISP[0-9]+',
                   'TDIM[0-9]+',
                   'THEAP']


def is_column_keyword(keyword):
    for c in COLUMN_KEYWORDS:
        if re.match(c, keyword) is not None:
            return True
    return False


def is_fits(origin, filepath, fileobj, *args, **kwargs):
    """
    Determine whether `origin` is a FITS file.

    Parameters
    ----------
    origin : str or readable file-like object
        Path or file object containing a potential FITS file.

    Returns
    -------
    is_fits : bool
        Returns `True` if the given file is a FITS file.
    """
    if fileobj is not None:
        pos = fileobj.tell()
        sig = fileobj.read(30)
        fileobj.seek(pos)
        return sig == FITS_SIGNATURE
    elif filepath is not None:
        if filepath.lower().endswith(('.fits', '.fits.gz', '.fit', '.fit.gz')):
            return True
    elif isinstance(args[0], (HDUList, TableHDU, BinTableHDU, GroupsHDU, PrimaryHDU, ImageHDU)):
        return True
    else:
        return False


def _put_header_in_meta(header, meta):

    for key, value, comment in header.cards:

        if key in ['COMMENT', 'HISTORY']:
            if key in meta:
                meta[key].append(value)
            else:
                meta[key] = [value]

        elif key in meta:  # key is duplicate

            if isinstance(meta[key], list):
                meta[key].append(value)
            else:
                meta[key] = [meta[key], value]

        elif (is_column_keyword(key.upper()) or
              key.upper() in REMOVE_KEYWORDS):

            pass

        else:

            meta[key] = value


def _put_meta_in_header(meta, header):

    for key, value in meta.items():

        if is_column_keyword(key.upper()) or key.upper() in REMOVE_KEYWORDS:

            warnings.warn(
                "Meta-data keyword {0} will be ignored since it conflicts "
                "with a FITS reserved keyword".format(key), AstropyUserWarning)

        if isinstance(value, list):
            for item in value:
                try:
                    header.append((key, item))
                except ValueError:
                    warnings.warn(
                        "Attribute `{0}` of type {1} cannot be written to "
                        "FITS files - skipping".format(key, type(value)),
                        AstropyUserWarning)
        else:
            try:
                header[key] = value
            except ValueError:
                warnings.warn(
                    "Attribute `{0}` of type {1} cannot be written to FITS "
                    "files - skipping".format(key, type(value)),
                    AstropyUserWarning)
