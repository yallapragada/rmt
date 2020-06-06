from gzip import open as gzopen
from bz2 import BZ2File as bzopen
from mimetypes import guess_type


def open_file_by_mimetype(filename, mode):

    """
    This function determines the compression MIME type of a file as gz, bz, or none, and returns
    an open file handle of the requested mode ('w', 'r', or 'a')
    """

    if mode != 'r' and mode != 'w' and mode != 'a':
        print ("please specific a valid mode:  w, r, a")
        return

    if guess_type(filename)[1] == 'gzip':
        try:
            fh = gzopen(filename, mode)
        except Exception as error:
            print ("Error opening file ", filename, ": ", error)
            return
    elif guess_type(filename) == 'bzip2':
        try:
            fh = bzopen(filename, mode)
        except Exception as error:
            print ("Error opening file ", filename, ": ", error)
            return
    else:
        try:
            fh = open(filename, mode)
        except Exception as error:
            print ("Error opening file ", filename, ": ", error)
            return

    return fh
