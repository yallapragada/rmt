def open_file_by_mimetype(fname, mode):
    '''This function determines the compression MIME type of a file as gz, bz, or none, and returns 
    an open file handle of the requested mode ('w', 'r', or 'a')'''
    from gzip import open as gzopen
    from bz2 import BZ2File as bzopen
    from mimetypes import guess_type

    if mode != 'r' and mode != 'w' and mode != 'a':
        print ("please specific a valid mode:  w, r, a")
        return

    if guess_type(fname)[1] == 'gzip': 
        try:
            fh = gzopen(fname, mode)
        except Exception as error:
            print ("Error opening file ", fname, ": ", error)
            return
    elif guess_type(fname) == 'bzip2':
        try:
            fh = bzopen(fname, mode)
        except Exception as error:
            print ("Error opening file ", fname, ": ", error)
            return
    else:
        try:
            fh = open(fname, mode)
        except Exception as error:
            print ("Error opening file ", fname, ": ", error)
            return

    return fh