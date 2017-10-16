#!env python
# Works with either Python 2 or 3.

from ctypes import (
    cdll, c_int, c_double, c_char_p,
    POINTER, Structure,
)

c_double_p = POINTER(c_double)


class its_struct(Structure):
    _fields_ = [
        ("nb", c_int),
        ("mcycle", c_int),
        ("fb_filename", c_char_p),
        ("norml_filename", c_char_p),
        ("gfsum", c_double),
        ("bgfsum", c_double),
        ("highT", c_double),
        ("lowT", c_double),
        ("beta0", c_double),
        ("fb", c_double_p),
        ("norml", c_double_p),
        ("normlold", c_double_p),
        ("mybeta", c_double_p),
        ("gf", c_double_p),
        ("ratio", c_double_p),
        ("pratio", c_double_p),
        ("rb", c_double_p),
        ("rbfb", c_double_p),
        ("bgf", c_double_p)
    ]
its_struct_p = POINTER(its_struct)


class Its(object):
    _lib = cdll.LoadLibrary("./its.so")
    _lib.its_init.restype = its_struct_p

    def __init__(self, objtemp, n=200, highT=400.0, lowT=280.0, fb_filename="fb.dat", norml_filename="norml.dat"):
        self.data_p = Its._lib.its_init(
            c_double(objtemp), n,
            c_double(highT), c_double(lowT),
            fb_filename, norml_filename
        )
        self.data = self.data_p[0]

    def __del__(self):
        Its._lib.its_free(self.data_p)

    def updatefb(self):
        Its._lib.its_updatefb(self.data_p)

    def force(self, v, vshift, step, f, natom):
        Its._lib.its_force(c_double(v), c_double(vshift), step, f.ctypes.data_as(c_double_p), natom)
