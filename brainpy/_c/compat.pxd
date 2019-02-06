from cpython cimport PyObject

cdef extern from "compat.h":
    char* PyStr_AsString(str string)
    char* PyStr_AsUTF8AndSize(str string, Py_ssize_t*)
    str PyStr_FromString(char* string)
    str PyStr_FromStringAndSize(char* string, Py_ssize_t)
    long PyInt_AsLong(object i)
    object PyInt_FromLong(long i)
    void PyStr_InternInPlace(PyObject** string)
    str PyStr_InternFromString(const char*)