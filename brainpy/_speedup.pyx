from cpython.list cimport PyList_GET_ITEM, PyList_GET_SIZE, PyList_Append
from cpython.int cimport PyInt_FromLong
from cpython.float cimport PyFloat_FromDouble, PyFloat_AsDouble


cpdef _update_elementary_symmetric_polynomial(list power_sum, list elementary_symmetric_polynomial, size_t order):
    cdef:
        size_t begin, end, k, j
        double el
        int sign
    begin = PyList_GET_SIZE(elementary_symmetric_polynomial)
    end = PyList_GET_SIZE(power_sum)
    for k in range(begin, end):
        if k == 0:
            PyList_Append(elementary_symmetric_polynomial, 1.0)
        elif k > order:
            PyList_Append(elementary_symmetric_polynomial, 0.)
        else:
            el = 0.
            for j in range(1, k + 1):
                sign = 1 if (j % 2) == 1 else -1
                el += sign * PyFloat_AsDouble(<object>PyList_GET_ITEM(power_sum, j)) * PyFloat_AsDouble(<object>PyList_GET_ITEM(elementary_symmetric_polynomial, k - j))
            el /= <double>(k)
            PyList_Append(elementary_symmetric_polynomial, el)

