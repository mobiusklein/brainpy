from brainpy._c.composition cimport (
    _ElementTable, Element, Isotope, Composition, ElementHashTable)

from brainpy._c.double_vector cimport DoubleVector

ctypedef DoubleVector dvec


cdef struct PolynomialParameters:
    dvec* elementary_symmetric_polynomial
    dvec* power_sum


cdef struct PhiConstants:
    int order
    Element* element
    PolynomialParameters* element_coefficients
    PolynomialParameters* mass_coefficients


cdef struct IsotopicConstants:
    int order
    PhiConstants** constants
    size_t size
    size_t used

cdef size_t DEFAULT_ISOTOPIC_CONSTANTS_SIZE = 7


cdef dvec* vietes(dvec* coefficients) nogil
cdef void _update_power_sum(dvec* ps_vec, dvec* esp_vec, int order) nogil
cdef void _update_elementary_symmetric_polynomial(dvec* ps_vec, dvec* esp_vec, int order) nogil
cdef void newton(dvec* ps_vec, dvec* esp_vec, int order) nogil
cdef dvec* compute_isotopic_coefficients(Element* element, bint with_mass, dvec* accumulator) nogil
cdef PolynomialParameters* make_polynomial_parameters(Element* element, bint with_mass, dvec* accumulator) nogil
cdef void print_polynomial_parameters(PolynomialParameters* params) nogil
cdef void free_polynomial_parameters(PolynomialParameters* params) nogil


cdef void print_phi_constants(PhiConstants* constants) nogil
cdef void free_phi_constants(PhiConstants* constants) nogil


cdef IsotopicConstants* make_isotopic_constants() nogil
cdef int isotopic_constants_resize(IsotopicConstants* ics) nogil
cdef void free_isotopic_constants(IsotopicConstants* isotopes) nogil

cdef void isotopic_constants_add_element(IsotopicConstants* isotopes, char* element_symbol) nogil
cdef int isotopic_constants_get(IsotopicConstants* isotopes, char* element_symbol, PhiConstants** out) nogil
cdef void isotopic_constants_update_coefficients(IsotopicConstants* isotopes) nogil

cdef double isotopic_constants_nth_element_power_sum(IsotopicConstants* isotopes, char* symbol, int order) nogil
cdef double isotopic_constants_nth_modified_element_power_sum(IsotopicConstants* isotopes, char* symbol, int order) nogil

cdef double isotopic_constants_nth_element_power_sum_by_index(IsotopicConstants* isotopes, size_t index, int order) nogil
cdef double isotopic_constants_nth_modified_element_power_sum_by_index(IsotopicConstants* isotopes, size_t index, int order) nogil

cdef void print_isotopic_constants(IsotopicConstants* isotopes) nogil
