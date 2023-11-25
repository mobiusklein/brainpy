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


cdef dvec* vietes(dvec* coefficients) noexcept nogil
cdef void _update_power_sum(dvec* ps_vec, dvec* esp_vec, int order) noexcept nogil
cdef void _update_elementary_symmetric_polynomial(dvec* ps_vec, dvec* esp_vec, int order) noexcept nogil
cdef void newton(dvec* ps_vec, dvec* esp_vec, int order) noexcept nogil
cdef dvec* compute_isotopic_coefficients(Element* element, bint with_mass, dvec* accumulator) noexcept nogil
cdef PolynomialParameters* make_polynomial_parameters(Element* element, bint with_mass, dvec* accumulator) noexcept nogil
cdef void print_polynomial_parameters(PolynomialParameters* params) noexcept nogil
cdef void free_polynomial_parameters(PolynomialParameters* params) noexcept nogil


cdef void print_phi_constants(PhiConstants* constants) noexcept nogil
cdef void free_phi_constants(PhiConstants* constants) noexcept nogil


cdef IsotopicConstants* make_isotopic_constants() noexcept nogil
cdef int isotopic_constants_resize(IsotopicConstants* ics) noexcept nogil
cdef void free_isotopic_constants(IsotopicConstants* isotopes) noexcept nogil

cdef void isotopic_constants_add_element(IsotopicConstants* isotopes, char* element_symbol) noexcept nogil
cdef int isotopic_constants_get(IsotopicConstants* isotopes, char* element_symbol, PhiConstants** out) noexcept nogil
cdef void isotopic_constants_update_coefficients(IsotopicConstants* isotopes) noexcept nogil

cdef double isotopic_constants_nth_element_power_sum(IsotopicConstants* isotopes, char* symbol, int order) noexcept nogil
cdef double isotopic_constants_nth_modified_element_power_sum(IsotopicConstants* isotopes, char* symbol, int order) noexcept nogil

cdef double isotopic_constants_nth_element_power_sum_by_index(IsotopicConstants* isotopes, size_t index, int order) noexcept nogil
cdef double isotopic_constants_nth_modified_element_power_sum_by_index(IsotopicConstants* isotopes, size_t index, int order) noexcept nogil

cdef void print_isotopic_constants(IsotopicConstants* isotopes) noexcept nogil
