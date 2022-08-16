# cython: embedsignature=True

cimport cython

from brainpy._c.composition cimport (
    Element, Isotope, Composition, ElementHashTable,
    element_max_neutron_shift, _parse_isotope_string,
    _ElementTable, element_hash_table_get, make_fixed_isotope_element,
    element_hash_table_put, make_element_hash_table, free_element_hash_table)

from brainpy._c.double_vector cimport(
    DoubleVector, make_double_vector, double_vector_append,
    free_double_vector, print_double_vector, double_vector_to_list,
    reset_double_vector, make_double_vector_with_size)

from libc.stdlib cimport malloc, free, realloc
from libc.string cimport strcmp
from libc cimport *

cdef extern from * nogil:
    int printf (const char *template, ...)

# -----------------------------------------------------------------------------
# PolynomialParameters Methods

@cython.cdivision
cdef dvec* vietes(dvec* coefficients) nogil:
    cdef:
        DoubleVector* elementary_symmetric_polynomial
        size_t i
        double tail
        int sign
        double el

    elementary_symmetric_polynomial = make_double_vector_with_size(
        coefficients.used)
    tail = coefficients.v[coefficients.used - 1]

    for i in range(coefficients.used):
        sign = 1 if (i % 2) == 0 else -1
        el = sign * coefficients.v[coefficients.used - i - 1] / tail
        double_vector_append(elementary_symmetric_polynomial, el)

    return elementary_symmetric_polynomial


cdef void _update_power_sum(dvec* ps_vec, dvec* esp_vec, int order) nogil:
    cdef:
        size_t begin, end, k, j
        int sign
        double temp_ps

    begin = ps_vec.used
    end = esp_vec.used

    for k in range(begin, end):
        if k == 0:
            double_vector_append(ps_vec, 0.0)
            continue
        temp_ps = 0.
        sign = -1

        for j in range(1, k):
            sign *= -1
            temp_ps += sign * (esp_vec.v[j]) * (ps_vec.v[k - j])
        sign *= -1
        temp_ps += sign * (esp_vec.v[k]) * k
        double_vector_append(ps_vec, temp_ps)


@cython.cdivision
cdef void _update_elementary_symmetric_polynomial(dvec* ps_vec, dvec* esp_vec, int order) nogil:
    cdef:
        size_t begin, end, k, j
        int sign
        double el

    begin = esp_vec.used
    end = ps_vec.used

    for k in range(begin, end):
        if k == 0:
            double_vector_append(esp_vec, 1.0)
        elif k > order:
            double_vector_append(esp_vec, 0.0)
        else:
            el = 0.
            for j in range(1, k + 1):
                sign = 1 if (j % 2) == 1 else -1
                el += sign * ps_vec.v[j] * esp_vec.v[k - j]
            el /= <double>k
            double_vector_append(esp_vec, el)


cdef void newton(dvec* ps_vec, dvec* esp_vec, int order) nogil:
    if ps_vec.used > esp_vec.used:
        _update_elementary_symmetric_polynomial(ps_vec, esp_vec, order)
    elif ps_vec.used < esp_vec.used:
        _update_power_sum(ps_vec, esp_vec, order)


cdef dvec* compute_isotopic_coefficients(Element* element, bint with_mass, dvec* accumulator) nogil:
    cdef:
        int max_isotope_number, current_order
        Isotope* isotope
        double coef
        size_t i, j, k

    max_isotope_number = element_max_neutron_shift(element)
    for i in range(element.isotopes.size):
        k = element.isotopes.size - i - 1
        isotope = &(element.isotopes.bins[k])
        current_order = max_isotope_number - isotope.neutron_shift
        if with_mass:
            coef = isotope.mass
        else:
            coef = 1.0
        if current_order > accumulator.used:
            for j in range(accumulator.used, current_order):
                double_vector_append(accumulator, 0.)
            double_vector_append(accumulator, isotope.abundance * coef)
        elif current_order == accumulator.used:
            double_vector_append(accumulator, isotope.abundance * coef)
        else:
            printf("Error, unordered isotopes for %s\n", element.symbol)
    return accumulator

cdef PolynomialParameters* make_polynomial_parameters(Element* element, bint with_mass, dvec* accumulator) nogil:
    cdef:
        dvec* elementary_symmetric_polynomial
        dvec* power_sum
        PolynomialParameters* result

    compute_isotopic_coefficients(element, with_mass, accumulator)

    elementary_symmetric_polynomial = vietes(accumulator)
    power_sum = make_double_vector_with_size(elementary_symmetric_polynomial.used + 4)
    newton(power_sum, elementary_symmetric_polynomial, accumulator.used - 1)
    result = <PolynomialParameters*>malloc(sizeof(PolynomialParameters))
    result.elementary_symmetric_polynomial = elementary_symmetric_polynomial
    result.power_sum = power_sum
    return result


cdef void print_polynomial_parameters(PolynomialParameters* params) nogil:
    printf("PolynomialParameters: %d\n", <int>params)
    printf("  ")
    print_double_vector(params.elementary_symmetric_polynomial)
    printf("  ")
    print_double_vector(params.power_sum)
    printf("\n")


cdef void free_polynomial_parameters(PolynomialParameters* params) nogil:
    free_double_vector(params.elementary_symmetric_polynomial)
    free_double_vector(params.power_sum)
    free(params)

# -----------------------------------------------------------------------------
# PhiConstants Methods


cdef void print_phi_constants(PhiConstants* constants) nogil:
    printf("PhiConstants: %d\n", <int>constants)
    printf("Element: %s, Order: %d\n", constants.element.symbol, constants.order)
    printf("Element Coefficients:\n")
    print_polynomial_parameters(constants.element_coefficients)
    printf("Mass Coefficients:\n")
    print_polynomial_parameters(constants.mass_coefficients)
    printf("\n")


cdef void free_phi_constants(PhiConstants* constants) nogil:
    free_polynomial_parameters(constants.element_coefficients)
    free_polynomial_parameters(constants.mass_coefficients)
    free(constants)


# -----------------------------------------------------------------------------
# IsotopicConstants Methods

cdef size_t DEFAULT_ISOTOPIC_CONSTANTS_SIZE = 7


cdef IsotopicConstants* make_isotopic_constants() nogil:
    cdef:
        IsotopicConstants* result
    result = <IsotopicConstants*>malloc(sizeof(IsotopicConstants))
    result.constants = <PhiConstants**>malloc(sizeof(PhiConstants*) * DEFAULT_ISOTOPIC_CONSTANTS_SIZE)
    result.size = DEFAULT_ISOTOPIC_CONSTANTS_SIZE
    result.used = 0
    return result


cdef int isotopic_constants_resize(IsotopicConstants* ics) nogil:
    ics.constants = <PhiConstants**>realloc(ics.constants, sizeof(PhiConstants*) * ics.size * 2)
    ics.size *= 2
    if ics.constants == NULL:
        return -1
    return 0


cdef void free_isotopic_constants(IsotopicConstants* isotopes) nogil:
    cdef:
        size_t i

    i = 0
    while i < isotopes.used:
        free_phi_constants(isotopes.constants[i])
        i += 1
    free(isotopes.constants)
    free(isotopes)


cdef void isotopic_constants_add_element(IsotopicConstants* isotopes, char* element_symbol) nogil:
    cdef:
        Element* element
        dvec* accumulator
        int order, status, isotope_number
        PolynomialParameters* element_parameters
        PolynomialParameters* mass_parameters
        PhiConstants* phi_constants
        char* element_buffer

    status = isotopic_constants_get(isotopes, element_symbol, &phi_constants)
    if status == 0:
        return

    phi_constants = NULL

    status = element_hash_table_get(_ElementTable, element_symbol, &element)
    if status == -1:
        printf("Could not resolve element_symbol %s\n", element_symbol)
        return
    accumulator = make_double_vector()
    order = element_max_neutron_shift(element)
    element_parameters = make_polynomial_parameters(element, False, accumulator)
    reset_double_vector(accumulator)
    mass_parameters = make_polynomial_parameters(element, True, accumulator)
    free_double_vector(accumulator)
    phi_constants = <PhiConstants*>malloc(sizeof(PhiConstants))
    phi_constants.order = order
    phi_constants.element = element
    phi_constants.element_coefficients = element_parameters
    phi_constants.mass_coefficients = mass_parameters
    if isotopes.used + 1 == isotopes.size:
        isotopic_constants_resize(isotopes)
    isotopes.constants[isotopes.used] = phi_constants
    isotopes.used += 1


cdef int isotopic_constants_get(IsotopicConstants* isotopes, char* element_symbol, PhiConstants** out) nogil:
    cdef:
        size_t i
        int status

    i = 0
    while i < isotopes.used:
        if strcmp(isotopes.constants[i].element.symbol, element_symbol) == 0:
            out[0] = (isotopes.constants[i])
            return 0
        i += 1
    return 1


cdef int isotopic_constants_get_by_index(IsotopicConstants* isotopes, size_t index, PhiConstants** out) nogil:
    if index > isotopes.used:
        return 1
    else:
        out[0] = isotopes.constants[index]
        return 0


cdef void isotopic_constants_update_coefficients(IsotopicConstants* isotopes) nogil:
    cdef:
        size_t i, j
        PhiConstants* constants

    for i in range(isotopes.used):
        constants = isotopes.constants[i]

        if isotopes.order < constants.order:
            continue

        for j in range(constants.order, isotopes.order + 1):
            double_vector_append(constants.element_coefficients.elementary_symmetric_polynomial, 0.)
            double_vector_append(constants.mass_coefficients.elementary_symmetric_polynomial, 0.)

        constants.order = constants.element_coefficients.elementary_symmetric_polynomial.used

        newton(constants.element_coefficients.power_sum,
               constants.element_coefficients.elementary_symmetric_polynomial,
               constants.order)
        newton(constants.mass_coefficients.power_sum,
               constants.mass_coefficients.elementary_symmetric_polynomial,
               constants.order)


cdef double isotopic_constants_nth_element_power_sum(IsotopicConstants* isotopes, char* symbol, int order) nogil:
    cdef:
        PhiConstants* constants
    isotopic_constants_get(isotopes, symbol, &constants)
    return constants.element_coefficients.power_sum.v[order]


cdef double isotopic_constants_nth_element_power_sum_by_index(IsotopicConstants* isotopes, size_t index, int order) nogil:
    cdef:
        PhiConstants* constants
    isotopic_constants_get_by_index(isotopes, index, &constants)
    return constants.element_coefficients.power_sum.v[order]


cdef double isotopic_constants_nth_modified_element_power_sum(IsotopicConstants* isotopes, char* symbol, int order) nogil:
    cdef:
        PhiConstants* constants
    isotopic_constants_get(isotopes, symbol, &constants)
    return constants.mass_coefficients.power_sum.v[order]


cdef double isotopic_constants_nth_modified_element_power_sum_by_index(IsotopicConstants* isotopes, size_t index, int order) nogil:
    cdef:
        PhiConstants* constants
    isotopic_constants_get_by_index(isotopes, index, &constants)
    return constants.mass_coefficients.power_sum.v[order]


cdef void print_isotopic_constants(IsotopicConstants* isotopes) nogil:
    cdef:
        size_t i

    i = 0
    for i in range(isotopes.used):
        printf("%d\n", i)
        print_phi_constants(isotopes.constants[i])


def main():
    cdef:
        IsotopicConstants* ic
        dvec* coefficients
        dvec* elementary_symmetric_polynomial
        PhiConstants* constant
        Element* elem
        char* sym

    sym = "O"
    ic = make_isotopic_constants()
    print(ic.used)
    isotopic_constants_add_element(ic, sym)
    print(ic.used)
    isotopic_constants_add_element(ic, "C")
    isotopic_constants_add_element(ic, "H")

    if isotopic_constants_get(ic, "O", &constant) == 0:
        print_phi_constants(constant)
    else:
        print("Nope")

    free_isotopic_constants(ic)
