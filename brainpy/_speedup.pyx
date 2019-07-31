# cython: embedsignature=True
# cython: profile=False

from cpython.list cimport PyList_GET_ITEM, PyList_GET_SIZE, PyList_Append
from cpython.int cimport PyInt_FromLong
from cpython.float cimport PyFloat_FromDouble, PyFloat_AsDouble
from cpython.dict cimport PyDict_GetItem, PyDict_SetItem, PyDict_Keys, PyDict_Next
from cpython.object cimport PyObject

from libc.math cimport log, exp

import operator
import re

mz_getter = operator.attrgetter("mz")

from mass_dict import nist_mass as _nist_mass
from brainpy._c.isotopic_distribution cimport TheoreticalPeak as Peak


nist_mass = _nist_mass
PROTON = nist_mass["H+"][0][0]


cdef double neutral_mass(double mz,  int z, double charge_carrier=PROTON):
    return (mz * abs(z)) - (z * charge_carrier)


cdef double mass_charge_ratio(double neutral_mass, int z, double charge_carrier=PROTON):
    return (neutral_mass + (z * charge_carrier)) / abs(z)


cdef void _update_elementary_symmetric_polynomial(list power_sum, list elementary_symmetric_polynomial, size_t order):
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

cdef void _update_power_sum(list ps_vec, list esp_vec, size_t order):
    cdef:
        size_t begin, end, k, j
        int sign
        double temp_ps
    begin = PyList_GET_SIZE(ps_vec)
    end = PyList_GET_SIZE(esp_vec)
    for k in range(begin, end):
        if k == 0:
            PyList_Append(ps_vec, 0.)
            continue
        temp_ps = 0.
        sign = -1
        for j in range(1, k):
            sign *= -1
            temp_ps += sign * PyFloat_AsDouble(<object>PyList_GET_ITEM(esp_vec, j)) * PyFloat_AsDouble(<object>PyList_GET_ITEM(ps_vec, k - j))
        sign *= -1
        temp_ps += sign * PyFloat_AsDouble(<object>PyList_GET_ITEM(esp_vec, k)) * k
        PyList_Append(ps_vec, temp_ps)

cdef void newton(list power_sum, list elementary_symmetric_polynomial, int order):
    if PyList_GET_SIZE(power_sum) > PyList_GET_SIZE(elementary_symmetric_polynomial):
        _update_elementary_symmetric_polynomial(power_sum, elementary_symmetric_polynomial, order)
    elif PyList_GET_SIZE(power_sum) < PyList_GET_SIZE(elementary_symmetric_polynomial):
        _update_power_sum(power_sum, elementary_symmetric_polynomial, order)

cdef list vietes(list coefficients):
    cdef:
        list elementary_symmetric_polynomial
        double tail, el
        size_t size, i
        int sign

    elementary_symmetric_polynomial = []
    tail = float(coefficients[-1])
    size = len(coefficients)

    for i in range(size):
        sign = 1 if (i % 2) == 0 else -1
        el = sign * coefficients[size - i - 1] / tail
        elementary_symmetric_polynomial.append(el)
    return elementary_symmetric_polynomial


cdef class PolynomialParameters(object):
    cdef:
        public list elementary_symmetric_polynomial
        public list power_sum
    def __init__(self, elementary_symmetric_polynomial, power_sum):
        self.elementary_symmetric_polynomial = elementary_symmetric_polynomial
        self.power_sum = power_sum

    def __iter__(self):
        yield self.power_sum
        yield self.elementary_symmetric_polynomial


cdef class PhiConstants(object):
    cdef:
        public int order
        public Element element
        public PolynomialParameters element_coefficients
        public PolynomialParameters mass_coefficients

    def __init__(self, order, element, element_coefficients, mass_coefficients):
        self.order = order
        self.element = element
        self.element_coefficients = element_coefficients
        self.mass_coefficients = mass_coefficients


cdef class Isotope(object):
    '''
    Isotope represents an elenent with an integer number of neutrons specified.

    Attributes
    ----------
    mass: float
        Measurable mass in Da
    abundance: float [0.0:1.0]
        The abundance of this isotope in nature. This number is between 0.0 and 1.0
    neutron_shift: int
        The number of neutrons different between this isotope and the "normal" form. May be 0
        if this represents that normal form.
    '''
    cdef:
        public double mass
        public double abundance
        public int neutron_shift
        public int neutrons

    def __init__(self, mass, abundance, neutron_shift, neutrons):
        self.mass = mass
        self.abundance = abundance
        self.neutron_shift = neutron_shift
        self.neutrons = neutrons

    def __repr__(self):
        return "Isotope(mass=%0.3f, abundance=%0.3f, neutron_shift=%d, neutrons=%d)" % (
            self.mass, self.abundance, self.neutron_shift, self.neutrons)



cdef int max_variants(dict composition):
    """Calculates the maximum number of isotopic variants that could be produced by a
    composition.

    Parameters
    ----------
    composition : Mapping
        Any Mapping type where keys are element symbols and values are integers

    Returns
    -------
        max_n_variants : int
    """
    max_n_variants = 0

    for element, count in composition.items():
        if element == "H+":
            continue
        try:
            max_n_variants += count * periodic_table[element].max_neutron_shift()
        except KeyError:
            pass

    return max_n_variants


cdef double calculate_mass(dict composition, dict mass_data=None):
    cdef:
        double mass
        object match
    mass = 0.0
    if mass_data is None:
        mass_data = nist_mass
    for element in composition:
            try:
                mass += (composition[element] * mass_data[element][0][0])
            except KeyError:
                match = re.search(r"(\S+)\[(\d+)\]", element)
                if match:
                    element_ = match.group(1)
                    isotope = int(match.group(2))
                    mass += composition[element] * mass_data[element_][isotope][0]
                else:
                    raise

    return mass


cdef str _make_isotope_string(str element, int isotope=0):
    if isotope == 0:
        return element
    else:
        return "%s[%d]" % (element, isotope)


cdef tuple _get_isotope(str element_string):
    cdef:
        object match
        str element_
        int isotope
    if "[" in element_string:
        match = re.search(r"(\S+)\[(\d+)\]", element_string)
        if match:
            element_ = match.group(1)
            isotope = int(match.group(2))
            return element_, isotope
    else:
        return element_string, 0


cpdef Element make_fixed_isotope_element(Element element, int neutrons):
    cdef:
        Isotope isotope
        Element el
    isotope = element.isotopes[neutrons - element.isotopes[0].neutrons]
    el = Element(element.symbol + ("[%d]" % neutrons), dict([
            (0, Isotope(isotope.mass, abundance=1.0, neutron_shift=0, neutrons=neutrons)),
        ]))
    return el


def _isotopes_of(element):
    freqs = dict()
    for i, mass_freqs in nist_mass[element].items():
        if i == 0:
            continue
        if mass_freqs[1] > 0:
            freqs[i] = mass_freqs
    if len(freqs) == 0:
        return dict()
    mono_neutrons = max(freqs.items(), key=lambda x: x[1][1])[0]
    freqs = list(sorted(
        [(k - mono_neutrons, Isotope(*v, neutron_shift=k - mono_neutrons, neutrons=k))
                                for k, v in freqs.items()], key=lambda x: x[0]))
    return dict(freqs)


cdef class Element(object):

    def __init__(self, str symbol, dict isotopes=None):
        if isotopes is None:
            isotopes = _isotopes_of(symbol)
        self.symbol = symbol
        self.isotopes = isotopes
        min_shift = 1000
        max_shift = 0
        for shift in self.isotopes:
            if shift > max_shift:
                max_shift = shift
            if shift < min_shift:
                min_shift = shift
        self._min_neutron_shift = min_shift
        self._max_neutron_shift = max_shift
        self._no_mass_elementary_symmetric_polynomial_cache = None
        self._no_mass_power_sum_cache = None
        self._mass_elementary_symmetric_polynomial_cache = None
        self._mass_power_sum_cache = None
        try:
            self._monoisotopic_mass = self.isotopes[0].mass
        except:
            self._monoisotopic_mass = nist_mass[self.symbol][0][0]

    def __iter__(self):
        for key in sorted(self.isotopes.keys()):
            yield self.isotopes[key]

    def max_neutron_shift(self):
        return self._max_neutron_shift

    def min_neutron_shift(self):
        return self._min_neutron_shift

    cpdef double monoisotopic_mass(self):
        return self._monoisotopic_mass


periodic_table = periodic_table = {k: Element(k) for k in nist_mass}


cdef class IsotopicConstants(dict):
    cdef:
        public long _order

    def __init__(self, order):
        self._order = 0
        self.order = order

    property order:
        def __get__(self):
            return self._order

        def __set__(self, value):
            self._order = value
            self.update_coefficients()

    cdef PolynomialParameters coefficients(self, Element element, bint with_mass=False):
        cdef:
            int max_isotope_number, current_order
            list accumulator, isotope_keys
            Isotope isotope
            double coef
            size_t isotope_iter, i

        if with_mass:
            if element._mass_elementary_symmetric_polynomial_cache is not None:
                return PolynomialParameters(
                    list(element._mass_elementary_symmetric_polynomial_cache),
                    list(element._mass_power_sum_cache))
        else:
            if element._no_mass_elementary_symmetric_polynomial_cache is not None:
                return PolynomialParameters(
                    list(element._no_mass_elementary_symmetric_polynomial_cache),
                    list(element._no_mass_power_sum_cache))

        max_isotope_number = element._max_neutron_shift
        isotope_keys = sorted(element.isotopes, reverse=True)
        accumulator = []
        for isotope_iter in range(PyList_GET_SIZE(isotope_keys)):
            isotope = <Isotope>PyDict_GetItem(element.isotopes, <object>PyList_GET_ITEM(isotope_keys, isotope_iter))
            current_order = max_isotope_number - isotope.neutron_shift
            if with_mass:
                coef = isotope.mass
            else:
                coef = 1.

            if current_order > len(accumulator):
                for i in range(len(accumulator), current_order):
                    PyList_Append(accumulator, 0.)
                PyList_Append(accumulator, isotope.abundance * coef)
            elif current_order == len(accumulator):
                PyList_Append(accumulator, isotope.abundance * coef)
            else:
                raise Exception("The list of neutron shifts is not ordered.")

        elementary_symmetric_polynomial = vietes(accumulator)
        power_sum = []
        newton(power_sum, elementary_symmetric_polynomial, len(accumulator) - 1)
        if with_mass:
            if element._mass_elementary_symmetric_polynomial_cache is None:
                element._mass_elementary_symmetric_polynomial_cache = list(elementary_symmetric_polynomial)
                element._mass_power_sum_cache = list(power_sum)
        else:
            if element._no_mass_elementary_symmetric_polynomial_cache is None:
                element._no_mass_elementary_symmetric_polynomial_cache = list(elementary_symmetric_polynomial)
                element._no_mass_power_sum_cache = list(power_sum)
        return PolynomialParameters(elementary_symmetric_polynomial, power_sum)

    cdef void add_element(self, str symbol):
        cdef:
            Element element
            int order
            int isotope
            str symbol_parsed
            PolynomialParameters element_parameters, mass_parameters

        if symbol in self:
            return
        try:
            element = periodic_table[symbol]
        except KeyError:
            symbol_parsed, isotope = _get_isotope(symbol)
            if isotope == 0:
                raise KeyError(symbol)
            element = make_fixed_isotope_element(periodic_table[symbol_parsed], isotope)
        order = element.max_neutron_shift()
        element_parameters = self.coefficients(element)
        mass_parameters = self.coefficients(element, True)
        self[symbol] = PhiConstants(order, element, element_parameters, mass_parameters)

    cdef void update_coefficients(self):
        cdef:
            str symbol
            PhiConstants phi_constants
            size_t i

        for symbol, phi_constants in self.items():
            if self.order < phi_constants.order:
                continue

            for i in range(phi_constants.order, self.order + 1):
                phi_constants.element_coefficients.elementary_symmetric_polynomial.append(0.)
                phi_constants.mass_coefficients.elementary_symmetric_polynomial.append(0.)

            phi_constants.order = len(phi_constants.element_coefficients.elementary_symmetric_polynomial)
            newton(phi_constants.element_coefficients.power_sum,
                   phi_constants.element_coefficients.elementary_symmetric_polynomial,
                   phi_constants.order)
            newton(phi_constants.mass_coefficients.power_sum,
                   phi_constants.mass_coefficients.elementary_symmetric_polynomial,
                   phi_constants.order)

    cdef double nth_element_power_sum(self, str symbol, int order):
        cdef:
            PhiConstants constants
        constants = <PhiConstants>self[symbol]
        return constants.element_coefficients.power_sum[order]

    cdef double nth_modified_element_power_sum(self, str symbol, int order):
        cdef:
            PhiConstants constants
        constants = <PhiConstants>self[symbol]
        return constants.mass_coefficients.power_sum[order]


cdef class IsotopicDistribution(object):
    cdef:
        public dict composition
        public IsotopicConstants _isotopic_constants
        public int _order
        public double average_mass
        public Peak monoisotopic_peak

    def __init__(self, composition, order=-1):
        self.composition = dict(composition)
        self._isotopic_constants = IsotopicConstants(order)
        self._order = 0
        self.order = order
        self.average_mass = 0.
        self.monoisotopic_peak = self._create_monoisotopic_peak()

    property order:
        def __get__(self):
            return self._order

        def __set__(self, value):
            max_variant_count = max_variants(self.composition)
            if value == -1:
                self._order = max_variant_count
            else:
                self._order = min(value, max_variant_count)
            self._update_isotopic_constants()

    cpdef _update_isotopic_constants(self):
        cdef:
            str element
        for element in self.composition:
            self._isotopic_constants.add_element(element)
        self._isotopic_constants.order = self._order

    cdef Peak _create_monoisotopic_peak(self):
        mass = calculate_mass(self.composition)
        intensity = 0.
        for element in self.composition:
            if element == "H+":
                continue
            # intensity += log(periodic_table[element].isotopes[0].abundance)
            intensity += log(self._isotopic_constants[element].element.isotopes[0].abundance)
        intensity = exp(intensity)
        return Peak(mass, intensity, 0)

    cdef double _phi_value(self, int order):
        cdef:
            double phi
            str element
            double count
        phi = 0.
        for element, count in self.composition.items():
            if element == "H+":
                continue
            phi += self._isotopic_constants.nth_element_power_sum(element, order) * count
        return phi

    cdef double _modified_phi_value(self, str symbol, int order):
        cdef:
            double phi
            str element
            double count
            double coef
            Py_ssize_t pos
            PyObject* pk
            PyObject* pv

        phi = 0.
        pos = 0
        #for element, count in self.composition.items():
        while PyDict_Next(self.composition, &pos, &pk, &pv):
            element = <str>pk
            count = PyFloat_AsDouble(<object>pv)
            if element == "H+":
                continue
            # Count is one lower for this symbol because an isotope is present
            # accounted for in the call to `nth_modified_element_power_sum` at
            # the end?
            coef = (count if element != symbol else count - 1)
            phi += self._isotopic_constants.nth_element_power_sum(element, order) * coef

        phi += self._isotopic_constants.nth_modified_element_power_sum(symbol, order)
        return phi

    cpdef list phi_values(self):
        cdef:
            list power_sum
            size_t i
        power_sum = [0.]
        for i in range(1, self.order + 1):
            power_sum.append(self._phi_value(i))
        return power_sum

    cpdef list modified_phi_values(self, symbol):
        cdef:
            list power_sum
            size_t i
        power_sum = [0.]
        for i in range(1, self.order + 1):
            power_sum.append(self._modified_phi_value(symbol, i))
        return power_sum

    cpdef list probability(self):
        cdef:
            list phi_values
            int max_variant_count
            list probability_vector
            size_t i
            int sign

        phi_values = self.phi_values()
        max_variant_count = max_variants(self.composition)
        probability_vector = []
        newton(phi_values, probability_vector, max_variant_count)

        for i in range(0, len(probability_vector)):
            # The sign of each term in the probability vector (populated by
            # Newton's Identities by solving for the Elementary Symmetric Polynomial
            # given the Power Sums) alternates in the same order as `sign`.
            # This ensures that the probability vector is strictly positive.
            sign = 1 if i % 2 == 0 else -1
            # q(j) = q(0)  * e(j) * (-1)^j
            # intensity of the jth peak is |probability[j]| * the intensity of monoisotopic peak
            probability_vector[i] *= self.monoisotopic_peak.intensity * sign

        return probability_vector

    cpdef list center_mass(self, list probability_vector):
        cdef:
            list mass_vector
            list composition_elements
            int max_variant_count, sign
            dict ele_sym_poly_map
            str element
            list ele_sym_poly
            list power_sum
            size_t i, j
            Py_ssize_t k
            Element element_obj
            PhiConstants phi_obj
            double center, temp
            double _element_count, polynomial_term, _monoisotopic_mass
            double base_intensity
            PyObject* pk
            PyObject* pv

        mass_vector = []
        max_variant_count = max_variants(self.composition)

        base_intensity = self.monoisotopic_peak.intensity
        ele_sym_poly_map = dict()
        composition_elements = PyDict_Keys(self.composition)

        for j in range(PyList_GET_SIZE(composition_elements)):
            element = <str>PyList_GET_ITEM(composition_elements, j)
            if element == "H+":
                continue
            power_sum = self.modified_phi_values(element)
            ele_sym_poly = []
            newton(power_sum, ele_sym_poly, max_variant_count)
            ele_sym_poly_map[element] = ele_sym_poly
        for i in range(self._order + 1):
            sign = 1 if i % 2 == 0 else -1
            center = 0.0
            k = 0
            while(PyDict_Next(ele_sym_poly_map, &k, &pk, &pv)):

            #for element, ele_sym_poly in ele_sym_poly_map.items():
                element = <str>pk
                ele_sym_poly = <list>pv

                _element_count = PyFloat_AsDouble(<object>PyDict_GetItem(self.composition, element))

                polynomial_term = PyFloat_AsDouble(<object>PyList_GET_ITEM(ele_sym_poly, i))

                phi_obj = <PhiConstants>PyDict_GetItem(self._isotopic_constants, element)
                element_obj = phi_obj.element
                # element_obj = (<Element>PyDict_GetItem(periodic_table, element))
                _monoisotopic_mass = element_obj._monoisotopic_mass

                temp = _element_count
                temp *= sign * polynomial_term
                temp *= base_intensity * _monoisotopic_mass
                center += temp

            mass_vector.append((center / probability_vector[i]) if probability_vector[i] > 0 else 0)
        return mass_vector

    def aggregated_isotopic_variants(self, int charge=0, double charge_carrier=PROTON):
        '''
        Compute the m/z (or neutral mass when `charge` == 0) for each
        aggregated isotopic peak and their intensity relative to
        the monoisotopic peak.
        '''
        cdef:
            list probability_vector
            list center_mass_vector
            list peak_set
            double average_mass, adjusted_mz
            double total
            size_t i
            Peak peak
            double center_mass_i, intensity_i
        probability_vector = self.probability()
        center_mass_vector = self.center_mass(probability_vector)

        peak_set = []
        average_mass = 0.
        total = sum(probability_vector)

        for i in range(self.order + 1):
            center_mass_i = PyFloat_AsDouble(<object>PyList_GET_ITEM(center_mass_vector, i))
            if charge != 0:
                adjusted_mz = mass_charge_ratio(center_mass_i, charge, charge_carrier)
            else:
                adjusted_mz = center_mass_i

            intensity_i = PyFloat_AsDouble(<object>PyList_GET_ITEM(probability_vector, i))

            peak = Peak(adjusted_mz, intensity_i / total, charge)
            if peak.intensity < 0:
                continue
            peak_set.append(peak)
            average_mass += adjusted_mz * intensity_i

        average_mass /= total
        self.average_mass = average_mass
        peak_set.sort(key=mz_getter)
        return tuple(peak_set)
