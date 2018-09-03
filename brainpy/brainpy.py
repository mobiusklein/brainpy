from __future__ import absolute_import

import operator
import re

from collections import OrderedDict, Counter
from math import exp, log, sqrt
from sys import float_info

from brainpy.mass_dict import nist_mass
from brainpy.composition import (
    PyComposition,
    parse_formula,
    calculate_mass,
    _make_isotope_string,
    _get_isotope)

mz_getter = operator.attrgetter("mz")
PROTON = nist_mass["H+"][0][0]
MACHINE_EPSILON = float_info.epsilon
ZERO = 0.
ONE = 1.0


def neutral_mass(mz, z, charge_carrier=PROTON):
    return (mz * abs(z)) - (z * charge_carrier)


def mass_charge_ratio(neutral_mass, z, charge_carrier=PROTON):
    return (neutral_mass + (z * charge_carrier)) / abs(z)


def give_repr(cls):  # pragma: no cover
    r"""Patch a class to give it a generic __repr__ method
    that works by inspecting the instance dictionary.

    Parameters
    ----------
    cls: type
        The class to add a generic __repr__ to.

    Returns
    -------
    cls: type
        The passed class is returned
    """
    def reprer(self):
        attribs = ', '.join(["%s=%r" % (k, v) for k, v in self.__dict__.items() if not k.startswith("_")])
        wrap = "{self.__class__.__name__}({attribs})".format(self=self, attribs=attribs)
        return wrap
    cls.__repr__ = reprer
    return cls


@give_repr
class PolynomialParameters(object):
    def __init__(self, elementary_symmetric_polynomial, power_sum):
        self.elementary_symmetric_polynomial = elementary_symmetric_polynomial
        self.power_sum = power_sum

    def __iter__(self):
        yield self.power_sum
        yield self.elementary_symmetric_polynomial


@give_repr
class PhiConstants(object):
    def __init__(self, order, element, element_coefficients, mass_coefficients):
        self.order = order
        self.element = element
        self.element_coefficients = element_coefficients
        self.mass_coefficients = mass_coefficients


def newton(power_sum, elementary_symmetric_polynomial, order):
    r'''
    Given two lists of values, the first list being the `power sum`s of a
    polynomial, and the second being expressions of the roots of the
    polynomial as found by Viete's Formula, use information from the longer list to
    fill out the shorter list using Newton's Identities.

    .. note::
        Updates are done **in place**

    Parameters
    ----------
    power_sum: list of float
    elementary_symmetric_polynomial: list of float
    order: int
        The number of terms to expand to when updating `elementary_symmetric_polynomial`

    See Also
    --------
    https://en.wikipedia.org/wiki/Newton%27s_identities[https://en.wikipedia.org/wiki/Newton%27s_identities]
    '''
    if len(power_sum) > len(elementary_symmetric_polynomial):
        _update_elementary_symmetric_polynomial(power_sum, elementary_symmetric_polynomial, order)
    elif len(power_sum) < len(elementary_symmetric_polynomial):
        _update_power_sum(power_sum, elementary_symmetric_polynomial, order)


def _update_elementary_symmetric_polynomial(power_sum, elementary_symmetric_polynomial, order):
    begin = len(elementary_symmetric_polynomial)
    end = len(power_sum)
    for k in range(begin, end):
        if k == 0:
            elementary_symmetric_polynomial.append(1.0)
        elif k > order:
            elementary_symmetric_polynomial.append(0.)
        else:
            el = 0.
            for j in range(1, k + 1):
                sign = 1 if (j % 2) == 1 else -1
                el += sign * power_sum[j] * elementary_symmetric_polynomial[k - j]
            el /= float(k)
            elementary_symmetric_polynomial.append(el)


def _update_power_sum(ps_vec, esp_vec, order):
    begin = len(ps_vec)
    end = len(esp_vec)
    for k in range(begin, end):
        if k == 0:
            ps_vec.append(0.)
            continue
        temp_ps = 0.
        sign = -1
        for j in range(1, k):
            sign *= -1
            temp_ps += sign * esp_vec[j] * ps_vec[k - j]
        sign *= -1
        temp_ps += sign * esp_vec[k] * k
        ps_vec.append(temp_ps)


def vietes(coefficients):
    r'''
    Given the coefficients of a polynomial of a single variable,
    compute an elementary symmetric polynomial of the roots of the
    input polynomial by Viete's Formula:

    .. math::
        \sum_{1\le i_1<i_2<...<i_k \le n} x_{i_1}x_{i_2}...x_{i_k} = (-1)^k\frac{a_{n-k}}{a_n}

    Parameters
    ----------
    coefficients: list of float
        A list of coefficients of the input polynomial of arbitrary length (degree) `n`

    Returns
    -------
    list of float:
        A list containing the expression of the roots of the input polynomial of length `n`

    See Also
    --------
    https://en.wikipedia.org/wiki/Vieta%27s_formulas

    '''
    elementary_symmetric_polynomial = []
    tail = float(coefficients[-1])
    size = len(coefficients)

    for i in range(size):
        sign = 1 if (i % 2) == 0 else -1
        el = sign * coefficients[size - i - 1] / tail
        elementary_symmetric_polynomial.append(el)
    return elementary_symmetric_polynomial


@give_repr
class Isotope(object):
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
    def __init__(self, mass, abundance, neutron_shift, neutrons):
        self.mass = mass
        self.abundance = abundance
        self.neutrons = neutrons
        self.neutron_shift = neutron_shift


@give_repr
class Element(object):
    def __init__(self, symbol, isotopes=None):
        self.symbol = symbol
        self.isotopes = isotopes or _isotopes_of(symbol)

    def average_mass(self):  # pragma: no cover
        mass = 0.0
        for i, isotope in self.isotopes.items():
            mass += isotope.mass * isotope.abundance
        return mass

    def monoisotopic_mass(self):
        return self.isotopes[0].mass

    def __iter__(self):
        return iter(self.isotopes.values())

    def max_neutron_shift(self):
        if len(self.isotopes) == 0:
            return 0
        return list(self.isotopes.values())[-1].neutron_shift

    def min_neutron_shift(self):
        if len(self.isotopes) == 0:
            return -1
        return list(self.isotopes.values())[0].neutron_shift

    def __eq__(self, other):
        if other is None:
            return False
        if isinstance(other, str):
            return self.symbol == other
        else:
            return self.symbol == other.symbol

    def __ne__(self, other):
        return not self == other

    def __hash__(self):
        return hash(self.symbol)


def make_fixed_isotope_element(element, neutrons):
    isotope = element.isotopes[neutrons - element.isotopes[0].neutrons]
    el = Element(element.symbol + ("[%d]" % neutrons), OrderedDict([
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
        return OrderedDict()
    mono_neutrons = max(freqs.items(), key=lambda x: x[1][1])[0]
    freqs = OrderedDict(sorted(((k - mono_neutrons, Isotope(*v, neutron_shift=k - mono_neutrons, neutrons=k))
                                for k, v in freqs.items()), key=lambda x: x[0]))
    return freqs


periodic_table = {k: Element(k) for k in nist_mass}
# periodic_table = {k: e for k, e in periodic_table.items() if e.max_neutron_shift() != 0 and e.min_neutron_shift() >= 0}
periodic_table["H+"] = Element("H+")


class IsotopicConstants(dict):
    def __init__(self, order):
        self._order = 0
        self.order = order

    @property
    def order(self):
        return self._order

    @order.setter
    def order(self, value):
        self._order = value
        self.update_coefficients()

    def coefficients(self, element, with_mass=False):
        max_isotope_number = element.max_neutron_shift()
        accumulator = []
        for isotope in reversed(list(element)):
            current_order = max_isotope_number - isotope.neutron_shift
            if with_mass:
                coef = isotope.mass
            else:
                coef = 1.

            if current_order > len(accumulator):
                for i in range(len(accumulator), current_order):
                    accumulator.append(0.)
                accumulator.append(isotope.abundance * coef)
            elif current_order == len(accumulator):
                accumulator.append(isotope.abundance * coef)
            else:
                raise Exception("The list of neutron shifts is not ordered.")

        elementary_symmetric_polynomial = vietes(accumulator)
        power_sum = []
        newton(power_sum, elementary_symmetric_polynomial, len(accumulator) - 1)
        return PolynomialParameters(elementary_symmetric_polynomial, power_sum)

    def add_element(self, symbol):
        if symbol in self:
            return
        if symbol in periodic_table:
            element = periodic_table[symbol]
        else:
            symbol_parsed, isotope = _get_isotope(symbol)
            if isotope == 0:
                raise KeyError(symbol)
            element = make_fixed_isotope_element(periodic_table[symbol_parsed], isotope)
        order = element.max_neutron_shift()
        element_parameters = self.coefficients(element)
        mass_parameters = self.coefficients(element, True)
        self[symbol] = PhiConstants(order, element, element_parameters, mass_parameters)

    def update_coefficients(self):
        for symbol, phi_constants in self.items():
            if self.order < phi_constants.order:
                continue

            for i in range(phi_constants.order, self.order + 1):
                phi_constants.element_coefficients.elementary_symmetric_polynomial.append(0.)
                phi_constants.mass_coefficients.elementary_symmetric_polynomial.append(0.)

            phi_constants.order = len(phi_constants.element_coefficients.elementary_symmetric_polynomial)
            newton(*phi_constants.element_coefficients, order=phi_constants.order)
            newton(*phi_constants.mass_coefficients, order=phi_constants.order)

    def nth_element_power_sum(self, symbol, order):
        constants = self[symbol]
        return constants.element_coefficients.power_sum[order]

    def nth_modified_element_power_sum(self, symbol, order):
        constants = self[symbol]
        return constants.mass_coefficients.power_sum[order]


def max_variants(composition):
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


@give_repr
class Peak(object):
    """
    Represent a single theoretical peak centroid.

    Peaks are comparable, hashable, and can be copied by calling
    :meth:`clone`

    Attributes
    ----------
    charge : int
        The charge state of the peak
    intensity : float
        The height of the peak. Peaks created as part of a
        theoretical isotopic cluster will be have an intensity
        between 0 and 1.
    mz : float
        The mass-to-charge ratio of the peak
    """
    def __init__(self, mz, intensity, charge):
        self.mz = mz
        self.intensity = intensity
        self.charge = charge

    def __eq__(self, other):  # pragma: no cover
        equal = all(
            abs(self.mz - other.mz) < 1e-10,
            abs(self.intensity - other.intensity) < 1e-10,
            self.charge == other.charge)
        return equal

    def __ne__(self, other):  # pragma: no cover
        return not (self == other)

    def __hash__(self):  # pragma: no cover
        return hash(self.mz)

    def clone(self):  # pragma: no cover
        return self.__class__(self.mz, self.intensity, self.charge)


@give_repr
class IsotopicDistribution(object):
    """
    Constructs a theoretical isotopic distribution for a given composition
    out to a given number of peaks.

    Attributes
    ----------
    average_mass : float
        The average (weighted) mass of the resulting isotopic cluster
    composition : dict
        The composition to create the isotopic cluster for
    order : int
        The number of peaks to produce and the number of terms in the
        generating polynomial expression.
    """
    def __init__(self, composition, order=-1):
        self.composition = composition
        self._isotopic_constants = IsotopicConstants(order)
        self._order = 0
        self.order = order
        self.average_mass = 0.
        self.monoisotopic_peak = self._create_monoisotopic_peak()

    @property
    def order(self):
        return self._order

    @order.setter
    def order(self, value):
        max_variant_count = max_variants(self.composition)
        if value == -1:
            self._order = max_variant_count
        else:
            self._order = min(value, max_variant_count)
        self._update_isotopic_constants()

    def _update_isotopic_constants(self):
        for element in self.composition:
            self._isotopic_constants.add_element(element)
        self._isotopic_constants.order = self._order

    def _create_monoisotopic_peak(self):
        mass = calculate_mass(self.composition)
        intensity = 0.
        for element in self.composition:
            if element == "H+":
                continue
            # intensity += log(periodic_table[element].isotopes[0].abundance)
            intensity += log(self._isotopic_constants[element].element.isotopes[0].abundance)
        intensity = exp(intensity)
        return Peak(mass, intensity, 0)

    def _phi_value(self, order):
        phi = 0.
        for element, count in self.composition.items():
            if element == "H+":
                continue
            phi += self._isotopic_constants.nth_element_power_sum(element, order) * count
        return phi

    def _modified_phi_value(self, symbol, order):
        phi = 0.
        for element, count in self.composition.items():
            if element == "H+":
                continue
            # Count is one lower for this symbol because an isotope is present
            # accounted for in the call to `nth_modified_element_power_sum` at
            # the end?
            coef = count if element != symbol else count - 1
            phi += self._isotopic_constants.nth_element_power_sum(element, order) * coef

        phi += self._isotopic_constants.nth_modified_element_power_sum(symbol, order)
        return phi

    def phi_values(self):
        power_sum = [0.]
        for i in range(1, self.order + 1):
            power_sum.append(self._phi_value(i))
        return power_sum

    def modified_phi_values(self, symbol):
        power_sum = [0.]
        for i in range(1, self.order + 1):
            power_sum.append(self._modified_phi_value(symbol, i))
        return power_sum

    def probability(self):
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

    def center_mass(self, probability_vector):
        mass_vector = []
        max_variant_count = max_variants(self.composition)

        ele_sym_poly_map = dict()
        for element in self.composition:
            if element == "H+":
                continue
            power_sum = self.modified_phi_values(element)
            ele_sym_poly = []
            newton(power_sum, ele_sym_poly, max_variant_count)
            ele_sym_poly_map[element] = ele_sym_poly
        for i in range(self.order + 1):
            sign = 1 if i % 2 == 0 else -1
            center = 0.0
            for element, ele_sym_poly in ele_sym_poly_map.items():
                center += self.composition[element] * sign * ele_sym_poly[i] *\
                 self.monoisotopic_peak.intensity * self._isotopic_constants[element].element.monoisotopic_mass()
                 # self.monoisotopic_peak.intensity * periodic_table[element].monoisotopic_mass()
            mass_vector.append((center / probability_vector[i]) if probability_vector[i] > 0 else 0)
        return mass_vector

    def aggregated_isotopic_variants(self, charge=0, charge_carrier=PROTON):
        '''
        Compute the m/z (or neutral mass when `charge` == 0) for each
        aggregated isotopic peak and their intensity relative to
        the monoisotopic peak.

        Parameters
        ----------
        charge: int
            The charge state of the resulting theoretical isotopic cluster
        charge_carrier: float
            The mass added for each degree of charge

        Returns
        -------
        theoretical_isotopic_distribution: list
            A list of :class:`Peak` objects whose intensities are proportional to
            each other to reflect relative peak heights.

        '''
        probability_vector = self.probability()
        center_mass_vector = self.center_mass(probability_vector)

        peak_set = []
        average_mass = 0.
        total = sum(probability_vector)

        for i in range(self.order + 1):
            if charge != 0:
                adjusted_mz = mass_charge_ratio(center_mass_vector[i], charge, charge_carrier)
            else:
                adjusted_mz = center_mass_vector[i]
            if adjusted_mz < 1:
                continue
            peak = Peak(adjusted_mz, probability_vector[i] / total, charge)
            if peak.intensity < 0:
                continue
            peak_set.append(peak)
            average_mass += adjusted_mz * probability_vector[i]

        average_mass /= total
        self.average_mass = average_mass
        peak_set.sort(key=mz_getter)
        return tuple(peak_set)


def isotopic_variants(composition, npeaks=None, charge=0, charge_carrier=PROTON):
    '''
    Compute a peak list representing the theoretical isotopic cluster for `composition`.

    Parameters
    ----------
    composition : Mapping
        Any Mapping type where keys are element symbols and values are integers
    npeaks: int
        The number of peaks to include in the isotopic cluster, starting from the monoisotopic peak.
        If given a number below 1 or above the maximum number of isotopic variants, the maximum will
        be used. If `None`, a "reasonable" value is chosen by `int(sqrt(max_variants(composition)))`.
    charge: int
        The charge state of the isotopic cluster to produce. Defaults to 0, theoretical neutral mass.

    Returns
    -------
    list of Peaks

    See Also
    --------
    :class:`IsotopicDistribution`

    '''
    if npeaks is None:
        max_n_variants = max_variants(composition)
        npeaks = int(sqrt(max_n_variants) - 2)
        npeaks = max(npeaks, 3)
    else:
        # Monoisotopic Peak is not included
        npeaks -= 1
    return IsotopicDistribution(composition, npeaks).aggregated_isotopic_variants(
        charge, charge_carrier=charge_carrier)


try:
    _has_c = True
    _IsotopicDistribution = IsotopicDistribution
    from ._speedup import IsotopicDistribution
    from ._c.isotopic_distribution import pyisotopic_variants as isotopic_variants

except ImportError:
    _has_c = False
