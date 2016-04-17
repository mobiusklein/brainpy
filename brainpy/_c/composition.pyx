cimport cython
from cpython.list cimport PyList_New, PyList_SET_ITEM
from cpython.string cimport PyString_FromString, PyString_AsString
from cpython.int cimport PyInt_AsLong, PyInt_FromLong
from libc.string cimport strcmp, memcpy, strlen
from libc.stdlib cimport malloc, free, realloc, atoi
from libc.math cimport abs, fabs
from libc cimport *

cdef extern from * nogil:
    int printf (const char *template, ...)


from brainpy.mass_dict import nist_mass as __nist_mass

cdef dict nist_mass
nist_mass = __nist_mass


cdef double PROTON = nist_mass["H+"][0][0]

cdef double neutral_mass(double mz,  int z, double charge_carrier=PROTON) nogil:
    return (mz * fabs(z)) - (z * charge_carrier)

cdef double mass_charge_ratio(double neutral_mass, int z, double charge_carrier=PROTON) nogil:
    return (neutral_mass + (z * charge_carrier)) / fabs(z)

cdef IsotopeMap* make_isotope_map(list organized_isotope_data, size_t size):
    cdef:
        IsotopeMap* result
        size_t i
        int k
        tuple values

    result = <IsotopeMap*>malloc(sizeof(IsotopeMap))
    result.bins = <Isotope*>malloc(sizeof(Isotope) * size)
    result.size = size

    i = 0
    for k, values in organized_isotope_data:
        result.bins[i].mass = values[0]
        result.bins[i].abundance = values[1]
        result.bins[i].neutron_shift = values[2]
        result.bins[i].neutrons = values[3]
        i += 1
    return result

cdef Isotope* get_isotope_by_neutron_shift(IsotopeMap* isotopes, int neutron_shift) nogil:
    cdef:
        size_t i
        Isotope* isotope_item

    for i in range(isotopes.size):
        isotope_item = &(isotopes.bins[i])
        if isotope_item.neutron_shift == neutron_shift:
            return isotope_item
    return NULL

cdef Isotope* get_isotope_by_neutron_count(IsotopeMap* isotopes, int neutrons) nogil:
    cdef:
        size_t i
        Isotope* isotope_item

    for i in range(isotopes.size):
        isotope_item = &(isotopes.bins[i])
        if isotope_item.neutrons == neutrons:
            return isotope_item
    return NULL    

cdef void print_isotope_map(IsotopeMap* isotope_map):
    cdef:
        size_t i

    for i in range(isotope_map.size):
        printf("%f, %f, %d\n", isotope_map.bins[i].mass, isotope_map.bins[i].abundance, isotope_map.bins[i].neutron_shift)


cdef void free_isotope_map(IsotopeMap* isotopes) nogil:
    free(isotopes.bins)
    free(isotopes)

cdef struct Element:
    char* symbol
    IsotopeMap* isotopes

cdef double element_monoisotopic_mass(Element* element) nogil:
    return get_isotope_by_neutron_shift(element.isotopes, 0).mass

cdef double element_isotopic_mass(Element* element, int isotope_number) nogil:
    return get_isotope_by_neutron_count(element.isotopes, isotope_number).mass

cdef int element_min_neutron_shift(Element* element) nogil:
    return element.isotopes.bins[0].neutron_shift

cdef int element_max_neutron_shift(Element* element) nogil:
    return element.isotopes.bins[element.isotopes.size - 1].neutron_shift

def __select_1_1(object x):
    return x[1][1]

cdef void _isotopes_of(char* element_symbol, IsotopeMap** isotope_frequencies):
    cdef:
        dict freqs
        list freq_list
        int i, key
        tuple mass, value
        double abundance
        list bunches
        tuple mass_freqs
        int mono_neutrons
        object k, v

    freqs = dict()
    for i, mass_freqs in nist_mass[element_symbol].items():
        if i == 0:
            continue
        if mass_freqs[1] > 0:
            freqs[i] = mass_freqs

    if len(freqs) == 0:
        return

    freq_list = freqs.items() 
    freq_list.sort(key=__select_1_1)
    mono_neutrons = freq_list[-1][0]

    bunches = []
    for k, v in freqs.items():
        bunches.append((k - mono_neutrons, (v[0], v[1], k - mono_neutrons, k)))

    bunches.sort()

    # bunches = sorted([((k - mono_neutrons, (v[0], v[1], k - mono_neutrons, k))
    #                   for k, v in freqs.items())], key=lambda x: x[0])
    isotope_frequencies[0] = make_isotope_map(bunches, len(bunches))

cdef Element* make_element(char* symbol):
    cdef:
        Element* element
    element = <Element*>malloc(sizeof(Element))
    element.symbol = symbol    
    _isotopes_of(symbol, &element.isotopes)
    return element

cdef void make_element_in_place(char* symbol, Element* element):
    element.symbol = symbol
    _isotopes_of(symbol, &element.isotopes)

cdef void free_element(Element* element) nogil:
    free_isotope_map(element.isotopes)
    free(element)

cdef char* _parse_isotope_string(char* label, int* isotope_num, char* element_name) nogil:
    cdef:
        size_t i = 0
        bint in_bracket = False
        char current
        size_t name_end
        size_t num_start
        size_t num_end

        size_t size
        char[4] number_part

    name_end = 0
    num_start = 0
    num_end = 0
    size = strlen(label)
    for i in range(size):
        current = label[i]
        if in_bracket:
            if current == 93:
                break
            num_end += 1
        elif current == 91:
            in_bracket = True
            name_end = i
            num_start = i + 1
            num_end = num_start
        else:
            name_end += 1

    if num_start > 0:
        memcpy(number_part, label + num_start, num_end - num_start)
        number_part[num_end - num_start] = '\0'
        isotope_num[0] = atoi(number_part)
    else:
        isotope_num[0] = 0

    element_name[name_end] = '\0'
    memcpy(element_name, label, name_end)

    return element_name

cdef size_t hash_string(char *str) nogil:
    cdef:    
        size_t hash
        size_t i
        int c
    hash = 5381
    i = 0
    c = str[i]
    while (c):
        hash = ((hash << 5) + hash) + c
        i += 1
        c = str[i]
    return hash;

cdef struct PeriodicTable:
    Element** elements
    size_t size

cdef PeriodicTable* make_periodic_table():
    cdef: 
        Element* element 
        PeriodicTable* table 
        str pk
        char* k 
        set used 
        size_t i

    used = set()
    table = <PeriodicTable*>malloc(sizeof(PeriodicTable))
    table.elements = <Element**>malloc(len(nist_mass) * 2 * sizeof(Element*))
    table.size = len(nist_mass) * 2

    for pk in nist_mass:
        k = PyString_AsString(pk)
        i = hash_string(k) % table.size
        if i in used:
            while True:
                i += 1
                if i >= table.size:
                    i = 0
                if i not in used:
                    break

        table.elements[i] = make_element(k)
        used.add(i)
    return table

cdef int get_element_from_periodic_table(PeriodicTable* table, char* symbol, int* out) nogil:
    cdef:
        size_t i, j
    j = 0
    i = hash_string(symbol) % table.size
    while strcmp(table.elements[i].symbol, symbol) != 0:
        i += 1
        j += 1
        if i >= table.size:
            i = 0
        if j > table.size:
            return -1
    out[0] = i
    return 0

cdef int get_element_from_periodic_table2(PeriodicTable* table, char* symbol, Element** out) nogil:
    cdef:
        size_t i, j
    j = 0
    i = hash_string(symbol) % table.size
    while strcmp(table.elements[i].symbol, symbol) != 0:
        i += 1
        j += 1
        if i >= table.size:
            i = 0
        if j > table.size:
            return -1

    out[0] = table.elements[i]
    return 0

cdef PeriodicTable* _PeriodicTable
_PeriodicTable = make_periodic_table()

cdef struct Composition:
    char** elements
    double* counts
    size_t size
    size_t used

cdef Composition* make_composition() nogil:
    cdef:
        Composition* composition
    composition = <Composition*>malloc(sizeof(Composition))
    composition.elements = <char**>malloc(sizeof(char*) * 7)
    composition.counts = <count_type*>malloc(sizeof(count_type) * 7)
    composition.size = 7
    composition.used = 0
    return composition

cdef Composition* copy_composition(Composition* composition) nogil:
    cdef:
        Composition* result
        int status
        size_t i
    result = <Composition*>malloc(sizeof(Composition))
    result.elements = <char**>malloc(sizeof(char*) * composition.size)
    result.counts = <count_type*>malloc(sizeof(count_type) * composition.size)
    result.size = composition.size
    result.used = 0

    i = 0
    while i < composition.used:
        status = composition_set_element_count(result, composition.elements[i], composition.counts[i])
        i += 1
    return result

cdef void print_composition(Composition* composition) nogil:
    cdef:
        size_t i
        int status
    i = 0
    printf("Addr %d %d\n", <int>composition, <int>(composition.elements))
    printf("Size: %d, Used: %d\n", composition.size, composition.used)
    printf("{\n")
    while i < composition.used:
        printf("  %s: %d\n", composition.elements[i], composition.counts[i])
        i += 1
    printf("}\n\n")

cdef int composition_set_element_count(Composition* composition, char* element, count_type count) nogil:
    cdef:
        size_t i
        int status
        bint done
    done = False
    i = 0
    while (i < composition.used):
        if strcmp(element, composition.elements[i]) == 0:
            done = True
            composition.counts[i] = count
            return 0
        i += 1

    if not done:
        composition.used += 1
        if composition.used >= composition.size:
            status = composition_resize(composition)
            if status != 0:
                return -1
        composition.elements[i] = element
        composition.counts[i] = count
        return 0
    return 1

cdef int composition_get_element_count(Composition* composition, char* element, count_type* count) nogil:
    cdef:
        size_t i
        int status
        bint done
    done = False
    i = 0
    while (i < composition.used):
        if strcmp(element, composition.elements[i]) == 0:
            done = True
            count[0] = composition.counts[i]
            return 0
        i += 1

    if not done:
        count[0] = 0
    return 0

cdef int composition_inc_element_count(Composition* composition, char* element, count_type increment) nogil:
    cdef:
        size_t i
        int status
        bint done
    done = False
    i = 0
    while (i < composition.used):
        if strcmp(element, composition.elements[i]) == 0:
            done = True
            composition.counts[i] += increment
            return 0
        i += 1

    if not done:
        composition.used += 1
        if composition.used >= composition.size:
            status = composition_resize(composition)
            if status != 0:
                return -1
        composition.elements[i] = element
        composition.counts[i] = increment
        return 0
    return 1

cdef int composition_resize(Composition* composition) nogil:
    composition.elements = <char**>realloc(composition.elements, sizeof(char*) * composition.size * 2)
    composition.counts = <count_type*>realloc(composition.counts, sizeof(count_type) * composition.size * 2)
    composition.size *= 2
    if composition.counts == NULL:
        return -1
    return 0

cdef double composition_mass(Composition* composition) nogil:
    cdef:
        double mass
        Element* element
        char* element_label
        char[10] element_name
        int isotope_number
        int index
        size_t i
    i = 0
    while i < composition.used:
        element_label = composition.elements[i]
        _parse_isotope_string(element_label, &isotope_number, element_name)
        index = get_element_from_periodic_table2(_PeriodicTable, element_name, &element)
        if isotope_number == 0:
            mass += element_monoisotopic_mass(element) * composition.counts[i]
        else:
            mass += element_isotopic_mass(element, isotope_number) * composition.counts[i]
        i += 1

    return mass

cdef void free_composition(Composition* composition) nogil:
    free(composition.elements)
    free(composition.counts)
    free(composition)

cdef Composition* composition_add(Composition* composition_1, Composition* composition_2, int sign) nogil:
    cdef:
        Composition* result
        count_type value
        char* symbol
        int status
        size_t i

    i = 0
    result = copy_composition(composition_1)

    while i < composition_2.used:
        symbol = composition_2.elements[i]
        status = composition_get_element_count(composition_2, symbol, &value)
        i += 1
        if value == 0:
            continue
        status = composition_inc_element_count(result, symbol, sign * value)
    return result

cdef int composition_iadd(Composition* composition_1, Composition* composition_2, int sign) nogil:
    cdef:
        char* symbol
        count_type value
        int status
        size_t i
    i = 0
    status = 0
    while i < composition_2.used:
        symbol = composition_2.elements[i]
        status = composition_get_element_count(composition_2, symbol, &value)
        i += 1
        if value == 0:
            continue
        status = composition_inc_element_count(composition_1, symbol, sign * value)
        if status != 0:
            pass
    return status

cdef Composition* composition_mul(Composition* composition, int scale) nogil:
    cdef:
        Composition* result
        size_t i

    i = 0
    result = copy_composition(composition)
    while i < result.used:
        result.counts[i] *= scale
        i += 1
    return result

cdef void composition_imul(Composition* composition, int scale) nogil:
    cdef:
        size_t i

    i = 0
    while i < composition.used:
        composition.counts[i] *= scale
        i += 1

cdef dict composition_to_dict(Composition* composition):
    cdef:
        dict result
        str symbol
        size_t i
        char* symbol_c
        count_type value
    i = 0
    result = dict()
    while i < composition.used:
        symbol_c = composition.elements[i]
        composition_get_element_count(composition, symbol_c, &value)
        if value == 0:
            continue
        symbol = PyString_FromString(symbol_c)
        result[symbol] = value
        i += 1
    return result

cdef Composition* dict_to_composition(dict comp_dict):
    cdef:
        Composition* result
        str symbol
        size_t i
        char* symbol_c
        count_type value
    result = make_composition()
    for symbol, value in comp_dict.items():
        symbol_c = PyString_AsString(symbol)
        composition_set_element_count(result, symbol_c, value)
    return result

def main():
    cdef Element* elem
    cdef Composition* composition
    cdef Composition* sum_composition
    cdef int i, j
    cdef char *sym
    cdef char *k
    cdef char[10] element_name_buffer
    cdef dict copy

    sym = "O"    
    # i = get_element_from_periodic_table(_PeriodicTable, sym, &j)
    # if i == 0:
    #     elem = _PeriodicTable.elements[j]
    #     print element_monoisotopic_mass(elem), elem.symbol
    # else:
    #     print "Could not locate element"

    j = get_element_from_periodic_table2(_PeriodicTable, "C", &elem)
    # print elem.symbol, element_monoisotopic_mass(elem)

    composition = make_composition()
    composition_set_element_count(composition, "O[18]", 1)
    composition_set_element_count(composition, "H", 2)
    print_composition(composition)
    print "Massing"
    print composition_mass(composition)

    sum_composition = composition_add(composition, composition, 1)
    print "Summed"
    print composition_mass(sum_composition)
    print composition_mass(composition)
    composition_iadd(sum_composition, composition, 1)
    print composition_mass(sum_composition)

    # copy = composition_to_dict(sum_composition)
    # free_composition(sum_composition)
    # print(copy)
    # sum_composition = dict_to_composition(copy)
    # print_composition(composition_mul(sum_composition, 3))
    # k = "C[13]"
    # _parse_isotope_string(k, &j, element_name_buffer)
    # print element_name_buffer, j
    # print elem.symbol
    # print element_monoisotopic_mass(elem)
    # print element_isotopic_mass(elem, j)
    # print element_isotopic_mass(elem, 12)
    # print_isotope_map(elem.isotopes) 


cdef class PyComposition(object):
    def __cinit__(self, base=None, **kwargs):
        cdef:
            char* element
            count_type count
        self.impl = make_composition()
        self._clean = False
        self.cached_mass = 0.

        if base is not None and isinstance(base, dict):
            self.update(base)

        for element, count in kwargs.items():
            self[element] = count

    def __dealloc__(self):
        free_composition(self.impl)

    def __getitem__(self, str key):
        cdef:
            count_type count
            int status
            char* ckey
        ckey = PyString_AsString(key)
        status = composition_get_element_count(self.impl, ckey, &count)
        if status == 0:
            return count
        else:
            intern(key)
            composition_set_element_count(self.impl, ckey, 0)
            return 0.

    def __setitem__(self, str key, count_type count):
        cdef:
            int status
            char* ckey
        ckey = PyString_AsString(key)
        intern(key)
        status = composition_set_element_count(self.impl, ckey, count)
        self._clean = False

    def update(self, arg=None, **kwargs):
        if arg is not None:
            self.update(**arg)
        for key, value in kwargs.items():
            self[key] = value

    def keys(self):
        cdef:
            size_t i
            char* elem
            list keys
        i = 0
        keys = PyList_New(self.impl.used)
        while i < self.impl.used:
            elem = self.impl.elements[i]
            PyList_SET_ITEM(keys, i, PyString_FromString(elem))
            i += 1
        return keys

    def values(self):
        cdef:
            size_t i
            count_type value
            list counts
        counts = []
        i = 0
        while i < self.impl.used:
            value = self.impl.counts[i]
            counts.append(value)
            i += 1
        return counts

    def items(self):
        cdef:
            size_t i
            char* elem
            count_type value
            list items

        i = 0
        items = []
        while i < self.impl.used:
            elem = self.impl.elements[i]
            value = self.impl.counts[i]
            items.append((PyString_FromString(elem), value))
            i += 1
        return items

    def copy(self):
        inst = PyComposition()
        composition_iadd(inst.impl, self.impl, 1)
        return inst

    def __add__(self, other):
        cdef:
            PyComposition result
            PyComposition _other

        result = self.copy()

        if isinstance(other, dict):
            for key, value in other.items():
                result[key] += value
        elif isinstance(other, PyComposition):
            _other = other
            composition_iadd(result.impl, _other.impl, 1)
        else:
            return NotImplemented
        return result

    def __iadd__(self, other):
        cdef:
            PyComposition _other

        if isinstance(other, dict):
            for key, value in other.items():
                self[key] += value
        elif isinstance(other, PyComposition):
            _other = other
            composition_iadd(self.impl, _other.impl, 1)
        else:
            raise TypeError("Cannot add %s to PyComposition" % type(other))
        return self

    def __sub__(self, other):
        cdef:
            PyComposition result
            PyComposition _other

        result = self.copy()

        if isinstance(other, dict):
            for key, value in other.items():
                result[key] -= value
        elif isinstance(other, PyComposition):
            _other = other
            composition_iadd(result.impl, _other.impl, -1)
        else:
            return NotImplemented
        return result

    def __isub__(self, other):
        cdef:
            PyComposition _other

        if isinstance(other, dict):
            for key, value in other.items():
                self[key] -= value
        elif isinstance(other, PyComposition):
            _other = other
            composition_iadd(self.impl, _other.impl, -1)
        else:
            raise TypeError("Cannot subtract %s from PyComposition" % type(other))
        return self


    def __mul__(self, scale):
        cdef:
            PyComposition result
            PyComposition inst
            int scale_factor

        if isinstance(self, PyComposition):
            inst = self
            scale_factor = PyInt_AsLong(scale)
        else:
            inst = scale
            scale_factor = PyInt_AsLong(self)

        result = inst.copy()

        composition_imul(result.impl, scale_factor)

        return result

    def __imul__(self, int scale):
        composition_imul(self.impl, scale)
        return self

    def __iter__(self):
        return iter(self.keys())

    cpdef double mass(self):
        if self._clean:
            return self.cached_mass
        else:
            self.cached_mass = composition_mass(self.impl)
            self._clean = True
            return self.cached_mass

    def __len__(self):
        return self.impl.used

    def __repr__(self):
        return "{%s}" % ', '.join("%s: %d" % kv for kv in self.items())
