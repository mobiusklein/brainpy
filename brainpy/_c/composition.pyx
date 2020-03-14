# cython: embedsignature=True

cimport cython

# from cpython.string cimport PyStr_FromString, PyStr_AsString
# from cpython.int cimport PyInt_AsLong, PyInt_FromLong

from brainpy._c.compat cimport (
    PyStr_FromString, PyStr_AsString,
    PyInt_AsLong, PyInt_FromLong, PyStr_InternInPlace,
    PyStr_AsUTF8AndSize)

try:
    from collections import Mapping, MutableMapping
except ImportError:
    from collections.abc import Mapping, MutableMapping

from cpython.ref cimport Py_INCREF, Py_DECREF
from cpython.list cimport PyList_New, PyList_Append, PyList_Append
from cpython.dict cimport PyDict_SetItem, PyDict_GetItem, PyDict_Next, PyDict_SetItemString
from cpython.object cimport PyObject
from libc.string cimport strcmp, memcpy, strlen, strncpy
from libc.stdlib cimport malloc, free, realloc, atoi, calloc
from libc.math cimport abs, fabs
from libc cimport *

cdef extern from * nogil:
    int printf (const char *template, ...)
    int sprintf(char *, char *,...)


from brainpy.mass_dict import nist_mass as __nist_mass

cdef dict nist_mass
nist_mass = __nist_mass


cdef double PROTON = nist_mass["H+"][0][0]


cdef double neutral_mass(double mz,  int z, double charge_carrier=PROTON) nogil:
    return (mz * fabs(z)) - (z * charge_carrier)


@cython.cdivision
cdef double mass_charge_ratio(double neutral_mass, int z, double charge_carrier=PROTON) nogil:
    return (neutral_mass + (z * charge_carrier)) / fabs(z)


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
    memcpy(element_name, label, name_end)
    element_name[name_end] = '\0'
    return element_name

cdef char* _make_isotope_string(Element* element, Isotope* isotope, char* out) nogil:
    if isotope.neutron_shift == 0:
        sprintf(out, "%s", element.symbol)
        return out
    else:
        sprintf(out, "%s[%d]", element.symbol, isotope.neutrons)
        return out

cdef char* _make_fixed_isotope_string(Element* element, Isotope* isotope, char* out) nogil:
    sprintf(out, "%s[%d]", element.symbol, isotope.neutrons)
    return out

# -----------------------------------------------------------------------------
# Isotope and IsotopeMap Methods

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
        # if isotope_item.neutron_shift == neutron_shift:
        if isotopes.bins[i].neutron_shift == neutron_shift:
            return isotope_item
    return NULL


cdef int get_isotope_by_neutron_shift_index(IsotopeMap* isotopes, int neutron_shift) nogil:
    cdef:
        size_t i

    for i in range(isotopes.size):
        if isotopes.bins[i].neutron_shift == neutron_shift:
            return i
    return -1


cdef Isotope* get_isotope_by_neutron_count(IsotopeMap* isotopes, int neutrons) nogil:
    cdef:
        size_t i
        Isotope* isotope_item

    for i in range(isotopes.size):
        isotope_item = &(isotopes.bins[i])
        if isotope_item.neutrons == neutrons:
            return isotope_item
    return NULL


cdef void print_isotope_map(IsotopeMap* isotope_map) nogil:
    cdef:
        size_t i

    for i in range(isotope_map.size):
        printf("%f, %f, %d -> %d\n", isotope_map.bins[i].mass, isotope_map.bins[i].abundance,
                                     isotope_map.bins[i].neutron_shift, isotope_map.bins[i].neutrons)


cdef void free_isotope_map(IsotopeMap* isotopes) nogil:
    free(isotopes.bins)
    free(isotopes)


# -----------------------------------------------------------------------------
# Element Methods

cdef double element_monoisotopic_mass(Element* element) nogil:
    if element.monoisotopic_isotope_index < 0:
        return get_isotope_by_neutron_shift(element.isotopes, 0).mass
    else:
        return element.isotopes.bins[element.monoisotopic_isotope_index].mass

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
        dict freqs, element_data
        list freq_list
        int i, key
        tuple mass, value
        double abundance
        list bunches
        tuple mass_freqs
        int mono_neutrons
        object k, v
        str py_element_symbol

    freqs = dict()
    py_element_symbol = PyStr_FromString(element_symbol)

    try:
        element_data = nist_mass[py_element_symbol]
    except KeyError:
        element_data = nist_mass[py_element_symbol.encode('utf-8')]

    for i, mass_freqs in element_data.items():
        if i == 0:
            continue
        if mass_freqs[1] > 0:
            freqs[i] = mass_freqs

    if len(freqs) == 0:
        isotope_frequencies[0] = make_isotope_map([], 0)
        return

    freq_list = list(freqs.items())
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
    element.monoisotopic_isotope_index = -1
    element.monoisotopic_isotope_index = get_isotope_by_neutron_shift_index(
        element.isotopes, 0)
    return element

cdef void free_element(Element* element) nogil:
    free_isotope_map(element.isotopes)
    free(element)

cdef void print_element(Element* element) nogil:
    printf("Symbol: %s; Monoisotopic Index: %d\n", element.symbol,
                                                   element.monoisotopic_isotope_index)
    print_isotope_map(element.isotopes)


cdef bint ensure_fixed_isotope(char* string):
    cdef:
        char* element_name
        int isotope_count
        int found
        Element* elem
        Element* fixed_isotope_elem
    found = element_hash_table_get(_ElementTable, string, &elem)
    if found != 0:
        elem = NULL
    else:
        return True

    isotope_count = 0
    element_name = <char*>malloc(sizeof(char) * 10)
    _parse_isotope_string(string, &isotope_count, element_name)
    found = element_hash_table_get(_ElementTable, element_name, &elem)
    if found != 0:
        # with gil:
            raise KeyError(string)
    fixed_isotope_elem = make_fixed_isotope_element(elem, isotope_count)
    if fixed_isotope_elem == NULL:
        printf("Could not ensure fixed isotope %s\n", string)
        raise KeyError(string)
    element_hash_table_put(_ElementTable, fixed_isotope_elem)
    return False


cdef Element* make_fixed_isotope_element(Element* element, int neutron_count) nogil:
    cdef:
        int i
        size_t j
        Element* out
        IsotopeMap* isotope_map
        Isotope* reference
        char* symbol
    # print_element(element)
    out = <Element*>malloc(sizeof(Element))
    isotope_map = <IsotopeMap*>malloc(sizeof(IsotopeMap))
    isotope_map.bins = <Isotope*>malloc(sizeof(Isotope) * 1)
    isotope_map.size = 1
    out.isotopes = isotope_map
    reference = get_isotope_by_neutron_count(element.isotopes, neutron_count)
    if reference == NULL:
        printf("reference was NULL!\n")
        return NULL
    isotope_map.bins[0].mass = reference.mass
    isotope_map.bins[0].abundance = 1.0
    isotope_map.bins[0].neutron_shift = 0
    isotope_map.bins[0].neutrons = neutron_count
    out.symbol = <char*>malloc(sizeof(char) * 10)
    _make_fixed_isotope_string(element, reference, out.symbol)
    out.monoisotopic_isotope_index = 0
    return out


def test_make_fixed_isotope_element(str element, int neutron_count):
    cdef:
        Element* elem_obj
        Element* fixed_obj
        char* c_str_element
    c_str_element = PyStr_AsString(element)
    printf("Element String: %s\n", c_str_element)
    status = element_hash_table_get(_ElementTable, c_str_element, &elem_obj)
    print(status)
    printf("Symbol: %s; Monoisotopic Index: %d\n", elem_obj.symbol, elem_obj.monoisotopic_isotope_index)
    print_isotope_map(elem_obj.isotopes)
    print("Creating Fixed Isotope")
    fixed_obj = make_fixed_isotope_element(elem_obj, neutron_count)
    printf("Symbol: %s; Monoisotopic Index: %d\n", fixed_obj.symbol, fixed_obj.monoisotopic_isotope_index)
    print_isotope_map(fixed_obj.isotopes)


# -----------------------------------------------------------------------------
# ElementHashTable and ElementHashBucket Methods

cdef size_t hash_string(char *string_) nogil:
    cdef:
        size_t hash_value
        size_t i
        int c
    hash_value = 5381
    i = 0
    c = string_[i]
    while (c):
        hash_value = ((hash_value << 5) + hash_value) + c
        i += 1
        c = string_[i]
    return hash_value;


def test_hash_string(str string_):
    return hash_string(PyStr_AsString(string_))


cdef ElementHashTable* make_element_hash_table(size_t size) nogil:
    cdef:
        ElementHashTable* table
        size_t i

    table = <ElementHashTable*>malloc(sizeof(ElementHashTable))
    table.buckets = <ElementHashBucket*>malloc(sizeof(ElementHashBucket) * size)
    table.size = size

    for i in range(size):
        table.buckets[i].size = 6
        table.buckets[i].elements = <Element**>malloc(sizeof(Element*) * table.buckets[i].size)
        table.buckets[i].used = 0
    return table


cdef int element_hash_bucket_resize(ElementHashBucket* bucket) nogil:
    cdef:
        Element** elements
        size_t new_size

    new_size = bucket.size * 2
    elements = <Element**>realloc(bucket.elements, sizeof(Element*) * new_size)
    if elements == NULL:
        printf("element_hash_bucket_resize failed\n")
        return -1
    bucket.elements = elements
    return 0


cdef int element_hash_bucket_insert(ElementHashBucket* bucket, Element* element) nogil:
    cdef:
        size_t i
        int status
        Element* el
    if (bucket.used + 1) == bucket.size:
        status = element_hash_bucket_resize(bucket)
        if status != 0:
            printf("element_hash_bucket_insert failed with %s\n", element.symbol)
            return -1
    # printf("Inserting %s into bucket with total capacity %d\n", element.symbol, bucket.size)
    bucket.elements[bucket.used] = element
    # printf("Element %s put in slot %d\n", bucket.elements[bucket.used].symbol, bucket.used)
    bucket.used += 1
    return 0


cdef int element_hash_bucket_find(ElementHashBucket* bucket, char* symbol, Element** out) nogil:
    cdef:
        size_t i
    for i in range(bucket.used):
        if strcmp(bucket.elements[i].symbol, symbol) == 0:
            out[0] = bucket.elements[i]
            return 0
    return -1


cdef void free_element_hash_bucket(ElementHashBucket* bucket) nogil:
    free(bucket.elements)
    free(bucket)


@cython.cdivision(True)
cdef int element_hash_table_get(ElementHashTable* table, char* symbol, Element** out) nogil:
    cdef:
        size_t hash_value
        size_t position
        ElementHashBucket bucket
        int status

    hash_value = hash_string(symbol)
    position = hash_value % table.size
    bucket = (table.buckets[position])
    status = element_hash_bucket_find(&bucket, symbol, out)
    return status


@cython.cdivision(True)
cdef int element_hash_table_put(ElementHashTable* table, Element* element) nogil:
    cdef:
        size_t hash_value
        size_t position
        ElementHashBucket bucket
        int status

    hash_value = hash_string(element.symbol)
    position = hash_value % table.size
    bucket = (table.buckets[position])
    status = element_hash_bucket_insert(&bucket, element)
    table.buckets[position] = bucket
    return status


cdef size_t free_element_hash_table(ElementHashTable* table) nogil:
    cdef:
        size_t i

    for i in range(table.size):
        free_element_hash_bucket(&(table.buckets[i]))
    free(table)


cdef ElementHashTable* make_element_hash_table_populated(size_t size):
    cdef:
        Element* element
        Element* out_test
        ElementHashTable* table
        str pk  #!
        char* k
        size_t i
        int status

    table = make_element_hash_table(size)

    for pk in nist_mass:
        k = PyStr_AsString(pk)
        element = make_element(k)
        status = element_hash_table_put(table, element)
        if status != 0:
            printf("element_hash_table_put exited with status %d\n", status)
        status = element_hash_table_get(table, k, &out_test)
        if status != 0:
            printf("element_hash_table_get exited with status %d\n", status)
    return table


_ElementTable = make_element_hash_table_populated(256)


cdef ElementHashTable* get_system_element_hash_table() nogil:
    return _ElementTable


cdef int set_system_element_hash_table(ElementHashTable* table) nogil:
    _ElementTable[0] = table[0]


def show_element(str element):
    cdef:
        Element* elem_obj
        char* c_str_element
    c_str_element = PyStr_AsString(element)
    printf("Element String: %s\n", c_str_element)
    status = element_hash_table_get(_ElementTable, c_str_element, &elem_obj)
    print(status)
    printf("Symbol: %s; Monoisotopic Index: %d\n", elem_obj.symbol, elem_obj.monoisotopic_isotope_index)
    print_isotope_map(elem_obj.isotopes)


# -----------------------------------------------------------------------------
# Composition Methods

cdef Composition* make_composition() nogil:
    '''
    Create a new, empty Composition struct
    '''
    cdef:
        Composition* composition
    composition = <Composition*>malloc(sizeof(Composition))
    composition.elements = <char**>calloc(7, sizeof(char*))
    composition.counts = <count_type*>calloc(7, sizeof(count_type))
    composition.size = 7
    composition.used = 0
    composition.max_variants = 0
    return composition

cdef int composition_eq(Composition* composition_1, Composition* composition_2) nogil:
    '''
    Test two Composition instances for element-wise equality
    '''
    cdef:
        int status
        size_t i
        char* symbol
        count_type value_1, value_2

    if composition_1.used != composition_2.used:
        return 0
    i = 0
    while i < composition_2.used:
        symbol = composition_2.elements[i]
        status = composition_get_element_count(composition_2, symbol, &value_2)
        if status != 0:
            return 0
        status = composition_get_element_count(composition_1, symbol, &value_1)
        if status != 0:
            return 0
        if value_1 != value_2:
            return 0
        i += 1

    return 1

cdef Composition* copy_composition(Composition* composition) nogil:
    '''
    Create a new Composition instance whose element counts are copied from
    `composition`
    '''
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


cdef int fill_composition_from_composition(Composition* result, Composition* composition) nogil:
    cdef:
        int status
        size_t i

    if result.size < composition.used:
        composition_resize(result, )
    i = 0
    while i < composition.used:
        status = composition_set_element_count(result, composition.elements[i], composition.counts[i])
        composition.elements[i] = composition.elements[i]
        composition.counts[i] = composition.counts[i]
        i += 1
    return 0



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
    '''
    Set the count for `element` in `composition`

    Return Values:
    0: Success
    1: General Failure
    -1: Failure due to Out-of-Memory
    '''
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
    '''
    Get the count of `element` in `composition`. The count is 0 if `element` is not in `composition`

    Return Values:
    0: Success
    1: Not Found
    '''
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
    return 1

cdef int composition_inc_element_count(Composition* composition, char* element, count_type increment) nogil:
    '''
    Increase the count for `element` in `composition` by `increment`.

    Return Values:
    0: Success
    1: General Failure
    -1: Failure due to Out-of-Memory
    '''
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
    '''
    Increases the size of the parallel arrays in `composition`, doubling them in length

    Return Values:
    0: Success
    -1: Failure due to Out-of-Memory
    '''
    composition.elements = <char**>realloc(composition.elements, sizeof(char*) * composition.size * 2)
    composition.counts = <count_type*>realloc(composition.counts, sizeof(count_type) * composition.size * 2)
    composition.size *= 2
    if composition.counts == NULL:
        return -1
    return 0

cdef int composition_resize_to(Composition* composition, size_t size) nogil:
    '''
    Increases the size of the parallel arrays in `composition`, doubling them in length

    Return Values:
    0: Success
    -1: Failure due to Out-of-Memory
    '''
    composition.elements = <char**>realloc(composition.elements, sizeof(char*) * size)
    composition.counts = <count_type*>realloc(composition.counts, sizeof(count_type) * size)
    composition.size = size
    if composition.counts == NULL:
        return -1
    return 0

cdef double composition_mass(Composition* composition) nogil:
    '''
    Calculates the monoisotopic mass of `composition`
    '''
    cdef:
        double mass
        Element* element
        char* element_label
        char[10] element_name
        int isotope_number
        int status
        size_t i
    i = 0
    while i < composition.used:
        element_label = composition.elements[i]
        _parse_isotope_string(element_label, &isotope_number, element_name)
        status = element_hash_table_get(_ElementTable, element_name, &element)
        if status != 0:
            printf("Could not find element %s\n", element_label)
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

cdef Composition* composition_mul(Composition* composition, long scale) nogil:
    cdef:
        Composition* result
        size_t i

    i = 0
    result = copy_composition(composition)
    while i < result.used:
        result.counts[i] *= scale
        i += 1
    return result

cdef void composition_imul(Composition* composition, long scale) nogil:
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
        PyDict_SetItemString(result, symbol_c, PyInt_FromLong(value))
        i += 1
    return result

# appears to introduce memory corruption when called repeatedly.
# cdef Composition* dict_to_composition(dict comp_dict):
#     cdef:
#         Composition* result
#     result = make_composition()
#     fill_composition_from_dict(comp_dict, result)
#     return result


cdef Composition* dict_to_composition(dict comp_dict):
    cdef:
        Composition* result
        str symbol
        size_t i
        char* symbol_c
        count_type value
    result = make_composition()
    for symbol, value in comp_dict.items():
        symbol_c = PyStr_AsString(symbol)
        composition_set_element_count(result, symbol_c, value)
    return result


cdef int fill_composition_from_dict(dict comp_dict, Composition* composition) except 1:
    cdef:
        PyObject* pkey
        PyObject* pvalue
        Py_ssize_t pos
        str symbol
        size_t i
        char* symbol_c
        count_type value
    pos = 0
    while PyDict_Next(comp_dict, &pos, &pkey, &pvalue):
        if not isinstance(<object>pkey, str):
            raise TypeError("Composition keys must be strings!")
        PyStr_InternInPlace(&pkey)
        symbol = <str>pkey
        Py_INCREF(symbol)
        symbol_c = PyStr_AsString(symbol)
        i = strlen(symbol_c)
        if symbol_c[i - 1] == ']':
            ensure_fixed_isotope(symbol_c)
        value = PyInt_AsLong(int(<object>pvalue))
        composition_set_element_count(composition, symbol_c, value)
    return 0


cdef int composition_add_from_dict(Composition* composition, dict comp_dict, int sign) except 1:
    cdef:
        PyObject* pkey
        PyObject* pvalue
        Py_ssize_t pos
        str symbol
        size_t i
        char* symbol_c
        count_type value
    pos = 0
    while PyDict_Next(comp_dict, &pos, &pkey, &pvalue):
        if not isinstance(<object>pkey, str):
            raise TypeError("Composition keys must be strings!")
        PyStr_InternInPlace(&pkey)
        symbol = <str>pkey
        Py_INCREF(symbol)
        symbol_c = PyStr_AsString(symbol)
        i = strlen(symbol_c)
        if symbol_c[i - 1] == ']':
            ensure_fixed_isotope(symbol_c)
        value = PyInt_AsLong(int(<object>pvalue))
        composition_inc_element_count(composition, symbol_c, value * sign)
    return 0


# -----------------------------------------------------------------------------
# Prototyping and Testing

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

    j = element_hash_table_get(_ElementTable, "C", &elem)
    # print elem.symbol, element_monoisotopic_mass(elem)

    composition = make_composition()
    composition_set_element_count(composition, "O[18]", 1)
    composition_set_element_count(composition, "H", 2)
    print_composition(composition)
    print("Massing")
    print(composition_mass(composition))

    sum_composition = composition_add(composition, composition, 1)
    print("Summed")
    print(composition_mass(sum_composition))
    print(composition_mass(composition))
    composition_iadd(sum_composition, composition, 1)
    print(composition_mass(sum_composition))


def isotope_update_test(str element_symbol):
    cdef:
        char* c_element_symbol
        char* element
        char[10] isotope_string
        int isotope
        int status
        Element* elem_obj
        Element* fixed_isotope_element
        Isotope* isotope_obj
        str out

    element = <char*>malloc(sizeof(char) * 10)
    c_element_symbol = PyStr_AsString(element_symbol)
    printf("c_element_symbol %s\n", c_element_symbol)
    _parse_isotope_string(c_element_symbol, &isotope, element)
    printf("status %s\n", element)
    out = PyStr_FromString(element)
    printf("Element Symbol: %s, Isotope: %d\n", element, isotope)
    print("element_hash_table_get")
    status = element_hash_table_get(_ElementTable, element, &elem_obj)
    if status != 0:
        print("Error, ", element_symbol, element, isotope)
    print("Fetched element", elem_obj.symbol)
    isotope_obj = get_isotope_by_neutron_count(elem_obj.isotopes, isotope)
    print(isotope_obj.mass, isotope_obj.neutrons, isotope_obj.neutron_shift)
    _make_isotope_string(elem_obj, isotope_obj, isotope_string)
    print(isotope_string)
    fixed_isotope_element = make_fixed_isotope_element(elem_obj, isotope)
    printf("Fixed: %s, %f", fixed_isotope_element.symbol, element_monoisotopic_mass(fixed_isotope_element))


def isotope_parse_test(str element_symbol):
    cdef:
        char* c_element_symbol
        char* element
        char[10] isotope_string
        int isotope
    element = <char*>malloc(sizeof(char) * 10)
    c_element_symbol = PyStr_AsString(element_symbol)
    _parse_isotope_string(c_element_symbol, &isotope, element)
    print(element, isotope)
    free(element)


cdef class PyComposition(object):
    '''A mapping representing a chemical composition.

    Implements arithmetic operations, +/- is defined
    between a :class:`PyComposition` and a :class:`Mapping`-like
    object, and * is defined between a :class:`PyComposition` and
    an integer.
    '''

    @staticmethod
    cdef PyComposition _create(Composition* base):
        cdef:
            PyComposition inst

        inst = PyComposition.__new__(PyComposition)
        inst.impl = make_composition()
        if base != NULL:
            composition_iadd(inst.impl, base, 1)
        inst._clean = False
        return inst

    cdef void _set_impl(self, Composition* composition, bint free_existing=1):
        if free_existing:
            free_composition(self.impl)
        self.impl = composition
        self._clean = False

    cdef void _initialize_from_formula(self, str formula):
        cdef:
            ssize_t i, n
            char* cstr
            int state
        cstr = PyStr_AsString(formula)
        n = len(formula)
        state = initialize_composition_from_formula(cstr, n, self.impl)
        if state == 1:
            raise ValueError()
        elif state == 2:
            raise KeyError()

    def __init__(self, base=None, **kwargs):
        self.impl = make_composition()
        self._clean = False
        self.cached_mass = 0.

        if base is not None:
            if isinstance(base, (dict, PyComposition)):
                self.update(base)
            elif isinstance(base, str):
                self._initialize_from_formula(<str>base)

        if kwargs:
            self.update(<dict>kwargs)

    def __dealloc__(self):
        free_composition(self.impl)

    def __getitem__(self, str key):
        return self.getitem(key)

    def __setitem__(self, str key, count_type value):
        self.setitem(key, value)

    def __delitem__(self, str key):
        self.setitem(key, 0)

    def __contains__(self, str key):
        return self.getitem(key) != 0

    cdef count_type getitem(self, str key):
        cdef:
            count_type count
            int status
            PyObject* pkey
            char* ckey
        Py_INCREF(key)
        pkey = <PyObject*>key
        PyStr_InternInPlace(&pkey)
        Py_INCREF(<str>pkey)
        ckey = PyStr_AsString(<str>pkey)
        status = composition_get_element_count(self.impl, ckey, &count)
        if status == 0:
            return count
        else:
            return 0

    cdef void setitem(self, str key, count_type value):
        cdef:
            int status
            Py_ssize_t size
            PyObject* pkey
            char* ckey
        Py_INCREF(key)
        pkey = <PyObject*>key
        PyStr_InternInPlace(&pkey)
        Py_INCREF(<str>pkey)
        ckey = PyStr_AsUTF8AndSize(<str>pkey, &size)
        if ckey[size - 1] == ']':
            ensure_fixed_isotope(ckey)
        status = composition_set_element_count(self.impl, ckey, value)
        self._clean = False

    cdef void increment(self, str key, count_type value):
        cdef:
            int status
            Py_ssize_t size
            PyObject* pkey
            char* ckey
        Py_INCREF(key)
        pkey = <PyObject*>key
        PyStr_InternInPlace(&pkey)
        Py_INCREF(<str>pkey)
        ckey = PyStr_AsUTF8AndSize(<str>pkey, &size)
        if ckey[size - 1] == ']':
            ensure_fixed_isotope(ckey)
        status = composition_inc_element_count(self.impl, ckey, value)
        self._clean = False

    cpdef update(self, arg):
        cdef:
            Py_ssize_t i, n
            PyComposition composition
            dict d

        if isinstance(arg, PyComposition):
            composition = <PyComposition>arg
            for i in range(composition.impl.used):
                composition_set_element_count(
                    self.impl, composition.impl.elements[i],
                    composition.impl.counts[i])

        elif isinstance(arg, dict):
            d = <dict>arg
            fill_composition_from_dict(d, self.impl)
        else:
            for key, value in arg.items():
                self.setitem(key, value)

    cpdef list keys(self):
        cdef:
            size_t i
            int count
            char* elem
            list keys
            str key_str
        i = 0
        keys = []
        while i < self.impl.used:
            elem = self.impl.elements[i]
            count = self.impl.counts[i]
            i += 1
            if count == 0:
                continue
            key_str = PyStr_FromString(elem)
            PyList_Append(keys, key_str)
        return keys

    cpdef list values(self):
        cdef:
            size_t i
            count_type value
            list counts
        counts = []
        i = 0
        while i < self.impl.used:
            value = self.impl.counts[i]
            i += 1
            if value == 0:
                continue
            counts.append(value)
        return counts

    cpdef pop(self, str key, object default=None):
        value = self[key]
        self[key] = 0
        if value == 0:
            return default
        else:
            return value

    cpdef list items(self):
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
            i += 1
            if value == 0:
                continue
            items.append((PyStr_FromString(elem), value))
        return items

    cpdef PyComposition copy(self):
        inst = PyComposition._create(self.impl)
        return inst

    def __add__(self, other):
        cdef:
            PyComposition result
            PyComposition _other
            object temp

        if not isinstance(self, PyComposition):
            temp = self
            self = <PyComposition>other
            other = temp

        result = (<PyComposition>self).copy()

        if isinstance(other, PyComposition):
            _other = other
            composition_iadd(result.impl, _other.impl, 1)
        elif isinstance(other, dict):
            composition_add_from_dict(result.impl, <dict>other, 1)
        else:
            for key, value in other.items():
                if not isinstance(key, str):
                    raise TypeError("Composition keys must be strings!")
                result[key] += PyInt_AsLong(int(value))
        return result

    def __iadd__(self, other):
        cdef:
            PyComposition _other

        if isinstance(other, PyComposition):
            _other = other
            composition_iadd(self.impl, _other.impl, 1)
        elif isinstance(other, dict):
            composition_add_from_dict(self.impl, <dict>other, 1)
        else:
            for key, value in other.items():
                if not isinstance(key, str):
                    raise TypeError("Composition keys must be strings!")
                self[key] += PyInt_AsLong(int(value))
        return self

    def __sub__(self, other):
        cdef:
            PyComposition result
            PyComposition _other
            object temp

        if not isinstance(self, PyComposition):
            temp = self
            self = <PyComposition>other
            other = temp

        result = (<PyComposition>self).copy()

        if isinstance(other, PyComposition):
            _other = other
            composition_iadd(result.impl, _other.impl, -1)
        elif isinstance(other, dict):
            composition_add_from_dict(result.impl, <dict>other, -1)
        else:
            for key, value in other.items():
                if not isinstance(key, str):
                    raise TypeError("Composition keys must be strings!")
                result[key] -= PyInt_AsLong(int(value))
        return result

    def __isub__(self, other):
        cdef:
            PyComposition _other

        if isinstance(other, PyComposition):
            _other = other
            composition_iadd(self.impl, _other.impl, -1)
        elif isinstance(other, dict):
            composition_add_from_dict(self.impl, <dict>other, -1)
        else:
            for key, value in other.items():
                if not isinstance(key, str):
                    raise TypeError("Composition keys must be strings!")
                self[key] -= PyInt_AsLong(int(value))
        return self

    def __mul__(self, scale):
        cdef:
            PyComposition result
            PyComposition inst
            long scale_factor

        if isinstance(self, PyComposition):
            inst = <PyComposition>self
            if not isinstance(scale, int):
                raise TypeError("Cannot multiply a Composition by a non-integer!")
            scale_factor = PyInt_AsLong(scale)
        else:
            inst = <PyComposition>scale
            if not isinstance(self, int):
                raise TypeError("Cannot multiply a Composition by a non-integer!")
            scale_factor = PyInt_AsLong(self)

        result = inst.copy()

        composition_imul(result.impl, scale_factor)

        return result

    def __imul__(self, long scale):
        composition_imul(self.impl, scale)
        return self

    cdef void add_from(self, PyComposition other):
        composition_iadd(self.impl, other.impl, 1)

    cdef void subtract_from(self, PyComposition other):
        composition_iadd(self.impl, other.impl, -1)

    cdef void scale_by(self, long scale):
        composition_imul(self.impl, scale)

    def __iter__(self):
        return iter(self.keys())

    cpdef double mass(self):
        '''Calculate the monoisotopic mass of this chemical composition

        Returns
        -------
        float
        '''
        if self._clean:
            return self.cached_mass
        else:
            self.cached_mass = composition_mass(self.impl)
            self._clean = True
            return self.cached_mass

    def __len__(self):
        return self.impl.used

    def __repr__(self):
        return "PyComposition({%s})" % ', '.join("\"%s\": %d" % kv for kv in self.items())

    def __str__(self):
        return repr(self)

    cpdef bint __equality_pycomposition(self, PyComposition other):
        return composition_eq(self.impl, other.impl)

    cpdef bint __equality_dict(self, dict other):
        cdef:
            size_t n
            str key
            count_type value
            count_type my_value

        n = len(other)
        if n != self.impl.used:
            return False
        for key, value in other.items():
            my_value = self.getitem(key)
            if my_value != value:
                return False
        return True

    def __richcmp__(self, object other, int code):
        if not isinstance(self, PyComposition):
            temp = self
            self = <PyComposition>other
            other = temp
        if isinstance(other, PyComposition):
            if code == 2:
                return self.__equality_pycomposition(<PyComposition>other)
            elif code == 3:
                return not self.__equality_pycomposition(<PyComposition>other)
        elif isinstance(other, dict):
            if code == 2:
                return self.__equality_dict(<dict>other)
            elif code == 3:
                return not self.__equality_dict(<dict>other)
        else:
            other = dict(other)
            if code == 2:
                return self.__equality_dict(<dict>other)
            elif code == 3:
                return not self.__equality_dict(<dict>other)


Mapping.register(PyComposition)
MutableMapping.register(PyComposition)


cdef extern from "<stdlib.h>" nogil:
    int isdigit(int ch) nogil;
    int isuppser(int ch) nogil;


cdef enum:
      START
      ISOTOPE
      ELEMENT
      COUNT


cdef int initialize_composition_from_formula(char* formula, ssize_t n, Composition* composition) nogil:
    cdef:
        ssize_t i
        ssize_t elstart, elend
        ssize_t numstart, numend
        ssize_t isostart, isoend
        char* temp
        char a
        Element* elem
        Element* fixed_isotope_elem
        count_type count, prev_count
        int state, fixed_isotope, status, found

    elem = NULL
    prev_count = 0
    state = ELEMENT
    temp = <char*>malloc(sizeof(char) * 10)
    for i in range(10):
        temp[i] = "\0"
    elstart = isostart = isoend = 0
    for i in range(n):
        a = formula[i]
        if a == ']':
            if state == ISOTOPE:
                elend = i
                isoend = i
                numstart = i + 1
                state = COUNT
            else:
                return 1
        elif a == '[':
            if state == ELEMENT:
                elend = i
                isostart = i
                state = ISOTOPE
            else:
                return 1
        elif isdigit(a):
            if state == ISOTOPE or state == COUNT:
                pass
            elif state == ELEMENT:
                elend = i
                numstart = i
                state = COUNT
            else:
                return 1
        else:
            if state == ELEMENT:
                continue
            elif state == COUNT:
                numend = i
                strncpy(temp, formula + numstart, numend - numstart)
                temp[numend - numstart] = 0
                count = atoi(temp)
                strncpy(temp, formula + elstart, elend - elstart)
                temp[elend - elstart] = 0
                found = element_hash_table_get(_ElementTable, temp, &elem)
                if isostart != 0 and isoend != 0 and found != 0:
                    isostart += 1
                    strncpy(temp, formula + isostart, isoend - isostart)
                    temp[isoend - isostart] = 0
                    fixed_isotope = atoi(temp)
                    strncpy(temp, formula + elstart, isoend - elstart + 1)
                    temp[isoend - elstart + 2] = 0

                    found = element_hash_table_get(_ElementTable, temp, &fixed_isotope_elem)
                    if found != 0:
                        strncpy(temp, formula + elstart, isostart)
                        temp[isostart - elstart - 1] = 0
                        found = element_hash_table_get(_ElementTable, temp, &elem)
                        if found != 0:
                            return 2
                        fixed_isotope_elem = make_fixed_isotope_element(elem, fixed_isotope)
                        if fixed_isotope_elem == NULL:
                            return 3
                        element_hash_table_put(_ElementTable, fixed_isotope_elem)
                        elem = fixed_isotope_elem
                    else:
                        elem = fixed_isotope_elem
                prev_count = 0
                composition_inc_element_count(composition, elem.symbol, count)
                state = ELEMENT
                elstart = i
                isostart = 0
                isoend = 0

    if state == COUNT:
        numend = i + 1
        strncpy(temp, formula + numstart, numend - numstart)
        temp[numend - numstart] = 0
        count = atoi(temp)
        strncpy(temp, formula + elstart, elend - elstart)
        temp[elend - elstart] = 0
        found = element_hash_table_get(_ElementTable, temp, &elem)
        if isostart != 0 and isoend != 0 and found != 0:
            isostart += 1
            strncpy(temp, formula + isostart, isoend - isostart)
            temp[isoend - isostart] = 0
            fixed_isotope = atoi(temp)
            strncpy(temp, formula + elstart, isoend - elstart + 1)
            temp[isoend - elstart + 2] = 0
            found = element_hash_table_get(_ElementTable, temp, &fixed_isotope_elem)
            if found != 0:
                strncpy(temp, formula + elstart, isostart)
                temp[isostart - elstart - 1] = 0
                found = element_hash_table_get(_ElementTable, temp, &elem)
                if found != 0:
                    return 2
                fixed_isotope_elem = make_fixed_isotope_element(elem, fixed_isotope)
                if fixed_isotope_elem == NULL:
                    return 3
                element_hash_table_put(_ElementTable, fixed_isotope_elem)
                elem = fixed_isotope_elem
            else:
                elem = fixed_isotope_elem
        prev_count = 0
        composition_inc_element_count(composition, elem.symbol, count)
    free(temp)
    return 0


cpdef PyComposition parse_formula(str formula):
    """Parse a chemical formula and construct a :class:`PyComposition` object

    The formula must be made up of zero or more pieces following the pattern
    ``(<element>[A-Z][a-z]*)(<isotope>\[\d+\])?(<count>\d+)``. You cannot omit
    the `<count>` digits.

    Parameters
    ----------
    formula : :class:`str`

    Returns
    -------
    :class:`PyComposition`

    Examples
    --------
    >>> parse_formula("H2O1")
    PyComposition({"H": 2, "O": 1})
    >>> parse_formula("C34H53O15N7").mass()
    799.35996402671
    >>> parse_formula("C7H15C[13]1O6N[15]1")
    PyComposition({"C": 7, "H": 15, "C[13]": 1, "O": 6, "N[15]": 1})
    >>> parse_formula("C7H15C[13]1O6N[15]1").mass()
    223.09032693441

    Raises
    ------
    ValueError
        If the formula doesn't match the expected pattern
    """
    cdef:
        ssize_t i, n
        char* cstr
        int state
        PyComposition composition
    cstr = PyStr_AsString(formula)
    n = len(formula)
    composition = PyComposition()
    state = initialize_composition_from_formula(cstr, n, composition.impl)
    if state == 1:
        raise ValueError("Malformed Formula")
    elif state == 2:
        raise KeyError("Unrecognized Element")
    elif state == 3:
        raise KeyError("Unrecognized Isotope")
    return composition
