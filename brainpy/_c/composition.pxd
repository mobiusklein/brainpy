from brainpy.mass_dict import nist_mass as __nist_mass

ctypedef int count_type

cdef dict nist_mass
nist_mass = __nist_mass

cdef double PROTON = nist_mass["H+"][0][0]

cdef double neutral_mass(double mz,  int z, double charge_carrier=*) nogil
cdef double mass_charge_ratio(double neutral_mass, int z, double charge_carrier=*) nogil

cdef struct Isotope:
    double mass
    double abundance
    int neutrons
    int neutron_shift

cdef struct IsotopeMap:
    Isotope* bins
    size_t size

cdef struct Element:
    char* symbol
    IsotopeMap* isotopes

cdef struct PeriodicTable:
    Element** elements
    size_t size

cdef PeriodicTable* _PeriodicTable

cdef struct Composition:
    char** elements
    count_type* counts
    size_t size
    size_t used


cdef size_t hash_string(char *str) nogil


cdef IsotopeMap* make_isotope_map(list organized_isotope_data, size_t size)

cdef Isotope* get_isotope_by_neutron_shift(IsotopeMap* isotopes, int neutron_shift) nogil
cdef void free_isotope_map(IsotopeMap* isotopes) nogil



cdef void _isotopes_of(char* element_symbol, IsotopeMap** isotope_frequencies)
cdef Element* make_element(char* symbol)
cdef void make_element_in_place(char* symbol, Element* element)

cdef double element_monoisotopic_mass(Element* element) nogil
cdef int element_min_neutron_shift(Element* element) nogil
cdef int element_max_neutron_shift(Element* element) nogil
cdef void free_element(Element* element) nogil

cdef PeriodicTable* make_periodic_table()
cdef int get_element_from_periodic_table(PeriodicTable* table, char* symbol, int* out) nogil
cdef int get_element_from_periodic_table2(PeriodicTable* table, char* symbol, Element** out) nogil

cdef Composition* make_composition() nogil
cdef Composition* copy_composition(Composition* composition) nogil
cdef void print_composition(Composition* composition) nogil
cdef int composition_set_element_count(Composition* composition, char* element, count_type count) nogil
cdef int composition_get_element_count(Composition* composition, char* element, count_type* count) nogil
cdef int composition_inc_element_count(Composition* composition, char* element, count_type increment) nogil
cdef int composition_resize(Composition* composition) nogil
cdef double composition_mass(Composition* composition) nogil
cdef void free_composition(Composition* composition) nogil

cdef Composition* composition_add(Composition* composition_1, Composition* composition_2, int sign) nogil
cdef int composition_iadd(Composition* composition_1, Composition* composition_2, int sign) nogil
cdef Composition* composition_mul(Composition* composition, int scale) nogil
cdef void composition_imul(Composition* composition, int scale) nogil

cdef dict composition_to_dict(Composition* composition)
cdef Composition* dict_to_composition(dict comp_dict)

cdef class PyComposition(object):
    cdef:
        Composition* impl
        double cached_mass
        bint _clean
    cpdef double mass(self)