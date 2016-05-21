cdef dict nist_mass
cdef double PROTON
cdef dict periodic_table


cdef double neutral_mass(double mz,  int z, double charge_carrier=*)
cdef double mass_charge_ratio(double neutral_mass, int z, double charge_carrier=*)

cdef double calculate_mass(dict composition, dict mass_data=*)

cpdef Element make_fixed_isotope_element(Element element, int isotope_shift)

cdef class Element:
    cdef:
        public str symbol
        public dict isotopes
        double _monoisotopic_mass
        int _max_neutron_shift
        int _min_neutron_shift
        list _no_mass_elementary_symmetric_polynomial_cache
        list _no_mass_power_sum_cache
        list _mass_elementary_symmetric_polynomial_cache
        list _mass_power_sum_cache
    
    cpdef double monoisotopic_mass(self)