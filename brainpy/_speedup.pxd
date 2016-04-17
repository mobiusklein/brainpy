cdef dict nist_mass
cdef double PROTON
cdef dict periodic_table


cdef double neutral_mass(double mz,  int z, double charge_carrier=*)
cdef double mass_charge_ratio(double neutral_mass, int z, double charge_carrier=*)

cdef double calculate_mass(dict composition, dict mass_data=*)
