from cython cimport freelist
from brainpy._c.composition cimport Composition
from brainpy._c.isotopic_constants cimport IsotopicConstants
from brainpy._c.double_vector cimport DoubleVector


cdef struct Peak:
    double mz
    double intensity
    int charge


cdef struct PeakList:
    Peak* peaks
    size_t used
    size_t size


cdef PeakList* make_peak_list() nogil
cdef void free_peak_list(PeakList* peaklist) nogil

cdef Peak* make_peak(double mz, double intensity, int charge) nogil
cdef list peaklist_to_list(PeakList* peaklist)

cdef void print_peak(Peak* peak) nogil

cdef struct IsotopicDistribution:
    Composition* composition
    IsotopicConstants* _isotopic_constants
    int order
    double average_mass
    Peak* monoisotopic_peak


@freelist(100000)
cdef class TheoreticalPeak(object):
    cdef:
        public double mz
        public double intensity
        public int charge

    @staticmethod
    cdef TheoreticalPeak _create(double mz, double intensity, int charge)

    cpdef TheoreticalPeak clone(self)


cpdef list _isotopic_variants(object composition, object npeaks=*, int charge=*, double charge_carrier=*)


cdef PeakList* isotopic_variants(Composition* composition, int npeaks, int charge=*, double charge_carrier=*) nogil