from cython cimport freelist
from brainpy._c.composition cimport Composition, ElementHashTable, Element
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


cdef Peak* make_peak(double mz, double intensity, int charge) noexcept nogil
cdef void print_peak(Peak* peak) noexcept nogil

cdef PeakList* make_peak_list() noexcept nogil
cdef void free_peak_list(PeakList* peaklist) noexcept nogil
cdef int resize_peak_list(PeakList* peaklist) noexcept nogil
cdef void peak_list_append(PeakList* peaklist, Peak* peak) noexcept nogil
cdef void peak_list_reset(PeakList* peaklist) noexcept nogil

cdef PeakList* peak_list_ignore_below(PeakList* peaklist, double ignore_below, PeakList* result) noexcept nogil
cdef PeakList* peak_list_truncate_after(PeakList* peaklist, double truncate_after, PeakList* result) noexcept nogil
cdef void peak_list_shift(PeakList* peaklist, double shift) noexcept nogil

cdef list peaklist_to_list(PeakList* peaklist)


cdef struct ElementCache:
    Element** elements
    size_t used
    size_t size
    ElementHashTable* source


cdef ElementCache* make_element_cache(ElementHashTable* source) noexcept nogil
cdef void free_element_cache(ElementCache* cache) noexcept nogil
cdef int resize_element_cache(ElementCache* cache) noexcept nogil
cdef int element_cache_put(ElementCache* cache, Element** element) noexcept nogil
cdef int element_cache_get(ElementCache* cache, char* symbol, Element** out) noexcept nogil


cdef struct IsotopicDistribution:
    Composition* composition
    IsotopicConstants* _isotopic_constants
    int order
    double average_mass
    Peak* monoisotopic_peak
    ElementCache* cache


cdef IsotopicDistribution* make_isotopic_distribution(Composition* composition, int order, ElementCache* cache=*) noexcept nogil
cdef void free_isotopic_distribution(IsotopicDistribution* distribution) noexcept nogil


@freelist(100000)
cdef class TheoreticalPeak(object):
    cdef:
        public double mz
        public double intensity
        public int charge

    @staticmethod
    cdef TheoreticalPeak _create(double mz, double intensity, int charge)

    cpdef bint _eq(self, TheoreticalPeak other)
    cpdef TheoreticalPeak clone(self)


cdef int max_variants(Composition* composition, ElementCache* cache) noexcept nogil
cdef int guess_npeaks(Composition* composition_struct, size_t max_npeaks, ElementCache* cache=*) noexcept nogil

cpdef list _isotopic_variants(object composition, object npeaks=*, int charge=*, double charge_carrier=*)


cdef PeakList* isotopic_variants(Composition* composition, int npeaks, int charge=*, double charge_carrier=*) noexcept nogil


cdef size_t max_variants_approx(double mass, double lambda_factor=*, size_t maxiter=*, double threshold=*) noexcept nogil