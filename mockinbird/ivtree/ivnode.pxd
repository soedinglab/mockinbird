from cpython cimport bool


cdef class IVNode:
    cdef IVNode parent, left, right
    cdef long start, end, max_end
    cdef object data
    cdef bool black
