from cpython cimport bool


cdef class IVTreeNode:
    cdef IVTreeNode parent, left, right
    cdef long start, max_end, aug_end
    cdef IVNode list_head
    cdef bool black

    cdef insert(self, long end, object interval)
    cdef bool remove(self, object interval)
    cdef bool is_empty(self)


cdef class IVNode:
    cdef IVNode prev, next
    cdef long end
    cdef object interval
