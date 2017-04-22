from cpython cimport bool

from .ivtreenode cimport IVTreeNode as IVNode

cdef class IVTree:
    cdef IVNode sentinel
    cdef IVNode root

    cdef bool _overlaps(self, long a_start, long a_end,
                        long b_start, long b_end)
    cdef bool _update_aug(self, IVNode node)
    cdef _left_rotate(self, IVNode node)
    cdef _right_rotate(self, IVNode node)
    cdef _rb_transplant(self, IVNode node_dst, IVNode node_src)
    cdef _min_node(self, IVNode node)
    cdef _insert_rb_fixup(self, IVNode node)
    cdef _delete_rb_fixup(self, IVNode node)
