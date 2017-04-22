from cpython cimport bool

cdef class IVTreeNode:
    def __init__(self, long start, long end, bool black, object data):
        cdef IVNode node
        self.start = start
        self.max_end = end
        self.black = black
        node = IVNode()
        node.prev = None
        node.next = None
        node.end = end
        node.interval = data
        self.list_head = node

    cdef insert(self, long end, object interval):
        cdef IVNode prev, node, new_node
        new_node = IVNode()
        new_node.end = end
        new_node.interval = interval

        prev = None
        node = self.list_head
        while node is not None and node.end <= end:
            prev = node
            node = node.next
        if node is None:
            new_node.prev = prev
            new_node.next = None
            prev.next = new_node
            self.max_end = end
        else:
            new_node.next = node
            new_node.prev = node.prev
            if node.prev is not None:
                node.prev.next = new_node
            else:
                self.list_head = new_node
            node.prev = new_node
    
    cdef bool remove(self, object interval):
        cdef IVNode node
        node = self.list_head
        while node is not None and node.interval is not interval:
            node = node.next
        if node is None:
            return False
        if node.prev is None and node.next is None:
            self.list_head = None
        elif node.prev is None:
            node.next.prev = None
            self.list_head = node.next
        elif node.next is None:
            node.prev.next = None
            self.max_end = node.prev.end
        else:
            node.next.prev = node.prev
            node.prev.next = node.next
        return True

    def overlap(self, long start, long end):
        cdef IVNode node
        node = self.list_head
        while node is not None and node.end < start:
            node = node.next
        if node is not None and end >= self.start:
            while node is not None:
                yield node.interval
                node = node.next

    def get_all_intervals(self):
        cdef IVNode node
        node = self.list_head
        while node is not None:
            yield node.interval
            node = node.next

    cdef bool is_empty(self):
        return self.list_head is None


cdef class IVNode:
    pass
