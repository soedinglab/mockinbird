cdef class IVNode:
    def __init__(self, start, end, black, data):
        self.start = start
        self.end = end
        self.black = black
        self.data = data
