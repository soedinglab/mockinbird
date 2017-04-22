from cpython cimport bool
from .ivtreenode cimport IVTreeNode as IVNode

cdef class IVTree:
    """Interval tree implementation"""
    def __init__(self):
        # the sentinel is black, other attributes can be chosen arbitrarily
        self.sentinel = IVNode(0, 0, True, None)
        self.root = self.sentinel

    def query_overlap(self, interval):
        cdef IVNode sentinel, node
        cdef long start, end
        node = self.root
        sentinel = self.sentinel
        start = interval._int_start
        end = interval._int_end
        while node is not sentinel and not\
              self._overlaps(start, end, node.start, node.max_end):
            if node.left is not sentinel and node.left.aug_end >= start:
                node = node.left
            else:
                node = node.right
        if node is not sentinel:
            for iv in node.overlap(start, end):
                return iv
        else:
            return None

    def has_overlap(self, interval):
        return self.query_overlap(interval) is not None

    def query_all_overlaps(self, interval):
        cdef IVNode sentinel, node
        cdef long start, end
        cdef long stack_size
        sentinel = self.sentinel
        node = self.root
        start = interval._int_start
        end = interval._int_end
        stack = []
        stack_size = 0
        while node is not sentinel or stack_size > 0:
            if node is not sentinel:
                stack.append(node)
                stack_size += 1
                if node.left is not sentinel and\
                        node.left.aug_end >= start:
                    node = node.left
                else:
                    node = sentinel
            else:
                node = stack.pop()
                stack_size -= 1
                if self._overlaps(start, end, node.start, node.max_end):
                    yield from node.overlap(start, end)

                if node.right is not sentinel and\
                        node.start <= end:
                    node = node.right
                else:
                    node = sentinel

    def insert(self, interval):
        cdef IVNode sentinel, x, y, z, runner
        cdef long start, end
        start = interval._int_start
        end = interval._int_end
        sentinel = self.sentinel
        y = sentinel
        x = self.root
        while x is not sentinel:
            y = x
            if x.start > start:
                x = x.left
            elif x.start < start:
                x = x.right
            else:
                x.insert(end, interval)
                z = x
                break

        if x is sentinel:
            z = IVNode(start, end, False, interval)
            z.parent = y
            if y is sentinel:
                self.root = z
            elif start < y.start:
                y.left = z
            else:
                y.right = z
            z.left = sentinel
            z.right = sentinel
            z.black = False

        runner = z
        while runner is not sentinel:
            if not self._update_aug(runner):
                break
            runner = runner.parent

        if x is sentinel:
            self._insert_rb_fixup(z)

    def remove(self, interval):
        cdef long start, end
        cdef bool success, y_orig_black
        cdef IVNode sentinel, node, x, y, z
        start = interval._int_start
        end = interval._int_end
        sentinel = self.sentinel
        node = self.root
        while node is not sentinel:
            if node.start > start:
                node = node.left
            elif node.start < start:
                node = node.right
            else:
                break
        if node is sentinel:
            return False

        success = node.remove(interval)
        if success and not node.is_empty():
            while node is not sentinel:
                if not self._update_aug(node):
                    break
                node = node.parent
        elif success and node.is_empty():
            z = node
            y = z
            y_orig_black = y.black
            if z.left is sentinel:
                x = z.right
                self._rb_transplant(z, z.right)
            elif z.right is sentinel:
                x = z.left
                self._rb_transplant(z, z.left)
            else:
                y = self._min_node(z.right)
                y_orig_black = y.black
                x = y.right
                if y.parent is z:
                    x.parent = y
                else:
                    self._rb_transplant(y, y.right)
                    y.right = z.right
                    y.right.parent = y
                self._rb_transplant(z, y)
                y.left = z.left
                y.left.parent = y
                y.black = z.black

            node = x.parent
            while node is not y.parent:
                self._update_aug(node)
                node = node.parent
            node = y.parent
            while node is not sentinel:
                if not self._update_aug(node):
                    break
                node = node.parent
            if y_orig_black:
                self._delete_rb_fixup(x)
        return success

    def get_range(self):
        cdef IVNode node
        if self.root is self.sentinel:
            return None
        node = self._min_node(self.root)
        left_bound = node.start
        right_bound = self.root.aug_end
        return left_bound, right_bound

    def __iter__(self):
        cdef IVNode node, sentinel
        cdef long stack_size
        sentinel = self.sentinel
        node_stack = []
        stack_size = 0
        node = self.root
        while node is not sentinel or stack_size > 0:
            while node is not sentinel:
                node_stack.append(node)
                stack_size += 1
                node = node.left
            node = node_stack.pop()
            stack_size -= 1
            yield from node.get_all_intervals()
            node = node.right

    # helper functions
    cdef bool _overlaps(self, long a_start, long a_end,
                        long b_start, long b_end):
        return a_start <= b_end and b_start <= a_end

    cdef bool _update_aug(self, IVNode node):
        cdef bool changed
        cdef int aug_end
        aug_end = node.max_end
        if node.left.aug_end > aug_end:
            aug_end = node.left.aug_end
        if node.right.aug_end > aug_end:
            aug_end = node.right.aug_end
        changed = aug_end != node.aug_end
        node.aug_end = aug_end
        return changed

    cdef _left_rotate(self, IVNode node):
        cdef IVNode x,y
        x = node
        y = x.right
        x.right = y.left
        if y.left is not self.sentinel:
            y.left.parent = x
        y.parent = x.parent
        if x.parent is self.sentinel:
            self.root = y
        elif x is x.parent.left:
            x.parent.left = y
        else:
            x.parent.right = y
        y.left = x
        x.parent = y
        self._update_aug(x)
        self._update_aug(y)

    cdef _right_rotate(self, IVNode node):
        cdef IVNode x,y
        x = node
        y = x.left
        x.left = y.right
        if y.right is not self.sentinel:
            y.right.parent = x
        y.parent = x.parent
        if x.parent is self.sentinel:
            self.root = y
        elif x is x.parent.right:
            x.parent.right = y
        else:
            x.parent.left = y
        y.right = x
        x.parent = y
        self._update_aug(x)
        self._update_aug(y)

    cdef _rb_transplant(self, IVNode node_dst, IVNode node_src):
        if node_dst.parent is self.sentinel:
            self.root = node_src
        elif node_dst is node_dst.parent.left:
            node_dst.parent.left = node_src
        else:
            node_dst.parent.right = node_src
        node_src.parent = node_dst.parent

    cdef _min_node(self, IVNode node):
        while node.left is not self.sentinel:
            node = node.left
        return node

    cdef _insert_rb_fixup(self, IVNode node):
        cdef IVNode z, y
        z = node
        while not z.parent.black:
            if z.parent is z.parent.parent.left:
                y = z.parent.parent.right
                if not y.black:
                    z.parent.black = True
                    y.black = True
                    z.parent.parent.black = False
                    z = z.parent.parent
                else:
                    if z is z.parent.parent.right:
                        z = z.parent
                        self._left_rotate(z)
                    z.parent.black = True
                    z.parent.parent.black = False
                    self._right_rotate(z.parent.parent)
            else:
                y = z.parent.parent.left
                if not y.black:
                    z.parent.black = True
                    y.black = True
                    z.parent.parent.black = False
                    z = z.parent.parent
                else:
                    if z is z.parent.left:
                        z = z.parent
                        self._right_rotate(z)
                    z.parent.black = True
                    z.parent.parent.black = False
                    self._left_rotate(z.parent.parent)
        self.root.black = True

    cdef _delete_rb_fixup(self, IVNode node):
        cdef IVNode sentinel, root, w
        sentinel = self.sentinel
        root = self.root
        while node is not root and node.black:
            if node is node.parent.left:
                w = node.parent.right
                if not w.black:
                    w.black = True
                    node.parent.black = False
                    self._left_rotate(node.parent)
                    w = node.parent.right
                if w.left.black and w.right.black:
                    w.black = False
                    node = node.parent
                else:
                    if w.right.black:
                        w.left.black = True
                        w.black = False
                        self._right_rotate(w)
                        w = node.parent.right
                    w.black = node.parent.black
                    node.parent.black = True
                    w.right.black = True
                    self._left_rotate(node.parent)
                    node = root
            else:
                w = node.parent.left
                if not w.black:
                    w.black = True
                    node.parent.black = False
                    self._right_rotate(node.parent)
                    w = node.parent.left
                if w.right.black and w.left.black:
                    w.black = False
                    node = node.parent
                else:
                    if w.left.black:
                        w.right.black = True
                        w.black = False
                        self._left_rotate(w)
                        w = node.parent.left
                    w.black = node.parent.black
                    node.parent.black = True
                    w.left.black = True
                    self._right_rotate(node.parent)
                    node = root
        node.black = True
