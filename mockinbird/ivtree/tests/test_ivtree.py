import unittest
import random

from nose.tools import raises

from pyivtree import GenomicInterval, IVTree
from pyivtree.ivtreenode import IVTreeNode as IVNode
from pyivtree.tests import node_utils as nu
from pyivtree.tests import tree_utils as tu

random.seed(20140907)


class TreeStructTestCase(unittest.TestCase):

    def check_root_color(self, ivtree):
        root = tu.root(ivtree)
        self.assertTrue(nu.black(root))

    def check_sentinel_color(self, ivtree):
        sentinel = tu.sentinel(ivtree)
        self.assertTrue(nu.black(sentinel))

    def check_child_color(self, ivtree):
        stack = []
        sentinel = tu.sentinel(ivtree)
        node = tu.root(ivtree)
        while node is not sentinel or len(stack) > 0:
            while node is not sentinel:
                stack.append(node)
                node = nu.left(node)
            node = stack.pop()
            if not nu.black(node):
                self.assertTrue(nu.black(nu.left(node)))
                self.assertTrue(nu.black(nu.right(node)))
            node = nu.left(node)

    def check_black_height(self, ivtree):
        sentinel = tu.sentinel(ivtree)
        node = tu.root(ivtree)
        height_stack = []
        node_stack = []
        # -1 = left, +1 = right
        direction_stack = []
        while node is not sentinel or len(node_stack) != 0:
            while node is not sentinel:
                node_stack.append(node)
                direction_stack.append(-1)
                node = nu.left(node)
            height_stack.append(1)
            while direction_stack[-1] == 1:
                node = node_stack.pop()
                direction_stack.pop()
                black = nu.black(node)
                right_height = height_stack.pop()
                left_height = height_stack.pop()
                self.assertEqual(right_height, left_height)
                new_height = right_height + 1 if black else right_height
                if len(node_stack) == 0:
                    # terminating and returning the black height of the root
                    return new_height
                height_stack.append(new_height)
            node = nu.right(node_stack[-1])
            direction_stack[-1] = 1

    def check_augmentation(self, tree):
        sentinel = tu.sentinel(tree)
        node = tu.root(tree)
        node_stack = []
        while node is not sentinel or len(node_stack) > 0:
            while node is not sentinel:
                node_stack.append(node)
                node = nu.left(node)
            node = node_stack.pop()
            aug_end = nu.max_end(node)
            if nu.left(node) is not sentinel:
                left_aug_end = nu.aug_end(nu.left(node))
                if left_aug_end > aug_end:
                    aug_end = left_aug_end
            if nu.right(node) is not sentinel:
                right_aug_end = nu.aug_end(nu.right(node))
                if right_aug_end > aug_end:
                    aug_end = right_aug_end
            self.assertEqual(aug_end, nu.aug_end(node))
            node = nu.right(node)

    def check_structure(self, ivtree):
        self.check_root_color(ivtree)
        self.check_sentinel_color(ivtree)
        self.check_child_color(ivtree)
        self.check_black_height(ivtree)
        self.check_augmentation(ivtree)

    def create_custom_tree(self):
        ivtree = IVTree()
        sentinel = tu.sentinel(ivtree)

        iv1 = GenomicInterval(5, 6)
        iv2 = GenomicInterval(3, 4)
        iv3 = GenomicInterval(7, 8)
        iv4 = GenomicInterval(1, 2)
        iv5 = GenomicInterval(4, 5)

        node1 = IVNode(iv1.start, iv1.end, True, iv1)
        node2 = IVNode(iv2.start, iv2.end, False, iv2)
        node3 = IVNode(iv3.start, iv3.end, True, iv3)
        node4 = IVNode(iv4.start, iv4.end, True, iv4)
        node5 = IVNode(iv5.start, iv5.end, True, iv5)

        nu.set_left(node1, node2)
        nu.set_right(node1, node3)
        nu.set_left(node2, node4)
        nu.set_right(node2, node5)
        nu.set_left(node4, sentinel)
        nu.set_right(node4, sentinel)
        nu.set_left(node5, sentinel)
        nu.set_right(node5, sentinel)
        nu.set_left(node3, sentinel)
        nu.set_right(node3, sentinel)

        nu.set_aug_end(node1, 8)
        nu.set_aug_end(node2, 5)
        nu.set_aug_end(node3, 8)
        nu.set_aug_end(node4, 2)
        nu.set_aug_end(node5, 5)

        tu.set_root(ivtree, node1)
        return ivtree

    def create_intervals(self, iv_count, iv_start_min, iv_start_max,
                         iv_len_max):
        ivs = []
        for i in range(iv_count):
            iv_start = random.randrange(iv_start_min, iv_start_max)
            iv_end = iv_start + random.randrange(iv_len_max)
            ivs.append(GenomicInterval(iv_start, iv_end))
        return ivs


class IVTreeTest(TreeStructTestCase):

    def test_empty_tree(self):
        tree = IVTree()
        self.check_structure(tree)

    def test_custom_tree(self):
        tree = self.create_custom_tree()
        self.check_structure(tree)

    @raises(AssertionError)
    def test_bad_root_color(self):
        tree = self.create_custom_tree()
        root = tu.root(tree)
        nu.set_black(root, False)
        self.check_root_color(tree)

    @raises(AssertionError)
    def test_bad_sentinel_color(self):
        tree = self.create_custom_tree()
        sentinel = tu.sentinel(tree)
        nu.set_black(sentinel, False)
        self.check_sentinel_color(tree)

    @raises(AssertionError)
    def test_bad_child_color(self):
        tree = self.create_custom_tree()
        root = tu.root(tree)
        node4 = nu.left(nu.left(root))
        nu.set_black(node4, False)
        self.check_child_color(tree)

    @raises(AssertionError)
    def test_bad_black_height(self):
        tree = self.create_custom_tree()
        root = tu.root(tree)
        node3 = nu.right(root)
        nu.set_black(node3, False)
        self.check_black_height(tree)

    @raises(AssertionError)
    def test_bad_augmentation(self):
        tree = self.create_custom_tree()
        root = tu.root(tree)
        node2 = nu.left(root)
        nu.set_aug_end(node2, 4)
        self.check_augmentation(tree)

    def test_insertion(self):
        count, min_start, max_start, max_len = 200, 0, 200, 50
        tree = IVTree()
        for iv in self.create_intervals(count, min_start, max_start, max_len):
            tree.insert(iv)
            self.check_structure(tree)

    def test_deletion(self):
        count, min_start, max_start, max_len = 200, 0, 200, 50
        ivs = self.create_intervals(count, min_start, max_start, max_len)
        tree = IVTree()
        for iv in ivs:
            tree.insert(iv)
        self.check_structure(tree)
        random.shuffle(ivs)
        for iv in ivs:
            success = tree.remove(iv)
            self.assertTrue(success)
            self.check_structure(tree)

    def test_node_indel(self):
        iv1 = GenomicInterval(1, 2)
        iv2 = GenomicInterval(1, 3)
        node = IVNode(iv1._int_start, iv1._int_end, True, iv1)
        self.assertFalse(nu.is_empty(node))
        nu.insert(node, iv2)
        self.assertEqual(nu.max_end(node), iv2._int_end)
        self.assertTrue(nu.remove(node, iv2))
        self.assertEqual(nu.max_end(node), iv1._int_end)
        self.assertFalse(nu.is_empty(node))
        self.assertTrue(nu.remove(node, iv1))
        self.assertTrue(nu.is_empty(node))

    def test_node_overlap(self):
        iv1 = GenomicInterval(1, 5)
        iv2 = GenomicInterval(1, 10)
        node = IVNode(iv1._int_start, iv1._int_end, True, iv1)
        nu.insert(node, iv2)
        check_ivs = [
            (GenomicInterval(-5, 0), 0),
            (GenomicInterval(-2, 2), 2),
            (GenomicInterval(2, 4), 2),
            (GenomicInterval(4, 6), 2),
            (GenomicInterval(6, 8), 1),
            (GenomicInterval(10, 12), 1),
            (GenomicInterval(12, 14), 0),
        ]
        for iv, exp_ovl in check_ivs:
            obs_ovl = len(list(node.overlap(iv._int_start, iv._int_end)))
            self.assertEqual(obs_ovl, exp_ovl)

    def test_overlap_all(self):
        count, min_start, max_start, max_len = 200, 0, 500, 50
        ovl_iv = GenomicInterval(50, 75)
        tree = IVTree()
        ovl_count = 0
        for iv in self.create_intervals(count, min_start, max_start, max_len):
            tree.insert(iv)
            if iv.overlaps(ovl_iv):
                ovl_count += 1
        tree_ovl_count = len(list(tree.query_all_overlaps(ovl_iv)))
        self.assertEqual(ovl_count, tree_ovl_count)

    def test_range(self):
        tree = IVTree()
        self.assertIsNone(tree.get_range())
        iv1 = GenomicInterval(2, 5)
        iv2 = GenomicInterval(7, 8)
        iv3 = GenomicInterval(7, 10)
        iv4 = GenomicInterval(4, 8)
        for iv in iv1, iv2, iv3, iv4:
            tree.insert(iv)
        self.assertEqual(tree.get_range(), (2, 10))

    def test_iterate(self):
        iv1 = GenomicInterval(2, 5)
        iv2 = GenomicInterval(4, 6)
        iv3 = GenomicInterval(4, 10)
        iv4 = GenomicInterval(6, 8)
        ivs = [iv1, iv2, iv3, iv4]
        tree = IVTree()
        for iv in iv4, iv2, iv1, iv3:
            tree.insert(iv)
        for i, iv in enumerate(tree):
            self.assertIs(iv, ivs[i])

    def test_overlap_crosstest(self):
        count, min_start, max_start, max_len = 250, 1, 10000, 20
        step_size = 25
        ivs = self.create_intervals(count, min_start, max_start, max_len)
        loop_iv_start = min_start - step_size - 1
        tree = IVTree()
        for iv in ivs:
            tree.insert(iv)
        for i in range(loop_iv_start, max_start, step_size):
            test_iv = GenomicInterval(i, i + step_size)
            naive_ovl_count = 0
            for iv in ivs:
                if test_iv.overlaps(iv):
                    naive_ovl_count += 1
            ovl_all = list(tree.query_all_overlaps(test_iv))
            ovl_all_count = len(ovl_all)
            self.assertEqual(ovl_all_count, naive_ovl_count)
            ovl_iv = tree.query_overlap(test_iv)
            ovl = tree.has_overlap(test_iv)
            if ovl_all_count > 0:
                self.assertIsNotNone(ovl_iv)
                self.assertTrue(ovl)
            else:
                self.assertIsNone(ovl_iv)
                self.assertFalse(ovl)
