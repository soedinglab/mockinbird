from mockinbird.ivtree.ivtree cimport IVTree
from mockinbird.ivtree.ivtreenode cimport IVTreeNode as IVNode

def root(IVTree tree):
    return tree.root

def set_root(IVTree tree, IVNode node):
    tree.root = node

def sentinel(IVTree tree):
    return tree.sentinel
