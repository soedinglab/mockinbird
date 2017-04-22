from mockinbird.ivtree.ivtreenode cimport IVTreeNode as IVNode

def left(IVNode node):
    return node.left

def set_left(IVNode node, IVNode left):
    node.left = left

def right(IVNode node):
    return node.right

def set_right(IVNode node, IVNode right):
    node.right = right

def black(IVNode node):
    return node.black

def set_black(IVNode node, black):
    node.black = black

def max_end(IVNode node):
    return node.max_end

def data(IVNode node):
    return node.data

def aug_end(IVNode node):
    return node.aug_end

def set_aug_end(IVNode node, aug_end):
    node.aug_end = aug_end

def insert(IVNode node, iv):
    node.insert(iv._int_end, iv)

def remove(IVNode node, iv):
    return node.remove(iv)

def is_empty(IVNode node):
    return node.is_empty()
