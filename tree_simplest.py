#!/usr/bin/env python3

"""
tree_simplest.py
Probably the simplest implementation of a tree, k-d tree, kd-tree, 2-dimensional tree, or whatever name it has...

"""

__author__ = "Miroslav Broz (miroslav.broz@email.cz)"
__version__ = "Sep 14th 2023"

import math

def square_distance(a, b):
    s = 0.0
    for i in range(0, len(a)):
        s += (a[i]-b[i])**2
    return s

class Node(object):
    def __init__(self, p, left, right):
        self.p = p
        self.left = left
        self.right = right

class Tree(object):
    def __init__(self, k, points):
        self.k = k

        def build_tree(points, depth=0):
            if len(points) == 0:
                return None

            axis = depth % self.k
            points.sort(key=lambda x: x[axis])
            i = len(points) // 2
 
            return Node( \
                points[i], \
                build_tree(points[:i], depth+1), \
                build_tree(points[i+1:], depth+1), \
            )

        self.root = build_tree(points)

    def nearest_neighbor(self, p):
        best = [None, float('inf')]

        def recursive_search(node, depth=0):
            if node is None:
                return

            dist = square_distance(node.p, p)
            if dist < best[1]:
                best[:] = node.p, dist  # in-place!

            axis = depth % self.k
            diff = p[axis] - node.p[axis]
            near, away = (node.left, node.right) if diff <= 0.0 else (node.right, node.left)

            recursive_search(near, depth+1)
            if diff**2 < best[1]:
                recursive_search(away, depth+1)

        recursive_search(self.root)
        return best[0]

def main():
    points = [(2, 3), (5, 4), (9, 6), (4, 7), (8, 1), (7, 2)]

    tree = Tree(2, points)

    p = (8, 5)
    q = tree.nearest_neighbor(p)
    print('p = ', p)
    print('q = ', q)

if __name__ == "__main__":
    main()


