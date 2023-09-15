#!/usr/bin/env python3

"""
tree.py
A k-dimensional tree implementation with a nearest-neighbor search
and (full) multipole gravity computation.

References:

  Stadel, J.G.: Cosmological N-body simulations and their analysis.
  PhD Thesis, Univ. of Washington, 2001.

"""

__author__ = "Miroslav Broz (miroslav.broz@email.cz)"
__version__ = "Jan 20th 2016"

import math
import sys
import time
import random

def square_distance(a, b):
    """Squared distance between two points"""
    s = 0.0
    for i in range(0, len(a)):
        s += (a[i]-b[i])**2
    return s

def eval_multipole(mass, dx, dy, dz):
    """
    Evaluate gravitational multipole moments for a point-mass particles.
    see Stadel (2001), Eq. (2.5), M^_n_ = sum_{i in V} m_i x_i^_n_

    Beware, a "rank notation" is used here! It means simply:
    M^_0_ = M      ... monopole (there is NO index, i.e. the total mass)
    M^_1_ = M^i    ... dipole (i.e. zero in a centre-of-mass reference system, but we compute it anyway as a check)
    M^_2_ = M^ij   ... quadrupole, tensor of rank 2, of course
    M^_3_ = M^ijk  ... octupole
    M^_4_ = M^ijkl ... hexadecupole

    The sum is performed in an outer cycle.

    """
    M = mass

    Mx = M*dx
    My = M*dy
    Mz = M*dz

    # 2DO: precompute variables dx2, dx3, ...
    Mxx = M*dx*dx
    Myy = M*dy*dy
    Mzz = M*dz*dz
    Mxy = M*dx*dy
    Mxz = M*dx*dz
    Myz = M*dy*dz

    Mxxx = M*dx*dx*dx
    Myyy = M*dy*dy*dy
    Mzzz = M*dz*dz*dz
    Mxxy = M*dx*dx*dy
    Mxyy = M*dx*dy*dy
    Mxxz = M*dx*dx*dz
    Mxzz = M*dx*dz*dz
    Mxyz = M*dx*dy*dz
    Myyz = M*dy*dy*dz
    Myzz = M*dy*dz*dz

    Mxxxx = M*dx*dx*dx*dx
    Myyyy = M*dy*dy*dy*dy
    Mzzzz = M*dz*dz*dz*dz
    Mxxxy = M*dx*dx*dx*dy
    Mxxyy = M*dx*dx*dy*dy
    Mxyyy = M*dx*dy*dy*dy
    Mxxxz = M*dx*dx*dx*dz
    Mxxzz = M*dx*dx*dz*dz
    Mxzzz = M*dx*dz*dz*dz
    Myyyz = M*dy*dy*dy*dz
    Myyzz = M*dy*dy*dz*dz
    Myzzz = M*dy*dz*dz*dz
    Mxxyz = M*dx*dx*dy*dz
    Mxyyz = M*dx*dy*dy*dz
    Mxyzz = M*dx*dy*dz*dz

    return [ M, Mx, My, Mz, Mxx, Myy, Mzz, Mxy, Mxz, Myz, Mxxx, Myyy, Mzzz, Mxxy, Mxyy, Mxxz, Mxzz, Mxyz, Myyz, Myzz, Mxxxx, Myyyy, Mzzzz, Mxxxy, Mxxyy, Mxyyy, Mxxxz, Mxxzz, Mxzzz, Myyyz, Myyzz, Myzzz, Mxxyz, Mxyyz, Mxyzz ]

def parallel_axis_theorem(mom, dx, dy, dz):
    """
    Adjust moments for a shifted centre of mass (i.e. parallel-axis theorem)
    Stadel (2001), Eq. (2.34): M_d^_n_ = sum_m=0^n  M^(_n-m_ d^_m_)

    This one look simple, BUT it's written in a "rank" and "brace-permutation" notations!
    Everything enclosed in braces should be summed over all UNIQUE permutations.

    """
    M, Mx, My, Mz, Mxx, Myy, Mzz, Mxy, Mxz, Myz, Mxxx, Myyy, Mzzz, Mxxy, Mxyy, Mxxz, Mxzz, Mxyz, Myyz, Myzz, Mxxxx, Myyyy, Mzzzz, Mxxxy, Mxxyy, Mxyyy, Mxxxz, Mxxzz, Mxzzz, Myyyz, Myyzz, Myzzz, Mxxyz, Mxyyz, Mxyzz = mom

    # One has to begin with 4th-order --> and endup with 1st.
    # "How many times I can select e.g. dx from the indices of M?" <- see Eq. 2.9 and the example just below.
    # 2DO: check numerical factors!
    Mxxxx += 4*Mxxx*dx + 6*Mxx*dx*dx + Mx*dx*dx*dx + M*dx*dx*dx*dx
    Myyyy += 4*Myyy*dy + 6*Myy*dy*dy + My*dy*dy*dy + M*dy*dy*dy*dy
    Mzzzz += 4*Mzzz*dz + 6*Mzz*dz*dz + Mz*dz*dz*dz + M*dz*dz*dz*dz
    Mxxxy += Mxxx*dy + 3*Mxxy*dx + 3*Mxx*dx*dy + 3*Mxy*dx*dx + Mx*dx*dx*dy + M*dx*dx*dx*dy
    Mxxyy += 2*Mxxy*dy + 2*Mxyy*dx + Mxx*dy*dy + Myy*dx*dx + 4*Mxy*dx*dy + Mx*dx*dy*dy + M*dx*dx*dy*dy
    Mxyyy += 3*Mxyy*dy + Myyy*dx + 3*Mxy*dy*dy + 3*Myy*dx*dy + Mx*dy*dy*dy + M*dx*dy*dy*dy
    Mxxxz += Mxxx*dz + 3*Mxxz*dx + 3*Mxx*dx*dz + 3*Mxz*dx*dx + Mx*dx*dx*dz + M*dx*dx*dx*dz
    Mxxzz += 2*Mxxz*dz + 2*Mxzz*dx + Mxx*dz*dz + Mzz*dx*dx + 4*Mxz*dx*dz + Mx*dx*dz*dz + M*dx*dx*dz*dz
    Mxzzz += 3*Mxzz*dz + Mzzz*dx + 3*Mxz*dz*dz + 3*Mzz*dx*dz + Mx*dz*dz*dz + M*dx*dz*dz*dz
    Myyyz += Myyy*dz + 3*Myyz*dy + 3*Myy*dy*dz + 3*Myz*dy*dy + My*dy*dy*dz + M*dy*dy*dy*dz
    Myyzz += 2*Myyz*dz + 2*Myzz*dy + Myy*dz*dz + Mzz*dy*dy + 4*Myz*dy*dz + My*dy*dz*dz + M*dy*dy*dz*dz
    Myzzz += 3*Myzz*dz + Mzzz*dy + 3*Myz*dz*dz + 3*Mzz*dy*dz + My*dz*dz*dz + M*dy*dz*dz*dz
    Mxxyz += Mxxy*dz + Mxxz*dy + 2*Mxyz*dx + Mxx*dy*dz + 2*Mxy*dx*dz + 2*Mxz*dx*dy + Myz*dx*dx + Mx*dx*dy*dz + M*dx*dx*dy*dz
    Mxyyz += Mxyy*dz + 2*Mxyz*dy + Myyz*dx + 2*Mxy*dy*dz + Mxz*dy*dy + 2*Myz*dx*dy + Myy*dx*dz + Mx*dy*dy*dz + M*dx*dy*dy*dz
    Mxyzz += 2*Mxyz*dz + Mxzz*dy + Myzz*dx + Mxy*dz*dz + 2*Mxz*dy*dz + 2*Myz*dx*dz + Mzz*dx*dy + Mx*dy*dz*dz + M*dx*dy*dz*dz  # cf. permutations here!

    Mxxx += 3*Mxx*dx + Mx*dx*dx + M*dx*dx*dx
    Myyy += 3*Myy*dy + My*dy*dy + M*dy*dy*dy
    Mzzz += 3*Mzz*dz + Mz*dz*dz + M*dz*dz*dz
    Mxxy += Mxx*dy + 2*Mxy*dx + Mx*dx*dy + M*dx*dx*dy
    Mxyy += 2*Mxy*dy + Myy*dx + Mx*dy*dy + M*dx*dy*dy
    Mxxz += Mxx*dz + 2*Mxz*dx + Mx*dx*dz + M*dx*dx*dz
    Mxzz += 2*Mxz*dz + Mzz*dx + Mx*dz*dz + M*dx*dz*dz
    Mxyz += Mxy*dz + Mxz*dy + Myz*dx + Mx*dy*dz + M*dx*dy*dz  # cf. all relevant permutations here!
    Myyz += Myy*dz + 2*Myz*dy + My*dy*dz + M*dy*dy*dz
    Myzz += 2*Myz*dz + Mzz*dy + My*dz*dz + M*dy*dz*dz

    Mxx += Mx*dx + M*dx*dx
    Myy += My*dy + M*dy*dy
    Mzz += Mz*dz + M*dz*dz
    Mxy += Mx*dy + M*dx*dy
    Mxz += Mx*dz + M*dx*dz
    Myz += My*dz + M*dy*dz

    Mx += M*dx
    My += M*dy
    Mz += M*dz

    return [ M, Mx, My, Mz, Mxx, Myy, Mzz, Mxy, Mxz, Myz, Mxxx, Myyy, Mzzz, Mxxy, Mxyy, Mxxz, Mxzz, Mxyz, Myyz, Myzz, Mxxxx, Myyyy, Mzzzz, Mxxxy, Mxxyy, Mxyyy, Mxxxz, Mxxzz, Mxzzz, Myyyz, Myyzz, Myzzz, Mxxyz, Mxyyz, Mxyzz ]

def eval_gravity(mom, gam, dx, dy, dz, order=4):
    """
    Evaluate gravitational acceleration
    Adapted from ss_core/src/pkdgrav/meval.h

    """

    M, Mx, My, Mz, Mxx, Myy, Mzz, Mxy, Mxz, Myz, Mxxx, Myyy, Mzzz, Mxxy, Mxyy, Mxxz, Mxzz, Mxyz, Myyz, Myzz, Mxxxx, Myyyy, Mzzzz, Mxxxy, Mxxyy, Mxyyy, Mxxxz, Mxxzz, Mxzzz, Myyyz, Myyzz, Myzzz, Mxxyz, Mxyyz, Mxyzz = mom

    ax = 0.0
    ay = 0.0
    az = 0.0
    fPot = 0.0
    Qta = 0.0
    if order >= 4:
        Qmirx = (1.0/6.0)*(Mxzzz*dz*dz*dz + 3*Mxyzz*dy*dz*dz + 3*Mxyyz*dy*dy*dz + Mxyyy*dy*dy*dy + 3*Mxxzz*dx*dz*dz + 6*Mxxyz*dx*dy*dz + 3*Mxxyy*dx*dy*dy + 3*Mxxxz*dx*dx*dz + 3*Mxxxy*dx*dx*dy + Mxxxx*dx*dx*dx)
        Qmiry = (1.0/6.0)*(Myzzz*dz*dz*dz + 3*Mxyzz*dx*dz*dz + 3*Mxxyz*dx*dx*dz + Mxxxy*dx*dx*dx + 3*Myyzz*dy*dz*dz + 6*Mxyyz*dx*dy*dz + 3*Mxxyy*dx*dx*dy + 3*Myyyz*dy*dy*dz + 3*Mxyyy*dx*dy*dy + Myyyy*dy*dy*dy)
        Qmirz = (1.0/6.0)*(Myyyz*dy*dy*dy + 3*Mxyyz*dx*dy*dy + 3*Mxxyz*dx*dx*dy + Mxxxz*dx*dx*dx + 3*Myyzz*dy*dy*dz + 6*Mxyzz*dx*dy*dz + 3*Mxxzz*dx*dx*dz + 3*Myzzz*dy*dz*dz + 3*Mxzzz*dx*dz*dz + Mzzzz*dz*dz*dz)
        Qmir = (1.0/4.0)*(Qmirx*dx + Qmiry*dy + Qmirz*dz)
        Qxx = Mxxxx + Mxxyy + Mxxzz
        Qxy = Mxxxy + Mxyyy + Mxyzz
        Qxz = Mxxxz + Mxyyz + Mxzzz
        Qyy = Mxxyy + Myyyy + Myyzz
        Qyz = Mxxyz + Myyyz + Myzzz
        Qzz = Mxxzz + Myyzz + Mzzzz
        Qtr = (1.0/8.0)*(Qxx + Qyy + Qzz)
        Qhx = 0.5*(Qxx*dx + Qxy*dy + Qxz*dz)
        Qhy = 0.5*(Qxy*dx + Qyy*dy + Qyz*dz)
        Qhz = 0.5*(Qxz*dx + Qyz*dy + Qzz*dz)
        Qh = 0.5*(Qhx*dx + Qhy*dy + Qhz*dz)
        fPot -= gam[4]*Qmir - gam[3]*Qh + gam[2]*Qtr
        Qta += gam[5]*Qmir - gam[4]*Qh + gam[3]*Qtr
        ax += gam[4]*Qmirx - gam[3]*Qhx
        ay += gam[4]*Qmiry - gam[3]*Qhy
        az += gam[4]*Qmirz - gam[3]*Qhz
    if order >= 3:
        Qmirx = (1.0/2.0)*(Mxzz*dz*dz + 2*Mxyz*dy*dz + Mxyy*dy*dy + 2*Mxxz*dx*dz + 2*Mxxy*dx*dy + Mxxx*dx*dx)
        Qmiry = (1.0/2.0)*(Myzz*dz*dz + 2*Mxyz*dx*dz + Mxxy*dx*dx + 2*Myyz*dy*dz + 2*Mxyy*dx*dy + Myyy*dy*dy)
        Qmirz = (1.0/2.0)*(Myyz*dy*dy + 2*Mxyz*dx*dy + Mxxz*dx*dx + 2*Myzz*dy*dz + 2*Mxzz*dx*dz + Mzzz*dz*dz)
        Qmir = (1.0/3.0)*(Qmirx*dx + Qmiry*dy + Qmirz*dz)
        Qx = 0.5*(Mxxx + Mxyy + Mxzz)
        Qy = 0.5*(Mxxy + Myyy + Myzz)
        Qz = 0.5*(Mxxz + Myyz + Mzzz)
        Qtr = Qx*dx + Qy*dy + Qz*dz
        fPot -= gam[3]*Qmir - gam[2]*Qtr
        Qta += gam[4]*Qmir - gam[3]*Qtr
        ax += gam[3]*Qmirx - gam[2]*Qx
        ay += gam[3]*Qmiry - gam[2]*Qy
        az += gam[3]*Qmirz - gam[2]*Qz
    if order >= 2:
        Qmirx = (1.0/1.0)*(Mxz*dz + Mxy*dy + Mxx*dx)
        Qmiry = (1.0/1.0)*(Myz*dz + Mxy*dx + Myy*dy)
        Qmirz = (1.0/1.0)*(Myz*dy + Mxz*dx + Mzz*dz)
        Qmir = (1.0/2.0)*(Qmirx*dx + Qmiry*dy + Qmirz*dz)
        Qtr = 0.5*(Mxx + Myy + Mzz)
        fPot -= gam[2]*Qmir - gam[1]*Qtr
        Qta += gam[3]*Qmir - gam[2]*Qtr
        ax += gam[2]*Qmirx
        ay += gam[2]*Qmiry
        az += gam[2]*Qmirz
    if order >= 0:
        fPot -= gam[0]*M
        Qta += gam[1]*M
        ax -= dx*Qta
        ay -= dy*Qta
        az -= dz*Qta

    return [ax, ay, az]

class Particle(object):
    """Particle with all properties"""

    def __init__(self, r, m, i):
        """Define its radius vector, mass and index"""
        self.r = r
        self.m = m  # assuming GM units
        self.i = i

class Node(object):
    """One node of the k-dimensional tree"""

    NUMBER_OF_MOMENTS = 35

    def __init__(self, particle, left, right):
        """Define particle, left branch, right branch and other quantities (to be computed)"""
        self.particle = particle
        self.left = left
        self.right = right
        self.tm = None          # total mass
        self.cm = [None] * 3    # centre of mass
        self.sqsize = None      # squared size
        self.mom = [None] * self.NUMBER_OF_MOMENTS  # gravitational moments

class Tree(object):
    """A k-dimensional tree"""

    def __init__(self, k, particles):
        """Build a k-dimensional tree from a list of particles"""
        self.k = k
        self.particles = particles

        def build_tree(particles, depth=0):
            """Recursive build function"""
            if len(particles) == 0:
                return None
 
            axis = depth % self.k
            particles.sort(key=lambda particle: particle.r[axis])  # this is time sort with O(N log N)
            i = len(particles) // 2  # median index
 
            return Node( \
                particles[i], \
                build_tree(particles[:i], depth+1), \
                build_tree(particles[i+1:], depth+1), \
            )

        self.root = build_tree(particles)

    def nearest_neighbor(self, destination):
        """Nearest-neighbour search O(N log N)"""
        best = [None, float('inf')]  # particle & squared distance

        def recursive_search(node, depth=0):
            """Recursive search function"""
            if node is None:
                return

            particle, left, right = node.particle, node.left, node.right
            node_sqdist = square_distance(particle.r, destination)
            if node_sqdist < best[1]:
                best[:] = particle, node_sqdist  # in-place!

            axis = depth % self.k
            diff = destination[axis] - particle.r[axis]
            close, away = (left, right) if diff <= 0.0 else (right, left)

            recursive_search(close, depth+1)
            if diff**2 < best[1]:
                recursive_search(away, depth+1)

        recursive_search(self.root)

        return best[0], math.sqrt(best[1])

    def compute_mass(self):
        """Compute total masses of ALL tree nodes"""

        def recursive_mass(node):
            """Recursive mass function"""
            tm = node.particle.m
            for branch in node.left, node.right:
                if branch is not None:
                   tm += recursive_mass(branch)

            node.tm = tm
            return tm  # store in the tree and return it too!

        recursive_mass(self.root)
        
    def compute_centre(self):
        """
        Compute centres of masses for ALL tree nodes.
        The method self.compute_mass() has to be called first!

        """
        def recursive_centre(node):
            """Recursive centre-of-mass function"""
            cm = list(node.particle.r)

            for branch in node.left, node.right:
                if branch is not None:
                    centre = recursive_centre(branch)
                    for i in range(0, self.k):
                        cm[i] += branch.tm*centre[i]
            for i in range(0, self.k):
               cm[i] /= node.tm

            node.cm[:] = cm
            return cm  # store in the tree and return it too!

        recursive_centre(self.root)

    def compute_multipole(self, node=None):
        """
        Compute quadru-, octu- and hexadecapole moments for ALL tree nodes.
        The method self.compute_centre() has to be called first!

        """
        if node == None:
            node = self.root
        d = [0.0] * 3
        multipole_d = [0.0] * Node.NUMBER_OF_MOMENTS

        def recursive_multipole(node):
            """Recursive multipole function"""
            m = node.particle.m
            for i in range(0, self.k):
                d[i] = node.particle.r[i] - node.cm[i]

            mom = eval_multipole(m, d[0], d[1], d[2])

            for branch in node.left, node.right:
                if branch is not None:
                    multipole = recursive_multipole(branch)

                    # offset between the two centres of masses && parallel-axis theorem
                    for i in range(0, self.k):
                        d[i] = branch.cm[i] - node.cm[i]

                    multipole_offset = parallel_axis_theorem(multipole, d[0], d[1], d[2])

                    for i in range(0, len(mom)):
                        mom[i] += multipole_offset[i]

            node.mom[:] = mom  # in-place!
            return mom  # store in the tree and return it too!

        recursive_multipole(node)
        
    def compute_size(self):
        """Compute sizes corresponding to ALL tree nodes"""

        def recursive_size(node):
            minx = node.particle.r[0]
            miny = node.particle.r[1]
            maxx = minx
            maxy = miny
            for branch in node.left, node.right:
                if branch is not None:
                    x1, x2, y1, y2 = recursive_size(branch)
                    minx = min(minx, x1)
                    maxx = max(maxx, x2)
                    miny = min(miny, y1)
                    maxy = max(maxy, y2)

            sqsize = (maxx-minx)**2+(maxy-miny)**2  # Well, there are other estimates possible...
            node.sqsize = sqsize
            return minx, maxx, miny, maxy

        recursive_size(self.root)

    def compute_gravity(self, destination, phi=0.5):
        """
        Compute gravitational acceleration acting on destination ([x, y, z] vector);
        phi is the opening angle [in radians]

        """
        phi2 = phi*phi
        eps = 1.e-8

        def recursive_gravity(node):
            """Recursive gravity computation"""
            dist2 = square_distance(node.cm, destination)
            ag = [0.0] * self.k

            # do NOT open this node at all, because it's too small or far away
            if node.sqsize/(dist2+eps) < phi2:
                if dist2 > eps:
                    tmp = node.tm/(dist2*math.sqrt(dist2))
                    for i in range(0, self.k):
                        ag[i] = tmp*(node.cm[i]-destination[i])
                else:
                    ag = [0.0] * self.k  # no self-interactions

            # open this node and sum the corresponding particle and two branches
            else:
                dist2 = square_distance(node.particle.r, destination)
                if dist2 > eps:
                    tmp = node.particle.m/(dist2*math.sqrt(dist2))
                    for i in range(0, self.k):
                        ag[i] = tmp*(node.particle.r[i]-destination[i])
                else:
                    ag = [0.0] * self.k  # no self-interactions

                for branch in node.left, node.right:
                   if branch is not None:
                      accel = recursive_gravity(branch)
                      for i in range(0, self.k):
                         ag[i] += accel[i]
            return ag

        return recursive_gravity(self.root)

    def compute_gravity2(self, destination, phi=0.5, order=4):
        """
        Compute gravitational acceleration acting on destination ([x, y, z] vector);
        phi is the opening angle [in radians]

        """
        phi2 = phi*phi
        eps = 1.e-8

        def recursive_gravity(node):
            """Recursive gravity computation"""
            dist2 = square_distance(node.cm, destination)
            ag = [0.0] * self.k
            d = [0.0] * 3
            gam = [0.0] * 6

            # do NOT open this node at all, because it's too small or far away
            if node.sqsize/(dist2+eps) < phi2:

                gam[0] = 1.0/math.sqrt(dist2)
                dir2 = gam[0]*gam[0]
                gam[1] = gam[0]*dir2
                gam[2] = 3*gam[1]*dir2
                gam[3] = 5*gam[2]*dir2
                gam[4] = 7*gam[3]*dir2
                gam[5] = 9*gam[4]*dir2

                for i in range(0, self.k):
                    d[i] = destination[i] - node.cm[i]

                ag = eval_gravity(node.mom, gam, d[0], d[1], d[2], order=order)

            # open this node and sum the corresponding particle and two branches
            else:
                dist2 = square_distance(node.particle.r, destination)
                if dist2 > eps:
                    tmp = node.particle.m/(dist2*math.sqrt(dist2))
                    for i in range(0, self.k):
                        ag[i] = tmp*(node.particle.r[i]-destination[i])
                else:
                    ag = [0.0] * self.k  # no self-interactions

                for branch in node.left, node.right:
                   if branch is not None:
                      accel = recursive_gravity(branch)
                      for i in range(0, self.k):
                         ag[i] += accel[i]
            return ag

        return recursive_gravity(self.root)

    def brute_force(self, destination):
        """Compute gravity using a brute-force algorithm O(N^2)"""
        ag = [0.0] * self.k
        eps = 1.e-8
        for particle in self.particles:
            dist2 = square_distance(particle.r, destination)
            if dist2 > eps:
                tmp = particle.m/(dist2*math.sqrt(dist2))
                for i in range(0, self.k):
                    ag[i] += tmp*(particle.r[i]-destination[i])
        return ag

    def brute_mass(self, particles):
        """Sum the total mass"""
        tm = 0.0
        for particle in particles:
           tm += particle.m
        return tm

    def brute_centre(self, particles):
        """Compute centre of mass (and total mass too)"""
        tm = 0.0
        cm = [0.0] * self.k
        for particle in particles:
           tm += particle.m
           for i in range(0, self.k):
               cm[i] += particle.m*particle.r[i]
        for i in range(0, self.k):
           cm[i] /= tm
        return tm, cm

    def list_particles(self, node=None):
        """Get a list of particles under a given node"""
        if node == None:
            node = self.root
        particles = []

        def recursive_particles(node):
            particles.append(node.particle)
            for branch in node.left, node.right:
                if branch is not None:
                    recursive_particles(branch)

        recursive_particles(node)
        return particles

    def brute_multipole(self, node=None):
        """Compute gravitational moments for particles under a given node"""
        particles = self.list_particles(node)
        tm, cm = self.brute_centre(particles)
        d = [0.0] * 3
        mom = [0.0] * Node.NUMBER_OF_MOMENTS

        for particle in particles:
            m = particle.m
            for i in range(0, self.k):
                d[i] = particle.r[i] - cm[i]

            multipole = eval_multipole(m, d[0], d[1], d[2])

            for i in range(0, len(multipole)):
                mom[i] += multipole[i]

        return tm, cm, mom

    def prn(self, node=None, depth=0):
        """Print a tree with depths"""
        if node == None:
            node = self.root
        print("depth " + str(depth) + " ... " + str(node.particle.r) + "   total_mass = " + str(node.tm) + "   centre_of_mass = " + str(node.cm) + "   sqsize = " + str(node.sqsize))
        for branch in node.left, node.right:
            if branch is not None:
                self.prn(branch, depth+1)
    
    def prn_multipole(self, node=None, depth=0):
        """Print a tree with multipole moments"""
        if node == None:
            node = self.root

        for order in range(1,5):
            sys.stdout.write("depth " + str(depth) + " ... " + str(node.particle.r) + "   ")
            if order == 1:
                print("dipole       = " + str(node.mom[1:4]))
            elif order == 2:
                print("quadrupole   = " + str(node.mom[4:10]))
            elif order == 3:
                print("octupole     = " + str(node.mom[10:20]))
            elif order == 4:
                print("hexadecapole = " + str(node.mom[20:35]))
        print("")

        for branch in node.left, node.right:
            if branch is not None:
                self.prn_multipole(branch, depth+1)
    
    def plot(self):
        """Plot a 2D tree in a similar way as in Wikipedia"""
        import matplotlib.pyplot as plt
        points = []
        lines = []
        colors = []

        def recursive_plot(node, depth=0, x1=0, x2=10, y1=0, y2=10):
            """Recursive plot function"""
            if node is None:
                return
            r, left, right = node.particle.r, node.left, node.right
            axis = depth % self.k
            points.append(r)
            if axis == 0:
                lines.append([[r[axis], y1], [r[axis], y2]])
                colors.append('r')
                recursive_plot(left, depth+1, x1, r[axis], y1, y2)
                recursive_plot(right, depth+1, r[axis], x2, y1, y2)
            else:
                lines.append([[x1, r[axis]], [x2, r[axis]]])
                colors.append('b')
                recursive_plot(left, depth+1, x1, x2, y1, r[axis])
                recursive_plot(right, depth+1, x1, x2, r[axis], y2)
            
        recursive_plot(self.root)

        for i in range(0, len(lines)):
            plt.plot(*zip(*lines[i]), marker='', color=colors[i], ls='-')
        plt.plot(*zip(*points), marker='o', color='y', ls='')

        plt.axis([0,10,0,10])
        plt.xlabel('$x$', fontsize=10)
        plt.ylabel('$y$', fontsize=10)
        plt.axes().set_aspect('equal')
        plt.rcParams.update({'font.size': 18})
        plt.savefig('tree.png', bbox_inches='tight')

def main():
    """Example usage"""
    points = [(2, 3), (5, 4), (9, 6), (4, 7), (8, 1), (7, 2)]  # Wikipedia example
    particles = []
    n = len(points)
    for i in range(0,n):
        particles.append(Particle([points[i][0], points[i][1]], 1.0, i))

    tree = Tree(2, particles)
    tree.plot()

#    point = (8, 5, 0)  # inside
    point = (-10, 10, 0)  # outside
    particle, distance = tree.nearest_neighbor(point)
    print("particle.r = " + str(particle.r))
    print("distance = " + str(distance))

    tree.compute_mass()
    tree.compute_centre()
    tree.compute_size()
    tree.prn()
    print("")

    tm, cm, mom = tree.brute_multipole()
    print("# Gravitational moments computed by brute-force:")
    print("tm = " + str(tm))
    print("cm = " + str(cm))
    print("dipole       = " + str(mom[0:3]))
    print("quadrupole   = " + str(mom[3:9]))
    print("octupole     = " + str(mom[9:19]))
    print("hexadecapole = " + str(mom[19:34]))
    print("")

    tree.compute_multipole()
    tree.prn_multipole()

    ag0 = tree.brute_force(point)
    ag1 = tree.compute_gravity(point, phi=0.0)
    ag2 = tree.compute_gravity(point, phi=0.5)

    ag3 = tree.compute_gravity2(point, phi=0.5, order=0)
    ag4 = tree.compute_gravity2(point, phi=0.5, order=2)
    ag5 = tree.compute_gravity2(point, phi=0.5, order=3)
    ag6 = tree.compute_gravity2(point, phi=0.5, order=4)

    print("ag0 = " + str(ag0) + " <- brute-force algorithm")
    print("ag1 = " + str(ag1) + " <- kd-tree algorithm (phi = 0, i.e. a sum)")
    print("ag2 = " + str(ag2) + " <- kd-tree algorithm (monopoles only)")
    print("")
    print("ag3 = " + str(ag3) + " <- kd-tree algorithm (monopoles)")
    print("ag4 = " + str(ag4) + " <- kd-tree algorithm (quadrupoles)")
    print("ag5 = " + str(ag5) + " <- kd-tree algorithm (octupoles)")
    print("ag6 = " + str(ag6) + " <- kd-tree algorithm (hexadecapoles)")
    print("")

if __name__ == '__main__':
    main()


