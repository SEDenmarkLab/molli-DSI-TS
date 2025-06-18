from __future__ import annotations
import numpy as np
import re
from typing import *  # pylint: disable=unused-wildcard-import
from copy import deepcopy

from ..ftypes.xyz import parse_xyz

# pylint: disable=no-member

_HDR_PATTERN = re.compile(r"(?P<c>[0-9]+),(?P<r>[0-9]+)#(?P<g>.+)")

def find_centroid(coords):
    print("CENTROIIIID", coords)
    x=input()
    x_sum = sum(x for x, _, _ in coords)
    y_sum = sum(y for _, y, _ in coords)
    z_sum = sum(z for _, _, z in coords)
    count = len(coords)
    
    return float(x_sum/count, y_sum/count, z_sum/count)


def fit_plane(points):
    centroid = np.mean(points, axis=0)
    centered_points = points - centroid
    u, s, vh = np.linalg.svd(centered_points) # singular value decomposition generalizes eigenvalue/vectors to non-square matrices
    normal = vh[-1]
    return normal/np.linalg.norm(normal) # returns normal to plane
    
def angle_between_normals(normal1, normal2):

    cos_angle = np.dot(normal1, normal2)
    angle_rad = np.arccos(np.clip(cos_angle, -1.0, 1.0))
    angle_deg = np.degrees(angle_rad)
    return angle_deg
def deg_rad_conv(degrees):
    return float(degrees*0.0174533)

#rotate subgeometry (fragment)
def sub_rotation_matrix(axis, theta):

    #print("AXIS GOING IN", axis)
    axis = np.asarray(axis)
    axis = axis/np.sqrt(np.dot(axis,axis))
    a = np.cos(theta/2)
    #print(-axis)
    #print(np.sin(theta/2.0))
    b, c, d = -axis*np.sin(theta/2.0)
    aa, bb, cc, dd = a*a, b*b, c*c, d*d
    bc, ad, ac, ab, bd, cd = b * c, a * d, a * c, a * b, b * d, c * d
    
    return np.array([
        [aa + bb - cc - dd, 2 * (bc + ad), 2 * (bd - ac)],
        [2 * (bc - ad), aa + cc - bb - dd, 2 * (cd + ab)],
        [2 * (bd + ac), 2 * (cd - ab), aa + dd - bb - cc]
    ])

def sub_rotate_around_line(matrix, point1, point2, theta):
    """
    Rotate a matrix of points counterclockwise by theta radians around a line through
    point1 and point2. Each row in the matrix is a point.
    """
    # Define the axis of rotation
    axis = np.asarray(point2) - np.asarray(point1)
    
    # Apply the rotation
    #print(axis, theta)
    R = sub_rotation_matrix(axis, theta)
    rotated_matrix = np.dot(R, (matrix - point1).T).T + point1
    return rotated_matrix

def rotation_matrix(v1, v2, tol=1.0e-6):
    """
    Rotation Matrix (vector-to-vector definition)
    ---

    Computes a 3x3 rotation matrix that transforms v1/|v1| -> v2/|v2|
    tol detects a situation where dot(v1, v2) ~ -1.0
    and returns a diag(-1,-1,1) matrix instead. NOTE [Is this correct?!]
    This may be improved by figuring a more correct asymptotic behavior, but works for the present purposes.

    returns matrix [R], that satisfies the following equation:
    v1 @ [R] / |v1| == v2 / |v2|

    Inspired by
    https://en.wikipedia.org/wiki/Rotation_matrix
    https://math.stackexchange.com/questions/180418/calculate-rotation-matrix-to-align-vector-a-to-vector-b-in-3d

    """
    if not v1.size == v2.size == 3:
        raise ValueError("Vectors must have a size of 3")

    v1n = np.array(v1) / np.linalg.norm(v1)
    v2n = np.array(v2) / np.linalg.norm(v2)

    I = np.eye(3)

    Ux = np.outer(v1n, v2n) - np.outer(v2n, v1n)
    c = np.dot(v1n, v2n)

    if c <= -1 + tol:
        return np.diag([-1, -1, 1])

    else:
        return I + Ux + Ux @ Ux / (1 + c)


def distance(p1, p2):
    """
    Returns Euclidean distance between two points
    """
    _p1, _p2 = np.array(p1), np.array(p2)
    return np.sqrt(np.sum((_p1 - _p2) ** 2))


class CartesianGeometry:
    """
    Class that enables molecular geometry transformations.
    These are primary based on `numpy` matrix operations.
    Especially useful in rotations, translations, alignments etc

    *This class does not handle atom types.*
    """

    def __init__(
        self,
        coord: np.ndarray = None,
        precision=4,  # decimal places for geometry printouts
        unit="Angstrom",
        dtype=np.float32,
    ):
        """
        Docstring
        """
        self.precision = precision
        self.unit = unit
        if len(coord):
            _coord = np.array(coord, dtype=dtype)
            if _coord.shape[1] != 3:
                raise ValueError("")

            self.coord = _coord
            self.N = self.coord.shape[0]

    def center_geom(self):
        """
        Centers the molecule on it's geometrical center
        """
        centroid = np.average(self.coord, axis=(0,))
        self.translate(-1 * centroid)

    ### KSP MODIFIED
    def get_centroid(self):
        return np.average(self.coord, axis=(0,))
    ### \KSP MODIFIED
        
    def translate(self, vector: np.ndarray):
        """
        Translate the molecule by a given vector
        """
        assert vector.shape == (3,)
        for c in self.coord:
            c += vector

    def set_origin(self, i: int):
        """
        Translates the origin to a given coordinate
        """
        self.translate(-1 * self.coord[i])

    def scale(self, factor):
        """
        Simple multiplication of all coordinates by a factor.
        Useful for unit conversion.
        """
        self.coord *= factor

    def transform(self, matrix: np.ndarray):
        """
        Applies an in place transformation matrix to all coordinates.
        The only (!) check performed is that the matrix must be a 3x3.
        Normalization is NOT checked. Note that this transformation may
        result in inversion and/or inappropriate scaling.

        User assumes the responsibility to check the validity of transform matrix.

        NOTE: Add checks for trace and determinant?
        """
        if not matrix.shape == (3, 3):
            raise ValueError("Incorrect matrix received. Must be [3x3]")
        self.coord = self.coord @ matrix
        

    def randomize(self, std=0.1):
        """
        Each coordinate receives a random vector deviation.
        This should be very useful in optimization to eliminate flat hard-to-escape geometries
        """
        self.coord += np.random.normal(0, std, self.coord.shape)

    def delete(self, idx: int):
        """
        Delete a selected atom from the geometry
        """
        self.coord = np.delete(self.coord, idx, 0)
        self.N -= 1

    def get_distance(self, idx1: int, idx2: int):
        """
        Measure the Euclidean distance between two points
        """
        p1 = self.coord[idx1]
        p2 = self.coord[idx2]

        return distance(p1, p2)

    def get_angle(self, idx1: int, idx2: int, idx3: int):
        """
        Measure the angle
        """
        v1 = self.coord[idx1] - self.coord[idx2]
        v2 = self.coord[idx3] - self.coord[idx2]
        dt = np.dot(v1, v2) / (np.linalg.norm(v1) * np.linalg.norm(v2))
        return np.arccos(dt)

    def get_coord(self, idx: int):
        return self.coord[idx]

    def dumps(self):
        """
        Encode the geometry as string
        """
        gs = f"#{self.N},3:"

        for x, y, z in self.coord:
            p = self.precision
            gs += f"{x:0.{p}f},{y:0.{p}f},{z:0.{p}f};"

        return gs

    def bounding_box(self):
        """
        Return 6 coordinates that determine the rectangular space which fits all the atoms
        `return (xmin, ymin, zmin), (xmin, ymin, zmin)`
        """
        rmin = np.min(self.coord, axis=(0,))
        rmax = np.max(self.coord, axis=(0,))

        return rmin, rmax

    @classmethod
    def from_str(cls, s: str, u: str = "A"):
        """
        Decode geometry string
        """

        if u != "A":
            raise NotImplementedError

        m = re.match(r"#(?P<L>[0-9]+),(?P<D>[0-9]+):(?P<G>.+);", s)

        L = int(m.group("L"))
        D = int(m.group("D"))
        G = m.group("G")

        assert D == 3, "Only 3d coordinates supported for now"

        coord = []

        for xyz in G.split(";"):
            x, y, z = xyz.split(",")
            coord.append([float(x), float(y), float(z)])

        assert len(coord) == L, f"Expected {L} coords, found {len(coord)}"

        return cls(coord)

    @classmethod
    def from_xyz(cls, xyzs: str) -> (CartesianGeometry, List, str):
        """
        Parse an xyz list of lines
        if multixyz, returns a list of all individual geometries
        """
        parsed = parse_xyz(xyzs, single=False, assert_single=False)

        if len(parsed) == 1:
            coord, atoms, cmt = parsed[0]
            return [(cls(coord=coord), atoms, cmt)]

        elif len(parsed) < 1:
            raise SyntaxError

        else:
            res = []
            for coord, atoms, cmt in parsed:
                res.append((cls(coord=coord), atoms, cmt))
            return res

    def to_xyz(self, atoms: List, comment: str = None):
        """
        Create an xyz file string out of the geometry and atom list
        """
        N = self.coord.shape[0]
        assert N == len(atoms)

        res = f"{N}\n{comment if comment else 'Produced by molli'}\n"
        for i, xyz in enumerate(self.coord):
            x, y, z = xyz
            res += f"{atoms[i]} {x:>10.4f} {y:>10.4f} {z:>10.4f}\n"
        return res

    def clone_geometry(self, other: CartesianGeometry):
        self.coord = deepcopy(other.coord)
