import numpy as np
from numpy.linalg import norm
import math
from stl import Mesh, Mode
import argparse
from glob import glob
import os

import xml.etree.ElementTree as ET

el_list = ['H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al', 'Si', 'P',
           'S', 'Cl', 'Ar', 'K', 'Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn',
           'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr', 'Rb', 'Sr', 'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh',
           'Pd', 'Ag', 'Cd', 'In', 'Sn', 'Sb', 'Te', 'I', 'Xe', 'Cs', 'Ba', 'La', 'Ce', 'Pr', 'Nd',
           'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb', 'Lu', 'Hf', 'Ta', 'W', 'Re',
           'Os', 'Ir', 'Pt', 'Au', 'Hg', 'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn', 'Fr', 'Ra', 'Ac', 'Th',
           'Pa', 'U', 'Np', 'Pu', 'Am', 'Cm', 'Bk', 'Cf', 'Es', 'Fm', 'Md', 'No', 'Lr', 'Rf', 'Db',
           'Sg', 'Bh', 'Hs', 'Mt', 'Ds', 'Rg', 'Cn', 'Nh', 'Fl', 'Mc', 'Lv', 'Ts', 'Og']
#
# atomic_radia = [53, 120, 145, 105, 85, 70, 65, 60, 50, 160, 180, 150, 125, 110, 100, 100, 100, 71,
#                 220, 180, 160, 140, 135, 140, 140, 140, 135, 135, 135, 135, 130, 125, 115, 115, 115,
#                 None, 235, 200, 180, 155, 145, 145, 135, 130, 135, 140, 160, 155, 155, 145, 145, 140,
#                 140, None, 265, 215, 195, 185, 185, 185, 185, 185, 185, 180, 175, 175, 175, 175, 175,
#                 175, 175, 155, 145, 135, 135, 130, 135, 135, 135, 150, 190, 180, 160, 190, None, None,
#                 None, 215, 195, 180, 180, 175, 175, 175, 175, 176, None, None, None, None, None, None,
#                 None, None, None, None, None, None, None, None, None, None, None, None, None, None, None,
#                 None]

# from https://en.wikipedia.org/wiki/Atomic_radii_of_the_elements_(data_page)
# vdw for H was changed
vdw_radia = [90, 140, 182, 153, 192, 170, 155, 152, 147, 154, 227, 173, 184, 210, 180, 180, 175, 188, 275, 231,
             211, None, None, None, None, None, None, 163, 140, 139, 187, 211, 185, 190, 185, 202, 303, 249, None,
             None, None, None, None, None, None, 163, 172, 158, 193, 217, 206, 206, 198, 216, 343, 268, None, None,
             None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None,
             None, None, 175, 166, 155, 196, 202, 207, 197, 202, 220, 348, 283, None, None, None, 186, None, None,
             None, None, None, None, None, None, None, None, None, None, None,None, None, None, None,None, None, None,
             None, None, None, None, None, None]

el_radia_dict = dict(zip(el_list, vdw_radia))


class Atom(object):

    def __init__(self, id, element_type: str, position: np.ndarray):
        self.id = id
        self.element_type = element_type
        self.position = position

    def size(self):
        f = 0.003
        if self.element_type in el_radia_dict:
            if el_radia_dict[self.element_type] is None:
                return el_radia_dict['C'] * f

            return el_radia_dict[self.element_type] * f

        raise ValueError('Unknown atom')


class Bond(object):
    def __init__(self, atom1: int, atom2: int, order: int = 1):
        self.atom1 = atom1
        self.atom2 = atom2
        self.order = order


class Molecule(object):
    def __init__(self, atoms: list, bonds: list):
        self.atoms = atoms
        self.bonds = bonds


def parse_CML_file(filename):
    with open(filename) as f:
        xml_text = f.read()

    atoms = []
    bonds = []

    root = ET.fromstring(xml_text)

    for ar in root.iter('atomArray'):
        for atom in ar.iter('atom'):
            attr = atom.attrib
            new_atom = Atom(attr['id'],
                            attr['elementType'],
                            np.asarray([float(attr['x3']),  float(attr['y3']), float(attr['z3'])], dtype=np.float64))
            atoms.append(new_atom)

    for ar in root.iter('bondArray'):
        for bond in ar.iter('bond'):
            attr = bond.attrib
            at_refs = attr['atomRefs2']
            a1, a2 = at_refs.split(' ')
            atom1 = int(a1[1:]) - 1  # indexes in atom list
            atom2 = int(a2[1:]) - 1
            new_bond = Bond(atom1, atom2, int(attr['order']))
            bonds.append(new_bond)

    return Molecule(atoms, bonds)


def form_mesh(vertices: np.ndarray, faces: np.ndarray):
    obj = Mesh(np.zeros(faces.shape[0], dtype=Mesh.dtype))
    for i, f in enumerate(faces):
        obj.vectors[i][0] = vertices[f[0]]
        obj.vectors[i][1] = vertices[f[1]]
        obj.vectors[i][2] = vertices[f[2]]

    return obj


# from https://observablehq.com/@mourner/fast-icosphere-mesh
def subdivide(vertices: np.ndarray, faces: np.ndarray, order=1):
    if order == 0:
        return vertices, faces

    if order > 1:
        # recursively subdivide
        for i in range(order):
            vertices, faces = subdivide(vertices, faces, 1)

        return vertices, faces

    new_vertices = list(vertices)
    new_faces = list(faces)

    v = len(new_vertices)

    midcache = {}

    def add_midpoint(a: int, b: int):
        nonlocal v
        key = math.floor((a+b) * (a+b+1)) + (a if a < b else b)
        if key in midcache:
            return midcache.pop(key)
        midcache[key] = v
        new_vertices.append((new_vertices[a] + new_vertices[b]) / 2)
        idx = v
        v += 1
        return idx

    for face in faces:
        v1 = face[0]
        v2 = face[1]
        v3 = face[2]

        a = add_midpoint(v1, v2)
        b = add_midpoint(v2, v3)
        c = add_midpoint(v1, v3)

        new_faces.append([v1, a, c])
        new_faces.append([v2, b, a])
        new_faces.append([v3, c, b])
        new_faces.append([a, b, c])

    new_vertices = np.asarray(new_vertices)
    new_faces = np.asarray(new_faces)

    return new_vertices, new_faces


# inspiration from https://github.com/PyMesh/PyMesh/blob/main/python/pymesh/meshutils/generate_icosphere.py
# https://medium.com/@oscarsc/four-ways-to-create-a-mesh-for-a-sphere-d7956b825db4
def generate_icosphere(radius, center, refinement_order=3):
    """ Generate icosphere (subdivision surface of a regular `icosahedron`_).
    Args:
        radius (``float``): Radius of icosphere.
        # center (``numpy.ndarray``): Sphere center.
        refinement_order (``int``): (optional) Number of refinement.
    Returns:
        The (possibly refined) icosphere :py:class:`Mesh`.
    .. _icosahedron: http://mathworld.wolfram.com/Icosahedron.html
    """
    r = (1.0 + math.sqrt(5.0)) / 2.0
    vertices = np.asarray([
        [-1.0,   r, 0.0],
        [ 1.0,   r, 0.0],
        [-1.0,  -r, 0.0],
        [ 1.0,  -r, 0.0],
        [0.0, -1.0,   r],
        [0.0,  1.0,   r],
        [0.0, -1.0,  -r],
        [0.0,  1.0,  -r],
        [  r, 0.0, -1.0],
        [  r, 0.0,  1.0],
        [ -r, 0.0, -1.0],
        [ -r, 0.0,  1.0],
        ], dtype=np.float64)

    faces = np.asarray([
        [0, 11, 5],
        [0, 5, 1],
        [0, 1, 7],
        [0, 7, 10],
        [0, 10, 11],
        [1, 5, 9],
        [5, 11, 4],
        [11, 10, 2],
        [10, 7, 6],
        [7, 1, 8],
        [3, 9, 4],
        [3, 4, 2],
        [3, 2, 6],
        [3, 6, 8],
        [3, 8, 9],
        [5, 4, 9],
        [2, 4, 11],
        [6, 2, 10],
        [8, 6, 7],
        [9, 8, 1],
        ], dtype=int)

    vertices, faces = subdivide(vertices, faces, refinement_order)

    # normalize vertices
    vert_length = np.sqrt((vertices * vertices).sum(axis=1, keepdims=True))
    vertices = vertices / vert_length * radius + center

    return form_mesh(vertices, faces)


# from https://github.com/PyMesh/PyMesh/blob/main/python/pymesh/meshutils/generate_cylinder.py
def generate_cylinder(p0, p1, radius, num_segments=32):

    Z = np.asarray([0, 0, 1])

    p0 = np.array(p0, dtype=np.float64)
    p1 = np.array(p1, dtype=np.float64)

    angles = [2 * np.pi * i / float(num_segments) for i in range(num_segments)]
    rim = np.array([[math.cos(theta), math.sin(theta), 0.0] for theta in angles])

    axis = p1 - p0
    axis /= norm(axis)

    rot_axis = np.cross(Z, axis)  # rotational axis

    rot_matrix = Mesh.rotation_matrix(rot_axis, np.arccos(Z.dot(axis)))

    bottom_rim = rot_matrix.dot(rim.T).T * radius + p0
    top_rim = rot_matrix.dot(rim.T).T * radius + p1

    vertices = np.vstack([bottom_rim, top_rim])

    faces = np.array([
        [[i, (i + 1) % num_segments, i + num_segments],
         [i + num_segments, (i + 1) % num_segments, (i + 1) % num_segments + num_segments]]
        for i in range(num_segments)], dtype=int)
    faces = faces.reshape((-1, 3), order="C")

    return form_mesh(vertices, faces)


# def gen_cube(center=(0, 0, 0)):
#     vertices = np.array([
#         [-1, -1, -1],
#         [+1, -1, -1],
#         [+1, +1, -1],
#         [-1, +1, -1],
#         [-1, -1, +1],
#         [+1, -1, +1],
#         [+1, +1, +1],
#         [-1, +1, +1]], dtype=np.float64) + np.asarray(center)
#
#     # Define the 12 triangles composing the cube
#     faces = np.array([
#         [0, 3, 1],
#         [1, 3, 2],
#         [0, 4, 7],
#         [0, 7, 3],
#         [4, 5, 6],
#         [4, 6, 7],
#         [5, 1, 2],
#         [5, 2, 6],
#         [2, 3, 6],
#         [3, 7, 6],
#         [0, 1, 5],
#         [0, 5, 4]], dtype=int)
#
#     return form_mesh(vertices, faces)


def combine_objects(objs: list):
    return Mesh(np.concatenate([obj.data for obj in objs]))


def create_mol_model(mol: Molecule, atom_size=1, bond1_radius=1, bond2_radius=0.65, bond3_radius=0.45,
                     bond_num_segments=32, sphere_refinement_order=3, all_bonds_single=False, bd2=0.14, bd3=0.17,
                     **kwargs):
    """

    Creates the mesh for the whole molecule including bonds.


    :param mol: Molecule.
    :param atom_size: Relative factor that multiplies the van der Walls radius of each atom.
    :param bond1_radius: Relative factor that multiplies the radius of cylinder of single bond.
    :param bond2_radius: Relative factor that multiplies the radius of both cylinders in double bond.
    :param bond3_radius: Relative factor that multiplies the radius of all cylinders in triple bond.
    :param bond_num_segments: Number of segments the bond is composed of.
    :param sphere_refinement_order: The order of subdivisions of triangles for sphere generation.
    :param all_bonds_single: If True, all bonds will be single bonds regarding the actual order.
    :param bd2: Distance between cylinders in double bond.
    :param bd3: Distance between cylinders in triple bond.
    :return: Mesh model.
    """

    bond_factor = 0.17

    bond1_radius *= bond_factor
    bond2_radius *= bond_factor
    bond3_radius *= bond_factor

    objects = []
    for atom in mol.atoms:
        atom_obj = generate_icosphere(atom.size() * atom_size, atom.position, sphere_refinement_order)
        objects.append(atom_obj)

    for bond in mol.bonds:
        a1pos = mol.atoms[bond.atom1].position
        a2pos = mol.atoms[bond.atom2].position

        Z_ax = np.asarray([0, 0, 1])

        b_direct = generate_cylinder(a1pos, a2pos, bond1_radius, bond_num_segments)
        if not all_bonds_single:
            bond_ax = a2pos - a1pos
            perp_ax = np.cross(bond_ax, Z_ax)
            perp_ax /= norm(perp_ax)  # perpendicular axis to bond axis

            if bond.order == 1:
                objects.append(b_direct)
            elif bond.order == 2:
                v_diff = perp_ax * bd2
                b_up = generate_cylinder(a1pos + v_diff, a2pos + v_diff, bond2_radius, bond_num_segments)
                b_down = generate_cylinder(a1pos - v_diff, a2pos - v_diff, bond2_radius, bond_num_segments)
                objects += [b_up, b_down]
            elif bond.order == 3:
                v_diff = perp_ax * bd3
                b_direct = generate_cylinder(a1pos, a2pos, bond3_radius, bond_num_segments)
                b_up = generate_cylinder(a1pos + v_diff, a2pos + v_diff, bond3_radius, bond_num_segments)
                b_down = generate_cylinder(a1pos - v_diff, a2pos - v_diff, bond3_radius, bond_num_segments)
                objects += [b_direct, b_up, b_down]
        else:
            objects.append(b_direct)

    return combine_objects(objects)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("--atom_size", nargs="?", default=1.0, type=float,
                        help="Relative factor that multiplies the van der Walls radius of each atom. Default 1.")
    parser.add_argument("--bond1_radius", nargs="?", default=1.0, type=float,
                        help="Relative factor that multiplies the radius of cylinder of single bond. Default 1.")
    parser.add_argument("--bond2_radius", nargs="?", default=0.65, type=float,
                        help="Relative factor that multiplies the radius of both cylinders in double bond. Default 0.65.")
    parser.add_argument("--bond3_radius", nargs="?", default=0.45, type=float,
                        help="Relative factor that multiplies the radius of all cylinders in triple bond. Default 0.45.")
    parser.add_argument("--bond_num_segments", nargs="?", default=32, type=int,
                        help="Number of segments the bond is composed of. Default 32.")
    parser.add_argument("--sphere_refinement_order", nargs="?", default=3, type=int,
                        help="The order of subdivisions of triangles for sphere generation. Default 3.")
    parser.add_argument("--bd2", nargs="?", default=0.14, type=float,
                        help="Distance between cylinders in double bond. Default 0.14.")
    parser.add_argument("--bd3", nargs="?", default=0.17, type=float,
                        help="Distance between cylinders in triple bond. Default 0.17.")
    parser.add_argument("--all_bonds_single", action="store_true",
                        help="If True, all bonds will be single bonds regarding the actual order.")
    parser.add_argument("--stl_ascii_mode", action="store_true",
                        help="If True, stl model will be saved in ascii, otherwise in binary format.")

    parser.add_argument('files', nargs=argparse.ONE_OR_MORE)

    args, _ = parser.parse_known_args()

    # print(args.__dict__)

    fnames = []
    for fname in args.files:
        fnames += glob(fname)

    assert args.bond_num_segments > 2, "Number of bond segments must be > 2"

    for fpath in fnames:
        if not os.path.isfile(fpath):
            continue

        _dir, fname = os.path.split(fpath)  # get dir and filename
        fname, ext = os.path.splitext(fname)  # get filename without extension

        print(f'Processing \'{fname}{ext}\'...')

        molecule = parse_CML_file(fpath)
        model = create_mol_model(molecule, **args.__dict__)

        model.save(os.path.join(_dir, f'{fname}.stl'), mode=Mode.ASCII if args.stl_ascii_mode else Mode.BINARY)

        print(f'\'{fname}{ext}\' finished.\n')





