

import numpy as np
import math
from scipy.spatial import Delaunay

from stl import mesh
import stl
from numpy.linalg import norm

from cml_reader import parse_CML_file, Molecule


def form_mesh(vertices: np.ndarray, faces: np.ndarray):
    obj = mesh.Mesh(np.zeros(faces.shape[0], dtype=mesh.Mesh.dtype))
    for i, f in enumerate(faces):
        obj.vectors[i][0] = vertices[f[0]]
        obj.vectors[i][1] = vertices[f[1]]
        obj.vectors[i][2] = vertices[f[2]]

    return obj

# https://observablehq.com/@mourner/fast-icosphere-mesh
def subdivide(vertices: np.ndarray, faces: np.ndarray, order=1):
    if order == 0:
        return vertices, faces

    if order > 1:
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
        ])

    vertices, faces = subdivide(vertices, faces, refinement_order)

    # normalize vertices so the distance of each is the same
    vert_length = np.sqrt((vertices * vertices).sum(axis=1, keepdims=True))
    vertices = vertices / vert_length * radius + center

    return form_mesh(vertices, faces)


# from https://github.com/PyMesh/PyMesh/blob/main/python/pymesh/meshutils/generate_cylinder.py
def gen_cylinder(p0, p1, radius, num_segments=32):

    Z = np.asarray([0, 0, 1])

    p0 = np.array(p0, dtype=np.float64)
    p1 = np.array(p1, dtype=np.float64)

    angles = [2 * np.pi * i / float(num_segments) for i in range(num_segments)]
    rim = np.array([[math.cos(theta), math.sin(theta), 0.0] for theta in angles])

    axis = p1 - p0
    axis /= norm(axis)

    rot_axis = np.cross(Z, axis)  # rotational axis

    rot_matrix = stl.Mesh.rotation_matrix(rot_axis, np.arccos(Z.dot(axis)))

    bottom_rim = rot_matrix.dot(rim.T).T * radius + p0
    top_rim = rot_matrix.dot(rim.T).T * radius + p1

    vertices = np.vstack([bottom_rim, top_rim])

    faces = np.array([
        [[i, (i + 1) % num_segments, i + num_segments],
         [i + num_segments, (i + 1) % num_segments, (i + 1) % num_segments + num_segments]]
        for i in range(num_segments)], dtype=int)
    faces = faces.reshape((-1, 3), order="C")

    return form_mesh(vertices, faces)


def gen_cube(center=(0, 0, 0)):
    vertices = np.array([
        [-1, -1, -1],
        [+1, -1, -1],
        [+1, +1, -1],
        [-1, +1, -1],
        [-1, -1, +1],
        [+1, -1, +1],
        [+1, +1, +1],
        [-1, +1, +1]], dtype=np.float64) + np.asarray(center)

    # Define the 12 triangles composing the cube
    faces = np.array([
        [0, 3, 1],
        [1, 3, 2],
        [0, 4, 7],
        [0, 7, 3],
        [4, 5, 6],
        [4, 6, 7],
        [5, 1, 2],
        [5, 2, 6],
        [2, 3, 6],
        [3, 7, 6],
        [0, 1, 5],
        [0, 5, 4]], dtype=int)

    return form_mesh(vertices, faces)


def combine(objs: list):
    return mesh.Mesh(np.concatenate([obj.data for obj in objs]))


def random_point(limit=50):
    return (np.random.random(3) - 0.5) * 2 * limit



def create_mol_model(mol: Molecule, atom_size=0.3, bond_radius=0.15,
                     bond_num_segments=32, sphere_refinment_order=3,
                     create_bond_orders=False, bond_distance=0.13):

    objects = []
    for atom in mol.atoms:
        atom_obj = generate_icosphere(atom.size() * atom_size, atom.position, sphere_refinment_order)
        objects.append(atom_obj)

    for bond in mol.bonds:
        a1pos = mol.atoms[bond.atom1].position
        a2pos = mol.atoms[bond.atom2].position

        Z_ax = np.asarray([0, 0, 1])
        f = bond_distance

        b_direct = gen_cylinder(a1pos, a2pos, bond_radius, bond_num_segments)
        if create_bond_orders:
            bond_ax = a2pos - a1pos
            perp_ax = np.cross(bond_ax, Z_ax)
            perp_ax /= norm(perp_ax)  # perpendicular axis to bond axis
            b_up = gen_cylinder(a1pos + perp_ax * f, a2pos + perp_ax * f, bond_radius, bond_num_segments)
            b_down = gen_cylinder(a1pos - perp_ax * f, a2pos - perp_ax * f, bond_radius, bond_num_segments)

            if bond.order == 1:
                objects.append(b_direct)
            elif bond.order == 2:
                objects += [b_up, b_down]
            elif bond.order == 3:
                objects += [b_direct, b_up, b_down]
        else:
            objects.append(b_direct)

    return combine(objects)


if __name__ == "__main__":

    file = 'NHME2-PDP.cml'
    # file = 'bilirubin.cml'

    molecule = parse_CML_file(file)

    model = create_mol_model(molecule)

    model.save(f'{file}.stl', mode=stl.Mode.BINARY)

    #
    # objects = []
    #
    # for i in range(100):
    #     cube = gen_cube(random_point())
    #     cube.rotate(random_point(), theta=np.random.random(1)[0] * 2 * np.pi)
    #     objects.append(cube)
    #
    # for i in range(30):
    #     cyl = gen_cylinder(random_point(), random_point(), np.random.random(1)[0] * 2 + 0.1)
    #     objects.append(cyl)
    #
    # for i in range(50):
    #     rad = np.random.random(1)[0] * 4 + 0.2
    #     sphere = generate_icosphere(rad, random_point())
    #     objects.append(sphere)
    #
    # combined = combine(objects)
    # #
    # # # Write the mesh to file "cube.stl"
    # combined.save('objects.stl', mode=stl.Mode.BINARY)
    #
