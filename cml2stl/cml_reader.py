import xml.etree.ElementTree as ET
import numpy as np

el_list = ['H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al', 'Si', 'P',
           'S', 'Cl', 'Ar', 'K', 'Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn',
           'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr', 'Rb', 'Sr', 'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh',
           'Pd', 'Ag', 'Cd', 'In', 'Sn', 'Sb', 'Te', 'I', 'Xe', 'Cs', 'Ba', 'La', 'Ce', 'Pr', 'Nd',
           'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb', 'Lu', 'Hf', 'Ta', 'W', 'Re',
           'Os', 'Ir', 'Pt', 'Au', 'Hg', 'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn', 'Fr', 'Ra', 'Ac', 'Th',
           'Pa', 'U', 'Np', 'Pu', 'Am', 'Cm', 'Bk', 'Cf', 'Es', 'Fm', 'Md', 'No', 'Lr', 'Rf', 'Db',
           'Sg', 'Bh', 'Hs', 'Mt', 'Ds', 'Rg', 'Cn', 'Nh', 'Fl', 'Mc', 'Lv', 'Ts', 'Og']

atomic_radia = [53, 120, 145, 105, 85, 70, 65, 60, 50, 160, 180, 150, 125, 110, 100, 100, 100, 71,
                220, 180, 160, 140, 135, 140, 140, 140, 135, 135, 135, 135, 130, 125, 115, 115, 115,
                None, 235, 200, 180, 155, 145, 145, 135, 130, 135, 140, 160, 155, 155, 145, 145, 140,
                140, None, 265, 215, 195, 185, 185, 185, 185, 185, 185, 180, 175, 175, 175, 175, 175,
                175, 175, 155, 145, 135, 135, 130, 135, 135, 135, 150, 190, 180, 160, 190, None, None,
                None, 215, 195, 180, 180, 175, 175, 175, 175, 176, None, None, None, None, None, None,
                None, None, None, None, None, None, None, None, None, None, None, None, None, None, None,
                None]

vdw_radia = [90, 140, 182, 153, 192, 170, 155, 152, 147, 154, 227, 173, 184, 210, 180, 180, 175, 188, 275, 231,
211, None, None, None, None, None, None, 163, 140, 139, 187, 211, 185, 190, 185, 202, 303, 249, None,
None, None, None, None, None, None, 163, 172, 158, 193, 217, 206, 206, 198, 216, 343, 268,
None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None,
None, None, None, 175, 166, 155, 196, 202, 207, 197, 202, 220, 348, 283, None, None, None,
186, None, None, None, None, None, None, None, None, None, None, None, None, None,None, None, None, None,None, None, None,
 None, None, None, None, None, None]

# at_numbers = dict(zip(el_list, np.arange(1, len(el_list) + 1)))
el_radia_dict = dict(zip(el_list, vdw_radia))


class Atom(object):

    def __init__(self, id, element_type, position):
        self.id = id
        self.element_type = element_type
        self.position = position

    def size(self):
        f = 1/100
        if self.element_type in el_radia_dict:
            if el_radia_dict[self.element_type] is None:
                return el_radia_dict['C'] * f

            return el_radia_dict[self.element_type] * f

        raise ValueError('Unknown atom')


class Bond(object):
    def __init__(self, atom1: int, atom2: int, order=1):
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
            atom1 = int(a1[1:]) - 1  # indexes in array
            atom2 = int(a2[1:]) - 1
            new_bond = Bond(atom1, atom2, int(attr['order']))
            bonds.append(new_bond)

    return Molecule(atoms, bonds)


if __name__ == "__main__":

    file = 'NHME2-PDP.cml'

    mol = parse_CML_file(file)

    print(mol.atoms[20].size())





