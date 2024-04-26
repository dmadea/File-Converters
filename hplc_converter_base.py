
import numpy as np
import struct


def read_utf16(f):
    """modified string read method which works with Agilent files"""
    # determine length to read
    read_len, = struct.unpack('>B', f.read(1))
    # read values, decode, and strip
    data = f.read(2 * read_len)
    text = data.decode('utf16').strip()

    return text


def save_mat2csv(fname, matrix, times=None, wls=None, unit='', ):
    delimiter = ','
    times = np.arange(0, matrix.shape[0]) if times is None else times
    wls = np.arange(0, matrix.shape[1]) if wls is None else wls

    mat = np.hstack((times[:, None], matrix))
    buffer = f'unit: {unit} - Time | Wavelength->'
    buffer += delimiter + delimiter.join(f"{num:.6g}" for num in wls) + '\n'
    buffer += '\n'.join(delimiter.join(f"{num:.6g}" for num in row) for row in mat)

    with open(fname, 'w', encoding='utf8') as f:
        f.write(buffer)