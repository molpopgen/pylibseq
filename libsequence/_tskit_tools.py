import numpy as np
from ._libsequence import VariantMatrix
from ._libsequence import AlleleCountMatrix


def _process_temp_containers(temp, pos, ac):
    vm = VariantMatrix(np.array(temp), np.array(pos))
    vmac = AlleleCountMatrix(vm)
    if ac is None:
        ac = vmac
    else:
        ac = ac._merge(vmac)
    return ac


def AlleleCountMatrix_from_tree_sequence(ts, chunksize=50):
    ac = None
    temp = []
    pos = []
    for v in ts.variants():
        if len(temp) < chunksize:
            temp.append(np.array(v.genotypes, copy=True, dtype=np.int8))
            pos.append(v.position)
        else:
            ac = _process_temp_containers(temp, pos, ac)
            temp = [np.array(v.genotypes, copy=True, dtype=np.int8)]
            pos = [v.position]

    if len(temp) > 0:
        ac = _process_temp_containers(temp, pos, ac)

    return ac
