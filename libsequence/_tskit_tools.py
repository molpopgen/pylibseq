import numpy as np
from ._libsequence import VariantMatrix
from ._libsequence import AlleleCountMatrix


def _process_temp_containers(temp, pos, ac, max_allele_value):
    vm = VariantMatrix(np.array(temp), np.array(pos), max_allele_value)
    vmac = AlleleCountMatrix(vm)
    ac.extend(np.array(vmac, copy=False).flatten().tolist())
    return ac


def AlleleCountMatrix_from_tree_sequence(ts, chunksize=50, max_allele_value=1):
    ac = []
    temp = []
    pos = []
    for v in ts.variants():
        if len(temp) < chunksize:
            temp.append(np.array(v.genotypes, copy=True, dtype=np.int8))
            pos.append(v.position)
        else:
            ac = _process_temp_containers(
                temp, pos, ac,  max_allele_value)
            temp = [np.array(v.genotypes, copy=True, dtype=np.int8)]
            pos = [v.position]

    if len(temp) > 0:
        ac = _process_temp_containers(
            temp, pos, ac,  max_allele_value)
    return AlleleCountMatrix(ac, max_allele_value+1,
                             ts.num_mutations, ts.num_samples)
