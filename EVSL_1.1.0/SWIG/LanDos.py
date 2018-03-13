import numpy as np
import evsl
from scipy.sparse import csgraph
from scipy.sparse import csr_matrix
from scipy.sparse import coo_matrix
from scipy.sparse.sputils import get_index_dtype


class intv:
    lmin = 5
    lmax = 10
    a = 4
    b = 11


if __name__ == '__main__':
    msteps = 40
    npts = 200
    nvec = 100
    num = 20

    G = np.arange(num) * np.arange(num)[:, np.newaxis]
    csgraph.laplacian(G, normed=False)
    coomat = coo_matrix(G)
    M, N = coomat.shape
    idx_dtype = get_index_dtype((coomat.col, coomat.row), maxval=max(coomat.nnz,
                                                                    M))
    col = coomat.row.astype(idx_dtype, copy=False)
    row = coomat.row.astype(idx_dtype, copy=False)

    ecoo = evsl.cooMat()
    ecoo.nrows = N
    ecoo.ncols = M
    print(type(coomat.row), "\n")
    print(type(coomat.data), "\n")
    print(coomat.data, "\n")
    ecoo.vv = coomat.data
    ecoo.ir = coomat.row
    ecoo.jc = coomat.col
    ecsr = evsl.csrMat()

    csrmat = evsl.cooMat_to_csr(0, ecoo, ecsr)

    evsl.EVSLStart()
    evsl.SetStdEig()
    evsl.SetAMatrix(csrmat)

