import os
from cffi import FFI


from MasSpOT.optimal_transport import calc_distance, calc_steps

def spectral_distance(spec1, spec2, lam=1, eps=0.1, tol=1e-05, threshold=100,
                      max_iter=500, dist_upper_bound=100):

    mzs1, ints1 = zip(*[x for x in spec1.confs if x[1] > 0])
    mzs2, ints2 = zip(*[x for x in spec2.confs if x[1] > 0])

    res = calc_distance(
        mzs1, ints1, mzs2, ints2,
        lam, eps, tol, threshold, max_iter, dist_upper_bound)

    return res

def spectral_distance_moves(spec1, spec2, lam=1, eps=0.1, tol=1e-05, threshold=100,
                            max_iter=500, dist_upper_bound=100):

    mzs1, ints1 = zip(*[x for x in spec1.confs if x[1] > 0])
    mzs2, ints2 = zip(*[x for x in spec2.confs if x[1] > 0])

    res = calc_steps(
        mzs1, ints1, mzs2, ints2,
        lam, eps, tol, threshold, max_iter, dist_upper_bound)

    for id_left, id_right, transport in res:
        yield mzs1[id_left], mzs2[id_right], transport


from MasSpOT.optimal_transport_dense import calc_distance_dense, calc_steps_dense

# For test and legacy reasons
def spectral_distance_dense(m1, m2, dists=None, lam=1, eps=0.1, tol=1e-05,
                            threshold=1e+02, max_iter=500, method="TV"):
    if dists is not None:
        res = calc_distance_dense(m1, m2, dists, lam, eps, tol,
                                  threshold, max_iter, method)
    else:
        mzs1, ints1 = zip(*[x for x in m1.confs if x[1] > 0])
        mzs2, ints2 = zip(*[x for x in m2.confs if x[1] > 0])

        res = calc_distance_dense(mzs1, ints1, mzs2, ints2, lam, eps, tol,
                                  threshold, max_iter, method)

    return res

def spectral_distance_moves_dense(spec1, spec2, dists=None, lam=1, eps=0.1, tol=1e-05,
                                  threshold=100, max_iter=500, method="TV"):

    mzs1, ints1 = zip(*[x for x in spec1.confs if x[1] > 0])
    mzs2, ints2 = zip(*[x for x in spec2.confs if x[1] > 0])

    if dists is not None:
        res = calc_steps_dense(ints1, ints2, dists, lam, eps, tol, threshold,
                               max_iter, method)
    else:
        res = calc_steps_dense(mzs1, ints1, mzs2, ints2, lam, eps, tol,
                               threshold, max_iter, method)

    # TODO XXX remove these moves and return pure tp or something simpler
    # (sparse tp)
    # TODO XXX make interface same with distance
    for id_left, id_right, transport in res:
        yield mzs1[id_left], mzs2[id_right], transport
