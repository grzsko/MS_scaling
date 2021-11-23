from optimal_transport import calc_distance, calc_moves

def distance(ints1, ints2, mzs1, mzs2, lam=1, eps=0.1, tol=1e-05, threshold=100,
                      max_iter=500, dist_upper_bound=10):

    res = calc_distance(
        mzs1, ints1, mzs2, ints2,
        lam, eps, tol, threshold, max_iter, dist_upper_bound)

    return res

def distance_moves(ints1, ints2, mzs1, mzs2, lam=1, eps=0.1,
                            tol=1e-05, threshold=100, max_iter=500,
                            dist_upper_bound=10):

    res = calc_moves(
        mzs1, ints1, mzs2, ints2,
        lam, eps, tol, threshold, max_iter, dist_upper_bound)

    for id_left, id_right, transport in res:
        yield mzs1[id_left], mzs2[id_right], transport


from optimal_transport_dense import calc_moves_dense, calc_distance_dense

def distance_dense(ints1, ints2, mzs1=None, mzs2=None, dists=None, lam=1,
                   eps=0.1, tol=1e-05, threshold=1e+02, max_iter=500,
                   method="TV"):

    if dists is not None:
        return calc_distance_dense(ints1, ints2, dists, lam, eps, tol, threshold,
                                   max_iter, method)
    elif mzs1 is not None and mzs2 is not None:
        return calc_distance_dense(mzs1, ints1, mzs2, ints2, lam, eps, tol,
                                   threshold, max_iter, method)

    raise ValueError("Either mzs1 and mzs2 or dists must be not None.")

def distance_moves_dense(ints1, ints2, mzs1=None, mzs2=None, dists=None, lam=1,
                         eps=0.1, tol=1e-05, threshold=100, max_iter=500,
                         method="TV"):

    if dists is not None:
        return calc_moves_dense(ints1, ints2, dists, lam, eps, tol, threshold,
                               max_iter, method)
    elif mzs1 is not None and mzs2 is not None:
        return calc_moves_dense(mzs1, ints1, mzs2, ints2, lam, eps, tol,
                                threshold, max_iter, method)
