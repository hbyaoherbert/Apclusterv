import numbers
import warnings

import numpy as np

from sklearn.exceptions import ConvergenceWarning
from sklearn.base import BaseEstimator, ClusterMixin
from sklearn.utils import as_float_array, check_random_state
from sklearn.utils import check_scalar
from sklearn.utils.validation import check_is_fitted
from sklearn.metrics import euclidean_distances
from sklearn.metrics import pairwise_distances_argmin
def _equal_similarities_and_preferences(S, preference):
    def all_equal_preferences():
        return np.all(preference == preference.flat[0])

    def all_equal_similarities():
        # Create mask to ignore diagonal of S
        mask = np.ones(S.shape, dtype=bool)
        np.fill_diagonal(mask, 0)

        return np.all(S[mask].flat == S[mask].flat[0])

    return all_equal_preferences() and all_equal_similarities()
def adj_affinity_propagation(
    S,
    *,
    preference=None,
    convergence_iter=15,
    max_d_iter=10,
    max_p_iter = 10,
    pstep = 10,
    max_iter=200,
    damping=0.5,
    copy=True,
    verbose=False,
    return_n_iter=False,
    random_state=None,
):
    """Perform Affinity Propagation Clustering of data.

    Read more in the :ref:`User Guide <affinity_propagation>`.

    Parameters
    ----------

    S : array-like of shape (n_samples, n_samples)
        Matrix of similarities between points.

    preference : array-like of shape (n_samples,) or float, default=None
        Preferences for each point - points with larger values of
        preferences are more likely to be chosen as exemplars. The number of
        exemplars, i.e. of clusters, is influenced by the input preferences
        value. If the preferences are not passed as arguments, they will be
        set to the median of the input similarities (resulting in a moderate
        number of clusters). For a smaller amount of clusters, this can be set
        to the minimum value of the similarities.

    convergence_iter : int, default=15
        Number of iterations with no change in the number
        of estimated clusters that stops the convergence.

    max_iter : int, default=200
        Maximum number of iterations.
    max_d_iter: 
        Maximum super iterations for damping factor steps
    damping : float, default=0.5
        Damping factor between 0.5 and 1.

    copy : bool, default=True
        If copy is False, the affinity matrix is modified inplace by the
        algorithm, for memory efficiency.

    verbose : bool, default=False
        The verbosity level.

    return_n_iter : bool, default=False
        Whether or not to return the number of iterations.

    random_state : int, RandomState instance or None, default=None
        Pseudo-random number generator to control the starting state.
        Use an int for reproducible results across function calls.
        See the :term:`Glossary <random_state>`.

        .. versionadded:: 0.23
            this parameter was previously hardcoded as 0.

    Returns
    -------

    cluster_centers_indices : ndarray of shape (n_clusters,)
        Index of clusters centers.

    labels : ndarray of shape (n_samples,)
        Cluster labels for each point.

    n_iter : int
        Number of iterations run. Returned only if `return_n_iter` is
        set to True.

    Notes
    -----
    For an example, see :ref:`examples/cluster/plot_affinity_propagation.py
    <sphx_glr_auto_examples_cluster_plot_affinity_propagation.py>`.

    When the algorithm does not converge, it will still return a arrays of
    ``cluster_center_indices`` and labels if there are any exemplars/clusters,
    however they may be degenerate and should be used with caution.

    When all training samples have equal similarities and equal preferences,
    the assignment of cluster centers and labels depends on the preference.
    If the preference is smaller than the similarities, a single cluster center
    and label ``0`` for every sample will be returned. Otherwise, every
    training sample becomes its own cluster center and is assigned a unique
    label.

    References
    ----------
    Brendan J. Frey and Delbert Dueck, "Clustering by Passing Messages
    Between Data Points", Science Feb. 2007
    """
    S = as_float_array(S, copy=copy)
    n_samples = S.shape[0]
    
    if S.shape[0] != S.shape[1]:
        raise ValueError("S must be a square array (shape=%s)" % repr(S.shape))

    if preference is None:
        preference = np.median(S)
    
    preference = np.array(preference)

    if n_samples == 1 or _equal_similarities_and_preferences(S, preference):
        # It makes no sense to run the algorithm in this case, so return 1 or
        # n_samples clusters, depending on preferences
        warnings.warn(
            "All samples have mutually equal similarities. "
            "Returning arbitrary cluster center(s)."
        )

        if preference.flat[0] >= S.flat[n_samples - 1]:
            return (
                (np.arange(n_samples), np.arange(n_samples), 0)
                if return_n_iter
                else (np.arange(n_samples), np.arange(n_samples))
            )
        else:
            return (
                (np.array([0]), np.array([0] * n_samples), 0)
                if return_n_iter
                else (np.array([0]), np.array([0] * n_samples))
            )

    random_state = check_random_state(random_state)

    # Place preference on the diagonal of S
  
    #print(S[0,1],S[1,15])
    step_d_iter = max_d_iter
    num_clusters = 0
    unconvergable = False
    optimal = False
    for p_iter in range(max_p_iter):
        print("preference step p_iter %d" % p_iter)
        S.flat[:: (n_samples + 1)] = preference

        A = np.zeros((n_samples, n_samples))
        R = np.zeros((n_samples, n_samples))  # Initialize messages
        # Intermediate results
        tmp = np.zeros((n_samples, n_samples))

        # Remove degeneracies
        S += (
            np.finfo(S.dtype).eps * S + np.finfo(S.dtype).tiny * 100
        ) * random_state.standard_normal(size=(n_samples, n_samples))

        # Execute parallel affinity propagation updates
        e = np.zeros((n_samples, convergence_iter))

        ind = np.arange(n_samples)
        if unconvergable:
            break
        if optimal:
            break
        for d_iter in range(step_d_iter):
        
            if damping > 1:
                unconvergable = True
                break
            for it in range(max_iter):
                # tmp = A + S; compute responsibilities
                np.add(A, S, tmp)
                I = np.argmax(tmp, axis=1)
                Y = tmp[ind, I]  # np.max(A + S, axis=1)
                tmp[ind, I] = -np.inf
                Y2 = np.max(tmp, axis=1)

                # tmp = Rnew
                np.subtract(S, Y[:, None], tmp)
                tmp[ind, I] = S[ind, I] - Y2
                #print(tmp[0,1],tmp[1,15])
                # Damping
                tmp *= 1 - damping
                R *= damping
                R += tmp

                # tmp = Rp; compute availabilities
                np.maximum(R, 0, tmp)
                tmp.flat[:: n_samples + 1] = R.flat[:: n_samples + 1]

                # tmp = -Anew
                tmp -= np.sum(tmp, axis=0)
                dA = np.diag(tmp).copy()
                tmp.clip(0, np.inf, tmp)
                tmp.flat[:: n_samples + 1] = dA
                #print(tmp[0,1],tmp[1,15])
                # Damping
                tmp *= 1 - damping
                A *= damping
                A -= tmp

            # Check for convergence
                E = (np.diag(A) + np.diag(R)) > 0
                e[:, it % convergence_iter] = E
                K = np.sum(E, axis=0)

                if it >= convergence_iter:
                    se = np.sum(e, axis=1)
                    unconverged = np.sum((se == convergence_iter) + (se == 0)) != n_samples
                    if (not unconverged and (K > 0)) or (it == max_iter):
                        never_converged = False
                        if verbose:
                            print("Converged after %d iterations." % damping)
                        break
            else:
                never_converged = True
                if verbose:
                    print("Did not converge")

            I = np.flatnonzero(E)
            K = I.size  # Identify exemplars

            if K > 0:
                if never_converged:
                    warnings.warn(
                        "Affinity propagation d_iter did not converge, this model "
                        "may return degenerate cluster centers and labels.",
                        ConvergenceWarning,
                    )
                    print("daming %f not converged" % damping)
                c = np.argmax(S[:, I], axis=1)
                c[I] = np.arange(K)  # Identify clusters
                # Refine the final set of exemplars and clusters and return results
                for k in range(K):
                    ii = np.where(c == k)[0]
                    j = np.argmax(np.sum(S[ii[:, np.newaxis], ii], axis=0))
                    I[k] = ii[j]

                c = np.argmax(S[:, I], axis=1)
                c[I] = np.arange(K)
                labels = I[c]
                # Reduce labels to a sorted, gapless, list
                cluster_centers_indices = np.unique(labels)
                labels = np.searchsorted(cluster_centers_indices, labels)

                if not never_converged:
                    step_d_iter = 1
                    #print(preference,len(cluster_centers_indices))
                    preference += pstep
                    
                    if num_clusters >= len(cluster_centers_indices):
                        optimal = True
                    num_clusters = len(cluster_centers_indices)
                    break
                damping += 0.1
                if d_iter == max_d_iter - 1:
                    unconvergable = True
            else:
                warnings.warn(
                    "Affinity propagation did not converge and this model "
                    "will not have any cluster centers.",
                    ConvergenceWarning,
                )
                labels = np.array([-1] * n_samples)
                cluster_centers_indices = []
                if d_iter == max_d_iter - 1:
                    unconvergable = True
                damping += 0.1

    if return_n_iter:
        return cluster_centers_indices, labels, it + 1
    else:
        return cluster_centers_indices, labels