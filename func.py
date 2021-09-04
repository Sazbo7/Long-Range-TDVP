#Evaluating the polynomial decaying function (1/r^\alpha) as outlined by "Matrix product operator representations"
# B Pirvu et al., New Journal of Physics 12 (2010) 025012 doi:10.1088/1367-2630/12/2/025012
import numpy as np

def poly_decay(alpha, x):
    """
         Usage:
         y = poly_decay(alpha, x)

         :param alpha: Exponent of polynomial decay, alpha > 0.
         :param x: Real float or numpy.array() that is input of polynomial function.

        :return: Real float or numpy.array() output of polynomial function.
    """
    if alpha < 0:
        ValueError("Exponent must be infinite range or ")
    return x**-alpha;

def long_range_coeffs(function, alpha, Nsites, Napprox):
    r"""
    Evaluate k-order (Napprox) least squares fit of true "function" and sum_k exponentially decaying functions.
    $f(r) = 1/r^\alpha ~ \sum_k c_i (\lambda_i)^r$

    Requires QR decomposition of function matrix f(r) and subsequent least squares fit of Napprox eigenvalues.

    $ min(c_i, \lambda_i) \sum_i^{Nsites} |function(r_i) - \sum_k c_i \lambda_i^r_i |$

         Usage:
         kvals, lstsq = long_range_coeffs(poly_decay, alpha, Nsites, Napprox)

         :param function: True 1d function for which we are approximating
         :param alpha: Exponent argument for function
         :param Nsites: Number of sites over which to evaluate the function. Requires larger number of sites for smaller alpha
         :param Napprox: Number of unique exponentially decaying functions used to approximate "function". k = \infty is exact.

         :return eigs: Inputs of the exponential {\lambda_i}.
         :return coeffs: Coefficients for exponential i.
         :return error: np.array([Nsites]), absolute error between function and exponential approximation for each site.
    """

    num_rows = Nsites - Napprox + 1;
    num_cols = Napprox
    func_mat = np.zeros([num_rows, num_cols]);

    for i in range(num_rows):
        func_mat[i] = [function(alpha, j + i) for j in range(1, num_cols+1)];

    Q,R = np.linalg.qr(func_mat);
    Q,R = np.linalg.qr(func_mat);
    Q1 = Q[:num_rows - 1][:];
    Q2 = Q[1:][:];
    Q1_inv = np.linalg.pinv(Q1, rcond=1e-15, hermitian=False);
    L = Q1_inv @ Q2;

    eigs = np.linalg.eigvals(L);
    eigs.sort()
    eigs = eigs[::-1];

    Lmat = np.zeros([Nsites, Napprox]);
    fmat = np.zeros([Nsites]);

    for i in range(0, Nsites):
        Lmat[i] = eigs**(i+1);
        fmat[i] = function(alpha, i+1);
    coeffs = np.linalg.lstsq(Lmat, fmat, rcond=-1);

    sites = np.arange(0, Nsites, 1);
    error = function(alpha, sites) - appx_eval(eigs, coeffs, sites);

    return eigs, coeffs[0], error;

def appx_eval(eigs, coeffs, x):
    """
    Evaluate sum of exponentially decaying functions \sum_{|coeffs|} coeffs_i eigs_i^x.

    :param eigs: Set of exponentially decaying values
    :param coeffs: Coefficient for each exponential function
    :param x: Sites over which to evaluate approximate function decomposition

    :return est: Estimate function provided by exponential decay approximation.
    """

    est = 0;
    for i in range(len(eigs)):
        est += coeffs[i] * eigs[i]**x;
    return est;
