'''Repurposing example TDVP code provided by https://tenpy.readthedocs.io/en/latest/examples/advanced/mpo_exponential_decay.html for including sum on exponentially decaying MPOs to approximate long-range interactions.
'''

from tenpy.models.lattice import Chain
from tenpy.networks.mpo import MPO
from tenpy.networks.site import SpinHalfSite
from tenpy.models.model import MPOModel

class LongRangeHeisenberg(MPOModel):
    r"""Spin-1/2 Heisenberg Chain with long range decaying interactions approximated by sum of exponentially decaying MPOs.

    The Hamiltonian reads:

    .. math ::

        H = \sum_i \sum_{j>i} \sum_k \exp(-\frac{|j-i-1|}{\mathtt{xi[k]}}) (
                  \mathtt{Jxx[k]}/2 (S^{+}_i S^{-}_j + S^{-}_i S^{+}_j)
                + \mathtt{Jz[k]} S^z_i S^z_j                        ) \\
            - \sum_i \mathtt{hz} S^z_i

    All parameters are collected in a single dictionary `model_param` and read out with
    :func:`~tenpy.tools.params.get_parameter`.

    Parameters
    ----------
    L : int
        Length of the chain.
    Jxx, Jz, hz, xi: float
        Coupling parameters as defined for the Hamiltonian above and indexed by each exponentially decaying MPO.
    bc_MPS : {'finite' | 'infinte'}
        MPS boundary conditions.
    conserve : 'Sz' | 'parity' | None
        What should be conserved. See :class:`~tenpy.networks.Site.SpinHalfSite`.
    """
    def __init__(self, model_param):
        # model parameters
        L = get_parameter(model_param, 'L', 2, self.__class__)
        xi = get_parameter(model_param, 'xi', 0.5, self.__class__)
        Jxx = get_parameter(model_param, 'Jxx', 1., self.__class__)
        Jz = get_parameter(model_param, 'Jz', 1.5, self.__class__)
        hz = get_parameter(model_param, 'hz', 0., self.__class__)
        conserve = get_parameter(model_param, 'conserve', 'Sz', self.__class__)

        site = SpinHalfSite(conserve=conserve)
        lat = Chain(L, site, bc_MPS='finite', bc='open')
        Sz, Sp, Sm, Id = site.Sz, site.Sp, site.Sm, site.Id

        #Generate exponentially decaying MPO from W matrices and MPO.from_grids function. Generate MPO for each exponential decay parameter in array xi.
        MPO_List = [];
        for i in range(len(xi)):
            if xi[i] > 1.0:
                ValueError("Decay needs to be less than 1.0")

            grid = [[Id,   Sp,   Sm,   Sz,   -hz*Sz    ],
                    [None, xi[i]*Id, None, None, 0.5*Jxx[i]*Sm],
                    [None, None, xi[i]*Id, None, 0.5*Jxx[i]*Sp],
                    [None, None, None, xi[i]*Id, Jz[i]*Sz     ],
                    [None, None, None, None, Id        ]]
            grids = [grid] * L

            MPO_List.append(MPO.from_grids(lat.mps_sites(), grids, bc='finite', IdL=0, IdR=-1))

        H = MPO_List[0];

        #Perform naive sum over len(xi) exponentially decaying MPOs.
        for i in range(1, len(xi)):
            H = H.__add__(MPO_List[i]);
        MPOModel.__init__(self, lat, H)
