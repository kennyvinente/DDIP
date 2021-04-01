
from scipy.sparse import hstack, csr_matrix as sparse
from numpy import array, any, delete, unique, arange, nonzero, pi, \
    r_, ones, Inf, zeros
from numpy import flatnonzero as find


from idx_all import BUS_I, BUS_TYPE, REF, VA, VM, PD, GS, VMAX, VMIN
from idx_all import GEN_BUS, VG, PG, QG, PMAX, PMIN, QMAX, QMIN
from idx_all import RATE_A, F_BUS, T_BUS, BR_X, TAP, SHIFT, BR_STATUS, ANGMIN, ANGMAX
from idx_all import H_BUS


def makeBdc(baseMVA, bus, branch):
    """Builds the B matrices and phase shift injections for DC power flow.

    Returns the B matrices and phase shift injection vectors needed for a
    DC power flow.
    The bus real power injections are related to bus voltage angles by::
        P = Bbus * Va + PBusinj
    The real power flows at the from end the lines are related to the bus
    voltage angles by::
        Pf = Bf * Va + Pfinj
    Does appropriate conversions to p.u.

    @see: L{dcpf}

    @author: Carlos E. Murillo-Sanchez (PSERC Cornell & Universidad
    Autonoma de Manizales)
    @author: Ray Zimmerman (PSERC Cornell)
    """
    ## constants
    nb = bus.shape[0]          ## number of buses
    nl = branch.shape[0]       ## number of lines

    ## check that bus numbers are equal to indices to bus (one set of bus nums)
    if any(bus[:, BUS_I] != list(range(nb))):
        print('makeBdc: buses must be numbered consecutively in bus matrix\n')

    ## for each branch, compute the elements of the branch B matrix and the phase
    ## shift "quiescent" injections, where
    ##
    ##      | Pf |   | Bff  Bft |   | Vaf |   | Pfinj |
    ##      |    | = |          | * |     | + |       |
    ##      | Pt |   | Btf  Btt |   | Vat |   | Ptinj |
    ##
    stat = branch[:, BR_STATUS]               ## ones at in-service branches
    b = stat / branch[:, BR_X]                ## series susceptance
    tap = ones(nl)                            ## default tap ratio = 1
    i = find(branch[:, TAP])               ## indices of non-zero tap ratios
    tap[i] = branch[i, TAP]                   ## assign non-zero tap ratios
    b = b / tap

    ## build connection matrix Cft = Cf - Ct for line and from - to buses
    f = branch[:, F_BUS]                           ## list of "from" buses
    t = branch[:, T_BUS]                           ## list of "to" buses
    i = r_[range(nl), range(nl)]                   ## double set of row indices
    ## connection matrix
    Cft = sparse((r_[ones(nl), -ones(nl)], (i, r_[f, t])), (nl, nb))

    ## build Bf such that Bf * Va is the vector of real branch powers injected
    ## at each branch's "from" bus
    Bf = sparse((r_[b, -b], (i, r_[f, t])), shape = (nl, nb))## = spdiags(b, 0, nl, nl) * Cft

    ## build Bbus
    Bbus = Cft.T * Bf

    ## build phase shift injection vectors
    Pfinj = b * (-branch[:, SHIFT] * pi / 180)  ## injected at the from bus ...
    # Ptinj = -Pfinj                            ## and extracted at the to bus
    Pbusinj = Cft.T * Pfinj                ## Pbusinj = Cf * Pfinj + Ct * Ptinj

    return Bbus, Bf, Pbusinj, Pfinj

def makeAang(baseMVA, branch, nb):
    """Construct constraints for branch angle difference limits.

    Constructs the parameters for the following linear constraint limiting
    the voltage angle differences across branches, where C{Va} is the vector
    of bus voltage angles. C{nb} is the number of buses::

        lang <= Aang * Va <= uang

    C{iang} is the vector of indices of branches with angle difference limits.

    @author: Ray Zimmerman (PSERC Cornell)
    @author: Carlos E. Murillo-Sanchez (PSERC Cornell & Universidad
    Autonoma de Manizales)
    """

    iang = find(((branch[:, ANGMIN] != 0) & (branch[:, ANGMIN] > -360)) |
                ((branch[:, ANGMAX] != 0) & (branch[:, ANGMAX] <  360)))
    iangl = find(branch[iang, ANGMIN])
    iangh = find(branch[iang, ANGMAX])
    nang = len(iang)

    if nang > 0:
        ii = r_[arange(nang), arange(nang)]
        jj = r_[branch[iang, F_BUS], branch[iang, T_BUS]]
        Aang = sparse((r_[ones(nang), -ones(nang)],
                       (ii, jj)), (nang, nb))
        uang = Inf * ones(nang)
        lang = -uang
        lang[iangl] = branch[iang[iangl], ANGMIN] * pi / 180
        uang[iangh] = branch[iang[iangh], ANGMAX] * pi / 180
    else:
        Aang  = zeros((0, nb))
        lang  = array([])
        uang  = array([])

    return Aang, lang, uang, iang


def dcopf(ppc):

    ## data dimensions
    nb = ppc['bus'].shape[0]    ## number of buses
    nl = ppc['branch'].shape[0] ## number of branches
    ng = ppc['gen'].shape[0]    ## number of term gens
    nh = ppc['hidrogen'].shape[0] ## number of hydro gens

    ## create (read-only) copies of individual fields for convenience
    baseMVA, bus, gen, hidrogen, branch = ppc['baseMVA'], ppc['bus'], ppc['gen'], ppc['hidrogen'], ppc['branch']

    ## warn if there is more than one reference bus
    refs = find(bus[:, BUS_TYPE] == REF)

    ## set up initial variables and bounds
    gbus = gen[:, GEN_BUS].astype(int)
    hbus = hidrogen[:, H_BUS].astype(int)
    Va   = bus[:, VA] * (pi / 180.0)
    Vm   = bus[:, VM].copy()

    ## power mismatch constraints
    B, Bf, _, _ = makeBdc(baseMVA, bus, branch)
    neg_Cg = sparse((-ones(ng), (gen[:, GEN_BUS], arange(ng))), (nb, ng))   ## Pbus w.r.t. Pg
    neg_Ch = sparse((-ones(nh), (hidrogen[:, H_BUS], arange(nh))), (nb, nh))   ## Pbus w.r.t. Ph


    ## branch flow constraints
    il = find((branch[:, RATE_A] != 0) & (branch[:, RATE_A] < 1e10))
    nl2 = len(il)         ## number of constrained lines
    lpf = -Inf * ones(nl2)
    # upf = branch[il, RATE_A] / baseMVA - Pfinj[il]
    # upt = branch[il, RATE_A] / baseMVA + Pfinj[il]
    upf = branch[il, RATE_A] / baseMVA
    upt = branch[il, RATE_A] / baseMVA

    ## branch voltage angle difference limits
    Aang, lang, uang, iang  = makeAang(baseMVA, branch, nb)


    return B, Bf, neg_Cg, neg_Ch, Aang, lang, uang, iang
