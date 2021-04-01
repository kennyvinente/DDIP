# Copyright (c) 1996-2015 PSERC. All rights reserved.
# Use of this source code is governed by a BSD-style
# license that can be found in the LICENSE file.

"""Defines constants for named column indices to bus matrix.

Some examples of usage, after defining the constants using the line above,
are::

    Pd = bus[3, PD]     # get the real power demand at bus 4
    bus[:, VMIN] = 0.95 # set the min voltage magnitude to 0.95 at all buses

The index, name and meaning of each column of the bus matrix is given
below:

columns 0-12 must be included in input matrix (in case file)
    0.  C{BUS_I}       bus number (1 to 29997)
    1.  C{BUS_TYPE}    bus type (1 = PQ, 2 = PV, 3 = ref, 4 = isolated)
    2.  C{PD}          real power demand (MW)
    3.  C{QD}          reactive power demand (MVAr)
    4.  C{GS}          shunt conductance (MW at V = 1.0 p.u.)
    5.  C{BS}          shunt susceptance (MVAr at V = 1.0 p.u.)
    6.  C{BUS_AREA}    area number, 1-100
    7.  C{VM}          voltage magnitude (p.u.)
    8.  C{VA}          voltage angle (degrees)
    9.  C{BASE_KV}     base voltage (kV)
    10. C{ZONE}        loss zone (1-999)
    11. C{VMAX}        maximum voltage magnitude (p.u.)
    12. C{VMIN}        minimum voltage magnitude (p.u.)

columns 13-16 are added to matrix after OPF solution
they are typically not present in the input matrix

(assume OPF objective function has units, u)
    13. C{LAM_P}       Lagrange multiplier on real power mismatch (u/MW)
    14. C{LAM_Q}       Lagrange multiplier on reactive power mismatch (u/MVAr)
    15. C{MU_VMAX}     Kuhn-Tucker multiplier on upper voltage limit (u/p.u.)
    16. C{MU_VMIN}     Kuhn-Tucker multiplier on lower voltage limit (u/p.u.)

additional constants, used to assign/compare values in the C{BUS_TYPE} column
    1.  C{PQ}    PQ bus
    2.  C{PV}    PV bus
    3.  C{REF}   reference bus
    4.  C{NONE}  isolated bus

@author: Ray Zimmerman (PSERC Cornell)
@author: Richard Lincoln
"""

# define bus types
PQ      = 1
PV      = 2
REF     = 3
NONE    = 4

# define the indices
BUS_I       = 0    # bus number (1 to 29997)
BUS_TYPE    = 1    # bus type
PD          = 2    # Pd, real power demand (MW)
QD          = 3    # Qd, reactive power demand (MVAr)
GS          = 4    # Gs, shunt conductance (MW at V = 1.0 p.u.)
BS          = 5    # Bs, shunt susceptance (MVAr at V = 1.0 p.u.)
BUS_AREA    = 6    # area number, 1-100
VM          = 7    # Vm, voltage magnitude (p.u.)
VA          = 8    # Va, voltage angle (degrees)
BASE_KV     = 9    # baseKV, base voltage (kV)
ZONE        = 10   # zone, loss zone (1-999)
VMAX        = 11   # maxVm, maximum voltage magnitude (p.u.)
VMIN        = 12   # minVm, minimum voltage magnitude (p.u.)

# included in opf solution, not necessarily in input
# assume objective function has units, u
LAM_P       = 13   # Lagrange multiplier on real power mismatch (u/MW)
LAM_Q       = 14   # Lagrange multiplier on reactive power mismatch (u/MVAr)
MU_VMAX     = 15   # Kuhn-Tucker multiplier on upper voltage limit (u/p.u.)
MU_VMIN     = 16   # Kuhn-Tucker multiplier on lower voltage limit (u/p.u.)

#
# # Copyright (c) 1996-2015 PSERC. All rights reserved.
# # Use of this source code is governed by a BSD-style
# # license that can be found in the LICENSE file.
#
# """ Defines constants for named column indices to gen matrix.
#
# Some examples of usage, after defining the constants using the line above,
# are::
#
#     Pg = gen[3, PG]   # get the real power output of generator 4
#     gen[:, PMIN] = 0  # set to zero the minimum real power limit of all gens
#
# The index, name and meaning of each column of the gen matrix is given
# below:
#
# columns 0-20 must be included in input matrix (in case file)
#     0.  C{GEN_BUS}     bus number
#     1.  C{PG}          real power output (MW)
#     2.  C{QG}          reactive power output (MVAr)
#     3.  C{QMAX}        maximum reactive power output (MVAr)
#     4.  C{QMIN}        minimum reactive power output (MVAr)
#     5.  C{VG}          voltage magnitude setpoint (p.u.)
#     6.  C{MBASE}       total MVA base of machine, defaults to baseMVA
#     7.  C{GEN_STATUS}  1 - in service, 0 - out of service
#     8.  C{PMAX}        maximum real power output (MW)
#     9.  C{PMIN}        minimum real power output (MW)
#     10. C{PC1}         lower real power output of PQ capability curve (MW)
#     11. C{PC2}         upper real power output of PQ capability curve (MW)
#     12. C{QC1MIN}      minimum reactive power output at Pc1 (MVAr)
#     13. C{QC1MAX}      maximum reactive power output at Pc1 (MVAr)
#     14. C{QC2MIN}      minimum reactive power output at Pc2 (MVAr)
#     15. C{QC2MAX}      maximum reactive power output at Pc2 (MVAr)
#     16. C{RAMP_AGC}    ramp rate for load following/AGC (MW/min)
#     17. C{RAMP_10}     ramp rate for 10 minute reserves (MW)
#     18. C{RAMP_30}     ramp rate for 30 minute reserves (MW)
#     19. C{RAMP_Q}      ramp rate for reactive power (2 sec timescale) (MVAr/min)
#     20. C{APF}         area participation factor
#
# columns 21-24 are added to matrix after OPF solution
# they are typically not present in the input matrix
#
# (assume OPF objective function has units, u)
#     21. C{MU_PMAX}     Kuhn-Tucker multiplier on upper Pg limit (u/MW)
#     22. C{MU_PMIN}     Kuhn-Tucker multiplier on lower Pg limit (u/MW)
#     23. C{MU_QMAX}     Kuhn-Tucker multiplier on upper Qg limit (u/MVAr)
#     24. C{MU_QMIN}     Kuhn-Tucker multiplier on lower Qg limit (u/MVAr)
#
# @author: Ray Zimmerman (PSERC Cornell)
# @author: Richard Lincoln
# """

# define the indices
GEN_BUS     = 0    # bus number
PG          = 1    # Pg, real power output (MW)
QG          = 2    # Qg, reactive power output (MVAr)
QMAX        = 3    # Qmax, maximum reactive power output at Pmin (MVAr)
QMIN        = 4    # Qmin, minimum reactive power output at Pmin (MVAr)
VG          = 5    # Vg, voltage magnitude setpoint (p.u.)
MBASE       = 6    # mBase, total MVA base of this machine, defaults to baseMVA
GEN_STATUS  = 7    # status, 1 - machine in service, 0 - machine out of service
PMAX        = 8    # Pmax, maximum real power output (MW)
PMIN        = 9    # Pmin, minimum real power output (MW)
PC1         = 10   # Pc1, lower real power output of PQ capability curve (MW)
PC2         = 11   # Pc2, upper real power output of PQ capability curve (MW)
QC1MIN      = 12   # Qc1min, minimum reactive power output at Pc1 (MVAr)
QC1MAX      = 13   # Qc1max, maximum reactive power output at Pc1 (MVAr)
QC2MIN      = 14   # Qc2min, minimum reactive power output at Pc2 (MVAr)
QC2MAX      = 15   # Qc2max, maximum reactive power output at Pc2 (MVAr)
RAMP_AGC    = 16   # ramp rate for load following/AGC (MW/min)
RAMP_10     = 17   # ramp rate for 10 minute reserves (MW)
RAMP_30     = 18   # ramp rate for 30 minute reserves (MW)
RAMP_Q      = 19   # ramp rate for reactive power (2 sec timescale) (MVAr/min)
APF         = 20   # area participation factor

# included in opf solution, not necessarily in input
# assume objective function has units, u
MU_PMAX     = 21   # Kuhn-Tucker multiplier on upper Pg limit (u/MW)
MU_PMIN     = 22   # Kuhn-Tucker multiplier on lower Pg limit (u/MW)
MU_QMAX     = 23   # Kuhn-Tucker multiplier on upper Qg limit (u/MVAr)
MU_QMIN     = 24   # Kuhn-Tucker multiplier on lower Qg limit (u/MVAr)


##### Hidrogen data

H_BUS        = 0     # bus number
H_UPRIVER    = 1     # hydro located at upriver (only 1 hydro allowed atm)
H_DOWNRIVER  = 2     # hydro located at downriver (only 1 hydro allowed atm)
H_TRAVELTIME = 3     # traveltime to downriver (hours)
H_PMAX       = 4     # Pmax, maximum real power output (MW)
H_PMIN       = 5     # Pmin, minimum real power output (MW)
H_QMAX       = 6     # Qmax, maximum reactive power output at Pmin (MVAr)
H_QMIN       = 7     # Qmin, minimium reactive power output at Pmin (MVAr)
H_STATUS     = 8     # status, 1 - machine in service, 0 - machine out of service
H_VOLMIN     = 9     # Volmin, minimum volume allowed in reservoir (hm3)
H_VOLMAX     = 10    # Volmax, maximum volume allowed in reservoir (hm3)
H_SMIN       = 11    # Smin, minimum spillage (m3/s)
H_SMAX       = 12    # Smax, maximum spillage (m3/s)
H_qMIN       = 13    # qmin, minimum discharge (m3/s)
H_qMAX       = 14    # qmax, maximum discharge (m3/s)
H_V0         = 15    # initial volume (hm3)
H_Y0         = 16    # initial affluence (m3/s)
H_Q0         = 17    # initial discarge (m3/s)
H_S0         = 18    # initial spillage (m3/s)
H_HBMIN      = 19    # minimium Hb (m)
H_HBMAX      = 20    # maximum Hb (m)
H_NMAQ       = 21    # number of turbines
H_TYPE       = 22    # type of hydro: 0 - run of the river, 1 - reservoir
H_PROD         = 23  # productibility
H_MVAMAX     = 24    # Max power (MVA)
H_PG         = 25    # Current Pgen
H_QG         = 26    # Current Qgen

##### UC data
UC_TON      = 0     # number of hours that unit is on
UC_UPTIME   = 1     # minimum uptime (h)
UC_DOWNTIME = 2     # minimum downtime (h)
UC_RUP      = 3     # ramp up (MW)
UC_RDOWN    = 4     # ramp down (MW)
UC_P0       = 5     # initial power (P0 >= Pmin if on)


#
# """ Defines constants for named column indices to gencost matrix.
#
# Some examples of usage, after defining the constants using the line above,
# are::
#
#     start = gencost[3, STARTUP]       # get startup cost of generator 4
#     gencost[2, [MODEL, NCOST:COST+2]] = [POLYNOMIAL, 2, 30, 0]
#     # set the cost of generator 2 to a linear function COST = 30 * Pg
#
# The index, name and meaning of each column of the gencost matrix is given
# below:
#
# columns 1-5
#     1.  C{MODEL}       cost model, 1 - piecewise linear, 2 - polynomial
#     2.  C{STARTUP}     startup cost in US dollars
#     3.  C{SHUTDOWN}    shutdown cost in US dollars
#     4.  C{NCOST}       number of cost coefficients to follow for polynomial
#     cost function, or number of data points for piecewise linear
#     5.  C{COST}        1st column of cost parameters
#     cost data defining total cost function
#     For polynomial cost (highest order coeff first)::
#         e.g. cn, ..., c1, c0
#     where the polynomial is C{c0 + c1*P + ... + cn*P^n}
#     For piecewise linear cost::
#         x0, y0, x1, y1, x2, y2, ...
#     where C{x0 < x1 < x2 < ...} and the points C{(x0,y0), (x1,y1),
#     (x2,y2), ...} are the end- and break-points of the total cost function.
#
# additional constants, used to assign/compare values in the C{MODEL} column
#     1.  C{PW_LINEAR}   piecewise linear generator cost model
#     2.  C{POLYNOMIAL}  polynomial generator cost model
#
# @author: Ray Zimmerman (PSERC Cornell)
# @author: Richard Lincoln
# """
# define cost models
PW_LINEAR   = 1
POLYNOMIAL  = 2

# define the indices
MODEL       = 0
STARTUP     = 1
SHUTDOWN    = 2
NCOST       = 3
COST        = 4


#
# """Defines constants for named column indices to branch matrix.
#
# Some examples of usage, after defining the constants using the line above,
# are::
#
#     branch[3, BR_STATUS] = 0              # take branch 4 out of service
#     Ploss = branch[:, PF] + branch[:, PT] # compute real power loss vector
#
# The index, name and meaning of each column of the branch matrix is given
# below:
#
# columns 0-10 must be included in input matrix (in case file)
#     0.  C{F_BUS}       from bus number
#     1.  C{T_BUS}       to bus number
#     2.  C{BR_R}        resistance (p.u.)
#     3.  C{BR_X}        reactance (p.u.)
#     4.  C{BR_B}        total line charging susceptance (p.u.)
#     5.  C{RATE_A}      MVA rating A (long term rating)
#     6.  C{RATE_B}      MVA rating B (short term rating)
#     7.  C{RATE_C}      MVA rating C (emergency rating)
#     8.  C{TAP}         transformer off nominal turns ratio
#     9.  C{SHIFT}       transformer phase shift angle (degrees)
#     10. C{BR_STATUS}   initial branch status, 1 - in service, 0 - out of service
#     11. C{ANGMIN}      minimum angle difference, angle(Vf) - angle(Vt) (degrees)
#     12. C{ANGMAX}      maximum angle difference, angle(Vf) - angle(Vt) (degrees)
#
# columns 13-16 are added to matrix after power flow or OPF solution
# they are typically not present in the input matrix
#     13. C{PF}          real power injected at "from" bus end (MW)
#     14. C{QF}          reactive power injected at "from" bus end (MVAr)
#     15. C{PT}          real power injected at "to" bus end (MW)
#     16. C{QT}          reactive power injected at "to" bus end (MVAr)
#
# columns 17-18 are added to matrix after OPF solution
# they are typically not present in the input matrix
#
# (assume OPF objective function has units, C{u})
#     17. C{MU_SF}       Kuhn-Tucker multiplier on MVA limit at "from" bus (u/MVA)
#     18. C{MU_ST}       Kuhn-Tucker multiplier on MVA limit at "to" bus (u/MVA)
#
# columns 19-20 are added to matrix after SCOPF solution
# they are typically not present in the input matrix
#
# (assume OPF objective function has units, C{u})
#     19. C{MU_ANGMIN}   Kuhn-Tucker multiplier lower angle difference limit
#     20. C{MU_ANGMAX}   Kuhn-Tucker multiplier upper angle difference limit
#
# @author: Ray Zimmerman (PSERC Cornell)
# @author: Richard Lincoln
# """

# define the indices
F_BUS       = 0    # f, from bus number
T_BUS       = 1    # t, to bus number
BR_R        = 2    # r, resistance (p.u.)
BR_X        = 3    # x, reactance (p.u.)
BR_B        = 4    # b, total line charging susceptance (p.u.)
RATE_A      = 5    # rateA, MVA rating A (long term rating)
RATE_B      = 6    # rateB, MVA rating B (short term rating)
RATE_C      = 7    # rateC, MVA rating C (emergency rating)
TAP         = 8    # ratio, transformer off nominal turns ratio
SHIFT       = 9    # angle, transformer phase shift angle (degrees)
BR_STATUS   = 10   # initial branch status, 1 - in service, 0 - out of service
ANGMIN      = 11   # minimum angle difference, angle(Vf) - angle(Vt) (degrees)
ANGMAX      = 12   # maximum angle difference, angle(Vf) - angle(Vt) (degrees)

# included in power flow solution, not necessarily in input
PF          = 13   # real power injected at "from" bus end (MW)
QF          = 14   # reactive power injected at "from" bus end (MVAr)
PT          = 15   # real power injected at "to" bus end (MW)
QT          = 16   # reactive power injected at "to" bus end (MVAr)

# included in opf solution, not necessarily in input
# assume objective function has units, u
MU_SF       = 17   # Kuhn-Tucker multiplier on MVA limit at "from" bus (u/MVA)
MU_ST       = 18   # Kuhn-Tucker multiplier on MVA limit at "to" bus (u/MVA)
MU_ANGMIN   = 19   # Kuhn-Tucker multiplier lower angle difference limit
MU_ANGMAX   = 20   # Kuhn-Tucker multiplier upper angle difference limit

# """Defines constants for named column indices to areas matrix.
#
# The index, name and meaning of each column of the areas matrix is given below:
#
# columns 0-1
#     0.  C{AREA_I}           area number
#     1.  C{PRICE_REF_BUS}    price reference bus for this area
#
# @author: Ray Zimmerman (PSERC Cornell)
# @author: Richard Lincoln
# """

# define the indices
AREA_I          = 0    # area number
PRICE_REF_BUS   = 1    # price reference bus for this area

#
# """Defines constants for named column indices to dcline matrix.
#
# Some examples of usage, after defining the constants using the line above,
# are:
#
#   ppc.dcline(4, dcline['BR_STATUS']) = 0          take branch 4 out of service
#
# The index, name and meaning of each column of the branch matrix is given
# below:
#
# columns 1-17 must be included in input matrix (in case file)
#  1  F_BUS       f, "from" bus number
#  2  T_BUS       t,  "to"  bus number
#  3  BR_STATUS   initial branch status, 1 - in service, 0 - out of service
#  4  PF          MW flow at "from" bus ("from" -> "to")
#  5  PT          MW flow at  "to"  bus ("from" -> "to")
#  6  QF          MVAr injection at "from" bus ("from" -> "to")
#  7  QT          MVAr injection at  "to"  bus ("from" -> "to")
#  8  VF          voltage setpoint at "from" bus (p.u.)
#  9  VT          voltage setpoint at  "to"  bus (p.u.)
# 10  PMIN        lower limit on PF (MW flow at "from" end)
# 11  PMAX        upper limit on PF (MW flow at "from" end)
# 12  QMINF       lower limit on MVAr injection at "from" bus
# 13  QMAXF       upper limit on MVAr injection at "from" bus
# 14  QMINT       lower limit on MVAr injection at  "to"  bus
# 15  QMAXT       upper limit on MVAr injection at  "to"  bus
# 16  LOSS0       constant term of linear loss function (MW)
# 17  LOSS1       linear term of linear loss function (MW/MW)
#                 (loss = LOSS0 + LOSS1 * PF)
#
# columns 18-23 are added to matrix after OPF solution
# they are typically not present in the input matrix
#                 (assume OPF objective function has units, u)
# 18  MU_PMIN     Kuhn-Tucker multiplier on lower flow lim at "from" bus (u/MW)
# 19  MU_PMAX     Kuhn-Tucker multiplier on upper flow lim at "from" bus (u/MW)
# 20  MU_QMINF    Kuhn-Tucker multiplier on lower VAr lim at "from" bus (u/MVAr)
# 21  MU_QMAXF    Kuhn-Tucker multiplier on upper VAr lim at "from" bus (u/MVAr)
# 22  MU_QMINT    Kuhn-Tucker multiplier on lower VAr lim at  "to"  bus (u/MVAr)
# 23  MU_QMAXT    Kuhn-Tucker multiplier on upper VAr lim at  "to"  bus (u/MVAr)
#
# @see: L{toggle_dcline}
# """

## define the indices
dcline = {
    'F_BUS':     0,     ## f, "from" bus number
    'T_BUS':     1,     ## t,  "to"  bus number
    'BR_STATUS': 2,     ## initial branch status, 1 - in service, 0 - out of service
    'PF':        3,     ## MW flow at "from" bus ("from" -> "to")
    'PT':        4,     ## MW flow at  "to"  bus ("from" -> "to")
    'QF':        5,     ## MVAr injection at "from" bus ("from" -> "to")
    'QT':        6,     ## MVAr injection at  "to"  bus ("from" -> "to")
    'VF':        7,     ## voltage setpoint at "from" bus (p.u.)
    'VT':        8,     ## voltage setpoint at  "to"  bus (p.u.)
    'PMIN':      9,     ## lower limit on PF (MW flow at "from" end)
    'PMAX':     10,     ## upper limit on PF (MW flow at "from" end)
    'QMINF':    11,     ## lower limit on MVAr injection at "from" bus
    'QMAXF':    12,     ## upper limit on MVAr injection at "from" bus
    'QMINT':    13,     ## lower limit on MVAr injection at  "to"  bus
    'QMAXT':    14,     ## upper limit on MVAr injection at  "to"  bus
    'LOSS0':    15,     ## constant term of linear loss function (MW)
    'LOSS1':    16,     ## linear term of linear loss function (MW)
    'MU_PMIN':  17,     ## Kuhn-Tucker multiplier on lower flow lim at "from" bus (u/MW)
    'MU_PMAX':  18,     ## Kuhn-Tucker multiplier on upper flow lim at "from" bus (u/MW)
    'MU_QMINF': 19,     ## Kuhn-Tucker multiplier on lower VAr lim at "from" bus (u/MVAr)
    'MU_QMAXF': 20,     ## Kuhn-Tucker multiplier on upper VAr lim at "from" bus (u/MVAr)
    'MU_QMINT': 21,     ## Kuhn-Tucker multiplier on lower VAr lim at  "to"  bus (u/MVAr)
    'MU_QMAXT': 22      ## Kuhn-Tucker multiplier on upper VAr lim at  "to"  bus (u/MVAr)
}

