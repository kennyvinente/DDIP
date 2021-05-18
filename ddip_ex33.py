#### DDIP HTUC with no block aggregation

from numpy import array, sin, cos, zeros, linspace, append, random, max
import pandas as pd
from gurobipy import *


ntypes = 4
nperiods = 12

nblocks = 2

CVU = array([9.95,10.2,11,11])
CF = array([300,210,120,120])
CP = array([800,380,110,110])
cdef = 60

Pmin = array([70,40,30,30])
Pmax = array([500,250,150,150])

demand = array([600,800,700,950,867, 885, 672, 783, 771, 638, 876, 717])

upt = (array([2,2,2,2])/2).astype(int)
dnt = (array([2,2,4,4])/2).astype(int)

ro = 1
fcf0 = array([11.82,11.18])
fcf1 = array([208626,204594])

Vmax = 9000
V0 = Vmax
Y = zeros(nperiods)

def UC_model_1block(period,u0,w0,vol0,MILP_LP):

    m = Model()
    m.Params.OutputFlag = 0
    gt = m.addVars(ntypes, name="gt")
    if MILP_LP:
        v = m.addVars(ntypes, vtype='B', name="v")
        w = m.addVars(ntypes, vtype='B', name="w")
        u = m.addVars(ntypes, vtype='B', name="u")
        b = m.addVar(name="b_gh", vtype='B')
    else:
        v = m.addVars(ntypes, vtype='C', name="v", ub=1)
        w = m.addVars(ntypes, vtype='C', name="w", ub=1)
        u = m.addVars(ntypes, vtype='C', name="u", ub=1)
        b = m.addVar(name="b_gh", vtype='C', ub=1)

    gh = m.addVar(name="gh", ub=750)
    vol = m.addVar(name="vol", ub=Vmax)
    q = m.addVar(name="q", ub=750)
    s = m.addVar(name="s")
    y = m.addVar(name="vol*b")

    alpha = m.addVar(name='alpha')

    sl = m.addVar(name="def", ub=max(demand))

    Cost = m.addVar(name='Cost')
    Theta = m.addVar(name='Theta')

    z_u = m.addVars(ntypes,name='z_u')
    z_w = m.addVars((2,3),name='z_w')
    z_vol = m.addVar(name='z_vol')

    for i in range(ntypes):
        m.addConstr(z_u[i] == u0[i],name='aux_u%d'%i)
    for i in range(2,ntypes):
        m.addConstr(z_w[i] == w0[i-2],name='aux_w%d'%i) #these variables need to be updated to sum(w[k], k=T-1,...,DT)
    m.addConstr(z_vol == vol0, name='aux_vol')


    m.addConstrs(v[type] - w[type] == u[type] - z_u[type] for type in range(ntypes))
    m.addConstrs(v[type] <= u[type] for type in range(ntypes))
    for type in range(ntypes):
        if dnt[type] > 1 and period >= 1:
            m.addConstr(w[type] + z_w[type] <= 1 - u[type])
        else:
            m.addConstr(w[type] <= u[type])


    # m.addConstr(gh == ro*q + 0.005*vol*b)

    m.addConstr(y <= Vmax*b)
    m.addConstr(y <= vol)
    m.addConstr(y >= vol - Vmax*(1-b))

    m.addConstr(gh == ro*q + 0.005*y)


    m.addConstr(vol == z_vol - (q + s - Y[period-1]))

    if period == nperiods-1:
        m.addConstrs(alpha + fcf0[i]*vol >= fcf1[i] for i in range(2))

    m.addConstrs((gt[type] >= Pmin[type]*u[type]) for type in range(ntypes))

    m.addConstrs((gt[type] <= Pmax[type]*u[type]) for type in range(ntypes))

    m.addConstr(quicksum(gt[type] for type in range(ntypes)) + gh + sl == demand[period])

    active = quicksum(CF[type]*u[type] for type in range(ntypes))

    per_mw = quicksum(CVU[type]*gt[type] for type in range(ntypes))

    startup_obj = quicksum(CP[type]*v[type] for type in range(ntypes))

    deficit = cdef*sl

    m.addConstr(Cost == active + per_mw + startup_obj + deficit + alpha)

    m.setObjective(Cost + Theta)


    m._gt = gt
    m._v = v
    m._w = w
    m._u = u

    m._gh = gh
    m._vol = vol
    m._q = q
    m._s = s

    m._sl = sl
    m._alpha = alpha

    m._Cost = Cost
    m._Theta = Theta

    return m


maxiter = 5; # Maximum number of DDIP iterations
tolerance = 0.01; # Percentage tolerance stopping criteria

m_f = {}
m_b = {}
u = {}
w = {}
vol = {}
u[0,0] = [0,0,0,0]
w[0,0] = [0,0]
vol[0,0] = Vmax
Cost_f = {}
j=0
UB_Vector = array([])
best_UB_Vector = array([])
LB_Vector = array([])
best_LB_Vector = array([])
gap_Vector = array([])
best_gap_Vector = array([])
obj = {}
mu = {}


for period in range(1,nperiods+1): # Forward pass

    m_f[period] = UC_model_1block(period-1,u[j,period-1],w[j,period-1],vol[j,period-1],1);
    m_f[period].optimize();
    m_f[period].write('model_py.lp')

    u[j,period] = [m_f[period]._u[i].x for i in m_f[period]._u]
    w[j,period] = [m_f[period]._w[i].x for i in range(2,4)]
    vol[j,period] = m_f[period]._vol.x

    Cost_f[j,period] = m_f[period]._Cost.x
# End forward pass


GT,U,V,W = zeros((nperiods,ntypes)),zeros((nperiods,ntypes)),zeros((nperiods,ntypes)),zeros((nperiods,ntypes))
GH,VOL,Q,S = zeros(nperiods),zeros(nperiods),zeros(nperiods),zeros(nperiods)
SL = zeros(nperiods)
ALPHA = 0
for period in range(1,nperiods+1):
    GT[period-1] = array([m_f[period]._gt[i].x for i in range(ntypes)])
    U[period-1] = array([m_f[period]._u[i].x for i in range(ntypes)])
    V[period-1] = array([m_f[period]._v[i].x for i in range(ntypes)])
    W[period-1] = array([m_f[period]._w[i].x for i in range(ntypes)])
GT = GT.T
U = U.T
V = V.T
W = W.T

GH = array([m_f[i]._gh.x for i in range(1,nperiods+1)])
VOL = array([m_f[i]._vol.x for i in range(1,nperiods+1)])
Q = array([m_f[i]._q.x for i in range(1,nperiods+1)])
S = array([m_f[i]._s.x for i in range(1,nperiods+1)])
SL = array([m_f[i]._sl.x for i in range(1,nperiods+1)])
ALPHA = m_f[nperiods]._alpha.x

UB = sum(Cost_f[j,p] for p in range(1,nperiods+1));
UB_Vector = append(UB_Vector,UB)
best_UB = min(UB_Vector);
best_UB_Vector = append(best_UB_Vector, best_UB)


obj[j,nperiods+1] = 0; mu[j,nperiods+1] = 7*[0];

for period in range(nperiods,0,-1):
    m_b[period] = UC_model_1block(period-1,u[j,period-1],w[j,period-1],vol[j,period-1],0);

    if 1 <= period <= nperiods-1:
        m_f[period].addConstr(m_f[period]._Theta >= obj[j,period+1] + quicksum(mu[j,period+1][i]*(m_f[period]._u[i] - u[j,period][i]) for i in range(4)) +
                              quicksum(mu[j,period+1][4+i]*(m_f[period]._w[2+i] - w[j,period][i]) for i in range(2)) +
                              mu[j,period+1][6]*(m_f[period]._vol - vol[j,period] )) # Cut added to forward problem
        m_b[period].addConstr(m_b[period]._Theta >= obj[j,period+1] + quicksum(mu[j,period+1][i]*(m_b[period]._u[i] - u[j,period][i]) for i in range(4)) +
                              quicksum(mu[j,period+1][4+i]*(m_b[period]._w[2+i] - w[j,period][i]) for i in range(2)) +
                              mu[j,period+1][6]*(m_b[period]._vol - vol[j,period] )) # Cut added to forward problem)  # Cut added to backward problem

    m_b[period].optimize();
    m_b[period].write('model_py.lp')
    # mu[j,p] = m_b[p].pi[1];
    # mu[j,period] = m_b[period].pi[0:4];
    mu[j,period] = m_b[period].pi[0:7];
    obj[j,period] = m_b[period].objval
# Backward pass end

LB = m_b[1].objval
LB_Vector = append(LB_Vector,LB)
gap = abs(UB-LB)/UB*100; # Percentage gap
gap_Vector = append(gap_Vector,gap)
best_gap = abs(best_UB-LB)/best_UB*100; # Percentage gap
best_gap_Vector = append(best_gap_Vector,best_gap)
print("Iteration: %d, Current UB = %.2f, Best UB = %.2f, Current LB = %.2f, cur_gap = %.2f, best_gap = %.2f" %(j,UB,best_UB,LB,round(gap,3),round(best_gap,3)))


# Backward sweep for first iteration ended

maxiter = 25

#################  Iterations 2,3,... start  ###################
while best_gap >= tolerance and j < maxiter:
    flag = 0
    j += 1; # Iteration number

    u[j,0] = u[0,0]
    w[j,0] = w[0,0]
    vol[j,0] = Vmax

    # cuts added after each forward and backward
    # if j >= 1:
    #     for period in range(1,nperiods):
    #         m_f[period].addConstr(m_f[period]._Theta >= obj[j-1,period+1] + quicksum(mu[j-1,period+1][i]*(m_f[period]._u[i] - u[j-1,period][i]) for i in range(4)) +
    #                               quicksum(mu[j-1,period+1][4+i]*(m_f[period]._w[2+i] - w[j-1,period][i]) for i in range(2)) +
    #                               mu[j-1,period+1][6]*(m_f[period]._vol - vol[j-1,period] )) # Cut added to forward problem
    #         m_b[period].addConstr(m_b[period]._Theta >= obj[j-1,period+1] + quicksum(mu[j-1,period+1][i]*(m_b[period]._u[i] - u[j-1,period][i]) for i in range(4)) +
    #                               quicksum(mu[j-1,period+1][4+i]*(m_b[period]._w[2+i] - w[j-1,period][i]) for i in range(2)) +
    #                               mu[j-1,period+1][6]*(m_b[period]._vol - vol[j-1,period] )) # Cut added to forward problem)  # Cut added to backward problem

    for period in range(1,nperiods+1):
        for i in range(4):
            m_f[period].setAttr("RHS", m_f[period].getConstrByName("aux_u%d"%i), u[j,period-1][i])
        for i in range(2,4):
            m_f[period].setAttr("RHS", m_f[period].getConstrByName("aux_w%d"%i), w[j,period-1][i-2])
        m_f[period].setAttr("RHS", m_f[period].getConstrByName("aux_vol"), vol[j,period-1])
        m_f[period].optimize()
        m_f[period].write('model_py.lp')

        u[j, period] = [m_f[period]._u[i].x for i in m_f[period]._u]
        w[j, period] = [m_f[period]._w[i].x for i in range(2,4)]
        vol[j,period] = m_f[period]._vol.x
        Cost_f[j, period] = m_f[period]._Cost.x

        # if m_f[period].status == 3: #feas cut
        #     Feas[j,period-1] = {'u':u[j,period-1],'w':w[j,period-1]}
        #     m_f[period-1].addConstr(quicksum((1-u[j,period-1][i])*m_f[period-1]._u[i] for i in range(4)) + quicksum((u[j,period-1][i])*(1-m_f[period-1]._u[i]) for i in range(4)) +
        #                             quicksum((1-w[j,period-1][i])*m_f[period-1]._w[2+i] for i in range(2)) + quicksum((w[j,period-1][i])*(1-m_f[period-1]._w[2+i]) for i in range(2)) >= 1)
        #     m_b[period-1].addConstr(quicksum((1-u[j,period-1][i])*m_b[period-1]._u[i] for i in range(4)) + quicksum((u[j,period-1][i])*(1-m_b[period-1]._u[i]) for i in range(4)) +
        #                             quicksum((1-w[j,period-1][i])*m_b[period-1]._w[2+i] for i in range(2)) + quicksum((w[j,period-1][i])*(1-m_b[period-1]._w[2+i]) for i in range(2)) >= 1)
        #     j -=1;
        #     flag = 1
        #     break
        # else:
        #     u[j, period] = [m_f[period]._u[i].x for i in m_f[period]._u]
        #     w[j, period] = [m_f[period]._w[i].x for i in range(2,4)]
        #     Cost_f[j, period] = m_f[period]._Cost.x

    if flag == 1:
        continue

    UB = sum(Cost_f[j, p] for p in range(1, nperiods + 1));
    UB_Vector = append(UB_Vector, UB)
    best_UB = min(UB_Vector);
    best_UB_Vector = append(best_UB_Vector, best_UB)

    if UB-best_UB <= 1e-1:
        GT,U,V,W = zeros((nperiods,ntypes)),zeros((nperiods,ntypes)),zeros((nperiods,ntypes)),zeros((nperiods,ntypes))
        GH,VOL,Q,S = zeros(nperiods),zeros(nperiods),zeros(nperiods),zeros(nperiods)
        SL = zeros(nperiods)
        for period in range(1,nperiods+1):
            GT[period-1] = array([m_f[period]._gt[i].x for i in range(ntypes)])
            U[period-1] = array([m_f[period]._u[i].x for i in range(ntypes)])
            V[period-1] = array([m_f[period]._v[i].x for i in range(ntypes)])
            W[period-1] = array([m_f[period]._w[i].x for i in range(ntypes)])
        GT = GT.T
        U = U.T
        V = V.T
        W = W.T

        GH = array([m_f[i]._gh.x for i in range(1,nperiods+1)])
        VOL = array([m_f[i]._vol.x for i in range(1,nperiods+1)])
        Q = array([m_f[i]._q.x for i in range(1,nperiods+1)])
        S = array([m_f[i]._s.x for i in range(1,nperiods+1)])
        SL = array([m_f[i]._sl.x for i in range(1,nperiods+1)])

        ALPHA = m_f[nperiods]._alpha.x


    # Forward pass for iteration j ended
    # Backward pass starting for iteration j

    obj[j, nperiods + 1] = 0; mu[j, nperiods + 1] = 7 * [0];

    for period in range(nperiods, 0, -1):

        if 1 <= period <= nperiods-1:
            m_f[period].addConstr(m_f[period]._Theta >= obj[j,period+1] + quicksum(mu[j,period+1][i]*(m_f[period]._u[i] - u[j,period][i]) for i in range(4)) +
                                  quicksum(mu[j,period+1][4+i]*(m_f[period]._w[2+i] - w[j,period][i]) for i in range(2)) +
                                  mu[j,period+1][6]*(m_f[period]._vol - vol[j,period] )) # Cut added to forward problem
            m_b[period].addConstr(m_b[period]._Theta >= obj[j,period+1] + quicksum(mu[j,period+1][i]*(m_b[period]._u[i] - u[j,period][i]) for i in range(4)) +
                                  quicksum(mu[j,period+1][4+i]*(m_b[period]._w[2+i] - w[j,period][i]) for i in range(2)) +
                                  mu[j,period+1][6]*(m_b[period]._vol - vol[j,period] )) # Cut added to forward problem)  # Cut added to backward problem

        for i in range(4):
            m_b[period].setAttr("RHS", m_b[period].getConstrByName("aux_u%d"%i), u[j,period-1][i])
        for i in range(2, 4):
            m_b[period].setAttr("RHS", m_b[period].getConstrByName("aux_w%d" % i), w[j, period - 1][i - 2])
        m_b[period].setAttr("RHS", m_b[period].getConstrByName("aux_vol"), vol[j,period-1])

        m_b[period].optimize();
        m_b[period].write('model_py.lp')
        mu[j, period] = m_b[period].pi[0:7];
        obj[j, period] = m_b[period].objval

    LB = m_b[1].objval
    LB_Vector = append(LB_Vector, LB)
    gap = abs(UB - LB) / UB * 100;  # Percentage gap
    gap_Vector = append(gap_Vector, gap)
    best_gap = abs(best_UB - LB) / best_UB * 100;  # Percentage gap
    best_gap_Vector = append(best_gap_Vector, best_gap)
    print("Iteration: %d, Current UB = %.2f, Best UB = %.2f, Current LB = %.2f, cur_gap = %.2f, best_gap = %.2f" % (
    j, UB, best_UB, LB, round(gap, 3), round(best_gap, 3)))

    # Backward pass for j iteration ended
     # Iteration j ends
