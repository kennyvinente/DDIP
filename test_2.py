
from numpy import array, zeros, append
import pandas as pd
from gurobipy import *

ntypes = 4
nperiods = 5

CVU = array([9.95,10.2,11,11])
CF = array([300,210,120,120])
CP = array([800,380,110,110])
cdef = 60

Pmin = array([70,40,30,30])
Pmax = array([500,250,150,150])

demand = array([600,800,700,950,867, 885, 672, 783, 771, 638, 876, 717])


upt = (array([8,1,1,1])).astype(int)
dnt = (array([1,1,1,1])).astype(int)

ro = 1
fcf0 = array([11.82,11.18])
fcf1 = array([208626,204594])

Vmax = 9000
V0 = Vmax
Y = zeros(nperiods)


valida = 0

GT,V,W,U = zeros((ntypes,nperiods)),zeros((ntypes,nperiods)),zeros((ntypes,nperiods)),zeros((ntypes,nperiods))

def MILP_unique():

    model = Model('PowerGeneration')
    model.Params.OutputFlag = 0
    gt = model.addVars(ntypes, nperiods, name="gt")

    v = model.addVars(ntypes, nperiods, vtype='C', name="v", ub=1)
    w = model.addVars(ntypes, nperiods, vtype='C', name="w", ub=1)
    u = model.addVars(ntypes, nperiods, vtype='C', name="u", ub=1)
    b = model.addVars(nperiods, vtype="C", name="b_gh", ub=1)




    gh = model.addVars(nperiods, name="gh", ub=750)
    # vol = model.addVars(nperiods, name="vol", ub=Vmax)
    # q = model.addVars(nperiods, name="q", ub=750)
    # s = model.addVars(nperiods, name="s")

    # y = model.addVars(nperiods, name="vol*b")

    slp = model.addVars(nperiods, name='def+', ub=100)
    sln = model.addVars(nperiods, name='def-', ub=100)

    alpha = model.addVar(name='alpha')

    for period in range(nperiods):
        if period == 0:
            model.addConstrs(v[type,period] - w[type,period] == u[type,period] for type in range(ntypes))
        else:
            model.addConstrs(v[type,period] - w[type,period] == u[type,period] - u[type,period-1] for type in range(ntypes))


    for type in range(ntypes):
        for period in range(nperiods):
            model.addConstr(quicksum(v[type, l] for l in range(max(0,period+1-upt[type]),period + 1)) <= u[type, period], name='upt[%d,%d]'%(type,period+1))

    for type in range(ntypes):
        for period in range(nperiods):
            model.addConstr(quicksum(w[type, l] for l in range(max(0,period+1-dnt[type]),period + 1)) <= 1 - u[type, period])


    min_output = model.addConstrs((gt[type, period] >= Pmin[type]*u[type, period])
                                  for type in range(ntypes) for period in range(nperiods))

    max_output = model.addConstrs((gt[type, period] <= Pmax[type]*u[type, period])
                                  for type in range(ntypes) for period in range(nperiods))

    meet_demand = model.addConstrs(quicksum(gt[type, period] for type in range(ntypes))  + gh[period] + slp[period] - sln[period] == demand[period]
                                   for period in range(nperiods))

    # for period in range(nperiods):
    #     model.addConstr(y[period] <= Vmax*b[period])
    #     model.addConstr(y[period] <= vol[period])
    #     model.addConstr(y[period] >= vol[period] - Vmax*(1-b[period]))

    # hydro_generation = model.addConstrs(gh[period] == ro*q[period] + 0.005*vol[period]*b[period] for period in range(nperiods))
    # hydro_generation = model.addConstrs(gh[period] == ro*q[period] + 0.005*y[period] for period in range(nperiods))




    # for period in range(nperiods):
    #     if period == 0:
    #         model.addConstr(vol[period] == V0 - (q[period] + s[period] - Y[period]) )
    #     else:
    #         model.addConstr(vol[period] == vol[period-1] - (q[period] + s[period] - Y[period]))


    active = quicksum(CF[type]*u[type,period]
                        for type in range(ntypes) for period in range(nperiods))

    per_mw = quicksum(CVU[type]*gt[type,period]
                           for type in range(ntypes) for period in range(nperiods))

    startup_obj = quicksum(CP[type]*v[type,period]
                             for type in range(ntypes) for period in range(nperiods))

    deficit = quicksum(cdef*slp[period] + cdef*sln[period] for period in range(nperiods))


    # model.addConstr(alpha + 11.82*vol[nperiods-1] >= 208626)
    # model.addConstr(alpha + 11.18*vol[nperiods-1] >= 204594)


    model.setObjective(active + per_mw + startup_obj + deficit + alpha)


    model.write('teste_new.lp')
    model.optimize()


    return model.objval

MILP = MILP_unique()
print(MILP)
MILP_LP=0

def UC_model_1block(period,svars,MILP_LP):

    m = Model()
    m.Params.OutputFlag = 0
    gt = m.addVars(ntypes, name="gt")
    if MILP_LP:
        v = m.addVars(ntypes, vtype='B', name="v")
        w = m.addVars(ntypes, vtype='B', name="w")
        u = m.addVars(ntypes, vtype='B', name="u")
        # b = m.addVar(name="b_gh", vtype='B')
    else:
        v = m.addVars(ntypes, vtype='C', name="v", ub=1)
        w = m.addVars(ntypes, vtype='C', name="w", ub=1)
        u = m.addVars(ntypes, vtype='C', name="u", ub=1)
        # b = m.addVar(name="b_gh", vtype='C', ub=1)

    gh = m.addVar(name="gh", ub=750)
    # vol = m.addVar(name="vol", ub=Vmax)
    # q = m.addVar(name="q", ub=750)
    # s = m.addVar(name="s")
    # y = m.addVar(name="vol*b")

    alpha = m.addVar(name='alpha')

    slp = m.addVar(name="def", ub=100)
    sln = m.addVar(name="def", ub=100)

    Cost = m.addVar(name='Cost')
    Theta = m.addVar(name='Theta')

    aux_tu = m.addVars(ntypes,name='aux_u',ub=1)

    m.addConstrs(v[type] - w[type] == u[type] - aux_tu[type] for type in range(ntypes))

    # if period > 1:
    #     aux_tv = m.addVars(period-1,name='aux_v',ub=1)
    # if period == 2:
    #     m.addConstr(aux_tv[0] == svars['tv'][0][0,1], name='ztv%d%d'%(period,1))
    # if period == 3:
    #     m.addConstr(aux_tv[0] == svars['tv'][0][0,1], name='ztv%d%d'%(period,1))
    #     m.addConstr(aux_tv[1] == svars['tv'][0][0,2], name='ztv%d%d'%(period,2))
    # # m.addConstr(z_vol[0] == svars['vol'][0][0,period-1], name='zvol%d'%0)
    #

    #
    # if period == 1:
    #     m.addConstr(v[0] <= u[0])
    # else:
    #     m.addConstr(quicksum(aux_tv[l] for l in aux_tv) + v[0] <= u[0])

    aux_tv = m.addVars(period-1, name='aux_v', ub=1)

    for l in range(period-1):
        m.addConstr(aux_tv[l] == svars['tv'][ite][0,l+1], name='ztv%d%d'%(period,l))
    m.addConstr(quicksum(aux_tv[i] for i in aux_tv) + v[0] <= u[0])

    # aux_tv = m.addVars(int(upt[0]-1), name='aux_v', ub=1)
    # for l in range(upt[0]-1,0,-1):
    #     m.addConstr(aux_tv[l-1] == (svars['tv'][ite][0,l-upt[0]+period] if l-upt[0]+period >= 0 else 0), name='ztv%d%d'%(period,l-1))

    # m.addConstr(quicksum(aux_tv[i] for i in range(upt[0]-1)) + v[0] <= u[0])

    # z_vol = m.addVars(1,name='z_vol')

    # for type in range(ntypes):
    #     for period in range(nperiods):
    #         if period == 0:
    #             model.addConstr(aux_v[type,period] == 0, name='aux_v[%d,%d]'%(type,period+1))
    #             model.addConstr(aux_w[type,period] == 0, name='aux_w[%d,%d]'%(type,period+1))
    #         else:
    #             model.addConstr(aux_v[type,period] == quicksum(v[type,l] for l in range(max(0,period+1-upt[type]),period)), name='aux_v[%d,%d]'%(type,period+1))
    #             model.addConstr(aux_w[type,period] == quicksum(w[type,l] for l in range(max(0,period+1-dnt[type]),period)), name='aux_w[%d,%d]'%(type,period+1))
    #
    # for type in range(ntypes):
    #     for period in range(nperiods):
    #         # model.addConstr(quicksum(v[type, l] for l in range(max(0,period+1-upt[type]),period + 1)) <= u[type, period], name='upt[%d,%d]'%(type,period+1))
    #         model.addConstr(v[type,period] + aux_v[type,period] <= u[type, period], name='upt[%d,%d]'%(type,period+1))
    #
    # for type in range(ntypes):
    #     for period in range(nperiods):
    #         # model.addConstr(quicksum(w[type, l] for l in range(max(0,period+1-dnt[type]),period + 1)) <= 1 - u[type, period])
    #         model.addConstr(w[type,period] + aux_w[type,period] <= 1 - u[type, period])

    for type in range(ntypes):
        m.addConstr(aux_tu[type] == svars['tu'][0][type,period-1],name='ztu%d'%type)
        # for l in range(period+1-upt[type],period):
    # for l in range(upt[0]-1):
    #     if period+1-upt[0] + l < 0:
    #         m.addConstr(aux_tv[0,l] == 0, name='ztv%d%d'%(0,l))
    #     else:
    #         m.addConstr(aux_tv[0,l] == svars['tv'][0][0,period+1-upt[0] + l], name='ztv%d%d'%(0,l))



    for type in range(1,ntypes):
        m.addConstr(v[type] <= u[type])

    for type in range(ntypes):
        m.addConstr(w[type] <= 1 - u[type])


    # m.addConstr(y <= Vmax*b)
    # m.addConstr(y <= vol)
    # m.addConstr(y >= vol - Vmax*(1-b))
    #
    # m.addConstr(gh == ro*q + 0.005*y)


    # m.addConstr(vol == z_vol[0] - (q + s - Y[period-1]))

    # if period == nperiods:
    #     m.addConstrs(alpha + fcf0[i]*vol >= fcf1[i] for i in range(2))

    m.addConstrs((gt[type] >= Pmin[type]*u[type]) for type in range(ntypes))

    m.addConstrs((gt[type] <= Pmax[type]*u[type]) for type in range(ntypes))

    m.addConstr(quicksum(gt[type] for type in range(ntypes)) + gh + slp - sln == demand[period-1])

    active = quicksum(CF[type]*u[type] for type in range(ntypes))

    per_mw = quicksum(CVU[type]*gt[type] for type in range(ntypes))

    startup_obj = quicksum(CP[type]*v[type] for type in range(ntypes))

    deficit = cdef*(slp + sln)

    m.addConstr(Cost == active + per_mw + startup_obj + deficit + alpha)

    m.setObjective(Cost + Theta)


    m._gt = gt
    m._v = v
    m._w = w
    m._u = u

    m._aux_u = aux_tu
    if period > 1:
        m._aux_v = aux_tv

    m._gh = gh
    # m._vol = vol
    # m._q = q
    # m._s = s



    m._slp = slp
    m._sln = sln
    m._alpha = alpha

    m._Cost = Cost
    m._Theta = Theta

    return m


ite = 0

m_f = {}
m_b = {}
svars = {}
svars['vol'],svars['tu'],svars['tv'] = {}, {}, {}
svars['vol'][0] = zeros((1,nperiods+1))
svars['tu'][0] = zeros((ntypes,nperiods+1))
svars['tv'][0] = zeros((ntypes,nperiods+1))
svars['vol'][0][0,0] = 9000

Cost_f = {}
ite=0
UB_Vector = array([])
best_UB_Vector = array([])
LB_Vector = array([])
best_LB_Vector = array([])
gap_Vector = array([])
best_gap_Vector = array([])
obj = {}
mu = {}
period = 1
# et = time.time()

for period in range(1,nperiods+1): # Forward pass

    # m_f[period] = UC_model_1block(period,svars,MILP_LP)
    m_f[period] = UC_model_1block(period,svars,0);

    m_f[period].optimize();
    m_f[period].write('model_py_%d.lp'%(period))


    svars['tu'][ite][:,period] = array([m_f[period]._u[i].x for i in range(ntypes)])
    svars['tv'][ite][:,period] = array([m_f[period]._v[i].x for i in range(ntypes)])
    # svars['vol'][ite][:,period] = m_f[period]._vol.x

    Cost_f[ite,period] = m_f[period].objval - m_f[period]._Theta.x
# print('Elapsed time for building model: %.2f seconds\n'%(time.time() - et))
# End forward pass

UB = sum(Cost_f[ite,p] for p in range(1,nperiods+1));
UB_Vector = append(UB_Vector,UB)
best_UB = min(UB_Vector);
best_UB_Vector = append(best_UB_Vector, best_UB)


mu['vol'],mu['tu'],mu['tv']= {}, {}, {}
mu['vol'][0],mu['tu'][0],mu['tv'][0] = {}, {}, {}

obj[ite,nperiods+1] = 0; mu['tu'][ite][nperiods+1] = zeros(ntypes); mu['tv'][ite]={}
for l in range(ntypes):
    mu['tv'][ite][l] = zeros((upt[l],nperiods+1))
mu['vol'][ite][nperiods+1] = zeros(1);


for period in range(nperiods,0,-1):
    m_b[period] = UC_model_1block(period,svars,0);

    if 1 <= period <= nperiods-1:
        f1,f2 = 0,0
        for i in range(ntypes):
            f1 += mu['tu'][ite][period+1][i]*(m_f[period]._u[i] - svars['tu'][ite][i,period])
        for l in range(period):
            if l == period-1:
                f2 += mu['tv'][ite][0][l-1,period+1]*(m_f[period]._v[0] - svars['tv'][ite][0,l])
            else:
                f2 += mu['tv'][ite][0][l-1,period+1]*(m_f[period]._aux_v[l] - svars['tv'][ite][0,l])
        m_f[period].addConstr(m_f[period]._Theta >= obj[ite,period+1] + f1 + f2) # Cut added to forward problem
        f1,f2= 0,0
        for i in range(ntypes):
            f1 += mu['tu'][ite][period+1][i]*(m_b[period]._u[i] - svars['tu'][ite][i,period])
        for l in range(period):
            if l == period-1:
                f2 += mu['tv'][ite][0][l-1,period+1]*(m_b[period]._v[0] - svars['tv'][ite][0,l])
            else:
                f2 += mu['tv'][ite][0][l-1,period+1]*(m_b[period]._aux_v[l] - svars['tv'][ite][0,l])
        m_b[period].addConstr(m_b[period]._Theta >= obj[ite,period+1] + f1 + f2) # # Cut added to backward problem

    # if 1 <= period <= nperiods-1:
    #     f1,f2,f3 = 0,0,0
    #     for i in range(ntypes):
    #         f1 += mu['tu'][ite][period+1][i]*(m_f[period]._u[i] - svars['tu'][ite][i,period])
    #     for l in range(upt[0]-1,0,-1):
    #         if l-upt[0]+period < 0:
    #             break
    #         if l == upt[0]-1:
    #             f2 += mu['tv'][ite][0][l-1,period+1]*(m_f[period]._v[0] - svars['tv'][ite][0,period])
    #         else:
    #             f2 += mu['tv'][ite][0][l-1,period+1]*(m_f[period]._aux_v[l] - svars['tv'][ite][0,l-upt[0]+period])
    #     m_f[period].addConstr(m_f[period]._Theta >= obj[ite,period+1] + f1 + f2 + f3) # Cut added to forward problem
    #     f1,f2,f3 = 0,0,0
    #     for i in range(ntypes):
    #         f1 += mu['tu'][ite][period+1][i]*(m_b[period]._u[i] - svars['tu'][ite][i,period])
    #     for l in range(upt[0]-1,0,-1):
    #         if l-upt[0]+period < 0:
    #             break
    #         if l == upt[0]-1:
    #             f2 += mu['tv'][ite][0][l-1,period+1]*(m_b[period]._v[0] - svars['tv'][ite][0,period])
    #         else:
    #             f2 += mu['tv'][ite][0][l-1,period+1]*(m_b[period]._aux_v[l] - svars['tv'][ite][0,l-upt[0]+period])
    #     m_b[period].addConstr(m_b[period]._Theta >= obj[ite,period+1] + f1 + f2 + f3) # # Cut added to backward problem

    m_b[period].optimize();
    m_b[period].write('model_py_%d.lp'%(period))

    if period > 1:
        mu['tu'][ite][period] = array([m_b[period].getConstrByName("ztu%d"%i).pi for i in range(ntypes)])
        for l in range(period-1):
            mu['tv'][ite][0][l,period] = m_b[period].getConstrByName("ztv%d%d"%(period,l)).pi

    # if period > 1:
    #     mu['tu'][ite][period] = array([m_b[period].getConstrByName("ztu%d"%i).pi for i in range(ntypes)])
    #     for l in range(upt[0]-1,0,-1):
    #         mu['tv'][ite][0][l-1,period] = m_b[period].getConstrByName("ztv%d%d"%(period,l-1)).pi


    obj[ite,period] = m_b[period].objval


# for i in range(nperiods,0,-1):
#     print(mu['tv'][0][i])
#     print(mu['tu'][0][i])


LB = m_b[1].objval
LB_Vector = append(LB_Vector,LB)
gap = abs(UB-LB)/UB*100; # Percentage gap
gap_Vector = append(gap_Vector,gap)
best_gap = abs(best_UB-LB)/best_UB*100; # Percentage gap
best_gap_Vector = append(best_gap_Vector,best_gap)
print("Iteration: %d, Current UB = %.2f, Best UB = %.2f, Current LB = %.2f, cur_gap = %.2f, best_gap = %.2f" %(ite,UB,best_UB,LB,round(gap,3),round(best_gap,3)))


# Backward sweep for first iteration ended

maxiter = 0; # Maximum number of DDIP iterations
tolerance = 1e-8; # Percentage tolerance stopping criteria



#################  Iterations 2,3,... start  ###################
while best_gap >= tolerance and ite < maxiter:
    ite += 1; # Iteration number

    svars['tu'][ite] = zeros((ntypes,nperiods+1))
    svars['tv'][ite] = zeros((ntypes,nperiods+1))


    for period in range(1,nperiods+1):
        if period > 1:
            for i in range(ntypes): #### ajeitar
                m_f[period].setAttr("RHS", m_f[period].getConstrByName("ztu%d"%i), svars['tu'][ite][i,period-1])
            for l in range(upt[0]-1,0,-1):
                m_f[period].setAttr("RHS", m_f[period].getConstrByName("ztv%d%d"%(period,l-1)), (svars['tv'][ite][0,l-upt[0]+period] if l-upt[0]+period >= 0 else 0))

        m_f[period].optimize()
        m_f[period].write('model_py_%d.lp'%(period))

        svars['tu'][ite][:,period] = array([m_f[period]._u[i].x for i in range(ntypes)])
        svars['tv'][ite][:,period] = array([m_f[period]._v[i].x for i in range(ntypes)])
        # svars['vol'][ite][:,period] = m_f[period]._vol.x

        Cost_f[ite,period] = m_f[period].objval - m_f[period]._Theta.x

    UB = sum(Cost_f[ite,p] for p in range(1,nperiods+1));
    UB_Vector = append(UB_Vector,UB)
    best_UB = min(UB_Vector);
    best_UB_Vector = append(best_UB_Vector, best_UB)


    # Forward pass for iteration j ended
    # Backward pass starting for iteration j

    mu['vol'][ite],mu['tu'][ite],mu['tv'][ite] = {}, {}, {}

    obj[ite,nperiods+1] = 0; mu['tu'][ite][nperiods+1] = zeros(ntypes); mu['tv'][ite]={}
    for l in range(ntypes):
        mu['tv'][ite][l] = zeros((upt[l],nperiods+1))
    mu['vol'][ite][nperiods+1] = zeros(1);

    for period in range(nperiods,0,-1):
        if period > 1:
            for i in range(ntypes): #### ajeitar
                m_b[period].setAttr("RHS", m_b[period].getConstrByName("ztu%d"%i), svars['tu'][ite][i,period-1])
            for l in range(upt[0]-1,0,-1):
                m_b[period].setAttr("RHS", m_b[period].getConstrByName("ztv%d%d"%(period,l-1)), (svars['tv'][ite][0,l-upt[0]+period] if l-upt[0]+period >= 0 else 0))


        if 1 <= period <= nperiods-1:
            f1,f2,f3 = 0,0,0
            for i in range(ntypes):
                f1 += mu['tu'][ite][period+1][i]*(m_f[period]._u[i] - svars['tu'][ite][i,period])
            for l in range(upt[0]-1,0,-1):
                if l-upt[0]+period < 0:
                    break
                if l == upt[0]-1:
                    f2 += mu['tv'][ite][0][l-1,period+1]*(m_f[period]._v[0] - svars['tv'][ite][0,period])
                else:
                    f2 += mu['tv'][ite][0][l-1,period+1]*(m_f[period]._aux_v[l] - svars['tv'][ite][0,l-upt[0]+period])
            m_f[period].addConstr(m_f[period]._Theta >= obj[ite,period+1] + f1 + f2 + f3) # Cut added to forward problem
            f1,f2,f3 = 0,0,0
            for i in range(ntypes):
                f1 += mu['tu'][ite][period+1][i]*(m_b[period]._u[i] - svars['tu'][ite][i,period])
            for l in range(upt[0]-1,0,-1):
                if l-upt[0]+period < 0:
                    break
                if l == upt[0]-1:
                    f2 += mu['tv'][ite][0][l-1,period+1]*(m_b[period]._v[0] - svars['tv'][ite][0,period])
                else:
                    f2 += mu['tv'][ite][0][l-1,period+1]*(m_b[period]._aux_v[l] - svars['tv'][ite][0,l-upt[0]+period])
            m_b[period].addConstr(m_b[period]._Theta >= obj[ite,period+1] + f1 + f2 + f3) # # Cut added to backward problem

        m_b[period].optimize();
        m_b[period].write('model_py_%d.lp'%(period))

        if period > 1:
            mu['tu'][ite][period] = array([m_b[period].getConstrByName("ztu%d"%i).pi for i in range(ntypes)])
            for l in range(upt[0]-1,0,-1):
                mu['tv'][ite][0][l-1,period] = m_b[period].getConstrByName("ztv%d%d"%(period,l-1)).pi

        obj[ite,period] = m_b[period].objval


    LB = m_b[1].objval
    LB_Vector = append(LB_Vector,LB)
    gap = abs(UB-LB)/UB*100; # Percentage gap
    gap_Vector = append(gap_Vector,gap)
    best_gap = abs(best_UB-LB)/best_UB*100; # Percentage gap
    best_gap_Vector = append(best_gap_Vector,best_gap)
    print("Iteration: %d, Current UB = %.2f, Best UB = %.2f, Current LB = %.2f, cur_gap = %.2f, best_gap = %.2f" %(ite,UB,best_UB,LB,round(gap,3),round(best_gap,3)))

    # Backward pass for j iteration ended
     # Iteration j ends



