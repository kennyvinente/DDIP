#### DDIP HTUC with block aggregation


from ddip_formulations_blk_f import *
import formulations as form

NT = 24
nblocks = array([8,4,2,2,2,3,3]).astype(int)
nblocks = 4*ones(6).astype(int)

if sum(nblocks) != NT:
    print("The decomposition doesn't match the Time Horizon")

MIPGap = True

def UC_model_1block(nblocks,svars,ite,t0,MILP_LP):

    m = Model()
    m._NT = NT
    m.Params.OutputFlag = 0

    if MILP_LP:
        fph_model(m,nblocks,t0,1)
        term_model(m,nblocks,t0,1)
        if MIPGap:
            m.setParam('MIPGAP',0.01)
    else:
        fph_model(m,nblocks,t0,0)
        term_model(m,nblocks,t0,0)
    uct_model(m,svars,ite,nblocks,t0)
    network_model(m,nblocks,t0)
    reserv_model(m,svars,ite,nblocks,t0,NT)
    f1 = fobj(m)

    Theta = m.addMVar(1, name='Theta')
    m.setObjective(Theta+f1)

    m._Theta = Theta

    return m


# initial conditions, and dictionary with variables stored

m_f = {}
m_b = {}
svars = {}
svars['tu'],svars['tv'],svars['tw'],svars['gt'] = {}, {}, {}, {}
svars['vol'], svars['turb'],svars['vert'] = {},{},{}
svars['vol'][0] = zeros((NH,NT+1))
svars['turb'][0] = zeros((NH,NT+1))
svars['vert'][0] = zeros((NH,NT+1))
svars['tu'][0] = zeros((NG,NT+1))
svars['tv'][0] = zeros((NG,NT+1))
svars['tw'][0] = zeros((NG,NT+1))
svars['gt'][0] = zeros((NG,NT+1))
svars['vol'][0][:,0] = 0.01*df_hidr['V0']*(df_hidr['VMAX']-df_hidr['VMIN']) + df_hidr['VMIN']
svars['turb'][0][:,0] = df_hidr['Q0'].to_numpy().astype(float)
svars['vert'][0][:,0] = df_hidr['S0'].to_numpy().astype(float)
for i in range(NG):
    if df_term['TON'][i] > 0:
        svars['tu'][0][i][0] = 1
        if df_term['TON'][i] == 1:
            svars['tv'][0][i][0] = 1
            svars['tw'][0][i][0] = 0
        else:
            svars['tv'][0][i][0] = 0
            svars['tw'][0][i][0] = 0
    else:
        if df_term['TON'][i][0] == 1:
            svars['tv'][0][i][0] = 0
            svars['tw'][0][i][0] = 1
        else:
            svars['tv'][0][i][0] = 0
            svars['tw'][0][i][0] = 0
svars['gt'][0][:,0] = df_term['P0'].to_numpy()

Cost_f = {}
ite=0
UB_Vector = array([])
best_UB_Vector = array([])
LB_Vector = array([])
best_LB_Vector = array([])
gap_Vector = array([])
best_gap_Vector = array([])
period_vector = array([]).astype(int)
obj = {}
mu = {}

# et = time.time()
tet = 0

period = 1
# # # # # Forward pass

for blk in range(nblocks.shape[0]):
    if blk == 0:
        period = int(1)
    else:
        period += int(nblocks[blk-1])
    period_vector = append(period_vector,period)

    m_f[period] = UC_model_1block(int(nblocks[blk]),svars,ite,period,1);
    m_f[period].update()
    # m_f[period] = UC_model_1block(int(nblocks[blk]),svars,ite,period,0);

    # et = time.time()
    # m_f[period].optimize();
    # tet += time.time() - et
    #
    # m_f[period].write('model_py_%d.lp'%(period))
    #
    # for j in range(nblocks[blk]):
    #     svars['tu'][ite][:,period+j] = array([m_f[period]._tu[i,j].x[0] for i in range(NG)])
    #     svars['tv'][ite][:,period+j] = array([m_f[period]._tv[i,j].x[0] for i in range(NG)])
    #     svars['tw'][ite][:,period+j] = array([m_f[period]._tw[i,j].x[0] for i in range(NG)])
    #     svars['gt'][ite][:,period+j] = array([m_f[period]._gt[i,j].x[0] for i in range(NG)])
    #     svars['vol'][ite][:,period+j] = array([m_f[period]._vol[i,j].x[0] for i in range(NH)])
    #     svars['turb'][ite][:,period+j] = array([m_f[period]._turb[i,j].x[0] for i in range(NH)])
    #     svars['vert'][ite][:,period+j] = array([m_f[period]._vert[i,j].x[0] for i in range(NH)])
    #
    # Cost_f[ite,period] = m_f[period].objval - m_f[period]._Theta[0].x

# UB = sum(Cost_f[ite,p] for p in period_vector)
# UB_Vector = append(UB_Vector,UB)
# best_UB = min(UB_Vector);
# best_UB_Vector = append(best_UB_Vector, best_UB)

# # # # # End forward pass


# # # # # Pre-solving
model = Model()
model._FPH_MODEL = 'CH'
model._CQ = 3600*1e-6
model.setParam('TimeLimit', 600)
model.setParam('MIPGAP',0.01)

form.fph_model(model,NT,0)
form.term_model(model,NT,0)
form.uct_model(model, NT)
form.network_model(model,NT)
form.reserv_model(model, NT)
form.fobj(model,NT)

et = time.time()
model.optimize()
tet += time.time() - et


for j in range(NT):
    svars['tu'][ite][:,j+1] = array([model._tu[i,j].x[0] for i in range(NG)])
    svars['tv'][ite][:,j+1] = array([model._tv[i,j].x[0] for i in range(NG)])
    svars['tw'][ite][:,j+1] = array([model._tw[i,j].x[0] for i in range(NG)])
    svars['gt'][ite][:,j+1] = array([model._gt[i,j].x[0] for i in range(NG)])
    svars['vol'][ite][:,j+1] = array([model._vol[i,j].x[0] for i in range(NH)])
    svars['turb'][ite][:,j+1] = array([model._turb[i,j].x[0] for i in range(NH)])
    svars['vert'][ite][:,j+1] = array([model._vert[i,j].x[0] for i in range(NH)])

# m_presolve = {}
# for blk in range(nblocks.shape[0]):
#     if blk == 0:
#         period = int(1)
#     else:
#         period += int(nblocks[blk-1])
#     period_vector = append(period_vector,period)
#
#     m_presolve[period] = UC_model_1block(int(nblocks[blk]),svars,ite,period,0);
#
#     # et = time.time()
#     m_presolve[period].optimize();
#     # tet += time.time() - et
#
#
#     for j in range(nblocks[blk]):
#         svars['tu'][ite][:,period+j] = array([m_presolve[period]._tu[i,j].x[0] for i in range(NG)])
#         svars['tv'][ite][:,period+j] = array([m_presolve[period]._tv[i,j].x[0] for i in range(NG)])
#         svars['tw'][ite][:,period+j] = array([m_presolve[period]._tw[i,j].x[0] for i in range(NG)])
#         svars['gt'][ite][:,period+j] = array([m_presolve[period]._gt[i,j].x[0] for i in range(NG)])
#         svars['vol'][ite][:,period+j] = array([m_presolve[period]._vol[i,j].x[0] for i in range(NH)])
#         svars['turb'][ite][:,period+j] = array([m_presolve[period]._turb[i,j].x[0] for i in range(NH)])
#         svars['vert'][ite][:,period+j] = array([m_presolve[period]._vert[i,j].x[0] for i in range(NH)])


# # # # # End pre-solving




mu['vol'],mu['turb'],mu['vert'],mu['tu'],mu['tv'],mu['tw'],mu['gt'] = {}, {}, {}, {}, {},{},{}
mu['vol'][0],mu['turb'][0],mu['vert'][0],mu['tu'][0],mu['tv'][0],mu['tw'][0],mu['gt'][0] = {}, {}, {}, {}, {},{},{}
obj[ite,NT+1] = 0; mu['tu'][ite][NT+1] = zeros(NG);
mu['gt'][ite][NT+1] = zeros(NG);  mu['vol'][ite][NT+1] = zeros(NH);

for i in range(NG):
    if df_term['UPTIME'][i] > 1:
        mu['tv'][ite][i] = zeros((df_term['UPTIME'][i],NT+1))

for i in range(NG):
    if df_term['DOWNTIME'][i] > 1:
        mu['tw'][ite][i] = zeros((df_term['DOWNTIME'][i],NT+1))

for i in range(NH):
    if df_hidr['TRAVELTIME'][i] > 0:
        mu['turb'][ite][i] = zeros((df_hidr['TRAVELTIME'][i],NT+1))
        mu['vert'][ite][i] = zeros((df_hidr['TRAVELTIME'][i],NT+1))

# # # # #  Backward pass
for blk in range(nblocks.shape[0]):
    period = int(period_vector[-blk-1])
    m_b[period] = UC_model_1block(int(nblocks[-blk-1]),svars,ite,period,0)

    if blk > 0:
        f1,f2,f3,f4,f5,f6 = 0,0,0,0,0,0
        for i in range(NG):
            f1 += mu['tu'][ite][period_vector[-blk]][i]*(m_f[period]._tu[i,nblocks[-blk-1]-1] - svars['tu'][ite][i,period_vector[nblocks.shape[0]-blk]-1])

            if df_term['UPTIME'][i] > 1:
                c1,c2 = 0,1
                for l in range(df_term['UPTIME'][i]-1,0,-1):
                    if period_vector[nblocks.shape[0]-blk] - c1 - 1 <= 0:
                        break
                    if c1 < nblocks[-blk-1]:
                        f2 += mu['tv'][ite][i][l-1,period_vector[-blk]]*(m_f[period]._tv[i,nblocks[-blk-1] - c1 - 1] - svars['tv'][ite][i,period_vector[nblocks.shape[0]-blk] - c1 - 1])
                    else:
                        f2 += mu['tv'][ite][i][l-1,period_vector[-blk]]*(m_f[period]._aux_tv[i,period][df_term['UPTIME'][i]-1-c2] - svars['tv'][ite][i,period_vector[nblocks.shape[0]-blk] - c1 - 1])
                        c2+=1
                    c1 += 1

            if df_term['DOWNTIME'][i] > 1:
                c1,c2 = 0,1
                for l in range(df_term['DOWNTIME'][i]-1,0,-1):
                    if period_vector[nblocks.shape[0]-blk]-c1-1 <= 0:
                        break
                    if c1 < nblocks[-blk-1]:
                        f3 += mu['tw'][ite][i][l-1,period_vector[-blk]]*(m_f[period]._tw[i,nblocks[-blk-1]-c1-1] - svars['tw'][ite][i,period_vector[nblocks.shape[0]-blk]-c1-1])
                    else:
                        f3 += mu['tw'][ite][i][l-1,period_vector[-blk]]*(m_f[period]._aux_tw[i,period][df_term['DOWNTIME'][i]-1-c2] - svars['tw'][ite][i,period_vector[nblocks.shape[0]-blk]-c1-1])
                        c2+=1
                    c1 += 1

            f4 += mu['gt'][ite][period_vector[-blk]][i]*(m_f[period]._gt[i,nblocks[-blk-1]-1] - svars['gt'][ite][i,period_vector[nblocks.shape[0]-blk]-1])

        for i in range(NH):
            f5 += mu['vol'][ite][period_vector[-blk]][i]*(m_f[period]._vol[i,nblocks[-blk-1]-1] - svars['vol'][ite][i,period_vector[nblocks.shape[0]-blk]-1])

            if df_hidr['TRAVELTIME'][i] > 0:
                c1,c2 = 0,1
                for l in range(df_hidr['TRAVELTIME'][i],0,-1):
                    if period_vector[nblocks.shape[0]-blk]-c1-1 <= 0:
                        break
                    if c1 < nblocks[-blk-1]:
                        f6 += mu['turb'][ite][i][l-1,period_vector[-blk]]*(m_f[period]._turb[i,nblocks[-blk-1]-c1-1] - svars['turb'][ite][i,period_vector[nblocks.shape[0]-blk]-c1-1])
                        f6 += mu['vert'][ite][i][l-1,period_vector[-blk]]*(m_f[period]._vert[i,nblocks[-blk-1]-c1-1] - svars['vert'][ite][i,period_vector[nblocks.shape[0]-blk]-c1-1])
                    else:
                        f6 += mu['turb'][ite][i][l-1,period_vector[-blk]]*(m_f[period]._aux_turb[i,period][df_hidr['TRAVELTIME'][i]-1-c2] - svars['turb'][ite][i,period_vector[nblocks.shape[0]-blk]-c1-1])
                        f6 += mu['vert'][ite][i][l-1,period_vector[-blk]]*(m_f[period]._aux_vert[i,period][df_hidr['TRAVELTIME'][i]-1-c2] - svars['vert'][ite][i,period_vector[nblocks.shape[0]-blk]-c1-1])
                        c2 += 1
                    c1 += 1

        m_f[period].addConstr(m_f[period]._Theta[0] >= obj[ite,period_vector[-blk]] + f1 + f2 + f3 + f4 + f5 + f6) # Cut added to forward problem

        f1,f2,f3,f4,f5,f6 = 0,0,0,0,0,0
        for i in range(NG):
            f1 += mu['tu'][ite][period_vector[-blk]][i]*(m_b[period]._tu[i,nblocks[-blk-1]-1] - svars['tu'][ite][i,period_vector[nblocks.shape[0]-blk]-1])

            if df_term['UPTIME'][i] > 1:
                c1,c2 = 0,1
                for l in range(df_term['UPTIME'][i]-1,0,-1):
                    if period_vector[nblocks.shape[0]-blk] - c1 - 1 <= 0:
                        break
                    if c1 < nblocks[-blk-1]:
                        f2 += mu['tv'][ite][i][l-1,period_vector[-blk]]*(m_b[period]._tv[i,nblocks[-blk-1] - c1 - 1] - svars['tv'][ite][i,period_vector[nblocks.shape[0]-blk] - c1 - 1])
                    else:
                        f2 += mu['tv'][ite][i][l-1,period_vector[-blk]]*(m_b[period]._aux_tv[i,period][df_term['UPTIME'][i]-1-c2] - svars['tv'][ite][i,period_vector[nblocks.shape[0]-blk] - c1 - 1])
                        c2+=1
                    c1 += 1

            if df_term['DOWNTIME'][i] > 1:
                c1,c2 = 0,1
                for l in range(df_term['DOWNTIME'][i]-1,0,-1):
                    if period_vector[nblocks.shape[0]-blk]-c1-1 <= 0:
                        break
                    if c1 < nblocks[-blk-1]:
                        f3 += mu['tw'][ite][i][l-1,period_vector[-blk]]*(m_b[period]._tw[i,nblocks[-blk-1]-c1-1] - svars['tw'][ite][i,period_vector[nblocks.shape[0]-blk]-c1-1])
                    else:
                        f3 += mu['tw'][ite][i][l-1,period_vector[-blk]]*(m_b[period]._aux_tw[i,period][df_term['DOWNTIME'][i]-1-c2] - svars['tw'][ite][i,period_vector[nblocks.shape[0]-blk]-c1-1])
                        c2+=1
                    c1 += 1

            f4 += mu['gt'][ite][period_vector[-blk]][i]*(m_b[period]._gt[i,nblocks[-blk-1]-1] - svars['gt'][ite][i,period_vector[nblocks.shape[0]-blk]-1])

        for i in range(NH):
            f5 += mu['vol'][ite][period_vector[-blk]][i]*(m_b[period]._vol[i,nblocks[-blk-1]-1] - svars['vol'][ite][i,period_vector[nblocks.shape[0]-blk]-1])

            if df_hidr['TRAVELTIME'][i] > 0:
                c1,c2 = 0,1
                for l in range(df_hidr['TRAVELTIME'][i],0,-1):
                    if period_vector[nblocks.shape[0]-blk]-c1-1 <= 0:
                        break
                    if c1 < nblocks[-blk-1]:
                        f6 += mu['turb'][ite][i][l-1,period_vector[-blk]]*(m_b[period]._turb[i,nblocks[-blk-1]-c1-1] - svars['turb'][ite][i,period_vector[nblocks.shape[0]-blk]-c1-1])
                        f6 += mu['vert'][ite][i][l-1,period_vector[-blk]]*(m_b[period]._vert[i,nblocks[-blk-1]-c1-1] - svars['vert'][ite][i,period_vector[nblocks.shape[0]-blk]-c1-1])
                    else:
                        f6 += mu['turb'][ite][i][l-1,period_vector[-blk]]*(m_b[period]._aux_turb[i,period][df_hidr['TRAVELTIME'][i]-1-c2] - svars['turb'][ite][i,period_vector[nblocks.shape[0]-blk]-c1-1])
                        f6 += mu['vert'][ite][i][l-1,period_vector[-blk]]*(m_b[period]._aux_vert[i,period][df_hidr['TRAVELTIME'][i]-1-c2] - svars['vert'][ite][i,period_vector[nblocks.shape[0]-blk]-c1-1])
                        c2 += 1
                    c1 += 1


        m_b[period].addConstr(m_b[period]._Theta[0] >= obj[ite,period_vector[-blk]] + f1 + f2 + f3 + f4 + f5 + f6) # Cut added to backward problem


    et = time.time()
    m_b[period].optimize();
    tet += time.time() - et
    m_b[period].write('model_py_%d.lp'%(period))

    if period > 1:
        mu['tu'][ite][period] = array([m_b[period].getConstrByName("c_tu[%d]"%i).pi for i in range(NG)])
        mu['gt'][ite][period] = array([m_b[period].getConstrByName("c_gt[%d]"%i).pi for i in range(NG)])
        mu['vol'][ite][period] = array([m_b[period].getConstrByName("c_hvol[%d]"%i).pi for i in range(NH)])

        for i in range(NG):
            if df_term['UPTIME'][i] > 1:
                for l in range(df_term['UPTIME'][i]-1,0,-1):
                    mu['tv'][ite][i][l-1,period] = m_b[period].getConstrByName("c_tv[%d,%d]"%(i,l-1)).pi

        for i in range(NG):
            if df_term['DOWNTIME'][i] > 1:
                for l in range(df_term['DOWNTIME'][i]-1,0,-1):
                    mu['tw'][ite][i][l-1,period] = m_b[period].getConstrByName("c_tw[%d,%d]"%(i,l-1)).pi

        for i in range(NH):
            if df_hidr['TRAVELTIME'][i] > 0:
                for l in range(df_hidr['TRAVELTIME'][i],0,-1):
                    mu['turb'][ite][i][l-1,period] = m_b[period].getConstrByName("c_turb[%d,%d]"%(i,l-1)).pi
                    mu['vert'][ite][i][l-1,period] = m_b[period].getConstrByName("c_vert[%d,%d]"%(i,l-1)).pi

    obj[ite,period] = m_b[period].objval
# # # # # Backward pass end


target = 513950.03
print("Target: %.4f" %(513950.03))

# LB = m_b[1].objval
# LB_Vector = append(LB_Vector,LB)
# gap = abs(UB-LB)/UB*100; # Percentage gap
# gap_Vector = append(gap_Vector,gap)
# best_gap = abs(best_UB-LB)/best_UB*100; # Percentage gap
# best_gap_Vector = append(best_gap_Vector,best_gap)
# print("Iteration: %d, Current UB = %.2f, Best UB = %.2f, Current LB = %.2f, cur_gap = %.2f, best_gap = %.2f" %(ite,UB,best_UB,LB,round(gap,3),round(best_gap,3)))


# Backward sweep for first iteration ended

maxiter = 50; # Maximum number of DDIP iterations
tolerance = 1e-1; # Percentage tolerance stopping criteria

best_gap = 100
best_gap_Vector = append(best_gap_Vector,best_gap)
#################  Iterations 2,3,... start  ###################
while best_gap >= tolerance and ite < maxiter:
    ite += 1; # Iteration number

    # if ite > 5:
    #     if abs(best_gap_Vector[ite-1] - best_gap_Vector[ite-5]) < tolerance:
    #         break


    svars['vol'][ite] = zeros((NH,NT+1))
    svars['turb'][ite] = zeros((NH,NT+1))
    svars['vert'][ite] = zeros((NH,NT+1))
    svars['tu'][ite] = zeros((NG,NT+1))
    svars['tv'][ite] = zeros((NG,NT+1))
    svars['tw'][ite] = zeros((NG,NT+1))
    svars['gt'][ite] = zeros((NG,NT+1))
    svars['vol'][ite][:,0] = 0.01*df_hidr['V0']*(df_hidr['VMAX']-df_hidr['VMIN']) + df_hidr['VMIN']
    svars['turb'][ite][:,0] = df_hidr['Q0'].to_numpy().astype(float)
    svars['vert'][ite][:,0] = df_hidr['S0'].to_numpy().astype(float)

    for i in range(NG):
        if df_term['TON'][i] > 0:
            svars['tu'][ite][i][0] = 1
            if df_term['TON'][i] == 1:
                svars['tv'][ite][i][0] = 1
                svars['tw'][ite][i][0] = 0
            else:
                svars['tv'][ite][i][0] = 0
                svars['tw'][ite][i][0] = 0
        else:
            if df_term['TON'][i][0] == 1:
                svars['tv'][ite][i][0] = 0
                svars['tw'][ite][i][0] = 1
            else:
                svars['tv'][ite][i][0] = 0
                svars['tw'][ite][i][0] = 0
    svars['gt'][ite][:,0] = df_term['P0'].to_numpy()

    for blk in range(nblocks.shape[0]):
        period = period_vector[blk]

        if period > 1:
            for i in range(NG):
                m_f[period].setAttr("RHS", m_f[period].getConstrByName("c_tu[%d]"%i), svars['tu'][ite][i,period-1])

                if df_term['UPTIME'][i] > 1:
                    for l in range(df_term['UPTIME'][i]-1,0,-1):
                        m_f[period].setAttr("RHS", m_f[period].getConstrByName("c_tv[%d,%d]"%(i,l-1)), (svars['tv'][ite][i,l-df_term['UPTIME'][i]+period] if l-df_term['UPTIME'][i]+period >= 0 else 0))
                if df_term['DOWNTIME'][i] > 1:
                    for l in range(df_term['DOWNTIME'][i]-1,0,-1):
                        m_f[period].setAttr("RHS", m_f[period].getConstrByName("c_tw[%d,%d]"%(i,l-1)), (svars['tw'][ite][i,l-df_term['DOWNTIME'][i]+period] if l-df_term['DOWNTIME'][i]+period >= 0 else 0))
                m_f[period].setAttr("RHS", m_f[period].getConstrByName("c_gt[%d]"%i), svars['gt'][ite][i,period-1])
            for i in range(NH):
                m_f[period].setAttr("RHS", m_f[period].getConstrByName("c_hvol[%d]"%i), svars['vol'][ite][i,period-1])

                for l in range(int(df_hidr['TRAVELTIME'][i]),0,-1):
                    m_f[period].setAttr("RHS", m_f[period].getConstrByName("c_turb[%d,%d]"%(i,l-1)), (svars['turb'][ite][i,l-df_hidr['TRAVELTIME'][i]+period-1] if l-df_hidr['TRAVELTIME'][i]+period-1 >= 0 else 0))
                    m_f[period].setAttr("RHS", m_f[period].getConstrByName("c_vert[%d,%d]"%(i,l-1)), (svars['vert'][ite][i,l-df_hidr['TRAVELTIME'][i]+period-1] if l-df_hidr['TRAVELTIME'][i]+period-1 >= 0 else 0))


        et = time.time()
        m_f[period].optimize();
        tet += time.time() - et

        # m_f[period].write('model_py_%d.lp'%(period))

        for j in range(nblocks[blk]):
            svars['tu'][ite][:,period+j] = array([m_f[period]._tu[i,j].x[0] for i in range(NG)])
            svars['tv'][ite][:,period+j] = array([m_f[period]._tv[i,j].x[0] for i in range(NG)])
            svars['tw'][ite][:,period+j] = array([m_f[period]._tw[i,j].x[0] for i in range(NG)])
            svars['gt'][ite][:,period+j] = array([m_f[period]._gt[i,j].x[0] for i in range(NG)])
            svars['vol'][ite][:,period+j] = array([m_f[period]._vol[i,j].x[0] for i in range(NH)])
            svars['turb'][ite][:,period+j] = array([m_f[period]._turb[i,j].x[0] for i in range(NH)])
            svars['vert'][ite][:,period+j] = array([m_f[period]._vert[i,j].x[0] for i in range(NH)])

        Cost_f[ite,period] = m_f[period].objval - m_f[period]._Theta[0].x



    UB = sum(Cost_f[ite,p] for p in period_vector)
    UB_Vector = append(UB_Vector,UB)
    best_UB = min(UB_Vector);
    best_UB_Vector = append(best_UB_Vector, best_UB)



    # Forward pass for iteration j ended
    # Backward pass starting for iteration j


    mu['vol'][ite],mu['tu'][ite],mu['tv'][ite],mu['tw'][ite],mu['gt'][ite] = {}, {}, {}, {}, {}
    mu['turb'][ite],mu['vert'][ite] = {},{}
    obj[ite,NT+1] = 0; mu['tu'][ite][NT+1] = zeros(NG);
    mu['gt'][ite][NT+1] = zeros(NG);  mu['vol'][ite][NT+1] = zeros(NH);
    for i in range(NG):
        if df_term['UPTIME'][i] > 1:
            mu['tv'][ite][i] = zeros((df_term['UPTIME'][i],NT+1))
    for i in range(NG):
        if df_term['DOWNTIME'][i] > 1:
            mu['tw'][ite][i] = zeros((df_term['DOWNTIME'][i],NT+1))

    for i in range(NH):
        if df_hidr['TRAVELTIME'][i] > 0:
            mu['turb'][ite][i] = zeros((df_hidr['TRAVELTIME'][i],NT+1))
            mu['vert'][ite][i] = zeros((df_hidr['TRAVELTIME'][i],NT+1))


    for blk in range(nblocks.shape[0]):
        period = int(period_vector[-blk-1])
        if period > 1:
            for i in range(NG):
                m_b[period].setAttr("RHS", m_b[period].getConstrByName("c_tu[%d]"%i), svars['tu'][ite][i,period-1])

                if df_term['UPTIME'][i] > 1:
                    for l in range(df_term['UPTIME'][i]-1,0,-1):
                        m_b[period].setAttr("RHS", m_b[period].getConstrByName("c_tv[%d,%d]"%(i,l-1)), (svars['tv'][ite][i,l-df_term['UPTIME'][i]+period] if l-df_term['UPTIME'][i]+period >= 0 else 0))
                if df_term['DOWNTIME'][i] > 1:
                    for l in range(df_term['DOWNTIME'][i]-1,0,-1):
                        m_b[period].setAttr("RHS", m_b[period].getConstrByName("c_tw[%d,%d]"%(i,l-1)), (svars['tw'][ite][i,l-df_term['DOWNTIME'][i]+period] if l-df_term['DOWNTIME'][i]+period >= 0 else 0))
                m_b[period].setAttr("RHS", m_b[period].getConstrByName("c_gt[%d]"%i), svars['gt'][ite][i,period-1])
            for i in range(NH):
                m_b[period].setAttr("RHS", m_b[period].getConstrByName("c_hvol[%d]"%i), svars['vol'][ite][i,period-1])
                for l in range(int(df_hidr['TRAVELTIME'][i]),0,-1):
                    m_b[period].setAttr("RHS", m_b[period].getConstrByName("c_turb[%d,%d]"%(i,l-1)), (svars['turb'][ite][i,l-df_hidr['TRAVELTIME'][i]+period-1] if l-df_hidr['TRAVELTIME'][i]+period-1 >= 0 else 0))
                    m_b[period].setAttr("RHS", m_b[period].getConstrByName("c_vert[%d,%d]"%(i,l-1)), (svars['vert'][ite][i,l-df_hidr['TRAVELTIME'][i]+period-1] if l-df_hidr['TRAVELTIME'][i]+period-1 >= 0 else 0))

        # # # # # begin insert cuts
        if blk > 0:
            f1,f2,f3,f4,f5,f6 = 0,0,0,0,0,0
            for i in range(NG):
                f1 += mu['tu'][ite][period_vector[-blk]][i]*(m_f[period]._tu[i,nblocks[-blk-1]-1] - svars['tu'][ite][i,period_vector[nblocks.shape[0]-blk]-1])

                if df_term['UPTIME'][i] > 1:
                    c1,c2 = 0,1
                    for l in range(df_term['UPTIME'][i]-1,0,-1):
                        if period_vector[nblocks.shape[0]-blk] - c1 - 1 <= 0:
                            break
                        if c1 < nblocks[-blk-1]:
                            f2 += mu['tv'][ite][i][l-1,period_vector[-blk]]*(m_f[period]._tv[i,nblocks[-blk-1] - c1 - 1] - svars['tv'][ite][i,period_vector[nblocks.shape[0]-blk] - c1 - 1])
                        else:
                            f2 += mu['tv'][ite][i][l-1,period_vector[-blk]]*(m_f[period]._aux_tv[i,period][df_term['UPTIME'][i]-1-c2] - svars['tv'][ite][i,period_vector[nblocks.shape[0]-blk] - c1 - 1])
                            c2+=1
                        c1 += 1

                if df_term['DOWNTIME'][i] > 1:
                    c1,c2 = 0,1
                    for l in range(df_term['DOWNTIME'][i]-1,0,-1):
                        if period_vector[nblocks.shape[0]-blk]-c1-1 <= 0:
                            break
                        if c1 < nblocks[-blk-1]:
                            f3 += mu['tw'][ite][i][l-1,period_vector[-blk]]*(m_f[period]._tw[i,nblocks[-blk-1]-c1-1] - svars['tw'][ite][i,period_vector[nblocks.shape[0]-blk]-c1-1])
                        else:
                            f3 += mu['tw'][ite][i][l-1,period_vector[-blk]]*(m_f[period]._aux_tw[i,period][df_term['DOWNTIME'][i]-1-c2] - svars['tw'][ite][i,period_vector[nblocks.shape[0]-blk]-c1-1])
                            c2+=1
                        c1 += 1

                f4 += mu['gt'][ite][period_vector[-blk]][i]*(m_f[period]._gt[i,nblocks[-blk-1]-1] - svars['gt'][ite][i,period_vector[nblocks.shape[0]-blk]-1])

            for i in range(NH):
                f5 += mu['vol'][ite][period_vector[-blk]][i]*(m_f[period]._vol[i,nblocks[-blk-1]-1] - svars['vol'][ite][i,period_vector[nblocks.shape[0]-blk]-1])

                if df_hidr['TRAVELTIME'][i] > 0:
                    c1,c2 = 0,1
                    for l in range(df_hidr['TRAVELTIME'][i],0,-1):
                        if period_vector[nblocks.shape[0]-blk]-c1-1 <= 0:
                            break
                        if c1 < nblocks[-blk-1]:
                            f6 += mu['turb'][ite][i][l-1,period_vector[-blk]]*(m_f[period]._turb[i,nblocks[-blk-1]-c1-1] - svars['turb'][ite][i,period_vector[nblocks.shape[0]-blk]-c1-1])
                            f6 += mu['vert'][ite][i][l-1,period_vector[-blk]]*(m_f[period]._vert[i,nblocks[-blk-1]-c1-1] - svars['vert'][ite][i,period_vector[nblocks.shape[0]-blk]-c1-1])
                        else:
                            f6 += mu['turb'][ite][i][l-1,period_vector[-blk]]*(m_f[period]._aux_turb[i,period][df_hidr['TRAVELTIME'][i]-1-c2] - svars['turb'][ite][i,period_vector[nblocks.shape[0]-blk]-c1-1])
                            f6 += mu['vert'][ite][i][l-1,period_vector[-blk]]*(m_f[period]._aux_vert[i,period][df_hidr['TRAVELTIME'][i]-1-c2] - svars['vert'][ite][i,period_vector[nblocks.shape[0]-blk]-c1-1])
                            c2 += 1
                        c1 += 1

            m_f[period].addConstr(m_f[period]._Theta[0] >= obj[ite,period_vector[-blk]] + f1 + f2 + f3 + f4 + f5 + f6) # Cut added to forward problem

            f1,f2,f3,f4,f5,f6 = 0,0,0,0,0,0
            for i in range(NG):
                f1 += mu['tu'][ite][period_vector[-blk]][i]*(m_b[period]._tu[i,nblocks[-blk-1]-1] - svars['tu'][ite][i,period_vector[nblocks.shape[0]-blk]-1])

                if df_term['UPTIME'][i] > 1:
                    c1,c2 = 0,1
                    for l in range(df_term['UPTIME'][i]-1,0,-1):
                        if period_vector[nblocks.shape[0]-blk] - c1 - 1 <= 0:
                            break
                        if c1 < nblocks[-blk-1]:
                            f2 += mu['tv'][ite][i][l-1,period_vector[-blk]]*(m_b[period]._tv[i,nblocks[-blk-1] - c1 - 1] - svars['tv'][ite][i,period_vector[nblocks.shape[0]-blk] - c1 - 1])
                        else:
                            f2 += mu['tv'][ite][i][l-1,period_vector[-blk]]*(m_b[period]._aux_tv[i,period][df_term['UPTIME'][i]-1-c2] - svars['tv'][ite][i,period_vector[nblocks.shape[0]-blk] - c1 - 1])
                            c2+=1
                        c1 += 1

                if df_term['DOWNTIME'][i] > 1:
                    c1,c2 = 0,1
                    for l in range(df_term['DOWNTIME'][i]-1,0,-1):
                        if period_vector[nblocks.shape[0]-blk]-c1-1 <= 0:
                            break
                        if c1 < nblocks[-blk-1]:
                            f3 += mu['tw'][ite][i][l-1,period_vector[-blk]]*(m_b[period]._tw[i,nblocks[-blk-1]-c1-1] - svars['tw'][ite][i,period_vector[nblocks.shape[0]-blk]-c1-1])
                        else:
                            f3 += mu['tw'][ite][i][l-1,period_vector[-blk]]*(m_b[period]._aux_tw[i,period][df_term['DOWNTIME'][i]-1-c2] - svars['tw'][ite][i,period_vector[nblocks.shape[0]-blk]-c1-1])
                            c2+=1
                        c1 += 1

                f4 += mu['gt'][ite][period_vector[-blk]][i]*(m_b[period]._gt[i,nblocks[-blk-1]-1] - svars['gt'][ite][i,period_vector[nblocks.shape[0]-blk]-1])

            for i in range(NH):
                f5 += mu['vol'][ite][period_vector[-blk]][i]*(m_b[period]._vol[i,nblocks[-blk-1]-1] - svars['vol'][ite][i,period_vector[nblocks.shape[0]-blk]-1])

                if df_hidr['TRAVELTIME'][i] > 0:
                    c1,c2 = 0,1
                    for l in range(df_hidr['TRAVELTIME'][i],0,-1):
                        if period_vector[nblocks.shape[0]-blk]-c1-1 <= 0:
                            break
                        if c1 < nblocks[-blk-1]:
                            f6 += mu['turb'][ite][i][l-1,period_vector[-blk]]*(m_b[period]._turb[i,nblocks[-blk-1]-c1-1] - svars['turb'][ite][i,period_vector[nblocks.shape[0]-blk]-c1-1])
                            f6 += mu['vert'][ite][i][l-1,period_vector[-blk]]*(m_b[period]._vert[i,nblocks[-blk-1]-c1-1] - svars['vert'][ite][i,period_vector[nblocks.shape[0]-blk]-c1-1])
                        else:
                            f6 += mu['turb'][ite][i][l-1,period_vector[-blk]]*(m_b[period]._aux_turb[i,period][df_hidr['TRAVELTIME'][i]-1-c2] - svars['turb'][ite][i,period_vector[nblocks.shape[0]-blk]-c1-1])
                            f6 += mu['vert'][ite][i][l-1,period_vector[-blk]]*(m_b[period]._aux_vert[i,period][df_hidr['TRAVELTIME'][i]-1-c2] - svars['vert'][ite][i,period_vector[nblocks.shape[0]-blk]-c1-1])
                            c2 += 1
                        c1 += 1


            m_b[period].addConstr(m_b[period]._Theta[0] >= obj[ite,period_vector[-blk]] + f1 + f2 + f3 + f4 + f5 + f6) # Cut added to backward problem
        # # # # # end insert cuts

        et = time.time()
        m_b[period].optimize();
        tet += time.time() - et

        # m_b[period].write('model_py_%d.lp'%(period))

        if period > 1:
            mu['tu'][ite][period] = array([m_b[period].getConstrByName("c_tu[%d]"%i).pi for i in range(NG)])
            mu['gt'][ite][period] = array([m_b[period].getConstrByName("c_gt[%d]"%i).pi for i in range(NG)])
            mu['vol'][ite][period] = array([m_b[period].getConstrByName("c_hvol[%d]"%i).pi for i in range(NH)])

            for i in range(NG):
                if df_term['UPTIME'][i] > 1:
                    for l in range(df_term['UPTIME'][i]-1,0,-1):
                        mu['tv'][ite][i][l-1,period] = m_b[period].getConstrByName("c_tv[%d,%d]"%(i,l-1)).pi

            for i in range(NG):
                if df_term['DOWNTIME'][i] > 1:
                    for l in range(df_term['DOWNTIME'][i]-1,0,-1):
                        mu['tw'][ite][i][l-1,period] = m_b[period].getConstrByName("c_tw[%d,%d]"%(i,l-1)).pi

            for i in range(NH):
                if df_hidr['TRAVELTIME'][i] > 0:
                    for l in range(df_hidr['TRAVELTIME'][i],0,-1):
                        mu['turb'][ite][i][l-1,period] = m_b[period].getConstrByName("c_turb[%d,%d]"%(i,l-1)).pi
                        mu['vert'][ite][i][l-1,period] = m_b[period].getConstrByName("c_vert[%d,%d]"%(i,l-1)).pi

        obj[ite,period] = m_b[period].objval
    # Backward pass end


    LB = m_b[1].objval
    LB_Vector = append(LB_Vector, LB)
    gap = abs(UB - LB) / UB * 100;  # Percentage gap
    gap_Vector = append(gap_Vector, gap)
    best_gap = abs(best_UB - LB) / best_UB * 100;  # Percentage gap
    best_gap_Vector = append(best_gap_Vector, best_gap)
    print("Iteration: %d, Current UB = %.2f, Best UB = %.2f, Current LB = %.2f, cur_gap = %.2f, best_gap = %.2f" % (ite, UB, best_UB, LB, round(gap, 3), round(best_gap, 3)))

    # Backward pass for j iteration ended
     # Iteration j ends

print('total elapsed time: %f seconds'%tet)

import matplotlib.pyplot as plt

fig_title = 'MIPGAP: ' + str(1) + '%; Runtime: ' + str(round(tet,2)) + 's; CCB; Gap: ' + str(round(best_gap,2)) + '%'
# fig_title = 'Runtime: ' + str(round(tet,2)) + 's; CCB; Gap: ' + str(0) + '%'
t = linspace(0,LB_Vector.shape[0]-1,LB_Vector.shape[0])
plt.plot(t,best_UB_Vector,'g',label='UB')
plt.plot(t,LB_Vector,'r',label='LB')
plt.plot(t,ones(LB_Vector.shape[0])*513950.03,'k-.')
plt.xlabel('ite')
plt.ylabel('Cost')
plt.ylim([0,8*1e5])
plt.xlim([-1,LB_Vector.shape[0]])
plt.yticks([0,2e5,4e5,6e5,8e5])
plt.ticklabel_format(axis="y", style='sci', scilimits=(0,0))
# plt.legend()
plt.grid(True)
plt.title(fig_title)

figure = plt.gcf()
figure.set_size_inches(8, 6)

plt.show()

# plt.savefig('figures/fig61.png', dpi=100)
# plt.close()


# name_exp = "gap000001_CAB" #Cuts current backward
# savetxt(save_path+'/'+name_exp+'_bestUB', best_UB_Vector, fmt='%f')
# savetxt(save_path+'/'+name_exp+'_LB', LB_Vector, fmt='%f')
# savetxt(save_path+'/'+name_exp+'_UB', UB_Vector, fmt='%f')
# savetxt(save_path+'/'+name_exp+'_bestgap', best_gap_Vector, fmt='%f')


# def save_results():
#
#     name_exp = "gap01_CAB"
#     best_UB = loadtxt(save_path+'/'+name_exp+'_bestUB')
#     LB_Vector = loadtxt(save_path+'/'+name_exp+'_LB')
#     UB_Vector = loadtxt(save_path+'/'+name_exp+'_UB')
#     best_gap_Vector = loadtxt(save_path+'/'+name_exp+'_bestgap')
#
#     workbook = xlsxwriter.Workbook(save_path + '\\ddip_results1.xlsx')
#     worksheet = workbook.add_worksheet('exp_01_CAB')
#     worksheet.set_column('A:A', 20)
#
#     worksheet.write(1,1,"Best UB")
#     worksheet.write(1,2,"UB")
#     worksheet.write(1,3,"LB")
#     worksheet.write(1,4,"Best GAP")
#     for j in range(best_UB.shape[0]):
#         worksheet.write(j+2,1,best_UB[j])
#         worksheet.write(j+2,2,UB_Vector[j])
#         worksheet.write(j+2,3,LB_Vector[j])
#         worksheet.write(j+2,4,best_gap_Vector[j])
#
#     name_exp = "gap001_CAB"
#     best_UB = loadtxt(save_path+'/'+name_exp+'_bestUB')
#     LB_Vector = loadtxt(save_path+'/'+name_exp+'_LB')
#     UB_Vector = loadtxt(save_path+'/'+name_exp+'_UB')
#     best_gap_Vector = loadtxt(save_path+'/'+name_exp+'_bestgap')
#
#     worksheet = workbook.add_worksheet('exp_001_CAB')
#     worksheet.set_column('A:A', 20)
#
#     worksheet.write(1,1,"Best UB")
#     worksheet.write(1,2,"UB")
#     worksheet.write(1,3,"LB")
#     worksheet.write(1,4,"Best GAP")
#     for j in range(best_UB.shape[0]):
#         worksheet.write(j+2,1,best_UB[j])
#         worksheet.write(j+2,2,UB_Vector[j])
#         worksheet.write(j+2,3,LB_Vector[j])
#         worksheet.write(j+2,4,best_gap_Vector[j])
#
#     name_exp = "gap0001_CAB"
#     best_UB = loadtxt(save_path+'/'+name_exp+'_bestUB')
#     LB_Vector = loadtxt(save_path+'/'+name_exp+'_LB')
#     UB_Vector = loadtxt(save_path+'/'+name_exp+'_UB')
#     best_gap_Vector = loadtxt(save_path+'/'+name_exp+'_bestgap')
#
#     worksheet = workbook.add_worksheet('exp_0001_CAB')
#     worksheet.set_column('A:A', 20)
#
#     worksheet.write(1,1,"Best UB")
#     worksheet.write(1,2,"UB")
#     worksheet.write(1,3,"LB")
#     worksheet.write(1,4,"Best GAP")
#     for j in range(best_UB.shape[0]):
#         worksheet.write(j+2,1,best_UB[j])
#         worksheet.write(j+2,2,UB_Vector[j])
#         worksheet.write(j+2,3,LB_Vector[j])
#         worksheet.write(j+2,4,best_gap_Vector[j])
#
#     workbook.close()