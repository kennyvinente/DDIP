#### DDIP HTUC with no block aggregation

#target:

from ddip_formulations_blk import *

NT = 24
nblocks = 12

MIPGap = True

def UC_model_1block(nblocks,svars,ite,t0,MILP_LP):

    m = Model()
    m._NT = NT
    m.Params.OutputFlag = 0

    if MILP_LP:
        fph_model(m,nblocks,t0,1)
        term_model(m,nblocks,t0,1)
        if MIPGap:
            m.setParam('MIPGAP',0.1)
    else:
        fph_model(m,nblocks,t0,0)
        term_model(m,nblocks,t0,0)
    uct_model(m,svars,ite,nblocks,t0)
    network_model(m,nblocks,t0)
    # reserv_model(m,svars,ite,nblocks,t0,NT)
    f1 = fobj(m,nblocks)

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
obj = {}
mu = {}

# et = time.time()
tet = 0

period = 1
# # # # # Forward pass
for period in range(1,NT+1,nblocks):

    # m_f[period] = UC_model_1block(nblocks,svars,ite,period,1);
    m_f[period] = UC_model_1block(nblocks,svars,ite,period,0);
    et = time.time()
    m_f[period].optimize();
    tet += time.time() - et

    m_f[period].write('model_py_%d.lp'%(period))

    for j in range(nblocks):
        svars['tu'][ite][:,period+j] = array([m_f[period]._tu[i,j].x[0] for i in range(NG)])
        svars['tv'][ite][:,period+j] = array([m_f[period]._tv[i,j].x[0] for i in range(NG)])
        svars['tw'][ite][:,period+j] = array([m_f[period]._tw[i,j].x[0] for i in range(NG)])
        svars['gt'][ite][:,period+j] = array([m_f[period]._gt[i,j].x[0] for i in range(NG)])
        svars['vol'][ite][:,period+j] = array([m_f[period]._vol[i,j].x[0] for i in range(NH)])
        svars['turb'][ite][:,period+j] = array([m_f[period]._turb[i,j].x[0] for i in range(NH)])
        svars['vert'][ite][:,period+j] = array([m_f[period]._vert[i,j].x[0] for i in range(NH)])

    Cost_f[ite,period] = m_f[period].objval - m_f[period]._Theta[0].x
# # # # # End forward pass

UB = sum(Cost_f[ite,p] for p in range(1,NT+1,nblocks))
UB_Vector = append(UB_Vector,UB)
best_UB = min(UB_Vector);
best_UB_Vector = append(best_UB_Vector, best_UB)

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
for period in range(NT-nblocks+1,0,-nblocks):
    m_b[period] = UC_model_1block(nblocks,svars,ite,period,0)

    if 1 <= period <= NT-nblocks:
        f1,f2,f3,f4,f5,f6 = 0,0,0,0,0,0
        for i in range(NG):
            f1 += mu['tu'][ite][period+nblocks][i]*(m_f[period]._tu[i,nblocks-1] - svars['tu'][ite][i,period+nblocks-1])
            if df_term['UPTIME'][i] > 1:
                for l in range(df_term['UPTIME'][i]-1,0,-1):
                    if l - df_term['UPTIME'][i] + period + nblocks < 0:
                        break
                    if l == df_term['UPTIME'][i] - 1:
                        f2 += mu['tv'][ite][i][l-1,period+nblocks]*(m_f[period]._tv[i,period+nblocks -(df_term['UPTIME'][i]-l)-1] - svars['tv'][ite][i,period+nblocks-1])
                    else: #ver se, de acordo com o size do bloco, se é tv ou aux_tv
                        f2 += mu['tv'][ite][i][l-1,period+nblocks]*(m_f[period]._tv[i,period+nblocks -(df_term['UPTIME'][i]-l)-1] - svars['tv'][ite][i,period+nblocks -(df_term['UPTIME'][i]-l)])
                        # f2 += mu['tv'][ite][i][l-1,period+nblocks]*(m_f[period]._aux_tv[i,period][l] - svars['tv'][ite][i,period+nblocks -(df_term['UPTIME'][i]-l)-1])
        #     if df_term['DOWNTIME'][i] > 1:
        #         for l in range(df_term['DOWNTIME'][i]-1,0,-1):
        #             if l - df_term['DOWNTIME'][i] + period < 0:
        #                 break
        #             if l == df_term['DOWNTIME'][i] - 1:
        #                 f3 += mu['tw'][ite][i][l-1,period+1]*(m_f[period]._tw[i,0] - svars['tw'][ite][i,period])
        #             else:
        #                 f3 += mu['tw'][ite][i][l-1,period+1]*(m_f[period]._aux_tw[i,period][l] - svars['tw'][ite][i,l-df_term['DOWNTIME'][i]+period+1])
            f4 += mu['gt'][ite][period+nblocks][i]*(m_f[period]._gt[i,nblocks-1] - svars['gt'][ite][i,period+nblocks-1])
        # for i in range(NH):
        #     f5 += mu['vol'][ite][period+1][i]*(m_f[period]._vol[i,0] - svars['vol'][ite][i,period])
        #     if df_hidr['TRAVELTIME'][i] > 0:
        #         for l in range(df_hidr['TRAVELTIME'][i],0,-1):
        #             if l-df_hidr['TRAVELTIME'][i]+period-1 < 0:
        #                 break
        #             if l == df_hidr['TRAVELTIME'][i]:
        #                 f6 += mu['turb'][ite][i][l-1,period+1]*(m_f[period]._turb[i,0] - svars['turb'][ite][i,period])
        #                 f6 += mu['vert'][ite][i][l-1,period+1]*(m_f[period]._vert[i,0] - svars['vert'][ite][i,period])
        #             else:
        #                 f6 += mu['turb'][ite][i][l-1,period+1]*(m_f[period]._aux_turb[i,period][l] - svars['turb'][ite][i,l-df_hidr['TRAVELTIME'][i]+period])
        #                 f6 += mu['vert'][ite][i][l-1,period+1]*(m_f[period]._aux_vert[i,period][l] - svars['vert'][ite][i,l-df_hidr['TRAVELTIME'][i]+period])
        m_f[period].addConstr(m_f[period]._Theta[0] >= obj[ite,period+nblocks] + f1 + f2 + f3 + f4 + f5 + f6) # Cut added to forward problem

        f1,f2,f3,f4,f5,f6 = 0,0,0,0,0,0
        for i in range(NG):
            f1 += mu['tu'][ite][period+nblocks][i]*(m_b[period]._tu[i,nblocks-1] - svars['tu'][ite][i,period+nblocks-1])
            if df_term['UPTIME'][i] > 1:
                for l in range(df_term['UPTIME'][i]-1,0,-1):
                    if l - df_term['UPTIME'][i] + period + nblocks < 0:
                        break
                    if l == df_term['UPTIME'][i] - 1:
                        f2 += mu['tv'][ite][i][l-1,period+nblocks]*(m_b[period]._tv[i,period+nblocks -(df_term['UPTIME'][i]-l)-1] - svars['tv'][ite][i,period+nblocks-1])
                    else: #ver se, de acordo com o size do bloco, se é tv ou aux_tv
                        f2 += mu['tv'][ite][i][l-1,period+nblocks]*(m_b[period]._tv[i,period+nblocks -(df_term['UPTIME'][i]-l)-1] - svars['tv'][ite][i,period+nblocks -(df_term['UPTIME'][i]-l)])
                        # f2 += mu['tv'][ite][i][l-1,period+nblocks]*(m_b[period]._aux_tv[i,period][l] - svars['tv'][ite][i,period+nblocks -(df_term['UPTIME'][i]-l)-1])
        #     if df_term['DOWNTIME'][i] > 1:
        #         for l in range(df_term['DOWNTIME'][i]-1,0,-1):
        #             if l - df_term['DOWNTIME'][i] + period < 0:
        #                 break
        #             if l == df_term['DOWNTIME'][i] - 1:
        #                 f3 += mu['tw'][ite][i][l-1,period+1]*(m_b[period]._tw[i,0] - svars['tw'][ite][i,period])
        #             else:
        #                 f3 += mu['tw'][ite][i][l-1,period+1]*(m_b[period]._aux_tw[i,period][l] - svars['tw'][ite][i,l-df_term['DOWNTIME'][i]+period+1])
            f4 += mu['gt'][ite][period+nblocks][i]*(m_b[period]._gt[i,nblocks-1] - svars['gt'][ite][i,period+nblocks-1])
        # for i in range(NH):
        #     f5 += mu['vol'][ite][period+1][i]*(m_b[period]._vol[i,0] - svars['vol'][ite][i,period])
        #     if df_hidr['TRAVELTIME'][i] > 0:
        #         for l in range(df_hidr['TRAVELTIME'][i],0,-1):
        #             if l-df_hidr['TRAVELTIME'][i]+period-1 < 0:
        #                 break
        #             if l == df_hidr['TRAVELTIME'][i]:
        #                 f6 += mu['turb'][ite][i][l-1,period+1]*(m_b[period]._turb[i,0] - svars['turb'][ite][i,period])
        #                 f6 += mu['vert'][ite][i][l-1,period+1]*(m_b[period]._vert[i,0] - svars['vert'][ite][i,period])
        #             else:
        #                 f6 += mu['turb'][ite][i][l-1,period+1]*(m_b[period]._aux_turb[i,period][l] - svars['turb'][ite][i,l-df_hidr['TRAVELTIME'][i]+period])
        #                 f6 += mu['vert'][ite][i][l-1,period+1]*(m_b[period]._aux_vert[i,period][l] - svars['vert'][ite][i,l-df_hidr['TRAVELTIME'][i]+period])
        m_b[period].addConstr(m_b[period]._Theta[0] >= obj[ite,period+nblocks] + f1 + f2 + f3 + f4 + f5 + f6) # Cut added to backward problem


    et = time.time()
    m_b[period].optimize();
    tet += time.time() - et
    m_b[period].write('model_py_%d.lp'%(period))

    if period > 1:
        mu['tu'][ite][period] = array([m_b[period].getConstrByName("c_tu[%d]"%i).pi for i in range(NG)])
        mu['gt'][ite][period] = array([m_b[period].getConstrByName("c_gt[%d]"%i).pi for i in range(NG)])
    #     mu['vol'][ite][period] = array([m_b[period].getConstrByName("c_hvol[%d]"%i).pi for i in range(NH)])
    #
        for i in range(NG):
            if df_term['UPTIME'][i] > 1:
                for l in range(df_term['UPTIME'][i]-1,0,-1):
                    mu['tv'][ite][i][l-1,period] = m_b[period].getConstrByName("c_tv[%d,%d]"%(i,l-1)).pi
    #
    #     for i in range(NG):
    #         if df_term['DOWNTIME'][i] > 1:
    #             for l in range(df_term['DOWNTIME'][i]-1,0,-1):
    #                 mu['tw'][ite][i][l-1,period] = m_b[period].getConstrByName("c_tw[%d,%d]"%(i,l-1)).pi
    #
    #     for i in range(NH):
    #         if df_hidr['TRAVELTIME'][i] > 0:
    #             for l in range(df_hidr['TRAVELTIME'][i],0,-1):
    #                 mu['turb'][ite][i][l-1,period] = m_b[period].getConstrByName("c_turb[%d,%d]"%(i,l-1)).pi
    #                 mu['vert'][ite][i][l-1,period] = m_b[period].getConstrByName("c_vert[%d,%d]"%(i,l-1)).pi

    obj[ite,period] = m_b[period].objval
# # # # # Backward pass end

print("Target: %.4f" %(434373.558316281))

LB = m_b[1].objval
LB_Vector = append(LB_Vector,LB)
gap = abs(UB-LB)/UB*100; # Percentage gap
gap_Vector = append(gap_Vector,gap)
best_gap = abs(best_UB-LB)/best_UB*100; # Percentage gap
best_gap_Vector = append(best_gap_Vector,best_gap)
print("Iteration: %d, Current UB = %.2f, Best UB = %.2f, Current LB = %.2f, cur_gap = %.2f, best_gap = %.2f" %(ite,UB,best_UB,LB,round(gap,3),round(best_gap,3)))


# Backward sweep for first iteration ended

maxiter = 25; # Maximum number of DDIP iterations
tolerance = 1e-8; # Percentage tolerance stopping criteria


#################  Iterations 2,3,... start  ###################
while best_gap >= tolerance and ite < maxiter:
    ite += 1; # Iteration number

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

    # # # # # cuts added after each forward and backward
    # # # # # insert cuts here
    # for period in range(1,nblocks):
    #     f1,f2,f3,f4,f5,f6 = 0,0,0,0,0,0
    #     for i in range(NG):
    #         f1 += mu['tu'][ite-1][period+1][i]*(m_f[period]._tu[i,0] - svars['tu'][ite-1][i,period])
    #         if df_term['UPTIME'][i] > 1:
    #             for l in range(df_term['UPTIME'][i]-1,0,-1):
    #                 if l - df_term['UPTIME'][i] + period < 0:
    #                     break
    #                 if l == df_term['UPTIME'][i] - 1:
    #                     f2 += mu['tv'][ite-1][i][l-1,period+1]*(m_f[period]._tv[i,0] - svars['tv'][ite-1][i,period])
    #                 else:
    #                     f2 += mu['tv'][ite-1][i][l-1,period+1]*(m_f[period]._aux_tv[i,period][l] - svars['tv'][ite-1][i,l-df_term['UPTIME'][i]+period+1])
    #         if df_term['DOWNTIME'][i] > 1:
    #             for l in range(df_term['DOWNTIME'][i]-1,0,-1):
    #                 if l - df_term['DOWNTIME'][i] + period < 0:
    #                     break
    #                 if l == df_term['DOWNTIME'][i] - 1:
    #                     f3 += mu['tw'][ite-1][i][l-1,period+1]*(m_f[period]._tw[i,0] - svars['tw'][ite-1][i,period])
    #                 else:
    #                     f3 += mu['tw'][ite-1][i][l-1,period+1]*(m_f[period]._aux_tw[i,period][l] - svars['tw'][ite-1][i,l-df_term['DOWNTIME'][i]+period+1])
    #         f4 += mu['gt'][ite-1][period+1][i]*(m_f[period]._gt[i,0] - svars['gt'][ite-1][i,period])
    #     for i in range(NH):
    #         f5 += mu['vol'][ite-1][period+1][i]*(m_f[period]._vol[i,0] - svars['vol'][ite-1][i,period])
    #         if df_hidr['TRAVELTIME'][i] > 0:
    #             for l in range(df_hidr['TRAVELTIME'][i],0,-1):
    #                 if l-df_hidr['TRAVELTIME'][i]+period-1 < 0:
    #                     break
    #                 if l == df_hidr['TRAVELTIME'][i]:
    #                     f6 += mu['turb'][ite-1][i][l-1,period+1]*(m_f[period]._turb[i,0] - svars['turb'][ite-1][i,period])
    #                     f6 += mu['vert'][ite-1][i][l-1,period+1]*(m_f[period]._vert[i,0] - svars['vert'][ite-1][i,period])
    #                 else:
    #                     f6 += mu['turb'][ite-1][i][l-1,period+1]*(m_f[period]._aux_turb[i,period][l] - svars['turb'][ite-1][i,l-df_hidr['TRAVELTIME'][i]+period])
    #                     f6 += mu['vert'][ite-1][i][l-1,period+1]*(m_f[period]._aux_vert[i,period][l] - svars['vert'][ite-1][i,l-df_hidr['TRAVELTIME'][i]+period])
    #     m_f[period].addConstr(m_f[period]._Theta[0] >= obj[ite-1,period+1] + f1 + f2 + f3 + f4 + f5 + f6) # Cut added to forward problem
    #
    #     f1,f2,f3,f4,f5,f6 = 0,0,0,0,0,0
    #     for i in range(NG):
    #         f1 += mu['tu'][ite-1][period+1][i]*(m_b[period]._tu[i,0] - svars['tu'][ite-1][i,period])
    #         if df_term['UPTIME'][i] > 1:
    #             for l in range(df_term['UPTIME'][i]-1,0,-1):
    #                 if l - df_term['UPTIME'][i] + period < 0:
    #                     break
    #                 if l == df_term['UPTIME'][i] - 1:
    #                     f2 += mu['tv'][ite-1][i][l-1,period+1]*(m_b[period]._tv[i,0] - svars['tv'][ite-1][i,period])
    #                 else:
    #                     f2 += mu['tv'][ite-1][i][l-1,period+1]*(m_b[period]._aux_tv[i,period][l] - svars['tv'][ite-1][i,l-df_term['UPTIME'][i]+period+1])
    #         if df_term['DOWNTIME'][i] > 1:
    #             for l in range(df_term['DOWNTIME'][i]-1,0,-1):
    #                 if l - df_term['DOWNTIME'][i] + period < 0:
    #                     break
    #                 if l == df_term['DOWNTIME'][i] - 1:
    #                     f3 += mu['tw'][ite-1][i][l-1,period+1]*(m_b[period]._tw[i,0] - svars['tw'][ite-1][i,period])
    #                 else:
    #                     f3 += mu['tw'][ite-1][i][l-1,period+1]*(m_b[period]._aux_tw[i,period][l] - svars['tw'][ite-1][i,l-df_term['DOWNTIME'][i]+period+1])
    #         f4 += mu['gt'][ite-1][period+1][i]*(m_b[period]._gt[i,0] - svars['gt'][ite-1][i,period])
    #     for i in range(NH):
    #         f5 += mu['vol'][ite-1][period+1][i]*(m_b[period]._vol[i,0] - svars['vol'][ite-1][i,period])
    #         if df_hidr['TRAVELTIME'][i] > 0:
    #             for l in range(df_hidr['TRAVELTIME'][i],0,-1):
    #                 if l-df_hidr['TRAVELTIME'][i]+period-1 < 0:
    #                     break
    #                 if l == df_hidr['TRAVELTIME'][i]:
    #                     f6 += mu['turb'][ite-1][i][l-1,period+1]*(m_b[period]._turb[i,0] - svars['turb'][ite-1][i,period])
    #                     f6 += mu['vert'][ite-1][i][l-1,period+1]*(m_b[period]._vert[i,0] - svars['vert'][ite-1][i,period])
    #                 else:
    #                     f6 += mu['turb'][ite-1][i][l-1,period+1]*(m_b[period]._aux_turb[i,period][l] - svars['turb'][ite-1][i,l-df_hidr['TRAVELTIME'][i]+period])
    #                     f6 += mu['vert'][ite-1][i][l-1,period+1]*(m_b[period]._aux_vert[i,period][l] - svars['vert'][ite-1][i,l-df_hidr['TRAVELTIME'][i]+period])
    #     m_b[period].addConstr(m_b[period]._Theta[0] >= obj[ite-1,period+1] + f1 + f2 + f3 + f4 + f5 + f6) # Cut added to backward problem

    # # # # # end insert cuts

    for period in range(1,NT+1,nblocks):
        if period > 1:
            for i in range(NG): #### ajeitar
                m_f[period].setAttr("RHS", m_f[period].getConstrByName("c_tu[%d]"%i), svars['tu'][ite][i,period-1])

                if df_term['UPTIME'][i] > 1:
                    for l in range(df_term['UPTIME'][i]-1,0,-1):
                        m_f[period].setAttr("RHS", m_f[period].getConstrByName("c_tv[%d,%d]"%(i,l-1)), (svars['tv'][ite][i,l-df_term['UPTIME'][i]+period] if l-df_term['UPTIME'][i]+period >= 0 else 0))
            #     if df_term['DOWNTIME'][i] > 1:
            #         for l in range(df_term['DOWNTIME'][i]-1,0,-1):
            #             m_f[period].setAttr("RHS", m_f[period].getConstrByName("c_tw[%d,%d]"%(i,l-1)), (svars['tw'][ite][i,l-df_term['DOWNTIME'][i]+period] if l-df_term['DOWNTIME'][i]+period >= 0 else 0))
                m_f[period].setAttr("RHS", m_f[period].getConstrByName("c_gt[%d]"%i), svars['gt'][ite][i,period-1])
            # for i in range(NH):
            #     m_f[period].setAttr("RHS", m_f[period].getConstrByName("c_hvol[%d]"%i), svars['vol'][ite][i,period-1])
            #
            #     for l in range(int(df_hidr['TRAVELTIME'][i]),0,-1):
            #         m_f[period].setAttr("RHS", m_f[period].getConstrByName("c_turb[%d,%d]"%(i,l-1)), (svars['turb'][ite][i,l-df_hidr['TRAVELTIME'][i]+period-1] if l-df_hidr['TRAVELTIME'][i]+period-1 >= 0 else 0))
            #         m_f[period].setAttr("RHS", m_f[period].getConstrByName("c_vert[%d,%d]"%(i,l-1)), (svars['vert'][ite][i,l-df_hidr['TRAVELTIME'][i]+period-1] if l-df_hidr['TRAVELTIME'][i]+period-1 >= 0 else 0))


        et = time.time()
        m_f[period].optimize();
        tet += time.time() - et

        m_f[period].write('model_py_%d.lp'%(period))

        for j in range(nblocks):
            svars['tu'][ite][:,period+j] = array([m_f[period]._tu[i,j].x[0] for i in range(NG)])
            svars['tv'][ite][:,period+j] = array([m_f[period]._tv[i,j].x[0] for i in range(NG)])
            svars['tw'][ite][:,period+j] = array([m_f[period]._tw[i,j].x[0] for i in range(NG)])
            svars['gt'][ite][:,period+j] = array([m_f[period]._gt[i,j].x[0] for i in range(NG)])
            svars['vol'][ite][:,period+j] = array([m_f[period]._vol[i,j].x[0] for i in range(NH)])
            svars['turb'][ite][:,period+j] = array([m_f[period]._turb[i,j].x[0] for i in range(NH)])
            svars['vert'][ite][:,period+j] = array([m_f[period]._vert[i,j].x[0] for i in range(NH)])

        Cost_f[ite,period] = m_f[period].objval - m_f[period]._Theta[0].x



    UB = sum(Cost_f[ite,p] for p in range(1,NT+1,nblocks));
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


    for period in range(NT-nblocks+1,0,-nblocks):
        if period > 1:
            for i in range(NG): #### ajeitar
                m_b[period].setAttr("RHS", m_b[period].getConstrByName("c_tu[%d]"%i), svars['tu'][ite][i,period-1])

                if df_term['UPTIME'][i] > 1:
                    for l in range(df_term['UPTIME'][i]-1,0,-1):
                        m_b[period].setAttr("RHS", m_b[period].getConstrByName("c_tv[%d,%d]"%(i,l-1)), (svars['tv'][ite][i,l-df_term['UPTIME'][i]+period] if l-df_term['UPTIME'][i]+period >= 0 else 0))
            #     if df_term['DOWNTIME'][i] > 1:
            #         for l in range(df_term['DOWNTIME'][i]-1,0,-1):
            #             m_b[period].setAttr("RHS", m_b[period].getConstrByName("c_tw[%d,%d]"%(i,l-1)), (svars['tw'][ite][i,l-df_term['DOWNTIME'][i]+period] if l-df_term['DOWNTIME'][i]+period >= 0 else 0))
                m_b[period].setAttr("RHS", m_b[period].getConstrByName("c_gt[%d]"%i), svars['gt'][ite][i,period-1])
            # for i in range(NH):
            #     m_b[period].setAttr("RHS", m_b[period].getConstrByName("c_hvol[%d]"%i), svars['vol'][ite][i,period-1])
            #     for l in range(int(df_hidr['TRAVELTIME'][i]),0,-1):
            #         m_b[period].setAttr("RHS", m_b[period].getConstrByName("c_turb[%d,%d]"%(i,l-1)), (svars['turb'][ite][i,l-df_hidr['TRAVELTIME'][i]+period-1] if l-df_hidr['TRAVELTIME'][i]+period-1 >= 0 else 0))
            #         m_b[period].setAttr("RHS", m_b[period].getConstrByName("c_vert[%d,%d]"%(i,l-1)), (svars['vert'][ite][i,l-df_hidr['TRAVELTIME'][i]+period-1] if l-df_hidr['TRAVELTIME'][i]+period-1 >= 0 else 0))

        # # # # # begin insert cuts
        if 1 <= period <= NT-nblocks:
            f1,f2,f3,f4,f5,f6 = 0,0,0,0,0,0
            for i in range(NG):
                f1 += mu['tu'][ite][period+nblocks][i]*(m_f[period]._tu[i,nblocks-1] - svars['tu'][ite][i,period+nblocks-1])
                if df_term['UPTIME'][i] > 1:
                    for l in range(df_term['UPTIME'][i]-1,0,-1):
                        if l - df_term['UPTIME'][i] + period + nblocks < 0:
                            break
                        if l == df_term['UPTIME'][i] - 1:
                            f2 += mu['tv'][ite][i][l-1,period+nblocks]*(m_f[period]._tv[i,period+nblocks -(df_term['UPTIME'][i]-l)-1] - svars['tv'][ite][i,period+nblocks-1])
                        else: #ver se, de acordo com o size do bloco, se é tv ou aux_tv
                            f2 += mu['tv'][ite][i][l-1,period+nblocks]*(m_f[period]._tv[i,period+nblocks -(df_term['UPTIME'][i]-l)-1] - svars['tv'][ite][i,period+nblocks -(df_term['UPTIME'][i]-l)])
                            # f2 += mu['tv'][ite][i][l-1,period+nblocks]*(m_f[period]._aux_tv[i,period][l] - svars['tv'][ite][i,period+nblocks -(df_term['UPTIME'][i]-l)-1])
            #     if df_term['DOWNTIME'][i] > 1:
            #         for l in range(df_term['DOWNTIME'][i]-1,0,-1):
            #             if l - df_term['DOWNTIME'][i] + period < 0:
            #                 break
            #             if l == df_term['DOWNTIME'][i] - 1:
            #                 f3 += mu['tw'][ite][i][l-1,period+1]*(m_f[period]._tw[i,0] - svars['tw'][ite][i,period])
            #             else:
            #                 f3 += mu['tw'][ite][i][l-1,period+1]*(m_f[period]._aux_tw[i,period][l] - svars['tw'][ite][i,l-df_term['DOWNTIME'][i]+period+1])
                f4 += mu['gt'][ite][period+nblocks][i]*(m_f[period]._gt[i,nblocks-1] - svars['gt'][ite][i,period+nblocks-1])
            # for i in range(NH):
            #     f5 += mu['vol'][ite][period+1][i]*(m_f[period]._vol[i,0] - svars['vol'][ite][i,period])
            #     if df_hidr['TRAVELTIME'][i] > 0:
            #         for l in range(df_hidr['TRAVELTIME'][i],0,-1):
            #             if l-df_hidr['TRAVELTIME'][i]+period-1 < 0:
            #                 break
            #             if l == df_hidr['TRAVELTIME'][i]:
            #                 f6 += mu['turb'][ite][i][l-1,period+1]*(m_f[period]._turb[i,0] - svars['turb'][ite][i,period])
            #                 f6 += mu['vert'][ite][i][l-1,period+1]*(m_f[period]._vert[i,0] - svars['vert'][ite][i,period])
            #             else:
            #                 f6 += mu['turb'][ite][i][l-1,period+1]*(m_f[period]._aux_turb[i,period][l] - svars['turb'][ite][i,l-df_hidr['TRAVELTIME'][i]+period])
            #                 f6 += mu['vert'][ite][i][l-1,period+1]*(m_f[period]._aux_vert[i,period][l] - svars['vert'][ite][i,l-df_hidr['TRAVELTIME'][i]+period])
            m_f[period].addConstr(m_f[period]._Theta[0] >= obj[ite,period+nblocks] + f1 + f2 + f3 + f4 + f5 + f6) # Cut added to forward problem

            f1,f2,f3,f4,f5,f6 = 0,0,0,0,0,0
            for i in range(NG):
                f1 += mu['tu'][ite][period+nblocks][i]*(m_b[period]._tu[i,nblocks-1] - svars['tu'][ite][i,period+nblocks-1])
                if df_term['UPTIME'][i] > 1:
                    for l in range(df_term['UPTIME'][i]-1,0,-1):
                        if l - df_term['UPTIME'][i] + period + nblocks < 0:
                            break
                        if l == df_term['UPTIME'][i] - 1:
                            f2 += mu['tv'][ite][i][l-1,period+nblocks]*(m_b[period]._tv[i,period+nblocks -(df_term['UPTIME'][i]-l)-1] - svars['tv'][ite][i,period+nblocks-1])
                        else: #ver se, de acordo com o size do bloco, se é tv ou aux_tv
                            f2 += mu['tv'][ite][i][l-1,period+nblocks]*(m_b[period]._tv[i,period+nblocks -(df_term['UPTIME'][i]-l)-1] - svars['tv'][ite][i,period+nblocks -(df_term['UPTIME'][i]-l)])
                            # f2 += mu['tv'][ite][i][l-1,period+nblocks]*(m_b[period]._aux_tv[i,period][l] - svars['tv'][ite][i,period+nblocks -(df_term['UPTIME'][i]-l)-1])
            #     if df_term['DOWNTIME'][i] > 1:
            #         for l in range(df_term['DOWNTIME'][i]-1,0,-1):
            #             if l - df_term['DOWNTIME'][i] + period < 0:
            #                 break
            #             if l == df_term['DOWNTIME'][i] - 1:
            #                 f3 += mu['tw'][ite][i][l-1,period+1]*(m_b[period]._tw[i,0] - svars['tw'][ite][i,period])
            #             else:
            #                 f3 += mu['tw'][ite][i][l-1,period+1]*(m_b[period]._aux_tw[i,period][l] - svars['tw'][ite][i,l-df_term['DOWNTIME'][i]+period+1])
                f4 += mu['gt'][ite][period+nblocks][i]*(m_b[period]._gt[i,nblocks-1] - svars['gt'][ite][i,period+nblocks-1])
            # for i in range(NH):
            #     f5 += mu['vol'][ite][period+1][i]*(m_b[period]._vol[i,0] - svars['vol'][ite][i,period])
            #     if df_hidr['TRAVELTIME'][i] > 0:
            #         for l in range(df_hidr['TRAVELTIME'][i],0,-1):
            #             if l-df_hidr['TRAVELTIME'][i]+period-1 < 0:
            #                 break
            #             if l == df_hidr['TRAVELTIME'][i]:
            #                 f6 += mu['turb'][ite][i][l-1,period+1]*(m_b[period]._turb[i,0] - svars['turb'][ite][i,period])
            #                 f6 += mu['vert'][ite][i][l-1,period+1]*(m_b[period]._vert[i,0] - svars['vert'][ite][i,period])
            #             else:
            #                 f6 += mu['turb'][ite][i][l-1,period+1]*(m_b[period]._aux_turb[i,period][l] - svars['turb'][ite][i,l-df_hidr['TRAVELTIME'][i]+period])
            #                 f6 += mu['vert'][ite][i][l-1,period+1]*(m_b[period]._aux_vert[i,period][l] - svars['vert'][ite][i,l-df_hidr['TRAVELTIME'][i]+period])
            m_b[period].addConstr(m_b[period]._Theta[0] >= obj[ite,period+nblocks] + f1 + f2 + f3 + f4 + f5 + f6) # Cut added to backward problem
        # # # # # end insert cuts

        et = time.time()
        m_b[period].optimize();
        tet += time.time() - et

        m_b[period].write('model_py_%d.lp'%(period))

        if period > 1:
            mu['tu'][ite][period] = array([m_b[period].getConstrByName("c_tu[%d]"%i).pi for i in range(NG)])
            mu['gt'][ite][period] = array([m_b[period].getConstrByName("c_gt[%d]"%i).pi for i in range(NG)])
            # mu['vol'][ite][period] = array([m_b[period].getConstrByName("c_hvol[%d]"%i).pi for i in range(NH)])
            #
            for i in range(NG):
                if df_term['UPTIME'][i] > 1:
                    for l in range(df_term['UPTIME'][i]-1,0,-1):
                        mu['tv'][ite][i][l-1,period] = m_b[period].getConstrByName("c_tv[%d,%d]"%(i,l-1)).pi
            #
            # for i in range(NG):
            #     if df_term['DOWNTIME'][i] > 1:
            #         for l in range(df_term['DOWNTIME'][i]-1,0,-1):
            #             mu['tw'][ite][i][l-1,period] = m_b[period].getConstrByName("c_tw[%d,%d]"%(i,l-1)).pi
            #
            # for i in range(NH):
            #     if df_hidr['TRAVELTIME'][i] > 0:
            #         for l in range(df_hidr['TRAVELTIME'][i],0,-1):
            #             mu['turb'][ite][i][l-1,period] = m_b[period].getConstrByName("c_turb[%d,%d]"%(i,l-1)).pi
            #             mu['vert'][ite][i][l-1,period] = m_b[period].getConstrByName("c_vert[%d,%d]"%(i,l-1)).pi

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