# path = 'C:\\Users\\User\\Documents\\data118\\'
# path = 'C:\\Users\\User\\Documents\\dados_padrao\\'

path = 'C:\\Users\\User\\Documents\\data118en\\' #replace this with your current path

save_path = 'results'

import pandas as pd
from numpy import array, round, append, loadtxt, savetxt, unique, zeros, ones, insert, flatnonzero as find
from ext2int import ext2int
from dcopf import dcopf
from gurobipy import *
import time


baseMVA = 100
eps = 1e-4

df_hidr = pd.read_excel(path + 'hydrodata.xlsx')
NH = df_hidr.shape[0]

df_term = pd.read_excel(path + 'termdata.xlsx')
NG = df_term.shape[0]

df_bus = pd.read_excel(path + 'bus.xlsx')
NB = df_bus.shape[0]

df_branch = pd.read_excel(path + 'branch.xlsx')
NL = df_branch.shape[0]

df_load = pd.read_excel(path + 'load.xlsx')


def fph_model(model,NT):
    FPH_MODEL = model._FPH_MODEL

    gh = model.addMVar((NH,NT), name='gh')

    vol = model.addMVar((NH,NT), name='vol')
    turb = model.addMVar((NH,NT), name='turb')
    vert = model.addMVar((NH,NT), name='vert')

    for i in range(NH):
        for j in range(NT):
            vol[i,j].setAttr(GRB.Attr.VarName,"vol[%d,%d]"%(i,j+1))
            turb[i,j].setAttr(GRB.Attr.VarName,"turb[%d,%d]"%(i,j+1))
            vert[i,j].setAttr(GRB.Attr.VarName,"vert[%d,%d]"%(i,j+1))
            # print(gh[i,j].getAttr(GRB.Attr.VarName)) #example of how to print variables

    hu = model.addMVar((NH,NT),vtype='B', name='hon')
    pvu = model.addMVar((NH,NT),name='pvu')

    for i in range(NH):
        model.addConstr(vol[i,:] <= df_hidr['VMAX'][i])
        model.addConstr(turb[i,:] <= df_hidr['QMAX'][i]*df_hidr['NGU'][i]*hu[i,:])
        model.addConstr(turb[i,:] >= df_hidr['QMIN'][i]*hu[i,:])
        model.addConstr(vert[i,:] <= df_hidr['SMAX'][i])
        model.addConstr(vol[i,:] >= df_hidr['VMIN'][i])


    if FPH_MODEL == 'FREE':
        for i in range(NH):
            model.addConstr(gh[i,:] <= df_hidr['PMAX'][i]/baseMVA)

    if FPH_MODEL == 'CH':
        for i in range(NH):
            arr = loadtxt("HPF\\HPF_" + str(i + 1) + "_CH")
            model.addConstrs(baseMVA*gh[i,:] <= arr[k][0]*pvu[i,:] + arr[k][1]*turb[i,:] + arr[k][2]*vert[i,:] + arr[k][3]*hu[i,:] for k in range(arr.shape[0]))

        for i in range(NH):
            model.addConstr(pvu[i,:] <= df_hidr['VMAX'][i])
            model.addConstr(pvu[i,:] <= df_hidr['VMAX'][i]*hu[i,:])
            model.addConstr(pvu[i,:] >= df_hidr['VMIN'][i]*hu[i,:])
            model.addConstr(pvu[i,:] <= vol[i,:] - (1-hu[i,:])*df_hidr['VMIN'][i])
            model.addConstr(pvu[i,:] >= vol[i,:] - (1-hu[i,:])*df_hidr['VMAX'][i])

    model._gh = gh
    model._vol = vol
    model._turb = turb
    model._vert = vert

def term_model(model,NT):

    gt = model.addMVar((NG,NT),name="gt")

    tv = model.addMVar((NG,NT),name="tv", vtype='B')
    tw = model.addMVar((NG,NT),name="tw", vtype='B')
    tu = model.addMVar((NG,NT),name="tu", vtype='B')


    for i in range(NG):
        for j in range(NT):
            tv[i,j].setAttr(GRB.Attr.VarName,"tv[%d,%d]"%(i,j+1))
            tw[i,j].setAttr(GRB.Attr.VarName,"tw[%d,%d]"%(i,j+1))
            tu[i,j].setAttr(GRB.Attr.VarName,"tu[%d,%d]"%(i,j+1))
            gt[i,j].setAttr(GRB.Attr.VarName,"gt[%d,%d]"%(i,j+1))

    for i in range(NG):
        model.addConstr(gt[i,:] <= tu[i,:]*df_term['PMAX'][i]/baseMVA)
        model.addConstr(gt[i,:] >= tu[i,:]*df_term['PMIN'][i]/baseMVA)

    model._gt = gt
    model._tv = tv
    model._tw = tw
    model._tu = tu

def uct_model(model,NT):

    gt = model._gt
    tv = model._tv
    tw = model._tw
    tu = model._tu


    for i in range(NG):
        if df_term['TON'][i] > 0:
            model.addConstr(tv[i,0] - tw[i,0] == tu[i,0] - 1)
        else:
            model.addConstr(tv[i,0] - tw[i,0] == tu[i,0])
        model.addConstr(tv[i,1:] - tw[i,1:] == tu[i,1:] - tu[i,0:NT-1])


        for j in range(NT):
            if df_term['TON'][i] > 0:
                if (j+1) + df_term['TON'][i] <= df_term['UPTIME'][i]:
                    model.addConstr(tu[i,j] == 1)
                    model.addConstr(tv[i,j] == 0, name='upt[%d,%d]'%(i,j+1))
                    model.addConstr(tw[i,j] == 0)
                else:
                    model.addConstr(tv[i,max(0,j+1 - int(df_term['UPTIME'][i])):j+1].sum() <= tu[i,j], name='upt[%d,%d]'%(i,j+1))
            else:
                model.addConstr(tv[i,max(0,j+1 - int(df_term['UPTIME'][i])):j+1].sum() <= tu[i,j], name='upt[%d,%d]'%(i,j+1))


            if df_term['TON'][i] < 0:
                if (j+1) - df_term['TON'][i] <= df_term['DOWNTIME'][i]:
                    model.addConstr(tu[i,j] == 0)
                    model.addConstr(tv[i,j] == 0)
                    model.addConstr(tw[i,j] == 0, name='dwt[%d,%d]'%(i,j+1))
                else:
                    model.addConstr(tw[i,max(0,j+1 - int(df_term['DOWNTIME'][i])):j+1].sum() <=1 - tu[i,j], name='dwt[%d,%d]'%(i,j+1))

            else:
                model.addConstr(tw[i,max(0,j+1 - int(df_term['DOWNTIME'][i])):j+1].sum() <=1 - tu[i,j], name='dwt[%d,%d]'%(i,j+1))

        if df_term['TON'][i] > 0:
            model.addConstr(gt[i,0] <= min(df_term['RAMPUP'][i],1000)/baseMVA + df_term['P0'][i]/baseMVA)
        else:
            model.addConstr(gt[i,0] <= df_term['PMIN'][i]*tv[i,0]/baseMVA)
        model.addConstr(gt[i,1:] - gt[i,0:NT-1] <= min(df_term['RAMPUP'][i],1000)*tu[i,0:NT-1]/baseMVA + df_term['PMIN'][i]*tv[i,1:]/baseMVA)


        if df_term['TON'][i] > 0:
            model.addConstr(df_term['P0'][i]/baseMVA - gt[i, 0] <= min(df_term['RAMPDOWN'][i],1000)/baseMVA * tu[i, 0] + df_term['PMIN'][i] * tw[i, 0]/baseMVA)
        else:
            model.addConstr(-gt[i, 0] <= min(df_term['RAMPDOWN'][i],1000) * tu[i, 0]/baseMVA + df_term['PMIN'][i] * tw[i, 0]/baseMVA)

        model.addConstr(gt[i, 0:NT-1] - gt[i, 1:] <= min(df_term['RAMPDOWN'][i],1000) * tu[i, 1:]/baseMVA + df_term['PMIN'][i] * tw[i, 1:]/baseMVA)

def network_model(model,NT):

    ppc={}

    ppc['baseMVA'] = baseMVA

    # Bus, Pmax, Pmin, Qmax, Qmin, ..., status
    ppc['hidrogen'] = zeros((NH,11))
    for i in range(NH):
        ppc['hidrogen'][i] = array([df_hidr['BUS'][i],df_hidr['PMAX'][i],0,df_hidr['QMAX'][i],df_hidr['QMIN'][i],0,0,0,1,0,0])



    # bus, Pg, Qg, Qmax, Qmin, Vg, mBase, status, Pmax, Pmin
    ppc['gen'] = zeros((NG,10))
    for i in range(NG):
        ppc['gen'][i] = array([df_term['BUS'][i], 0, 0, df_term['QMAX'][i], df_term['QMIN'][i], 1.05, baseMVA, df_term['STATUS'][i], df_term['PMAX'][i], df_term['PMIN'][i]])

    # 1 startup shutdown n x1 y1 ... xn yn
    # 2 startup shutdown n c(n-1) ... c0
    ppc['gencost'] = zeros((NG,7))
    for i in range(NG):
        ppc['gencost'][i] = array([2, df_term['COST_START'][i], df_term['COST_SHUT'][i],3, df_term['COST_Q'][i], df_term['COST_L'][i], df_term['COST_F'][i]])


    # bus_i type Pd Qd Gs Bs area Vm Va baseKV zone Vmax Vmin
    ppc['bus'] = zeros((NB,13))
    for i in range(NB):
        ppc['bus'][i] = array([df_bus['NAME'][i], df_bus['TYPE'][i], df_bus['PD'][i], df_bus['QD'][i], df_bus['GS'][i], df_bus['BS'][i], df_bus['SUBSIST'][i],
                               df_bus['VM'][i], df_bus['VA'][i], df_bus['BASEKV'][i], df_bus['ZONE'][i], df_bus['VMAX'][i], df_bus['VMIN'][i]])


    # fbus, tbus, r, x, b, rateA, rateB, rateC, ratio, angle, status, angmin, angmax
    ppc['branch'] = zeros((NL,13))
    for i in range(NL):
        ppc['branch'][i] = df_branch.to_numpy()[i,1:14]


    Pload = zeros((NB,NT))
    for i in range(NB):
        for j in range(NT):
            Pload[i][j] = (df_bus['PD'][i]/sum(df_bus['PD']))*df_load['P_LOAD'][j]


    ppc = ext2int(ppc)


    nVa = model.addMVar((NB, NT), vtype='C', name='Va', lb=-GRB.INFINITY)

    nSlp = model.addMVar((NB, NT), vtype='C', name='nSl+', lb=0, ub=10)
    nSln = model.addMVar((NB, NT), vtype='C', name='nSL-', lb=0, ub=10)


    bus,gen,hidrogen,branch = ppc['bus'],ppc['gen'],ppc['hidrogen'],ppc['branch']


    B, Bf, neg_Cg, neg_Ch, Aang, lang, uang, iang = dcopf(ppc)
    model.addConstrs( (B @ nVa[:,j] + neg_Cg @ model._gt[:,j] + neg_Ch @ model._gh[:,j] + Pload[:,j]/baseMVA - nSlp[:,j] + nSln[:,j] == zeros(NB) for j in range(NT)), name='eq_DC')


    # # # DC Line Flow Limits
    model.addConstrs( (Bf @ nVa[:,j] <= branch[:,5]/baseMVA for j in range(NT)), name='eq_FL+')
    model.addConstrs( (Bf @ nVa[:,j] >= -branch[:,5]/baseMVA for j in range(NT)), name='eq_FL-')
    model.addConstrs( (Aang @ nVa[:,j] <= branch[:,12]*math.pi/180 for j in range(NT)), name='eq_theta+')
    model.addConstrs( (Aang @ nVa[:,j] >= branch[:,11]*math.pi/180 for j in range(NT)), name='eq_theta-')


    xxx = find(bus[:,1] == 3)
    if xxx.shape[0] > 0:
        for i in xxx:
            for j in range(NT):
                model.addConstr(nVa[i,j] == 0, name='eq_Vref[%d,%d]'%(i,j))

    model._nSlp = nSlp
    model._nSln = nSln

def reserv_model(model,NT):
    CTE = 3600*1e-6
    vmeta = 0 #% of acceptable outflow (hm3) on large reservoirs
    alpha = model.addVar(vtype='C', name='alpha')

    tviag = df_hidr['TRAVELTIME'].to_numpy()

    afl = df_hidr['Y0'].to_numpy()

    # # VOL META
    for i in range(NH):
        if df_hidr['TYPE'][i] == 1:
            model.addConstr(model._vol[i,NT-1] >= (1-vmeta/100)*(0.01*df_hidr['V0'][i]*(df_hidr['VMAX'][i]-df_hidr['VMIN'][i])+df_hidr['VMIN'][i]), name='vol_target_'+df_hidr['NAME'][i])


    # water balance
    for i in range(NH):
        temp = find(df_hidr['DOWNSTREAM'].to_numpy().astype(int) == i+1)
        V0 = 0.01*df_hidr['V0'][i]*(df_hidr['VMAX'][i]-df_hidr['VMIN'][i])+df_hidr['VMIN'][i]

        if temp.shape[0] > 0:
            for j in range(NT):
                if j == 0:
                    model.addConstr(model._vol[i,j] == V0 - CTE*(model._turb[i,j]+model._vert[i,j]-afl[i]-sum(df_hidr['Q0'][list(temp)])-sum(df_hidr['S0'][list(temp)])), name='wb_'+df_hidr['NAME'][i]+str(j))
                else:
                    xx1 = find(tviag[list(temp)] >= j+1)
                    if xx1.shape[0] > 0:
                        xx2 = 0
                        for k in range(temp.shape[0]):
                            if tviag[temp[k]] >= j + 1:
                                xx2 += df_hidr['Q0'][temp[k]] + df_hidr['S0'][temp[k]]
                            else:
                                xx2 += model._turb[temp[k],j-tviag[temp[k]]] + model._vert[temp[k],j-tviag[temp[k]]]
                        model.addConstr(model._vol[i,j] == model._vol[i,j-1]-CTE*(model._turb[i,j]+model._vert[i,j]-afl[i]-xx2), name='wb_'+df_hidr['NAME'][i]+str(j))
                    else:
                        xx2 = 0
                        for k in range(temp.shape[0]):
                            if df_hidr['TRAVELTIME'][temp[k]] > 0:
                                xx2 += model._turb[temp[k],j-tviag[temp[k]]] + model._vert[temp[k],j-tviag[temp[k]]]
                            else:
                                xx2 += model._turb[temp[k],j-tviag[temp[k]]] + model._vert[temp[k],j-tviag[temp[k]]]
                        model.addConstr(model._vol[i,j] == model._vol[i,j-1]- CTE*(model._turb[i,j]+model._vert[i,j]-afl[i]-xx2), name='wb_'+df_hidr['NAME'][i]+str(j))
        else:
            model.addConstr(model._vol[i,0] == V0 - CTE*(model._turb[i,0] + model._vert[i,0] - afl[i]), name='wb_'+df_hidr['NAME'][i])
            model.addConstr(model._vol[i,1:] == model._vol[i,0:NT-1] - CTE*(model._turb[i,1:] + model._vert[i,1:] - afl[i]), name='wb_'+df_hidr['NAME'][i])

    model._alpha = alpha

def fobj(model,NT):


    cdef = round(3*max(df_term['COST_L'].to_numpy()),2)

    fobj = 0
    for i in range(NG):
        fobj += baseMVA*df_term['COST_L'][i]*model._gt[i,:].sum()
        fobj += df_term['COST_F'][i]*model._tu[i,:].sum()
        fobj += df_term['COST_START'][i]*model._tv[i,:].sum()
        fobj += df_term['COST_SHUT'][i]*model._tw[i,:].sum()


    for i in range(NB):
        fobj += baseMVA*cdef*model._nSlp[i,:].sum()
        fobj += baseMVA*cdef*model._nSln[i,:].sum()



    model.setObjective(fobj)
