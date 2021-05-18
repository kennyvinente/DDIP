
path = 'C:\\Users\\User\\Documents\\data118en\\'

save_path = 'results'

import pandas as pd
from numpy import array, round, append, loadtxt, savetxt, unique, zeros, ones, insert, linspace, flatnonzero as find
from ext2int import ext2int
from dcopf import dcopf
import time
import xlsxwriter
from gurobipy import *
from options import *

import scipy.sparse as sp

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


def fph_model(model,nblocks,t0,MILP_LP):
    FPH_MODEL = 'CH'
    NT = model._NT

    gh = model.addMVar((NH,nblocks), name='gh')
    vol = model.addMVar((NH,nblocks), name='vol')
    turb = model.addMVar((NH,nblocks), name='turb')
    vert = model.addMVar((NH,nblocks), name='vert')



    for i in range(NH):
        for j in range(nblocks):
            vol[i,j].setAttr(GRB.Attr.VarName,"vol[%d,%d]"%(i,j+t0))
            turb[i,j].setAttr(GRB.Attr.VarName,"turb[%d,%d]"%(i,j+t0))
            vert[i,j].setAttr(GRB.Attr.VarName,"vert[%d,%d]"%(i,j+t0))
            gh[i,j].setAttr(GRB.Attr.VarName,"gh[%d,%d]"%(i,j+t0))
            # print(gh[i,j].getAttr(GRB.Attr.VarName)) #exemplo de imprimir variaveis

    if MILP_LP:
        hu = model.addMVar((NH,nblocks),vtype='B', name='hon')
    else:
        hu = model.addMVar((NH,nblocks),vtype='C', name='hon',ub=1)
    pvu = model.addMVar((NH,nblocks),name='pvu')

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
            # for j in range(NT):
            #     model.addConstrs(baseMVA*gh[i,j] <= (round(arr[k][0],5) if abs(arr[k][0] >= 1e-4) else 0)*vol[i,j] + round(arr[k][1],5)*turb[i,j] + round(arr[k][2],5) for k in range(arr.shape[0]))

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

def term_model(model,nblocks,t0,MILP_LP):

    gt = model.addMVar((NG,nblocks),name="gt")

    if MILP_LP:
        tv = model.addMVar((NG,nblocks),name="tv", vtype='B')
        tw = model.addMVar((NG,nblocks),name="tw", vtype='B')
        tu = model.addMVar((NG,nblocks),name="tu", vtype='B')
    else:
        tv = model.addMVar((NG,nblocks),name="tv", vtype='C', ub=1)
        tw = model.addMVar((NG,nblocks),name="tw", vtype='C', ub=1)
        tu = model.addMVar((NG,nblocks),name="tu", vtype='C', ub=1)

    for i in range(NG):
        model.addConstr(gt[i,:] <= tu[i,:]*df_term['PMAX'][i]/baseMVA)
        model.addConstr(gt[i,:] >= tu[i,:]*df_term['PMIN'][i]/baseMVA)

    for i in range(NG):
        for j in range(nblocks):
            tv[i,j].setAttr(GRB.Attr.VarName,"tv[%d,%d]"%(i,j+t0))
            tw[i,j].setAttr(GRB.Attr.VarName,"tw[%d,%d]"%(i,j+t0))
            tu[i,j].setAttr(GRB.Attr.VarName,"tu[%d,%d]"%(i,j+t0))
            gt[i,j].setAttr(GRB.Attr.VarName,"gt[%d,%d]"%(i,j+t0))

    model._gt = gt
    model._tv = tv
    model._tw = tw
    model._tu = tu

def uct_model(model,svars,ite,nblocks,t0):

    gt = model._gt
    tv = model._tv
    tw = model._tw
    tu = model._tu

    aux_tu = model.addMVar(NG,name="aux_tu")
    aux_tv = {}
    aux_tw = {}
    aux_gt = model.addMVar(NG,name="aux_gt")

    model.addConstr(aux_tu == svars['tu'][ite][:,t0-1], name='c_tu')
    model.addConstr(aux_gt == svars['gt'][ite][:,t0-1], name='c_gt')


    model.addConstr(tv[:,0]-tw[:,0] == tu[:,0] - aux_tu)
    for j in range(1,nblocks):
        model.addConstr(tv[:,j]-tw[:,j] == tu[:,j] - tu[:,j-1])


    # # # UPTIME, DOWNTIME, RAMP
    for i in range(NG):
        aux_tv[i,t0] = model.addMVar(int(df_term['UPTIME'][i]-1), name='aux_tv[%d,%d]'%(i,t0), ub=1)
        for l in range(df_term['UPTIME'][i]-1,0,-1):
            model.addConstr(aux_tv[i,t0][l-1] == (svars['tv'][ite][i,l-df_term['UPTIME'][i]+t0] if l-df_term['UPTIME'][i]+t0 >= 0 else 0), name='c_tv[%d,%d]'%(i,l-1))

        for j in range(t0,t0+nblocks):
            if df_term['TON'][i] > 0:
                if j + df_term['TON'][i] <= df_term['UPTIME'][i]:
                    model.addConstr(tu[i,j-t0] == 1)
                    model.addConstr(tv[i,j-t0] == 0)
                    model.addConstr(tw[i,j-t0] == 0)
                else:
                    model.addConstr(aux_tv[i,t0][0:max(df_term['UPTIME'][i]-(j-t0) - 1,0)].sum() + tv[i,max(j-t0-df_term['UPTIME'][i]+1,0):j-t0+1].sum() <= tu[i,j-t0])

            else:
                model.addConstr(aux_tv[i,t0][0:max(df_term['UPTIME'][i]-(j-t0) - 1,0)].sum() + tv[i,max(j-t0-df_term['UPTIME'][i]+1,0):j-t0+1].sum() <= tu[i,j-t0])

    # arr = ones(13)
    # arr_ = linspace(0,12,13)
    # j=21
    # print(max(df_term['UPTIME'][i]-(j-t0) - 1,0))
    # print(max(j-t0-df_term['UPTIME'][i]+1,0))
    # print(j-t0+1)

    #
    #     aux_tw[i,t0] = model.addMVar(int(df_term['DOWNTIME'][i]-1), name='aux_tw[%d,%d]'%(i,t0), ub=1)
    #     for l in range(df_term['DOWNTIME'][i]-1,0,-1):
    #         model.addConstr(aux_tw[i,t0][l-1] == (svars['tw'][ite][i,l-df_term['DOWNTIME'][i]+t0] if l-df_term['DOWNTIME'][i]+t0 >= 0 else 0), name='c_tw[%d,%d]'%(i,l-1))
    #
    #     for j in range(t0,t0+nblocks):
    #         if df_term['TON'][i] < 0:
    #             if j + df_term['TON'][i] <= df_term['DOWNTIME'][i]:
    #                 model.addConstr(tu[i,j-t0] == 0)
    #                 model.addConstr(tv[i,j-t0] == 0)
    #                 model.addConstr(tw[i,j-t0] == 0)
    #             else:
    #                 model.addConstr(aux_tw[i,t0].sum() + tw[i,j-t0] <= 1 - tu[i,j-t0])
    #
    #         else:
    #             model.addConstr(aux_tw[i,t0].sum() + tw[i,j-t0] <= 1 - tu[i,j-t0])
    #
    #

        if t0 == 1:
            if df_term['TON'][i] > 0:
                model.addConstr(gt[i,0] <= min(df_term['RAMPUP'][i],1000)/baseMVA + df_term['P0'][i]/baseMVA)
            else:
                model.addConstr(gt[i,0] <= df_term['PMIN'][i]*tv[i,0]/baseMVA)
        else:
            model.addConstr(gt[i,0] - aux_gt[i] <= min(df_term['RAMPUP'][i],1000)*aux_tu[i]/baseMVA + df_term['PMIN'][i]*tv[i,0]/baseMVA)

        for j in range(1,nblocks):
            model.addConstr(gt[i,j] - gt[i,j-1] <= min(df_term['RAMPUP'][i],1000)*tu[i,j-1]/baseMVA + df_term['PMIN'][i]*tv[i,j]/baseMVA)

        if t0 == 1:
            if df_term['TON'][i] > 0:
                model.addConstr(df_term['P0'][i]/baseMVA - gt[i,0] <= min(df_term['RAMPDOWN'][i],1000)/baseMVA*tu[i,0] + df_term['PMIN'][i]*tw[i,0]/baseMVA)
            else:
                model.addConstr(-gt[i,0] <= min(df_term['RAMPDOWN'][i],1000)*tu[i,0]/baseMVA + df_term['PMIN'][i]*tw[i,0]/baseMVA)
        else:
            model.addConstr(aux_gt[i] - gt[i,0] <= min(df_term['RAMPDOWN'][i],1000)*tu[i,0]/baseMVA + df_term['PMIN'][i]*tw[i,0]/baseMVA)

        for j in range(1,nblocks):
            model.addConstr(gt[i,j-1] - gt[i,j] <= min(df_term['RAMPDOWN'][i],1000)*tu[i,j]/baseMVA + df_term['PMIN'][i]*tw[i,j]/baseMVA)

    model._aux_tv = aux_tv
    model._aux_tw = aux_tw
    model._aux_tu = aux_tu
    model._aux_gt = aux_gt


def network_model(model,nblocks,t0):
    # if LOAD_MODEL == 2:
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


    Pload = zeros((NB,nblocks))
    for i in range(NB):
        for j in range(t0,t0+nblocks):
            Pload[i][j-t0] = (df_bus['PD'][i]/sum(df_bus['PD']))*df_load['P_LOAD'][j-1]


    ppc = ext2int(ppc)


    nVa = model.addMVar((NB, nblocks), vtype='C', name='Va', lb=-GRB.INFINITY)

    nSlp = model.addMVar((NB, nblocks), vtype='C', name='nSl+', lb=0, ub=10)
    nSln = model.addMVar((NB, nblocks), vtype='C', name='nSL-', lb=0, ub=10)


    bus,gen,hidrogen,branch = ppc['bus'],ppc['gen'],ppc['hidrogen'],ppc['branch']


    B, Bf, neg_Cg, neg_Ch, Aang, lang, uang, iang = dcopf(ppc)
    model.addConstrs( (B @ nVa[:,j] + neg_Cg @ model._gt[:,j] + neg_Ch @ model._gh[:,j] + Pload[:,j]/baseMVA - nSlp[:,j] + nSln[:,j] == zeros(NB) for j in range(nblocks)), name='eq_DC')


    # # # DC Line Flow Limits
    model.addConstrs( (Bf @ nVa[:,j] <= branch[:,5]/baseMVA for j in range(nblocks)), name='eq_FL+')
    model.addConstrs( (Bf @ nVa[:,j] >= -branch[:,5]/baseMVA for j in range(nblocks)), name='eq_FL-')
    model.addConstrs( (Aang @ nVa[:,j] <= branch[:,12]*math.pi/180 for j in range(nblocks)), name='eq_theta+')
    model.addConstrs( (Aang @ nVa[:,j] >= branch[:,11]*math.pi/180 for j in range(nblocks)), name='eq_theta-')


    xxx = find(bus[:,1] == 3)
    if xxx.shape[0] > 0:
        for i in xxx:
            for j in range(nblocks):
                model.addConstr(nVa[i,j] == 0, name='eq_Vref[%d,%d]'%(i,j))

    model._nSlp = nSlp
    model._nSln = nSln

def reserv_model(model,svars,ite,nblocks,t0,NT):
    CTE = 3600*1e-6
    vmeta = 0 #% de esvaziamento aceitavel nos reservatorios de acumulacao
    # alpha = model.addVar(vtype='C', name='alpha')

    tviag = df_hidr['TRAVELTIME'].to_numpy()

    afl = df_hidr['Y0'].to_numpy()

    nwbp = model.addMVar((NH, nblocks), vtype='C', name='nwb+', lb=0, ub=10000)
    nwbn = model.addMVar((NH, nblocks), vtype='C', name='nwb-', lb=0, ub=10000)


    aux_vol = model.addMVar(NH,name="aux_vol")
    aux_turb = {}
    aux_vert = {}
    model.addConstr(aux_vol == svars['vol'][ite][:,t0-1], name='c_hvol')

    # # VOL META

    if t0 == NT:
        for i in range(NH):
            if df_hidr['TYPE'][i] == 1:
                model.addConstr(model._vol[i,:] >= (1-vmeta/100)*(0.01*df_hidr['V0'][i]*(df_hidr['VMAX'][i]-df_hidr['VMIN'][i])+df_hidr['VMIN'][i]), name='vol_target_'+df_hidr['NAME'][i])


    # # water balance
    for i in range(NH):
        aux_turb[i,t0] = model.addMVar(int(df_hidr['TRAVELTIME'][i]), name='aux_turb[%d,%d]'%(i,t0), ub=df_hidr['QMAX'][i]*df_hidr['NGU'][i])
        aux_vert[i,t0] = model.addMVar(int(df_hidr['TRAVELTIME'][i]), name='aux_vert[%d,%d]'%(i,t0), ub= df_hidr['SMAX'][i])
        for l in range(int(df_hidr['TRAVELTIME'][i]),0,-1):
            model.addConstr(aux_turb[i,t0][l-1] == (svars['turb'][ite][i,l-df_hidr['TRAVELTIME'][i]+t0-1] if l-df_hidr['TRAVELTIME'][i]+t0-1 >= 0 else 0), name='c_turb[%d,%d]'%(i,l-1))
            model.addConstr(aux_vert[i,t0][l-1] == (svars['vert'][ite][i,l-df_hidr['TRAVELTIME'][i]+t0-1] if l-df_hidr['TRAVELTIME'][i]+t0-1 >= 0 else 0), name='c_vert[%d,%d]'%(i,l-1))

    for i in range(NH):
        V0 = 0.01*df_hidr['V0'][i]*(df_hidr['VMAX'][i]-df_hidr['VMIN'][i])+df_hidr['VMIN'][i]
        temp = find(df_hidr['DOWNSTREAM'].to_numpy().astype(int) == i+1)

        if temp.shape[0] > 0:
            if t0 == 1:
                model.addConstr(model._vol[i,:] + nwbp[i,0] - nwbn[i,0] == V0 - CTE*(model._turb[i,:]+model._vert[i,:]-afl[i]-sum(df_hidr['Q0'][list(temp)])-sum(df_hidr['S0'][list(temp)])))
            else:
                xx1 = find(tviag[list(temp)] >= t0)
                if xx1.shape[0] > 0:
                    xx2 = 0
                    for k in range(temp.shape[0]):
                        if tviag[temp[k]] >= t0:
                            xx2 += df_hidr['Q0'][temp[k]] + df_hidr['S0'][temp[k]]
                        else:
                            xx2 += aux_turb[temp[k],t0][(df_hidr['TRAVELTIME'][k]-1)-(tviag[temp[k]]-1)] + aux_vert[temp[k],t0][(df_hidr['TRAVELTIME'][k]-1)-(tviag[temp[k]]-1)]
                    model.addConstr(model._vol[i,:] + nwbp[i,0] - nwbn[i,0] == aux_vol[i]-CTE*(model._turb[i,:] + model._vert[i,:] - afl[i]-xx2))
                else:
                    xx2 = 0
                    for k in range(temp.shape[0]):
                        if df_hidr['TRAVELTIME'][temp[k]] > 0: #esta é a parte com problema
                            xx2 += aux_turb[temp[k],t0][(df_hidr['TRAVELTIME'][temp[k]]-1)-(tviag[temp[k]]-1)] + aux_vert[temp[k],t0][(df_hidr['TRAVELTIME'][temp[k]]-1)-(tviag[temp[k]]-1)]
                        else:
                            xx2 += model._turb[temp[k],:] + model._vert[temp[k],:]
                    model.addConstr(model._vol[i,:] + nwbp[i,0] - nwbn[i,0] == aux_vol[i]- CTE*(model._turb[i,:] + model._vert[i,:] - afl[i]-xx2))

        else:
            if t0 == 1:
                model.addConstr(model._vol[i,:] + nwbp[i,0] - nwbn[i,0] == V0 - CTE*(model._turb[i,:] + model._vert[i,:] - afl[i]))
            else:
                model.addConstr(model._vol[i,:] + nwbp[i,0] - nwbn[i,0] == aux_vol[i] - CTE*(model._turb[i,:] + model._vert[i,:] - afl[i]))

    # model._alpha = alpha
    model._nwbp = nwbp
    model._nwbn = nwbn
    model._aux_turb = aux_turb
    model._aux_vert = aux_vert


def fobj(model,nblocks):

    cdef = round(3*max(df_term['COST_L'].to_numpy()),2)

    f1 = 0
    for i in range(NG):
        f1 += baseMVA*df_term['COST_L'][i]*model._gt[i,:].sum()
        f1 += df_term['COST_F'][i]*model._tu[i,:].sum()
        f1 += df_term['COST_START'][i]*model._tv[i,:].sum()
        f1 += df_term['COST_SHUT'][i]*model._tw[i,:].sum()

    for i in range(NB):
        f1 += baseMVA*cdef*model._nSlp[i,:].sum()
        f1 += baseMVA*cdef*model._nSln[i,:].sum()

    # for i in range(NH):
    #     f1 += baseMVA*cdef*model._nwbp[i,:].sum()
    #     f1 += baseMVA*cdef*model._nwbn[i,:].sum()

    return f1

