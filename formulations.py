# path = 'C:\\Users\\User\\Documents\\data118\\'
# path = 'C:\\Users\\User\\Documents\\dados_padrao\\'

path = 'C:\\Users\\User\\Documents\\data118en\\'

save_path = 'results'

import pandas as pd
from numpy import array, sum, round, append, loadtxt, savetxt, unique, zeros, ones, insert, flatnonzero as find
from ext2int import ext2int
from dcopf import dcopf
import os
import time
import xlsxwriter
from gurobipy import *
from options import *

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

def piecewise1d(model, x1, x2, z, fx, intervs, method="MC", dim=1, params=None):
    vmin = params['vmin']
    vmax = params['vmax']
    vref = 0.5*(vmax + vmin)

    if method == "MC":
        n = len(intervs)
        if dim == 1:
            y = model.addVars(n, vtype='B', name='y')
            xhat = model.addVars(n, vtype='C', name='xhat')
            zhat = model.addVars(n, vtype='C', name='zhat')

            model.addConstr(quicksum(y[i] for i in range(n)) <= 1)
            model.addConstr(quicksum(xhat[i] for i in range(n)) == x1)
            model.addConstr(quicksum(zhat[i] for i in range(n)) == baseMVA*z)

            for i in range(n):
                model.addConstr(xhat[i] >= intervs[i][0]*y[i])
            for i in range(n):
                model.addConstr(xhat[i] <= intervs[i][1]*y[i])
            for i in range(n):
                model.addConstr(zhat[i] == fx(vref,intervs[i][0])*y[i] + (fx(vref,intervs[i][1])-fx(vref,intervs[i][0]))/(intervs[i][1]-intervs[i][0])*(xhat[i]-intervs[i][0]*y[i]))
        else:
            y = model.addVars(n, vtype='B', name='y')
            xhat = model.addVars(n, vtype='C', name='xhat')
            zhat = model.addVars(n, vtype='C', name='zhat')

            model.addConstr(quicksum(y[i] for i in range(n)) <= 1)
            model.addConstr(quicksum(xhat[i] for i in range(n)) == x1)
            model.addConstr(quicksum(zhat[i] for i in range(n)) == baseMVA*z)

            for i in range(n):
                model.addConstr(xhat[i] >= intervs[i][0]*y[i])
            for i in range(n):
                model.addConstr(xhat[i] <= intervs[i][1]*y[i])

            pvy = model.addVars(n, vtype='C', name='pvy')
            for i in range(n):
                # beta_k = 0.5*((fx(vmax,intervs[i][1])-fx(vmin,intervs[i][1]))/(vmax-vmin))
                beta_p = ((fx(vmax, intervs[i][0]) - fx(vref, intervs[i][0])) / (vmax - vref))
                # model.addConstr(pvy[i] <= vmax)
                # model.addConstr(vmin*y[i] <= pvy[i])
                # model.addConstr(pvy[i] <= vmax*y[i])
                # model.addConstr(x2 - (1-y[i])*vmax <= pvy[i])
                # model.addConstr(pvy[i] <= x2 - (1-y[i])*vmin)
                # model.addConstr(pvy[i] <= x2 + (1-y[i])*vmax)

                # model.addConstr(zhat[i] == fx(Vfix,intervs[i][0])*y[i] +
                #                 (fx(Vfix,intervs[i][1])-fx(Vfix,intervs[i][0]))/(intervs[i][1]-intervs[i][0])*(xhat[i]-intervs[i][0]*y[i])
                #                 + beta_k*pvy[i] - beta_k*0.5*(vmax+vmin)*y[i] ) #+ beta_k*(x2 - 0.5*(vmax-vmin))*y[i]
                # model.addConstr(zhat[i] == fx(Vfix,intervs[i][0])*y[i] + (fx(Vfix,intervs[i][1])-fx(Vfix,intervs[i][0]))/(intervs[i][1]-intervs[i][0])*(xhat[i]-intervs[i][0]*y[i])) #+ beta_k*(x2 - 0.5*(vmax+vmin))*y[i]
                # model.addConstr(zhat[i] == fx(Vfix,intervs[i][0])*y[i] + (fx(Vfix,intervs[i][1])-fx(Vfix,intervs[i][0]))/(intervs[i][1]-intervs[i][0])*(xhat[i]-intervs[i][0]*y[i]) + beta_k*pvy[i] - beta_k*Vfix*y[i] ) #+ beta_k*(x2 - 0.5*(vmax-vmin))*y[i]
                model.addConstr(zhat[i] == fx(vref,intervs[i][0])*y[i] + (fx(vref,intervs[i][1])-fx(vref,intervs[i][0]))/(intervs[i][1]-intervs[i][0])*(xhat[i]-intervs[i][0]*y[i]) + beta_p*(x2 - vref)*y[i] ) #+ beta_k*(x2 - 0.5*(vmax-vmin))*y[i]

    else:
        print('erro, formulações PWL disponiveis: MC')
        exit(0)



    return 0 #z

def find_optzones(id_plant,df):
    nmaq = int(df['NMAQ'][id_plant - 1])
    qmax = df['QMAX'][id_plant - 1]
    qmin = df['QMIN'][id_plant - 1]

    opt_zones = array([[qmin, qmax]])
    for i in range(2, nmaq + 1):
        if i * qmin <= (i - 1) * qmax:
            break
        else:
            opt_zones = append(opt_zones, array([[i * qmin, i * qmax]]), axis=0)
    opt_zones[-1][1] = qmax * nmaq

    return opt_zones

def fph_model(model,NT):
    FPH_MODEL = model._FPH_MODEL

    CQ = model._CQ


    afl = df_hidr['Y0'].to_numpy()

    gh = model.addMVar((NH,NT), name='gh')
    vol = model.addMVar((NH,NT), name='vol')
    turb = model.addMVar((NH,NT), name='turb')
    vert = model.addMVar((NH,NT), name='vert')

    hu = model.addMVar((NH,NT),vtype='B', name='hon')
    pvu = model.addMVar((NH,NT),name='pvu')

    for i in range(NH):
        model.addConstr(vol[i,:] <= df_hidr['VMAX'][i])
        model.addConstr(turb[i,:] <= df_hidr['QMAX'][i]*df_hidr['NGU'][i]*hu[i,:])
        model.addConstr(turb[i,:] >= df_hidr['QMIN'][i]*hu[i,:])
        model.addConstr(vert[i,:] <= df_hidr['SMAX'][i])
        model.addConstr(vol[i,:] >= df_hidr['VMIN'][i])


    # for i in range(NH):
    #     for j in range(NT):
    #         gh[i,j] = model.addVar(name="GH_%s_%d"%(df_hidr['NAME'][i],j+1), lb=0)
    #         vol[i,j] = model.addVar(name="VOL_%s_%d"%(df_hidr['NAME'][i],j+1), lb=df_hidr['VMIN'][i], ub=df_hidr['VMAX'][i])
    #         turb[i,j] = model.addVar(name="TURB_%s_%d"%(df_hidr['NAME'][i],j+1), lb=0, ub=df_hidr['QMAX'][i]*df_hidr['NGU'][i])
    #         vert[i,j] = model.addVar(name="VERT_%s_%d"%(df_hidr['NAME'][i],j+1), lb=0, ub=df_hidr['SMAX'][i])

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


    if FPH_MODEL == 'UCH':
        hu = {}
        for i in range(NH):
            # arr = loadtxt(path + "FPH\\UHE_" + str(i + 1) + "\\CH_def")
            arr = loadtxt(path + "FPH\\UHE_" + str(i + 1) + "\\CH_" + VREF)
            PMIN = min(arr.dot([df['VMIN'][i],df['QMIN'][i],1]))
            for j in range(NT):
                hu[i,j] = model.addVar(name="LIG_%s_%d"%(df['NOME'][i],j+1), vtype='B')

                model.addConstr(baseMVA*gh[i,j] <= hu[i,j]*df['PMAX'][i])
                model.addConstr(baseMVA*gh[i,j] >= hu[i,j]*PMIN)
                model.addConstr(turb[i,j] <= hu[i,j]*df['QMAX'][i]*df['NMAQ'][i])
                model.addConstr(turb[i,j] >= hu[i,j]*df['QMIN'][i])

                model.addConstrs(baseMVA*gh[i,j] <= (round(arr[k][0],5) if abs(arr[k][0] >= 1e-4) else 0)*vol[i,j] + round(arr[k][1],5)*turb[i,j] + round(arr[k][2],5) for k in range(arr.shape[0]))

    if FPH_MODEL == 'PCH':
        y,x1hat,zhat = {},{},{}
        for i in range(NH):
            zones = find_optzones(i+1,df)
            for j in range(NT):
                for k in range(zones.shape[0]):
                    y[i,j,k] = model.addVar(vtype='B')
                    x1hat[i,j,k] = model.addVar()
                    zhat[i,j,k] = model.addVar()

                model.addConstr(quicksum(y[i,j,k] for k in range(zones.shape[0])) <= 1)
                model.addConstr(quicksum(x1hat[i,j,k] for k in range(zones.shape[0])) == turb[i,j])
                model.addConstr(quicksum(zhat[i,j,k] for k in range(zones.shape[0])) == gh[i,j])

                for k in range(zones.shape[0]):
                    model.addConstr(x1hat[i,j,k] <= y[i,j,k]*zones[k][1])
                    model.addConstr(x1hat[i,j,k] >= y[i,j,k]*zones[k][0])
                    arr = loadtxt(path + "FPH\\UHE_" + str(i + 1) + "\\PCH_" + str(gamma) + "_fx" + str(k) + "_inf")
                    model.addConstrs(baseMVA*zhat[i,j,k] >= (round(arr[l][0],5) if abs(arr[l][0] >= 1e-4) else 0)*vol[i,j]*y[i,j,k] + round(arr[l][1],5)*x1hat[i,j,k] + round(arr[l][2],5)*y[i,j,k] for l in range(arr.shape[0]))
                    arr = loadtxt(path + "FPH\\UHE_" + str(i + 1) + "\\PCH_" + str(gamma) + "_fx" + str(k) + "_sup")
                    model.addConstrs(baseMVA*zhat[i,j,k] <= (round(arr[l][0],5) if abs(arr[l][0] >= 1e-4) else 0)*vol[i,j]*y[i,j,k] + round(arr[l][1],5)*x1hat[i,j,k] + round(arr[l][2],5)*y[i,j,k] for l in range(arr.shape[0]))


    if FPH_MODEL == 'PWL':
        PWL_MODEL = array([2,1,1,1,1,1,1,1,1,1,1,1,2,2,1])
        # PWL_MODEL = ones(15)

        for ind_h in range(NH):
            def fph_N(V, Q):
                GH = 0
                S = 0
                for n in range(nmaq):
                    N = n + 1
                    if Q > qmax:
                        if (Q / N) > qmax + eps or (Q / N) < qmin - eps:
                            continue
                    else:
                        if Q < qmin:
                            continue
                        else:
                            N = 1
                    hb = F[0] + F[1] * V + F[2] * V ** 2 + F[3] * V ** 3 + F[4] * V ** 4 - \
                         (G[0] + G[1] * (Q + S) + G[2] * (Q + S) ** 2 + G[3] * (Q + S) ** 3 + G[4] * (Q + S) ** 4)
                    if hloss_flag == 1:
                        hl = hb * (1 - hloss / 100)
                    elif hloss_flag == 2:
                        hl = hb - hloss
                    else:
                        hl = hb - hloss * (Q / N) ** 2

                    rend = I[0] + I[1] * (Q / N) + I[2] * hl + I[3] * (Q / N) * hl + I[4] * (Q / N) ** 2 + I[5] * hl ** 2
                    rend = min(max(rend, 0), 1)
                    GH = max(GH, 0.00981 * rend * hl * Q)
                return GH

            id_plant = ind_h + 1
            arr = loadtxt(path + 'FPH\\UHE_' + str(id_plant) + '\\pwl_zones')

            F = zeros(5)
            G = zeros(5)
            I = zeros(6)

            for i in range(5):
                pos = 'F' + str(i)
                F[i] = df[pos][id_plant - 1]
                pos = 'G' + str(i)
                G[i] = df[pos][id_plant - 1]
                pos = 'I' + str(i)
                I[i] = df[pos][id_plant - 1]

            I[5] = df['I5'][id_plant - 1]
            hloss = df['H0'][id_plant - 1]
            hloss_flag = df['H1'][id_plant - 1]

            qmax = df['QMAX'][id_plant - 1]
            qmin = df['QMIN'][id_plant - 1]
            nmaq = int(df['NMAQ'][id_plant - 1])
            vmax = df['VMAX'][id_plant - 1]
            vmin = df['VMIN'][id_plant - 1]

            v0 = 0.01 * gamma * (vmax - vmin) + vmin
            if df['TYPE'][id_plant - 1] == 1:
                v1 = max(vmin, v0 - NT * CQ * qmax * nmaq)
            else:
                v1 = max(vmin, v0 - NT * CQ * (qmax * nmaq + df['SMAX'][id_plant - 1]))
            if find(df['JUSANTE'] == id_plant).shape[0] > 0:
                Qups = df['QMAX'][find(df['JUSANTE'] == id_plant)[0]] * df['NMAQ'][find(df['JUSANTE'] == id_plant)[0]]
            else:
                Qups = 0
            v2 = min(vmax, v0 + NT * CQ * (afl[id_plant - 1] + Qups))

            vref = 0.5 * (v1 + v2)

            intervs = []
            for i in range(arr.shape[0]):
                temp = round(arr[i][1:3], 5)
                for j in range(temp.shape[0] - 1):
                    intervs.append((temp[j], temp[j + 1]))

            params = {}
            params['vref'] = vref
            params['vmin'] = v1
            params['vmax'] = v2

            for j in range(NT):
                piecewise1d(model, turb[ind_h,j], vol[ind_h,j], gh[ind_h,j], fph_N, intervs, method='MC', dim = PWL_MODEL[ind_h], params=params)



    model._gh = gh
    model._vol = vol
    model._turb = turb
    model._vert = vert

def term_model(model,NT):

    gt = model.addMVar((NG,NT),name="gt")
    tv = model.addMVar((NG,NT),name="tstart", vtype='B')
    tw = model.addMVar((NG,NT),name="tshut", vtype='B')
    tu = model.addMVar((NG,NT),name="ton", vtype='B')

    # gt, tv, tw, tu = {},{},{},{}
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
                    model.addConstr(tv[i,j] == 0)
                    model.addConstr(tw[i,j] == 0)
                else:
                    if j+1 < df_term['UPTIME'][i]:
                        model.addConstr(tv[i,0:j+1].sum() <= tu[i,j])
                    else:
                        model.addConstr(tv[i,j+1 - int(df_term['UPTIME'][i]):j+1].sum() <= tu[i,j])
            else:
                if (j + 1) < df_term['UPTIME'][i]:
                    model.addConstr(tv[i,0:j+1] <= tu[i,j])
                else:
                    model.addConstr(tv[i,j+1 - int(df_term['UPTIME'][i]):j+1]  <= tu[i,j])

            if df_term['TON'][i] < 0:
                if (j+1) - df_term['TON'][i] <= df_term['DOWNTIME'][i]:
                    model.addConstr(tu[i,j] == 0)
                    model.addConstr(tv[i,j] == 0)
                    model.addConstr(tw[i,j] == 0)
                elif ((j + 1) -df_term['TON'][i] >= df_term['DOWNTIME'][i]) and ((j + 1) < df_term['DOWNTIME'][i]):
                    model.addConstr(tw[i,0:j+1].sum() <= 1 - tu[i,j])
                else:
                    model.addConstr(tw[i,j + 1 - int(df_term['DOWNTIME'][i]):j+1].sum()  <= 1 - tu[i,j])
            else:
                if (j + 1) < df_term['DOWNTIME'][i]:
                    model.addConstr(tw[i, 0:j+1].sum() <= 1 - tu[i, j])
                else:
                    model.addConstr(tw[i, j + 1 - int(df_term['DOWNTIME'][i]):j+1].sum() <= 1 - tu[i,j])

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
    # if LOAD_MODEL == 1:
    #     df = pd.read_excel(path + 'subs.xlsx')
    #     arr_ = zeros((df.shape[0],3))
    #     for i in range(df.shape[0]):
    #         arr_[i][0] = df['SUBS_DE'][i]
    #         arr_[i][1] = df['SUBS_PARA'][i]
    #         arr_[i][2] = df['PMAX'][i]
    #
    #
    #     df = pd.read_excel(path + 'bus.xlsx')
    #     dem_ref = df['PD'].to_numpy().astype(float)
    #     A_bs = (df['SUBSIST'].to_numpy()).astype(int) #vetor incidencia barra x subsistema
    #     NSUBS = unique(A_bs).shape[0] #numero de subsistemas
    #     NB = df.shape[0]
    #
    #
    #     df = pd.read_excel(path + 'hidrodata.xlsx')
    #
    #     A_ghb = (df['BARRA'].to_numpy()).astype(int)
    #
    #     Aghb = zeros((NSUBS,NH))  #matriz incidencia geradores hidraulicos x barra
    #     for i in range(NSUBS):
    #         for j in range(NH):
    #             if A_bs[A_ghb[j]-1] == i+1:
    #                 Aghb[i][j] = 1
    #
    #     df = pd.read_excel(path + 'termdata.xlsx')
    #     NG = df.shape[0]
    #     A_gtb = (df['BARRA'].to_numpy()).astype(int)
    #     Agtb = zeros((NSUBS,NG))  #matriz incidencia geradores termeletricos x barra
    #     for i in range(NSUBS):
    #         for j in range(NG):
    #             if A_bs[A_gtb[j]-1] == i+1:
    #                 Agtb[i][j] = 1
    #
    #     df = pd.read_excel(path + 'load.xlsx')
    #     total_demP = df['P_LOAD'].to_numpy().astype(float)
    #
    #     load = zeros((NT,NB)) #NT x NBARRA
    #     for j in range(NT):
    #         load[j] = round(total_demP[j]*dem_ref/sum(dem_ref),2)
    #
    #     INT = {}
    #     for k in range(arr_.shape[0]):
    #         for j in range(NT):
    #             INT[k,j] = model.addVar(name='INT[%d,%d]'%(k,j+1), lb=-arr_[k][2]/baseMVA,ub=arr_[k][2]/baseMVA)
    #     for k in range(NSUBS):
    #         for j in range(NT):
    #             temp = 0
    #             xxx = find(arr_[:,0].astype(int) == k+1)
    #             if xxx.shape[0] > 0:
    #                 for i in xxx:
    #                     temp += INT[i,j]
    #             xxx = find(arr_[:,1].astype(int) == k+1)
    #             if xxx.shape[0] > 0:
    #                 for i in xxx:
    #                     temp -= INT[i,j]
    #             PL = 0
    #             for i in range(NB):
    #                 PL += load[j][i] if A_bs[i] == k+1 else 0
    #             model.addConstr( quicksum(Agtb[k][i]*model._gt[i,j] for i in range(NG)) + quicksum(Aghb[k][i]*model._gh[i,j] for i in range(NH)) + temp == PL/baseMVA )
    #
    #
    #     # # # somente valido para arquivos DESSEM
    #     # df = pd.read_excel(path + 'load.xlsx')
    #     # NB = df.shape[0]
    #     # FC = zeros((NSUBS,NT)) #fator de carga por subsistema
    #     # for i in range(NSUBS):
    #     #     FC[i] = df['FC_S'+str(i+1)].to_numpy()
    #     #
    #     # load = zeros((NT,NB)) #NT x NBARRA
    #     # for i in range(NT):
    #     #     if df['PATAMAR'][i] == 'LEVE':
    #     #         load[i] = df['P_LEVE'].to_numpy()
    #     #     elif df['PATAMAR'][i] == 'MEDIA':
    #     #         load[i] = df['P_MEDIA'].to_numpy()
    #     #     else:
    #     #         load[i] = df['P_PESADA'].to_numpy()
    #     #
    #     # # otimizar esta parte do codigo
    #     # for i in range(NB):
    #     #     for k in range(NSUBS):
    #     #         if A_bs[i] == k+1:
    #     #             for j in range(NT):
    #     #                 load[j][i] *= FC[k][j]/100
    #     #
    #     # INT = {}
    #     # for k in range(arr_.shape[0]):
    #     #     for j in range(NT):
    #     #         INT[k,j] = model.addVar(name='INT[%d,%d]'%(k,j+1), lb=-arr_[k][2],ub=arr_[k][2])
    #     # for k in range(NSUBS):
    #     #     for j in range(NT):
    #     #         temp = 0
    #     #         xxx = find(int(arr_[:,0]) == k+1)
    #     #         if xxx.shape[0] > 0:
    #     #             for i in xxx:
    #     #                 temp += INT[i,j]
    #     #         xxx = find(int(arr_[:,1]) == k+1)
    #     #         if xxx.shape[0] > 0:
    #     #             for i in xxx:
    #     #                 temp -= INT[i,j]
    #     #         PL = 0
    #     #         for i in range(NB):
    #     #             PL += load[j][i] if A_bs[i] == k+1 else 0
    #     #         model.addConstr( quicksum(Agtb[k][i]*model._gt[i,j] for i in range(NG)) + quicksum(Aghb[k][i]*model._gh[i,j] for i in range(NH)) + temp == PL )
    #
    #     model._INT = INT

    if LOAD_MODEL == 2:

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
    vmeta = 0 #% de esvaziamento aceitavel nos reservatorios de acumulacao
    alpha = model.addVar(vtype='C', name='alpha')

    tviag = df_hidr['TRAVELTIME'].to_numpy()

    afl = df_hidr['Y0'].to_numpy()

    # # FCF
    # for l in range(fcf_coef.shape[0]):
    #     model.addConstr(alpha + quicksum(fcf_coef[l][i]*hvol[i,NT-1] for i in range(NH)) >= fcf_val[l], name='eq_FCF')

    # for ind_h in range(NH):
    #     if hidrogen[ind_h][H_TYPE] == 1:
    #         for j in range(NT):
    #             model.addConstr(hS[ind_h,j] == 0)


    # # VOL META

    for i in range(NH):
        if df_hidr['TYPE'][i] == 1:
            model.addConstr(model._vol[i,NT-1] >= (1-vmeta/100)*(0.01*df_hidr['V0'][i]*(df_hidr['VMAX'][i]-df_hidr['VMIN'][i])+df_hidr['VMIN'][i]), name='vol_target_'+df_hidr['NAME'][i])


    # water balance
    for i in range(NH):
        temp = find(df_hidr['DOWNSTREAM'].to_numpy().astype(int) == i+1)
        if temp.shape[0] > 0:
            for j in range(NT):
                if j == 0:
                    model.addConstr(model._vol[i,j] == 0.01*df_hidr['V0'][i]*(df_hidr['VMAX'][i]-df_hidr['VMIN'][i])+df_hidr['VMIN'][i]-
                                    CTE*(model._turb[i,j]+model._vert[i,j]-afl[i]-sum(df_hidr['Q0'][list(temp)])-sum(df_hidr['S0'][list(temp)])))
                else:
                    xx1 = find(tviag[list(temp)] >= j+1)
                    if xx1.shape[0] > 0:
                        xx2 = 0
                        for k in range(temp.shape[0]):
                            if tviag[temp[k]] >= j + 1:
                                xx2 += df_hidr['Q0'][temp[k]] + df_hidr['S0'][temp[k]]
                            else:
                                xx2 += model._turb[temp[k],j-tviag[temp[k]]] + model._vert[temp[k],j-tviag[temp[k]]]
                        model.addConstr(model._vol[i,j] == model._vol[i,j-1]-CTE*(model._turb[i,j]+model._vert[i,j]-afl[i]-xx2))
                    else:
                        xx2 = 0
                        for k in range(temp.shape[0]):
                            xx2 += model._turb[temp[k],j-tviag[temp[k]]] + model._vert[temp[k],j-tviag[temp[k]]]
                        model.addConstr(model._vol[i,j] == model._vol[i,j-1]- CTE*(model._turb[i,j]+model._vert[i,j]-afl[i]-xx2))
        else:
            model.addConstr(model._vol[i,0] == 0.01*df_hidr['V0'][i]*(df_hidr['VMAX'][i]-df_hidr['VMIN'][i])+df_hidr['VMIN'][i] - CTE*(model._turb[i,0] + model._vert[i,0] - afl[i]))
            model.addConstr(model._vol[i,1:] == model._vol[i,0:NT-1] - CTE*(model._turb[i,1:] + model._vert[i,1:] - afl[i]))

    model._alpha = alpha

def fobj(model,NT):


    cdef = round(3*max(df_term['COST_L'].to_numpy()),2)

    fobj = 0
    for i in range(NG):
        fobj += baseMVA*df_term['COST_L'][i]*model._gt[i,:].sum()
        fobj += df_term['COST_F'][i]*model._tu[i,:].sum()
        fobj += df_term['COST_START'][i]*model._tv[i,:].sum()
        fobj += df_term['COST_SHUT'][i]*model._tw[i,:].sum()
    # for i in range(NG):
    #     for j in range(NT):
    #         fobj += baseMVA*df_term['COST_L'][i]*model._gt[i,j] + df_term['COST_F'][i]*model._tu[i,j] + df_term['COST_START'][i]*model._tv[i,j] + df_term['COST_SHUT'][i]*model._tw[i,j]

    for i in range(NB):
        fobj += baseMVA*cdef*model._nSlp[i,:].sum()
        fobj += baseMVA*cdef*model._nSln[i,:].sum()

    # for i in range(NH):
    #     for j in range(NT):
    #         fobj += df['TYPE'][i]*cdef*model._vert[i,j]

    model.setObjective(fobj)

def uch_model(model,NT):
    Dp = 0.25 # Delta Potencia

    df = pd.read_excel(path + 'hidrodata.xlsx')
    NH = df.shape[0]
    for i in range(NH):
        for j in range(1,NT):
            model.addConstr(model._gh[i,j] >= (1-Dp)*model._gh[i,j-1])
            model.addConstr(model._gh[i,j] <= (1+Dp)*model._gh[i,j-1])

### post-processing
def save_results(model, NT, name_exp):
    df = pd.read_excel(path + 'termdata.xlsx')
    NG = df.shape[0]
    df = pd.read_excel(path + 'hydrodata.xlsx')
    NH = df.shape[0]


    tpgval = array([[round(model._gt[i,j].x[0],2) for j in range(NT)] for i in range(NG)])
    hpgval = array([[round(model._gh[i,j].x[0],2) for j in range(NT)] for i in range(NH)])
    hqval = array([[round(model._turb[i,j].x[0],2) for j in range(NT)] for i in range(NH)])
    hSval = array([[round(model._vert[i,j].x[0],2) for j in range(NT)] for i in range(NH)])
    hvolval = array([[round(model._vol[i,j].x[0],2) for j in range(NT)] for i in range(NH)])
    hvolval = insert(hvolval, 0, round((0.01*df['V0']*(df['VMAX']-df['VMIN'])+df['VMIN']).to_numpy(),2), axis=1)
    # nflval = np.array([[100 * nfl[i, j].x for j in range(NT)] for i in range(NL)])
    # slack = np.array([[nSl[i, j].x for j in range(NT)] for i in range(NB)])

    nSlp = array([[round(model._nSlp[i,j].x[0],2) for j in range(NT)] for i in range(NB)])
    nSln = array([[round(model._nSln[i,j].x[0],2) for j in range(NT)] for i in range(NB)])

    # savetxt(save_path+'/GT', tpgval, fmt='%f')
    # savetxt(save_path+'/GH', hpgval, fmt='%f')
    # savetxt(save_path+'/TURB', hqval, fmt='%f')
    # savetxt(save_path+'/VERT', hSval, fmt='%f')
    # savetxt(save_path+'/VOL', hvolval, fmt='%f')

    df = pd.read_excel(path + 'hydrodata.xlsx')
    NH = df.shape[0]
    workbook = xlsxwriter.Workbook(save_path + '\\results_' + name_exp + '_new.xlsx')
    worksheet = workbook.add_worksheet('gen_hydro')
    worksheet.set_column('A:A', 20)
    for j in range(NH):
        worksheet.write(j+2,0,df['NAME'][j])
        for i in range(NT):
            worksheet.write(j+2,i+1,hpgval[j][i])

    worksheet = workbook.add_worksheet('qtur')
    worksheet.set_column('A:A', 20)
    for j in range(NH):
        worksheet.write(j+2,0,df['NAME'][j])
        for i in range(NT):
            worksheet.write(j+2,i+1,hqval[j][i])

    worksheet = workbook.add_worksheet('vol')
    worksheet.set_column('A:A', 20)
    for j in range(NH):
        worksheet.write(j+2,0,df['NAME'][j])
        for i in range(NT+1):
            worksheet.write(j+2,i+1,hvolval[j][i])

    worksheet = workbook.add_worksheet('vert')
    worksheet.set_column('A:A', 20)
    for j in range(NH):
        worksheet.write(j+2,0,df['NAME'][j])
        for i in range(NT):
            worksheet.write(j+2,i+1,hSval[j][i])

    df = pd.read_excel(path + 'termdata.xlsx')
    NG = df.shape[0]
    worksheet = workbook.add_worksheet('gen_term')
    worksheet.set_column('A:A', 20)
    for j in range(NG):
        worksheet.write(j+2,0,'GT_' + str(df['NAME'][j]))
        for i in range(NT):
            worksheet.write(j+2,i+1,tpgval[j][i])

    worksheet = workbook.add_worksheet('outras_infos')
    worksheet.set_column('A:A', 20)
    worksheet.write(1,0,'EXP_TYPE')
    worksheet.write(1,1,name_exp.split("_")[0])
    worksheet.write(2,0,'Total_GH')
    worksheet.write(2,1,baseMVA*sum(hpgval))
    worksheet.write(3,0,'Total_GT')
    worksheet.write(3,1,baseMVA*sum(tpgval))
    worksheet.write(4,0,'Fobj')
    worksheet.write(4,1,model.objval)
    worksheet.write(5,0,'Runtime (s)')
    worksheet.write(5,1,model.Runtime)
    worksheet.write(6,0,'Nbv')
    worksheet.write(6,1,model.NumBinVars)
    worksheet.write(7,0,'Ncons')
    worksheet.write(7,1,model.NumConstrs)
    worksheet.write(8,0,'gap (%)')
    worksheet.write(8,1,round(model.MIPGAP*100,2))
    worksheet.write(9,0,'Slacks rede')
    worksheet.write(9,1,round(baseMVA*sum(nSlp),2)+round(baseMVA*sum(nSln),2))

    workbook.close()


def measure_fpherror(model,name_exp):
    eps = 1e-3
    name_exp = model._name_exp
    def fph_N(V, Q, S):
        GH = 0
        for n in range(nmaq):
            N = n + 1
            if Q > qmax:
                if (Q / N) > qmax + 0.0001 or (Q / N) < qmin - 0.0001:
                    continue
            else:
                if Q < qmin:
                    continue
                else:
                    N = 1
            hb = F[0] + F[1] * V + F[2] * V ** 2 + F[3] * V ** 3 + F[4] * V ** 4 - (
                    G[0] + G[1] * (Q + S) + G[2] * (Q + S) ** 2 + G[3] * (Q + S) ** 3 + G[4] * (Q + S) ** 4)
            hl = hb*(1-hloss/100) if hloss_flag == 1 else hb - hloss
            rend = I[0] + I[1]*(Q/N) + I[2]*hl + I[3]*(Q/N)*hl + I[4]*(Q/N)**2 + I[5]*hl**2
            rend = min(max(rend, 0), 1)
            GH = max(GH, 0.00981 * rend * hl * Q)
        return GH

    def find_optzones(id_plant):
        nmaq = int(df['NMAQ'][id_plant])
        qmax = df['QMAX'][id_plant]
        qmin = df['QMIN'][id_plant]

        opt_zones = array([[qmin, qmax]])
        for i in range(2, nmaq + 1):
            if i * qmin <= (i - 1) * qmax:
                break
            else:
                opt_zones = append(opt_zones, array([[i * qmin, i * qmax]]), axis=0)
        opt_zones[-1][1] = qmax * nmaq

        return opt_zones

    GT = loadtxt(save_path+'/GT')

    df = pd.read_excel(path + 'hidrodata.xlsx')
    GH = loadtxt(save_path+'/GH')
    GH *= baseMVA
    TURB = loadtxt(save_path+'/TURB')
    VERT = loadtxt(save_path+'/VERT')
    VOL = loadtxt(save_path+'/VOL')

    NH = GH.shape[0]
    NT = GH.shape[1]
    GH_violado = zeros(GH.shape[0])
    GH_erro = zeros(GH.shape[0])

    for i in range(NH):
        qmin = df['QMIN'][i]
        qmax = df['QMAX'][i]
        nmaq = int(df['NMAQ'][i])
        F = zeros(5)
        G = zeros(5)
        I = zeros(6)
        for j in range(5):
            pos = 'F' + str(j)
            F[j] = df[pos][i]
            pos = 'G' + str(j)
            G[j] = df[pos][i]
            pos = 'I' + str(j)
            I[j] = df[pos][i]
        I[5] = df['I5'][i]
        hloss = df['H0'][i]
        hloss_flag = df['H1'][i]

        zones = find_optzones(i)
        for j in range(NT):
            flag = 0
            for k in range(zones.shape[0]):
                if TURB[i][j] <= eps:
                    flag = 1
                    break
                if TURB[i][j] >= zones[k][0] -eps and TURB[i][j] <= zones[k][1] + eps:
                    flag = 1
                    break
            if flag == 0:
                GH_violado[i] += GH[i][j]
            else:
                GH_erro[i] += abs(fph_N(VOL[i][j],TURB[i][j],VERT[i][j]) - GH[i][j])


    workbook = xlsxwriter.Workbook(save_path + '\\uhe_results_' + name_exp + '.xlsx')
    worksheet = workbook.add_worksheet('erros_usinas')

    worksheet.set_column('A:A', 20)
    bold = workbook.add_format({'bold': True})
    for j in range(NH):
        worksheet.write(j+2,0,df['NOME'][j])
    worksheet.write('A2', 'UHE', bold)
    worksheet.write(1,1,'Total GH(MW)')
    for j in range(NH):
        worksheet.write(j+2,1,sum(GH[j]))
    worksheet.write(1,2,'GH vio(MW)')
    for j in range(NH):
        worksheet.write(j+2,2,GH_violado[j])
    worksheet.write(1,3,'GH erro(MW)')
    for j in range(NH):
        worksheet.write(j+2,3,GH_erro[j])
    worksheet.write(1,4,'GH vio(%)')
    for j in range(NH):
        worksheet.write(j+2,4,100*GH_violado[j]/sum(GH[j]) if sum(GH[j]) >= eps else ('Inf' if GH_violado[j] >= eps else 0) )
    worksheet.write(1,5,'GH erro(%)')
    for j in range(NH):
        worksheet.write(j+2,5,100*GH_erro[j]/sum(GH[j]) if sum(GH[j]) >= eps else ('Inf' if GH_erro[j] >= eps else 0) )
    worksheet.write(1,6,'GH/GHmax(%)U')
    for j in range(NH):
        worksheet.write(j+2,6,100*sum(GH[j])/(NT*df['PMAX'][j]))
    worksheet.write(1,7,'GH/GHmax(%)S')
    for j in range(NH):
        worksheet.write(j+2,7,100*sum(GH[j])/(sum(GH)))
    # worksheet.write(1,8,'erro-GH2')
    # for j in range(nT):
    #     worksheet.write(j+2,8,erro_gh[1,j])
    #

    worksheet = workbook.add_worksheet('outras_infos')
    worksheet.set_column('A:A', 20)
    worksheet.write(1,0,'EXP_TYPE')
    worksheet.write(1,1,name_exp.split("_")[0])
    worksheet.write(2,0,'gamma')
    worksheet.write(2,1,int(name_exp.split("_")[1]))
    worksheet.write(3,0,'afl')
    worksheet.write(3,1,name_exp.split("_")[2])
    worksheet.write(4,0,'Total_GH')
    worksheet.write(4,1,sum(GH))
    worksheet.write(5,0,'Total_GT')
    worksheet.write(5,1,baseMVA*sum(GT))
    worksheet.write(6,0,'Fobj')
    worksheet.write(6,1,model.objval)
    worksheet.write(7,0,'Runtime (s)')
    worksheet.write(7,1,model.Runtime)
    worksheet.write(8,0,'Nbv')
    worksheet.write(8,1,model.NumBinVars)
    worksheet.write(9,0,'Ncons')
    worksheet.write(9,1,model.NumConstrs)
    worksheet.write(10,0,'gap (%)')
    worksheet.write(10,1,round(model.MIPGAP*100,2))

    workbook.close()