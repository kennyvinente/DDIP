from formulations import *


if __name__ =='__main__':

    NT = 24
    FPH_MODEL = 'CH'
    name_exp = FPH_MODEL

    model = Model('UC')
    model._name_exp = name_exp
    model._FPH_MODEL = FPH_MODEL
    model._CQ = 3600*1e-6
    model.setParam('TimeLimit', 600)
    model.setParam('MIPGAP',0)

    et = time.time()
    fph_model(model,NT)
    term_model(model,NT)
    uct_model(model, NT)
    network_model(model,NT)
    reserv_model(model, NT)
    # uch_model(m, NT)
    fobj(model,NT)


    print('Elapsed time for building model: %.2f seconds\n'%(time.time() - et))
    model.write('teste_new.lp')

    # model.Params.NonConvex = 2
    model.optimize()


    if model.Status == GRB.OPTIMAL:
        save_results(model, NT, name_exp)
        # measure_fpherror(model, name_exp)