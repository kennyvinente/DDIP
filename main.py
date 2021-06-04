
from formulations import *

NT = 24
FPH_MODEL = 'CH'
name_exp = FPH_MODEL

model = Model()
model._name_exp = name_exp
model._FPH_MODEL = FPH_MODEL
model._CQ = 3600*1e-6
model.setParam('TimeLimit', 600)
model.setParam('MIPGAP',0.01)

et = time.time()
fph_model(model,NT)
term_model(model,NT)
uct_model(model, NT)
network_model(model,NT)
reserv_model(model, NT)
fobj(model,NT)

print('Elapsed time for building model: %.2f seconds\n'%(time.time() - et))
model.write('model_or.lp')

model.Params.NonConvex = 2

model.optimize()

print("Fobj: %.2f\n"%model.objval)
print("Runtime: %.2f\n"%model.Runtime)
print("MIPGap: %.2f\n"%float(model.MIPGap*100))


# p = model.relax()
# p.update()
# # p.write('model_or.lp')
# p.optimize()
# print("Fobj: %.2f\n"%p.objval)
# print("Runtime: %.2f\n"%p.Runtime)

