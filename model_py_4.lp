\ LP format - for model browsing. Use MPS format to capture full model detail.
Minimize
  Cost + Theta
Subject To
 R0: v[0] - w[0] - u[0] + aux_u[0] = 0
 R1: v[1] - w[1] - u[1] + aux_u[1] = 0
 R2: v[2] - w[2] - u[2] + aux_u[2] = 0
 R3: v[3] - w[3] - u[3] + aux_u[3] = 0
 ztv46: aux_v[6] = 0
 ztv45: aux_v[5] = 0
 ztv44: aux_v[4] = 0
 ztv43: aux_v[3] = 0
 ztv42: aux_v[2] = 0
 ztv41: aux_v[1] = 0
 ztv40: aux_v[0] = 0
 R11: v[0] - u[0] + aux_v[0] + aux_v[1] + aux_v[2] + aux_v[3] + aux_v[4]
   + aux_v[5] + aux_v[6] <= 0
 ztu0: aux_u[0] = 0
 ztu1: aux_u[1] = 0
 ztu2: aux_u[2] = 0
 ztu3: aux_u[3] = 0
 R16: v[1] - u[1] <= 0
 R17: v[2] - u[2] <= 0
 R18: v[3] - u[3] <= 0
 R19: w[0] + u[0] <= 1
 R20: w[1] + u[1] <= 1
 R21: w[2] + u[2] <= 1
 R22: w[3] + u[3] <= 1
 R23: gt[0] - 70 u[0] >= 0
 R24: gt[1] - 40 u[1] >= 0
 R25: gt[2] - 30 u[2] >= 0
 R26: gt[3] - 30 u[3] >= 0
 R27: gt[0] - 500 u[0] <= 0
 R28: gt[1] - 250 u[1] <= 0
 R29: gt[2] - 150 u[2] <= 0
 R30: gt[3] - 150 u[3] <= 0
 R31: gt[0] + gt[1] + gt[2] + gt[3] + gh + def - def = 950
 R32: - 9.95 gt[0] - 10.2 gt[1] - 11 gt[2] - 11 gt[3] - 800 v[0] - 380 v[1]
   - 110 v[2] - 110 v[3] - 300 u[0] - 210 u[1] - 120 u[2] - 120 u[3]
   - alpha - 60 def - 60 def + Cost = 0
Bounds
 v[0] <= 1
 v[1] <= 1
 v[2] <= 1
 v[3] <= 1
 w[0] <= 1
 w[1] <= 1
 w[2] <= 1
 w[3] <= 1
 u[0] <= 1
 u[1] <= 1
 u[2] <= 1
 u[3] <= 1
 gh <= 750
 def <= 100
 def <= 100
 aux_u[0] <= 1
 aux_u[1] <= 1
 aux_u[2] <= 1
 aux_u[3] <= 1
 aux_v[0] <= 1
 aux_v[1] <= 1
 aux_v[2] <= 1
 aux_v[3] <= 1
 aux_v[4] <= 1
 aux_v[5] <= 1
 aux_v[6] <= 1
End
