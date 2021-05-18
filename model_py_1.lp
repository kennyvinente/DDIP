\ LP format - for model browsing. Use MPS format to capture full model detail.
Minimize
  0 s + Cost + Theta
Subject To
 ztu0: aux_u[0] = 0
 ztu1: aux_u[1] = 0
 ztu2: aux_u[2] = 0
 ztu3: aux_u[3] = 0
 ztv00: aux_v[0,0] = 0
 ztv01: aux_v[0,1] = 0
 ztv02: aux_v[0,2] = 0
 ztv03: aux_v[0,3] = 0
 ztv04: aux_v[0,4] = 0
 ztv05: aux_v[0,5] = 0
 ztv06: aux_v[0,6] = 0
 ztv07: - v[0] + aux_v[0,7] = 0
 zvol0: z_vol[0] = 9000
 R13: v[0] - w[0] - u[0] + aux_u[0] = 0
 R14: v[1] - w[1] - u[1] + aux_u[1] = 0
 R15: v[2] - w[2] - u[2] + aux_u[2] = 0
 R16: v[3] - w[3] - u[3] + aux_u[3] = 0
 R17: - u[0] + aux_v[0,0] + aux_v[0,1] + aux_v[0,2] + aux_v[0,3]
   + aux_v[0,4] + aux_v[0,5] + aux_v[0,6] + aux_v[0,7] <= 0
 R18: v[1] - u[1] <= 0
 R19: v[2] - u[2] <= 0
 R20: v[3] - u[3] <= 0
 R21: w[0] + u[0] <= 1
 R22: w[1] + u[1] <= 1
 R23: w[2] + u[2] <= 1
 R24: w[3] + u[3] <= 1
 R25: - 9000 b_gh + vol*b <= 0
 R26: - vol + vol*b <= 0
 R27: - 9000 b_gh - vol + vol*b >= -9000
 R28: gh - q - 0.005 vol*b = 0
 R29: gt[0] - 70 u[0] >= 0
 R30: gt[1] - 40 u[1] >= 0
 R31: gt[2] - 30 u[2] >= 0
 R32: gt[3] - 30 u[3] >= 0
 R33: gt[0] - 500 u[0] <= 0
 R34: gt[1] - 250 u[1] <= 0
 R35: gt[2] - 150 u[2] <= 0
 R36: gt[3] - 150 u[3] <= 0
 R37: gt[0] + gt[1] + gt[2] + gt[3] + gh + def - def = 600
 R38: - 9.95 gt[0] - 10.2 gt[1] - 11 gt[2] - 11 gt[3] - 800 v[0] - 380 v[1]
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
 b_gh <= 1
 gh <= 750
 vol <= 9000
 q <= 750
 def <= 100
 def <= 100
 aux_u[0] <= 1
 aux_u[1] <= 1
 aux_u[2] <= 1
 aux_u[3] <= 1
 aux_v[0,0] <= 1
 aux_v[0,1] <= 1
 aux_v[0,2] <= 1
 aux_v[0,3] <= 1
 aux_v[0,4] <= 1
 aux_v[0,5] <= 1
 aux_v[0,6] <= 1
 aux_v[0,7] <= 1
End
