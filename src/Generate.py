import numpy as np
import ZMP

q = np.zeros((851, 12))
H=0.90
x_start=0
theta=0

gait = ZMP.generate_gait(x_start,theta)
for i in range(len(gait)):
    pc = [gait[i, 1], gait[i, 2], H]
    pfl = gait[i, 3:6]
    pfr = gait[i, 6:9]
    L_Hip_Yaw, L_Hip_Roll, L_Hip_Pitch, L_knee, L_Ankle_Pitch, L_Ankle_Roll = ZMP.ik(pc, pfl, 1)
    R_Hip_Yaw, R_Hip_Roll, R_Hip_Pitch, R_knee, R_Ankle_Pitch, R_Ankle_Roll = ZMP.ik(pc, pfr, 2)
    q[i-1, 0:12] = [L_Hip_Yaw, L_Hip_Roll, L_Hip_Pitch, L_knee, L_Ankle_Pitch, L_Ankle_Roll,
                          R_Hip_Yaw, R_Hip_Roll, R_Hip_Pitch, R_knee, R_Ankle_Pitch, R_Ankle_Roll]