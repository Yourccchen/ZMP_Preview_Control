import numpy as np
import matplotlib.pyplot as plt
from scipy.linalg import solve_discrete_are
from numpy.linalg import pinv

# Preview_Control
def Preview_Control_Improvement(pref, zc, dt):
    # Input parameters: pref is the ZMP reference position (desired step point, 1D data varying over time)
    # zc is the center of mass height, dt is the control time step
    g = 9.8
    T = (len(pref) - 1) * dt
    t = np.arange(0, T + dt, dt)

    N = 5000
    pref = np.concatenate((pref, np.ones(N) * pref[-1]))

    Q = 10000
    R = 1
    A = np.array([[1, dt, dt**2 / 2],
                  [0, 1, dt],
                  [0, 0, 1]])
    B = np.array([[dt**3 / 6], [dt**2 / 2], [dt]])
    C = np.array([[1, 0, -zc / g]])
    Ap = np.block([[1, C @ A], [np.zeros((3, 1)), A]])
    Bp = np.vstack((C @ B, B))
    Cp = np.array([[1, 0, 0, 0]])
    QQ = Cp.T * Q * Cp
    Pp = solve_discrete_are(Ap, Bp, QQ, R)
    Kp = pinv(R + Bp.T @ Pp @ Bp) @ (Bp.T @ Pp @ Ap)

    fp = np.zeros(N)
    R = np.array([R])
    Q = np.array([Q])
    for i in range(N):
        fp[i] = pinv(R + Bp.T @ Pp @ Bp) @ Bp.T @ (np.linalg.matrix_power((Ap - Bp @ Kp).T, i) @ Cp.T @ Q)

    u = np.zeros(len(t))
    du = np.zeros(len(t))
    x = np.zeros((3, len(t)))  # Center of mass position, velocity, and acceleration
    dx = np.zeros((3, len(t)))
    p = np.zeros(len(t))  # ZMP Position
    x[:, 0] = pref[0]
    p[:]= pref[0]

    xs = np.vstack((p, dx))

    for i in range(len(t)-1):
        du[i] = -Kp @ xs[:, i] + fp @ pref[i + 1:i + N + 1]
        xs[:, i + 1] = Ap @ xs[:, i] + Bp @ np.array([du[i]])
        dx[:, i + 1] = xs[1:4, i + 1]
        x[:, i + 1] = x[:, i] + dx[:, i + 1]
        p[i] = Cp @ xs[:, i]

    plt.plot(t, pref[:len(t)], label='Expected ZMP')
    plt.plot(t, p, label='Real ZMP')
    plt.plot(t, x[0, :len(t)], label='CoM')
    plt.grid()
    plt.legend()
    plt.show()

    gait = np.column_stack((t, pref[:len(t)], p, x[0, :len(t)]))
    return gait

# generate_gait
def generate_gait(x_start,theta):
    # Walking mode: start-walk-stop, get the change in the center of mass and foot positions over time
    # x_start: the starting point for generating the gait
    # theta: the counterclockwise angle relative to the positive direction of the x-axis.
    W0 = 0.2  # Hip width
    zc = 1  # Center of mass height
    Tstep = 0.5  # Half step cycle
    double_support_ratio = 0.1  # Proportion of double support time
    Lstep = 0.5  # Step length
    Nstep = 15  # Number of half steps
    dt = 0.01  # Time step
    H = 0.1  # Foot lifting height

    Ttot = Tstep * (Nstep + 2) # one cycle before and after for stabilization
    t = np.arange(0, Ttot + dt, dt)

    # ZMP
    pref = np.zeros((len(t), 3))   # The first column is time, the second column is ZMP x, the third column is ZMP y
    pref[:, 0] = t
    for i in range(2, Nstep + 2):
        if i % 2 == 0:  # Odd steps, right leg support
            pref[(i - 1) * round(Tstep / dt) + 1:i * round(Tstep / dt), 2] = -W0 / 2
        else:
            pref[(i - 1) * round(Tstep / dt) + 1:i * round(Tstep / dt), 2] = W0 / 2
        pref[(i - 1) * round(Tstep / dt) + 1:i * round(Tstep / dt) + 1, 1] = x_start + (i - 2) * Lstep
    pref[:round(Tstep / dt)+1, 1] = x_start
    pref[(Nstep + 1) * round(Tstep / dt) + 1:, 1] = x_start + (i - 2) * Lstep

    x = pref[:, 1]
    y = pref[:, 2]

    rotation_matrix = np.array([
        [np.cos(theta), -np.sin(theta)],
        [np.sin(theta), np.cos(theta)]
    ])

    rotated_xy = np.dot(rotation_matrix, np.vstack((x, y)))

    pref[:, 1] = rotated_xy[0, :]
    pref[:, 2] = rotated_xy[1, :]

    # CoM position calculation
    gaitx = Preview_Control_Improvement(pref[:, 1], zc, dt)
    gaity = Preview_Control_Improvement(pref[:, 2], zc, dt)


    Lfoot = np.zeros((len(t), 4))  #
    Lfoot[:, 0] = t
    Lfoot[:, 1] = 10**10
    Lfoot[:, 2] = W0 / 2
    Rfoot = np.zeros((len(t), 4))
    Rfoot[:, 0] = t
    Rfoot[:, 1] = 10**10
    Rfoot[:, 2] = -W0 / 2
    doubleN = round(double_support_ratio * (Tstep / dt) / 2)

    Lfoot[:round(Tstep / dt) + doubleN, 1] = pref[0, 1]
    Rfoot[:round(Tstep / dt) + doubleN, 1] = pref[0, 1]
    for i in range(len(t) - 1):
        if pref[i, 2] == 0 and pref[i + 1, 2] > 0:  # Left foot touch down moment
            if i + 1 - doubleN < 0:
                Lfoot[:, 1] = pref[i + 1, 1]
            else:
                Lfoot[i + 1 - doubleN:, 1] = pref[i + 1, 1]
        if pref[i, 2] > 0 and pref[i + 1, 2] == 0:  # Left foot lift off moment
            if i + 1 + doubleN < len(t):
                Lfoot[i + 1 + doubleN:, 1] = 10 ** 10
        if pref[i, 2] == 0 and pref[i + 1, 2] < 0:  # Right foot touch down moment
            if i + 1 - doubleN < 0:
                Rfoot[:, 1] = pref[i + 1, 1]
            else:
                Rfoot[i + 1 - doubleN:, 1] = pref[i + 1, 1]
        if pref[i, 2] < 0 and pref[i + 1, 2] == 0:  # Right foot lift off moment
            if i + 1 + doubleN < len(t):
                Rfoot[i + 1 + doubleN:, 1] = 10 ** 10


    Lfoot[len(t) - round(Tstep / dt) - doubleN:, 1] = pref[-1, 1]
    Rfoot[len(t) - round(Tstep / dt) - doubleN:, 1] = pref[-1, 1]

    # Interpolation for lift-off periods
    # Left leg
    while True:
        a0 = None
        for i in range(len(t) - 1):
            if Lfoot[i, 1] < 10**10 - 1 and Lfoot[i + 1, 1] > 10**10 - 1:  # 左脚离地
                a0 = Lfoot[i, 1]
                t0 = i
            if Lfoot[i, 1] > 10**10 - 1 and Lfoot[i + 1, 1] < 10**10 - 1:  # 左脚着地
                a1 = Lfoot[i + 1, 1]
                t1 = i + 1
        if a0 is None:
            break
        else:
            q, qd, qdd = Quintic_polynomial_interpolation_general(a0, a1, 0, 0, 0, 0, t1 - t0, np.arange(t0, t1) - t0)
            Lfoot[t0:t1, 1] = q
            q1, qd, qdd = Quintic_polynomial_interpolation_general(0, H, 0, 0, 0, 0, round((t1 - t0) / 2), np.arange(t0, t0 + round((t1 - t0) / 2)) - t0)
            q2, qd, qdd = Quintic_polynomial_interpolation_general(H, 0, 0, 0, 0, 0, t1 - (t0 + round((t1 - t0) / 2) + 1), np.arange(t0 + round((t1 - t0) / 2) + 1, t1) - (t0 + round((t1 - t0) / 2) + 1))
            Lfoot[t0:t1-1, 3] = np.concatenate((q1, q2))

    # Right leg
    while True:
        a0 = None
        for i in range(len(t) - 1):
            if Rfoot[i, 1] < 10**10 - 1 and Rfoot[i + 1, 1] > 10**10 - 1:  # 右脚离地
                a0 = Rfoot[i, 1]
                t0 = i
            if Rfoot[i, 1] > 10**10 - 1 and Rfoot[i + 1, 1] < 10**10 - 1:  # 右脚着地
                a1 = Rfoot[i + 1, 1]
                t1 = i + 1
        if a0 is None:
            break
        else:
            q, qd, qdd = Quintic_polynomial_interpolation_general(a0, a1, 0, 0, 0, 0, t1 - t0, np.arange(t0, t1) - t0)
            Rfoot[t0:t1, 1] = q
            q1, qd, qdd = Quintic_polynomial_interpolation_general(0, H, 0, 0, 0, 0, round((t1 - t0) / 2), np.arange(t0, t0 + round((t1 - t0) / 2)) - t0)
            q2, qd, qdd = Quintic_polynomial_interpolation_general(H, 0, 0, 0, 0, 0, t1 - (t0 + round((t1 - t0) / 2) + 1), np.arange(t0 + round((t1 - t0) / 2) + 1, t1) - (t0 + round((t1 - t0) / 2) + 1))
            Rfoot[t0:t1-1, 3] = np.concatenate((q1, q2))

    gait = np.column_stack((t, gaitx[:, 3], gaity[:, 3], Lfoot[:, 1:4], Rfoot[:, 1:4]))
    return gait

def ik(pc, pf, lr):
    # Foot inverse dynamics
    # Given the position of the center of mass (pc) and the foot center position (pf), calculate the foot joint angles
    # The output is the joint angles from top to bottom, with the default knee joint being bent forward
    if lr == 1:  # left
        phip = pc + np.array([0, 0.1, -0.175])  # 髋关节位置
    else:  # right
        phip = pc + np.array([0, -0.1, -0.175])  # 髋关节位置

    thigh_length = 0.42
    shin_length = 0.4
    pankle = pf + np.array([0, 0, 0.05])  # 蹽关节位置
    Lleg = np.linalg.norm(phip - pankle)  # 脚长

    # Knee angle
    L_knee = np.pi - 2 * np.arcsin(Lleg)
    # L_knee = np.arccos((thigh_length ** 2 + shin_length ** 2 - Lleg ** 2) / (2 * thigh_length * shin_length))

    # Side sway
    L_Hip_Roll = np.arcsin((pankle[1] - phip[1]) / np.linalg.norm(phip[1:3] - pankle[1:3]))
    L_Ankle_Roll = -L_Hip_Roll

    # Forward sway
    L_Hip_Pitch = -np.arcsin((pankle[0] - phip[0]) / Lleg) - L_knee / 2
    L_Ankle_Pitch = np.arcsin((pankle[0] - phip[0]) / Lleg) - L_knee / 2

    L_Hip_Yaw = 0

    return L_Hip_Yaw, L_Hip_Roll, L_Hip_Pitch, L_knee, L_Ankle_Pitch, L_Ankle_Roll

# Quintic polynomial interpolation
def Quintic_polynomial_interpolation_general(q0, q1, qd0, qd1, qdd0, qdd1, T, t):
    # Given initial and final positions, velocities, accelerations, and total time, return the position, velocity, and acceleration
    # after quintic interpolation at time t (0<t<T). t can also be a time sequence.
    A = np.array([[0, 0, 0, 0, 0, 1],
                  [T**5, T**4, T**3, T**2, T, 1],
                  [0, 0, 0, 0, 1, 0],
                  [5 * T**4, 4 * T**3, 3 * T**2, 2 * T, 1, 0],
                  [0, 0, 0, 2, 0, 0],
                  [20 * T**3, 12 * T**2, 6 * T, 2, 0, 0]])
    B = np.array([q0, q1, qd0, qd1, qdd0, qdd1])
    p = np.linalg.solve(A, B)
    q = p[0] * t**5 + p[1] * t**4 + p[2] * t**3 + p[3] * t**2 + p[4] * t + p[5]
    qd = 5 * p[0] * t**4 + 4 * p[1] * t**3 + 3 * p[2] * t**2 + 2 * p[3] * t + p[4]
    qdd = 20 * p[0] * t**3 + 12 * p[1] * t**2 + 6 * p[2] * t + 2 * p[3]
    return q, qd, qdd