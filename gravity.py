import numpy as np
from numba import njit
import matplotlib.pyplot as plt
from matplotlib import animation
from time import time


@njit
def gravity0(m, P0, V0, T, DELTA, ll):

    n = int(T/DELTA)
    G = 6.673e-11
    N = len(P0)
    l1 = np.minimum(ll, n)

    Pt = P0
    Vt = V0

    Rx = Ry = Rz = Sx = Sy = Sz = L = 0.0
    k = 0
    i = j = 0
    D = n//l1

    P = np.zeros((N, 3, n//D+1), dtype=np.double)
    V = np.zeros((N, 3, n//D+1), dtype=np.double)
    P[:, :, 0] = P0
    V[:, :, 0] = V0
    for k in range(n):
        for j in range(N):
            Sx = Sy = Sz = 0.0
            for i in range(N):
                Rx = Pt[i, 0]-Pt[j, 0]
                Ry = Pt[i, 1]-Pt[j, 1]
                Rz = Pt[i, 2]-Pt[j, 2]

                L = (Rx**2 + Ry**2 + Rz**2)**1.5 + (i == j)
                Sx += m[i]*Rx/L
                Sy += m[i]*Ry/L
                Sz += m[i]*Rz/L

            Pt[j, 0] += DELTA*Vt[j, 0] + DELTA*DELTA*G*Sx/2
            Pt[j, 1] += DELTA*Vt[j, 1] + DELTA*DELTA*G*Sy/2
            Pt[j, 2] += DELTA*Vt[j, 2] + DELTA*DELTA*G*Sz/2

            Vt[j, 0] += DELTA*G*Sx
            Vt[j, 1] += DELTA*G*Sy
            Vt[j, 2] += DELTA*G*Sz

            if (k+1) % D == 0:
                P[j, 0, (k+1)//D] = Pt[j, 0]
                P[j, 1, (k+1)//D] = Pt[j, 1]
                P[j, 2, (k+1)//D] = Pt[j, 2]
                V[j, 0, (k+1)//D] = Vt[j, 0]
                V[j, 1, (k+1)//D] = Vt[j, 1]
                V[j, 2, (k+1)//D] = Vt[j, 2]
    return P, V, T, m


def gravity(m, P0, V0, T, DELTA, ll=1000):
    m = np.array(m, dtype=np.double)
    P0 = np.array(P0, dtype=np.double)
    V0 = np.array(V0, dtype=np.double)
    t = time()
    res = gravity0(m, P0, V0, T, DELTA, ll)
    print(time()-t)
    P = res[0]
    m = res[3]
    end = res[2]
    b = P.shape[0]
    cols = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd',
            '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf']

    sm = 1.5*np.log(np.abs(m))
    sm1 = 4+sm-np.min(sm)
    pointsize = 5*sm1/np.min(sm1)
    SHAPE = P.shape
    N = SHAPE[2]
    lD = len(P)
    xmin = np.min(P[:, 0, :])
    xmax = np.max(P[:, 0, :])
    ymin = np.min(P[:, 1, :])
    ymax = np.max(P[:, 1, :])
    zmin = np.min(P[:, 2, :])
    zmax = np.max(P[:, 2, :])
    xl = xmax-xmin
    yl = ymax-ymin
    zl = zmax-zmin

    Ml = np.max([xl, yl, zl])
    if xl == 0:
        xl = 20*Ml
    if yl == 0:
        yl = 20*Ml
    if zl == 0:
        zl = 20*Ml

    fig = plt.figure(figsize=(16, 10))
    ax = fig.add_subplot(projection='3d')
    time_text = ax.text(xmin-0.015*xl, ymax+0.01*yl, zmax +
                        0.01*zl, '', fontsize=11, color=cols[0])
    dt = float((end/SHAPE[2])/(60*60*24))

    def update(num, P, line):
        for j in range(b):
            line[j].set_data(P[j, :2, :num])
            line[j].set_3d_properties(P[j, 2, :num])
        for j in range(b, 2*b):
            line[j].set_data(P[j-b, :2, num])
            line[j].set_3d_properties(P[j-b, 2, num])
        for i in range(2*lD, 3*lD):
            line[i].set_data([[P[i-2*lD, 0, num], P[i-2*lD, 0, num]],
                              [P[i-2*lD, 1, num], P[i-2*lD, 1, num]]])
            line[i].set_3d_properties([zmin-0.01*zl, P[i-2*lD, 2, num]])
        for i in range(3*lD, 4*lD):
            line[i].set_data(P[i-3*lD, :2, :num])
            line[i].set_3d_properties(zmin-0.01*zl)
            time_text.set_text('time = %.1f days' %
                               (num*dt))

    line = []
    for j in range(P.shape[0]):
        line += ax.plot(P[j, 0, 0:1], P[j, 1, 0:1],
                        P[j, 2, 0:1], linewidth=1)
    for i in range(P.shape[0]):
        line += ax.plot([], [], [], marker='o',
                        color=cols[i % 10], markersize=10)
    for i in range(lD):
        line += ax.plot([], [], [], linestyle='--',
                        linewidth=0.9, color=cols[i % 10])
    for i in range(lD):
        line += ax.plot([], [], [], linestyle='--',
                        linewidth=0.8, color=cols[i % 10])

    ax.set_xlim3d([xmin-0.01*xl, xmax+0.01*xl])
    ax.set_xlabel('X')

    ax.set_ylim3d([ymin-0.01*yl, ymax+0.01*yl])
    ax.set_ylabel('Y')

    ax.set_zlim3d([zmin-0.01*zl, zmax+0.01*zl])
    ax.set_zlabel('Z')
    # ax.set_facecolor('black')
    # plt.axis('off')
    plt.tight_layout()
    ani = animation.FuncAnimation(fig, update, N, fargs=(
        P, line), interval=10, blit=False)
    # ani.save('matplot003.gif', writer='imagemagick')
    plt.show()
