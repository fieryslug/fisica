import math

def update(P, Q, Hmltn, dt, dQ=0.0001, dP=0.0001):
    N = len(P)
    N_ = len(Q)
    if N != N_:
        raise ValueError("P and Q must have the same length!")

    P_ = [0 for i in range(N)]
    Q_ = [0 for i in range(N)]
    H = Hmltn(P, Q)
    for i in range(N):
        P[i] += dP
        dHdP = (Hmltn(P, Q) - H) / dP
        P[i] -= dP
        Q[i] += dQ
        dHdQ = (Hmltn(P, Q) - H) / dQ
        Q[i] -= dQ
        P_[i] = P[i] + (-dHdQ * dt)
        Q_[i] = Q[i] + dHdP * dt

    return (P_, Q_,)
        


def Hamilton(P, Q):
    pass
