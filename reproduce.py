"""Independent 2-track ideal sliding-point calculation and frozen checks."""
from __future__ import annotations

import json
import math
from pathlib import Path

from scipy.integrate import quad, solve_ivp
from scipy.optimize import brentq

R = json.loads((Path(__file__).parent/"reference"/"result.json").read_text())


def close(label, actual, expected, tol=3e-7):
    if not math.isfinite(actual) or abs(actual-expected)>tol:
        raise AssertionError(f"{label}: {actual} != {expected}")


def arrival(acceleration, length):
    def finish(_t, state):
        return state[0]-length
    finish.direction, finish.terminal = 1, True
    sol=solve_ivp(lambda _t,y:(y[1],acceleration(y[0])),
                  (0,2),(0,0),events=finish,method="DOP853",
                  rtol=1e-12,atol=1e-13,max_step=.004)
    if not sol.success or len(sol.t_events[0])!=1:
        raise AssertionError("first-arrival event missing")
    return float(sol.t_events[0][0]),float(sol.y_events[0][0][1])


def run():
    X,H,G = R["horizontal_m"],R["drop_m"],R["gravity_m_s2"]
    close("registered x",X,1.0,1e-12)
    close("registered drop",H,.7,1e-12)
    close("registered g",G,9.81,1e-12)
    L=math.hypot(X,H)
    theta=brentq(lambda u:(1-math.cos(u))/(u-math.sin(u))-H/X,2.1,math.pi,
                 xtol=1e-14)
    a=X/(theta-math.sin(theta))
    close("cycloid x endpoint",a*(theta-math.sin(theta)),X,1e-11)
    close("cycloid y endpoint",a*(1-math.cos(theta)),H,1e-11)
    assert 0<theta<math.pi, "curve must remain downhill"
    C=4*a*(1-math.cos(theta/2))
    Tstraight=math.sqrt(2*L*L/(G*H))
    Tcycloid=theta*math.sqrt(a/G)
    speed=math.sqrt(2*G*H)
    # First independent reference is ds/v quadrature. Substitute u=sqrt(s)
    # to remove the integrable rest singularity on the straight segment.
    straight_quad=quad(lambda u:2*math.sqrt(L/(2*G*H)),0,math.sqrt(L),
                       epsabs=1e-12,epsrel=1e-12)[0]
    close("straight quadrature",straight_quad,Tstraight,1e-10)
    curve_quad=quad(lambda u:math.sqrt(a/G),0,theta,
                    epsabs=1e-12,epsrel=1e-12)[0]
    close("cycloid quadrature",curve_quad,Tcycloid,1e-10)
    arc_quad=quad(lambda u:2*a*math.sin(u/2),0,theta,
                  epsabs=1e-12,epsrel=1e-12)[0]
    close("cycloid arc quadrature",arc_quad,C,1e-10)
    # Independent arc-length initial-value solves and first crossing events.
    ts,vs=arrival(lambda _s:G*H/L,L)
    tc,vc=arrival(lambda s:G*(1-s/(4*a)),C)
    close("straight IVP",ts,Tstraight,2e-8)
    close("cycloid IVP",tc,Tcycloid,2e-8)
    close("straight energy",vs,speed,2e-8)
    close("cycloid energy",vc,speed,2e-8)
    values={"straight_length_m":L,"cycloid_length_m":C,
            "straight_time_s":Tstraight,"cycloid_time_s":Tcycloid,
            "length_extra_percent":100*(C/L-1),
            "time_saved_percent":100*(1-Tcycloid/Tstraight),
            "final_speed_m_s":speed}
    for name,value in values.items():
        close(name,value,R[name])
    assert C>L and Tcycloid<Tstraight
    wrong_straight,wrong_curve=L/speed,C/speed
    assert wrong_curve>wrong_straight, "negative control must reverse winner"
    assert abs(wrong_curve-Tcycloid)>.1, "wrong shortcut must be rejected"
    print("PASS: frozen values, independent quadrature/IVP, energy, negative control")
    print(json.dumps(values,indent=2))


if __name__=="__main__":
    run()
