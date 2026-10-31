"""
MBAM core for a synchronous-generator DAE  (Saric, Saric, Transtrum, Stankovic, IEEE TPWRS 2021, Sections IV-V)

Model  (two-axis, sub-transient; paper (A1)-(A2) extended with e''q, e''d, Milano 15.19 without gamma/T_AA terms)
    d(delta)/dt = Ob (w - 1)
    d(w)/dt     = (Pm - Pg - D (w - 1)) / 2H
    d(e'q)/dt   = (vf - e'q - (xd - x'd) id) / T'd0
    d(e'd)/dt   = (-e'd + (xq - x'q) iq) / T'q0
    d(e''q)/dt  = (-e''q + e'q - (x'd - x''d) id) / T''d0
    d(e''d)/dt  = (-e''d + e'd + (x'q - x''q) iq) / T''q0
    id = (e''q - V cos(delta - th)) / x''d ,  iq = (V sin(delta - th) - e''d) / x''q
    Pg = V sin(delta-th) id + V cos(delta-th) iq ,  Qg = V cos(delta-th) id - V sin(delta-th) iq
Inputs are the measured terminal V(t), th(t), vf(t), Pm(t).  Outputs compared with data: Pg, Qg, w, delta - th.

MBAM:  r(x) = (y_model - y_data)/sigma ;  J = dr/dx ; g = J^T J (Fisher information, x = log of parameters)
       geodesic  x'' = - g^{-1} J^T (d^2 r / dx dx)[x', x']   started along the sloppiest eigenvector of g,
       followed until the manifold boundary (g singular); the parameters that run to 0 / inf define the limiting model.
"""
import numpy as np
import scipy.linalg as sla
from scipy.optimize import least_squares

try:
    from numba import njit
    HAVE_NUMBA = True
except ImportError:                                   # slow fallback, same results
    HAVE_NUMBA = False
    def njit(*a, **k):
        if len(a) == 1 and callable(a[0]): return a[0]
        return lambda f: f

NAMES = ['H2', 'D', 'xd', 'xq', 'xdp', 'xqp', 'xdpp', 'xqpp', 'Td0p', 'Tq0p', 'Td0pp', 'Tq0pp', 'kf', 'km', 'c', 'eq0', 'ed0']
IS_LOG = np.array([True] * 14 + [False] * 3)
IDX = {n: i for i, n in enumerate(NAMES)}
BOUNDS = dict(H2=(0.5, 300.0), D=(1e-3, 200.0), xd=(0.05, 6.0), xq=(0.05, 6.0), xdp=(0.02, 3.0), xqp=(0.02, 3.0), xdpp=(0.02, 3.0),
              xqpp=(0.02, 3.0), Td0p=(0.2, 60.0), Tq0p=(0.02, 10.0), Td0pp=(0.004, 1.0), Tq0pp=(0.004, 1.0), kf=(1e-6, 1.0), km=(0.3, 3.0))
DESC = {'H2': '2H inertia [s]', 'D': 'damping [pu]', 'xd': 'xd [pu]', 'xq': 'xq [pu]', "xdp": "x'd [pu]", "xqp": "x'q [pu]",
        'xdpp': "x''d [pu]", 'xqpp': "x''q [pu]", 'Td0p': "T'd0 [s]", 'Tq0p': "T'q0 [s]", 'Td0pp': "T''d0 [s]", 'Tq0pp': "T''q0 [s]",
        'kf': 'field-voltage scale [pu/raw]', 'km': 'Pm scale [-]', 'c': 'angle offset [rad]', 'eq0': "e'q(t0) [pu]", 'ed0': "e'd(t0) [pu]"}


@njit(cache=True)
def _sim(Phi, Psi, a1, a2, a3, a4, xdpp, xqpp, Tv, H2, D, Ob, ns, h, Vf, thf, vff, Pmf, d0, w0, e0, outP, outQ, outW, outD):
    M = len(Vf)
    d = d0; w = w0
    e = e0.copy()
    for k in range(M):
        if k % ns == 0:
            n = k // ns
            dd = d - thf[k]
            Vd = Vf[k] * np.sin(dd); Vq = Vf[k] * np.cos(dd)
            idc = (e[2] - Vq) / xdpp; iqc = (Vd - e[3]) / xqpp
            outP[n] = Vd * idc + Vq * iqc; outQ[n] = Vq * idc - Vd * iqc
            outW[n] = w; outD[n] = dd
        if k == M - 1:
            break
        # --- step k -> k+1
        dd = d - thf[k]
        Vd = Vf[k] * np.sin(dd); Vq = Vf[k] * np.cos(dd)
        idc = (e[2] - Vq) / xdpp; iqc = (Vd - e[3]) / xqpp
        Pg = Vd * idc + Vq * iqc
        b = np.empty(4)
        b[0] = (vff[k] + a1 * Vq) / Tv[0]; b[1] = a2 * Vd / Tv[1]; b[2] = a3 * Vq / Tv[2]; b[3] = a4 * Vd / Tv[3]
        f_w = (Pmf[k] - Pg - D * (w - 1.0)) / H2
        dp = d + h * Ob * (w - 1.0); wp = w + h * f_w
        ep = Phi @ e + Psi @ b
        dd1 = dp - thf[k + 1]
        Vd1 = Vf[k + 1] * np.sin(dd1); Vq1 = Vf[k + 1] * np.cos(dd1)
        id1 = (ep[2] - Vq1) / xdpp; iq1 = (Vd1 - ep[3]) / xqpp
        Pg1 = Vd1 * id1 + Vq1 * iq1
        b1 = np.empty(4)
        b1[0] = (vff[k + 1] + a1 * Vq1) / Tv[0]; b1[1] = a2 * Vd1 / Tv[1]; b1[2] = a3 * Vq1 / Tv[2]; b1[3] = a4 * Vd1 / Tv[3]
        f_w1 = (Pmf[k + 1] - Pg1 - D * (wp - 1.0)) / H2
        d = d + 0.5 * h * Ob * ((w - 1.0) + (wp - 1.0))
        w = w + 0.5 * h * (f_w + f_w1)
        e = Phi @ e + Psi @ (0.5 * (b + b1))


class SGProblem:
    """One generator: data window + model structure + parameter handling."""

    def __init__(self, t, V, th, vf_raw, Pm_raw, P, Q, w, dmt, nominal, Ob, salient=False, ns=4, sig_frac=0.005,
                 delta0_meas=None, fit_mask=None):
        self.t = t; self.N = len(t); self.ns = ns; self.Ob = Ob; self.salient = salient
        self.fit_mask = np.ones(self.N, dtype=bool) if fit_mask is None else np.asarray(fit_mask, dtype=bool)
        if self.fit_mask.shape != (self.N,) or not self.fit_mask.any():
            raise ValueError("fit_mask must be a nonempty boolean vector with one entry per sample")
        self.h = (t[1] - t[0]) / ns
        tf = np.linspace(t[0], t[-1], (self.N - 1) * ns + 1)
        self.Vf = np.interp(tf, t, V); self.thf = np.interp(tf, t, th)
        self.vraw = np.interp(tf, t, vf_raw); self.Pmraw = np.interp(tf, t, Pm_raw)
        self.y = dict(P=P, Q=Q, w=w, dmt=dmt)           # measured outputs (dmt = delta - theta_terminal)
        self.sig = {
            k: max(sig_frac * (np.max(v[self.fit_mask]) - np.min(v[self.fit_mask])), 1e-9)
            for k, v in self.y.items()
        }
        self.w0 = w[0]; self.th0 = th[0]; self.V0 = V[0]
        self.p_base = np.array([nominal[n] for n in NAMES], float)
        self.fixed = {}                                  # name -> fixed physical value
        if salient:
            self.fixed['ed0'] = 0.0; self.fixed['Tq0p'] = 1.0
        self.outs = [np.zeros(self.N) for _ in range(4)]

    # ---------- parameter <-> coordinate mapping (free parameters only)
    def free_names(self):
        return [n for n in NAMES if n not in self.fixed]

    def full_p(self, x, free=None):
        free = free or self.free_names()
        p = self.p_base.copy()
        for n, v in self.fixed.items(): p[IDX[n]] = v
        for n, xi in zip(free, x):
            i = IDX[n]; p[i] = np.exp(np.clip(xi, -25.0, 25.0)) if IS_LOG[i] else xi
        if self.salient: p[IDX['xqp']] = p[IDX['xq']]
        return p

    def to_x(self, p, free=None):
        free = free or self.free_names()
        return np.array([np.log(p[IDX[n]]) if IS_LOG[IDX[n]] else p[IDX[n]] for n in free])

    # ---------- simulation
    def simulate(self, p):
        g = lambda n: p[IDX[n]]
        xd, xq, xdp, xqp, xdpp, xqpp = g('xd'), g('xq'), g('xdp'), g('xqp'), g('xdpp'), g('xqpp')
        Tv = np.array([g('Td0p'), g('Tq0p'), g('Td0pp'), g('Tq0pp')])
        a1, a2 = (xd - xdp) / xdpp, (xq - xqp) / xqpp
        a3, a4 = (xdp - xdpp) / xdpp, (xqp - xqpp) / xqpp
        A = np.zeros((4, 4))
        A[0, 0] = -1 / Tv[0]; A[0, 2] = -a1 / Tv[0]
        A[1, 1] = -1 / Tv[1]; A[1, 3] = -a2 / Tv[1]
        A[2, 2] = -(1 + a3) / Tv[2]; A[2, 0] = 1 / Tv[2]
        A[3, 3] = -(1 + a4) / Tv[3]; A[3, 1] = 1 / Tv[3]
        Phi = sla.expm(A * self.h)
        Psi = np.linalg.solve(A, Phi - np.eye(4))
        vff = g('kf') * self.vraw; Pmf = g('km') * self.Pmraw
        # initial state: e'q, e'd free; e''q, e''d quasi-steady with them
        dd0 = self.dmt0_model(g('c'))
        Vd0 = self.V0 * np.sin(dd0); Vq0 = self.V0 * np.cos(dd0)
        eq0, ed0 = g('eq0'), g('ed0')
        e0 = np.array([eq0, ed0, (eq0 + a3 * Vq0) / (1 + a3), (ed0 + a4 * Vd0) / (1 + a4)])
        d0 = self.th0 + dd0
        outP, outQ, outW, outD = [np.zeros(self.N) for _ in range(4)]
        _sim(Phi, Psi, a1, a2, a3, a4, xdpp, xqpp, Tv, g('H2'), g('D'), self.Ob, self.ns, self.h,
             self.Vf, self.thf, vff, Pmf, d0, self.w0, e0, outP, outQ, outW, outD)
        return outP, outQ, outW, outD

    def dmt0_model(self, c):
        return self.y['dmt'][0] + c                      # model load angle at t0 = measured + offset c

    # ---------- residuals
    def resid_p(self, p):
        P, Q, W, D = self.simulate(p)
        y = self.y
        # delta offset c enters the model through d0; the measured load angle is compared WITHOUT offset reapplied
        r = np.concatenate([(P - y['P']) / self.sig['P'], (Q - y['Q']) / self.sig['Q'],
                            (W - y['w']) / self.sig['w'], (D - (y['dmt'] + p[IDX['c']])) / self.sig['dmt']])
        return np.nan_to_num(r, nan=1e6, posinf=1e6, neginf=-1e6)

    def resid(self, x, free=None):
        r = self.resid_p(self.full_p(x, free)).reshape(4, self.N)
        return r[:, self.fit_mask].ravel()

    def jac(self, x, free=None, eps=1e-5):
        r0 = self.resid(x, free); J = np.empty((len(r0), len(x)))
        for i in range(len(x)):
            dx = np.zeros_like(x); dx[i] = eps
            J[:, i] = (self.resid(x + dx, free) - self.resid(x - dx, free)) / (2 * eps)
        return J

    def cost(self, x, free=None):
        r = self.resid(x, free); return 0.5 * float(r @ r)

    # ---------- least squares fit
    def fit(self, p_start, free=None, max_nfev=300, verbose=0, xtol=1e-12, ftol=1e-12, gtol=1e-10):
        free = free or self.free_names()
        x0 = self.to_x(p_start, free)
        lo, hi = [], []
        for n in free:
            i = IDX[n]
            if IS_LOG[i]:
                b = BOUNDS.get(n, (1e-6, 1e6)); lo.append(np.log(b[0])); hi.append(np.log(b[1]))
            else: lo.append(-20.0); hi.append(20.0)
        lo, hi = np.array(lo), np.array(hi)
        x0 = np.clip(x0, lo + 1e-9, hi - 1e-9)
        sol = least_squares(lambda x: self.resid(x, free), x0, bounds=(lo, hi), x_scale='jac', max_nfev=max_nfev,
                            xtol=xtol, ftol=ftol, gtol=gtol, verbose=verbose)
        return self.full_p(sol.x, free), sol.cost, sol


# ====================================================================== MBAM ======================================
def fisher(prob, x, free):
    J = prob.jac(x, free)
    g = J.T @ J
    w, v = np.linalg.eigh(g)
    return J, g, w, v


def geodesic(prob, x0, free, v0, dtau=0.05, max_steps=400, max_excursion=9.0, h_fd=0.05, lam_stop=1e-9, verbose=False):
    """Integrate x'' = -g^{-1} J^T A_vv  (Transtrum et al.) with RK2 (midpoint); stop at the manifold boundary."""
    free_log = np.array([IS_LOG[IDX[n]] for n in free])
    def accel(x, v):
        J = prob.jac(x, free)
        g = J.T @ J + 1e-12 * np.eye(len(x))
        r0 = prob.resid(x, free)
        Avv = (prob.resid(x + h_fd * v, free) + prob.resid(x - h_fd * v, free) - 2 * r0) / h_fd ** 2
        return -np.linalg.solve(g, J.T @ Avv), np.linalg.eigvalsh(J.T @ J)
    x = x0.copy(); v = v0.copy()
    path = [(0.0, x.copy(), v.copy(), None, prob.cost(x, free))]
    reason = 'max steps'
    for k in range(max_steps):
        try:
            a, ev = accel(x, v)
            xm = x + 0.5 * dtau * v; vm = v + 0.5 * dtau * a
            am, ev_m = accel(xm, vm)
            x_new = x + dtau * vm; v_new = v + dtau * am
            _, ev = accel(x_new, v_new)
        except (ZeroDivisionError, FloatingPointError, np.linalg.LinAlgError):
            reason = 'numerical failure (boundary)'; break
        x, v = x_new, v_new
        path.append(((k + 1) * dtau, x.copy(), v.copy(), ev.copy(), prob.cost(x, free)))
        exc = np.max(np.abs(x - x0)[free_log]) if free_log.any() else 0
        if verbose and k % 10 == 0: print(f"   tau={(k+1)*dtau:5.2f} lam_min={ev[0]:.2e} lam_max={ev[-1]:.2e} max|dlogp|={exc:.2f}")
        if ev[0] < lam_stop * max(ev[-1], 1e-30): reason = 'metric singular (boundary)'; break
        if exc > max_excursion: reason = f'parameter excursion > {max_excursion}'; break
        if not np.all(np.isfinite(x)): reason = 'non-finite'; break
        if np.linalg.norm(v) > 1e3: reason = 'velocity blow-up (boundary)'; break
    return path, reason
