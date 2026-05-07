import numpy as np

class BeamSolver:
    def __init__(self, L, EI_coeffs, AsG_coeffs, N):
        self.L = L
        self.N = N
        self.EI = np.zeros(N + 1)
        self.AsG = np.zeros(N + 1)
        self.EI[:len(EI_coeffs)] = EI_coeffs
        self.AsG[:len(AsG_coeffs)] = AsG_coeffs

    def _compute_recurrence(self, vi, ti, FYi, Mi):
        c_theta = np.zeros(self.N + 1)
        c_v = np.zeros(self.N + 1)
        c_theta[0], c_theta[1] = ti, -Mi / self.EI[0]
        
        if self.N >= 1:
            c_theta[2] = (-self.EI[1] * c_theta[1] + FYi) / (2 * self.EI[0])
            
        for k in range(2, self.N):
            sum_val = sum((k1 + 1) * c_theta[k1 + 1] * self.EI[k - k1] for k1 in range(k))
            c_theta[k + 1] = -sum_val / ((k + 1) * self.EI[0])

        c_v[0], c_v[1] = vi, c_theta[0] - FYi / self.AsG[0]

        for k in range(1, self.N):
            sub_sum = sum(((k1 + 1) * c_v[k1 + 1] - c_theta[k1]) * self.AsG[k - k1] for k1 in range(k))
            c_v[k + 1] = (c_theta[k] / (k + 1)) - sub_sum / ((k + 1) * self.AsG[0])

        return c_v, c_theta

    def solve(self):
        v_f, t_f = self._compute_recurrence(0, 0, 1, 0)
        v_m, t_m = self._compute_recurrence(0, 0, 0, 1)
        L_pow = self.L ** np.arange(self.N + 1)
        A = np.array([[np.sum(v_f * L_pow), np.sum(v_m * L_pow)], 
                      [np.sum(t_f * L_pow), np.sum(t_m * L_pow)]])

        cases = {"psi2": (1, 0, 0, 0), "psi3": (0, 1, 0, 0), "psi5": (0, 0, 1, 0), "psi6": (0, 0, 0, 1)}
        results = {}
        for name, (vi, ti, vj, tj) in cases.items():
            v_init, t_init = self._compute_recurrence(vi, ti, 0, 0)
            B = np.array([vj - np.sum(v_init * L_pow), tj - np.sum(t_init * L_pow)])
            sol = np.linalg.solve(A, B)
            results[name + "_v"] = v_init + sol[0] * v_f + sol[1] * v_m
            results[name + "_theta"] = t_init + sol[0] * t_f + sol[1] * t_m
            results[name + "_FYi"] = sol[0]
            results[name + "_Mi"] = sol[1]
        return results

    def get_stiffness_matrix(self, res):
        # K_ij = Reacción en i debido a desplazamiento unitario en j
        K = np.zeros((4, 4))
        # Grados de libertad: 0:v_i, 1:theta_i, 2:v_j, 3:theta_j
        names = ["psi2", "psi3", "psi5", "psi6"]
        for j, name in enumerate(names):
            # Reacciones en i (nodo inicial x=0)
            K[0, j] = -res[name + "_FYi"]
            K[1, j] = -res[name + "_Mi"]
            
            # Reacciones en j (nodo final x=L)
            # V(L) = AsG(L) * (v'(L) - theta(L))
            v_poly = res[name + "_v"]
            t_poly = res[name + "_theta"]
            dv_L = np.sum(np.polyder(v_poly[::-1]) * (self.L ** np.arange(self.N)))
            t_L = np.sum(t_poly * (self.L ** np.arange(self.N + 1)))
            EI_L = np.sum(self.EI * (self.L ** np.arange(self.N + 1)))
            AsG_L = np.sum(self.AsG * (self.L ** np.arange(self.N + 1)))
            
            K[2, j] = AsG_L * (dv_L - t_L) # Cortante en L
            # M(L) = EI(L) * theta'(L)
            dt_L = np.sum(np.polyder(t_poly[::-1]) * (self.L ** np.arange(self.N)))
            K[3, j] = EI_L * dt_L
        return K