import numpy as np

class BeamSolver:
    def __init__(self, L, EI_coeffs, AsG_coeffs, N):
        self.L = L
        self.N = N
        # Coeficientes de los polinomios (aseguramos longitud N+1)
        self.EI = np.zeros(N + 1)
        self.AsG = np.zeros(N + 1)
        self.EI[:len(EI_coeffs)] = EI_coeffs
        self.AsG[:len(AsG_coeffs)] = AsG_coeffs

    def _compute_recurrence(self, vi, ti, FYi, Mi):
        """Calcula los coeficientes de v y theta mediante recurrencia numérica"""
        c_theta = np.zeros(self.N + 1)
        c_v = np.zeros(self.N + 1)

        c_theta[0] = ti
        c_theta[1] = -Mi / self.EI[0]
        
        if self.N >= 1:
            c_theta[2] = (-self.EI[1] * c_theta[1] + FYi) / (2 * self.EI[0])
            
        for k in range(2, self.N):
            sum_val = 0
            for k1 in range(k):
                sum_val += (k1 + 1) * c_theta[k1 + 1] * self.EI[k - k1]
            c_theta[k + 1] = -sum_val / ((k + 1) * self.EI[0])

        c_v[0] = vi
        c_v[1] = c_theta[0] - FYi / self.AsG[0]

        for k in range(1, self.N):
            main_v = c_theta[k] / (k + 1)
            sub_sum = 0
            for k1 in range(k):
                sub_sum += ((k1 + 1) * c_v[k1 + 1] - c_theta[k1]) * self.AsG[k - k1]
            c_v[k + 1] = main_v - sub_sum / ((k + 1) * self.AsG[0])

        return c_v, c_theta

    def get_shape_functions(self):
        """Calcula las 4 funciones de forma (psi2, psi3, psi5, psi6)"""
        # Para resolver FYi y Mi, evaluamos la respuesta a FYi=1 y Mi=1
        # v(L) = v_hom(L) + FYi*v_f(L) + Mi*v_m(L)
        
        # 1. Respuesta a condiciones iniciales nulas pero fuerzas unitarias
        v_f, t_f = self._compute_recurrence(0, 0, 1, 0) # FYi=1
        v_m, t_m = self._compute_recurrence(0, 0, 0, 1) # Mi=1
        
        # Valores en L (usando Horner para evaluar polinomios)
        L_pow = self.L ** np.arange(self.N + 1)
        vf_L, tf_L = np.sum(v_f * L_pow), np.sum(t_f * L_pow)
        vm_L, tm_L = np.sum(v_m * L_pow), np.sum(t_m * L_pow)

        # Sistema de matrices para resolver FYi y Mi dado vi, ti, vj, tj
        A = np.array([[vf_L, vm_L], [tf_L, tm_L]])

        results = {}
        # Casos de carga unitaria para funciones de forma
        cases = {
            "psi2": (1, 0, 0, 0), # vi=1
            "psi3": (0, 1, 0, 0), # ti=1
            "psi5": (0, 0, 1, 0), # vj=1
            "psi6": (0, 0, 0, 1)  # tj=1
        }

        for name, (vi, ti, vj, tj) in cases.items():
            # Respuesta debida solo a vi, ti
            v_init, t_init = self._compute_recurrence(vi, ti, 0, 0)
            vi_L, ti_L = np.sum(v_init * L_pow), np.sum(t_init * L_pow)
            
            # B = [vj - v_init(L), tj - t_init(L)]
            B = np.array([vj - vi_L, tj - ti_L])
            sol = np.linalg.solve(A, B)
            
            # Coeficientes finales: base + (FYi_sol * base_FY) + (Mi_sol * base_Mi)
            final_v = v_init + sol[0] * v_f + sol[1] * v_m
            final_t = t_init + sol[0] * t_f + sol[1] * t_m
            results[name + "_v"] = final_v
            results[name + "_theta"] = final_t
            
        return results