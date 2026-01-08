import numpy as np
from scipy.integrate import quad
from scipy.special import erfc
from colossus.cosmology import cosmology

class HaloMassFunction:
    def __init__(self, omega_m=0.3089, z=0, mdef='m200b'):
        self.omega_m = omega_m
        self.z = z
        self.rho_crit = 277536627245.708  # M_sun / (h Mpc)^3
        self.rho_m = omega_m * self.rho_crit
        self.cosmo = cosmology.setCosmology('planck15')
        print(self.cosmo)
        self.D0 = self.D_unnormalized(0.0)
        self.pk_table_path = 'power_uchuu_log10.dat'
        self.ps_args = dict(model='uchuu_table', path=self.pk_table_path) if self.pk_table_path else dict(model='camb')
        self.mdef = mdef
        
        if self.mdef == 'm200b':
            self.aa, self.bb, self.DD, self.EE, self.FF = 1.089, 0.652, 1.0, 0.17, 0.087

        else:
            raise ValueError(f"mdef '{self.mdef}' not valid. Usa 'm200b'")
    

    def RtoM(self, M):
        return (3 * M / (4 * np.pi * self.omega_m * self.rho_crit))**(1/3)

    def E(self, z):
        return np.sqrt(self.omega_m * (1 + z)**3 + (1 - self.omega_m))

    def D_unnormalized(self, z):
        integral, _ = quad(lambda zp: (1 + zp) / (self.E(zp)**3), z, np.inf)
        return (5 * self.omega_m * self.E(z) / 2) * integral

    def sigma(self, M):
        M = np.atleast_1d(M)
        R = self.RtoM(M)
        sigma_std = self.cosmo.sigma(R, self.z, ps_args=self.ps_args)
        x = sigma_std / 1.676
        U2 = (-0.00221 * x**3 + 0.03835 * x**2 + 0.17810 * x - 0.01507)**2
        sigma_mod = np.sqrt(sigma_std**2 + U2)
        return sigma_mod[0] if np.isscalar(M) else sigma_mod

    def b(self, m_val):
        m = np.array([1e16, 1e15, 1e14, 6.5e10, 1e10, 1e9, 1e8, 1e7, 1e6])
        b = np.array([0.5259, 0.415, 0.328, 0.1764, 0.1552, 0.1308, 0.1179, 0.1045, 0.094])
        coeffs = np.polyfit(np.log10(m), np.log10(b), 4)
        return 10**np.polyval(coeffs, np.log10(m_val))

    def c(self, m_val):
        m = np.array([3e15, 3e14, 3e13, 3e12, 3e11, 3e10, 3e9, 3e8, 3e7, 3e6, 1e10, 1e9, 1e8, 1e7, 1e6])
        b = np.array([0.613, 0.474, 0.373, 0.301, 0.249, 0.209, 0.1794, 0.1560, 0.1355, 0.1223, 0.1942, 0.168, 0.1466, 0.1298, 0.1161])
        coeffs = np.polyfit(np.log10(m), np.log10(b), 4)
        return 10**np.polyval(coeffs, np.log10(m_val))

    def F(self, m_array):
        m_array = np.atleast_1d(m_array)  
        R = self.RtoM(m_array)
        b_val = self.b(m_array)
        sig = self.sigma(m_array)
        x = self.cosmo.sigma(R, self.z, ps_args=self.ps_args) / 1.676

        term1 = (1 + 0.845*x - 0.04*x**2 + 0.0025*x**3)**self.bb
        term2 = self.aa * 1.365 * (1 + self.EE * b_val - self.FF * b_val**2)**self.DD
        delta_c = term1 * term2

        c_m = self.c(m_array)
        cte = delta_c / (np.sqrt(2) * sig)

        xi = np.linspace(0, 1, 1000) 
        xi2 = xi**2
        xi_mat = xi[np.newaxis, :]
        c_m_mat = c_m[:, np.newaxis]
        cte_mat = cte[:, np.newaxis]
        integrand = erfc(cte_mat * np.sqrt((1 - np.exp(-c_m_mat * xi_mat**2)) / (1 + np.exp(-c_m_mat * xi_mat**2)))) * xi2
        I = np.trapz(integrand, xi, axis=1)
        V = 3 * I

        F_val = erfc(0.98 * cte) / V
        return F_val if F_val.size > 1 else F_val[0]

    def n0(self, m): ##this gives dn/dlnM
        s = 0.01
        Fm = self.F(m)
        Fm_s = self.F((1 + s) * m)
        der = (Fm - Fm_s) / s
        return der * self.rho_m / (m * (1 + s/2))
