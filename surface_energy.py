import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt

import pandas as pd

import subprocess

import utils
import numpy as np
from scipy.optimize import curve_fit
from lmfit import Model

from ctypes import *

from io import StringIO, BytesIO
import base64

#subprocess.call(['./setup.sh'])

# lib = cdll.LoadLibrary("numerical_integral.so")
# lib.numerical_integral.argtypes = [
#     c_double, c_double, c_double, c_double, c_double, c_double, c_double,
#     c_double, c_double, c_double, c_double, c_double, c_double, c_double,
#     c_double, c_double,
#     POINTER(c_double),
#     POINTER(c_double),
#     POINTER(c_double),
#     POINTER(c_double)
# ]


def theta_func(x, b, theta_c):
    # # Eqn 20.
    # c_numer = b * np.sqrt(x)
    # c_denom = 1. + c_numer
    # chem = (theta_c * c_numer) / (c_denom)
    # return chem
    return theta_c * b * x**(0.5) / (1 + b * x**(0.5))


def theta_two_func(x, b, theta_c, c, theta_p):
    c_numer = b * np.sqrt(x)
    c_denom = 1. + c_numer
    chem = (theta_c * c_numer) / (c_denom)
    p_numer = c * x
    p_denom = (1 - x) * (1 + (c - 1) * x)
    phys = (theta_p * p_numer) / p_denom
    theta = chem + phys
    return theta


def dH_func(x, a, b):
    return a * x + b


def dH_one_func(x, D, d, E, e, f):
    # Eqn 22. x is theta
    a = D * np.exp(-x / d)
    b = E * (f * x - x**2) * np.exp(-x / e)
    return a + b


def dH_two_func(x, a, b):
    return a * np.exp(-x / b)


class Simulator:
    # exp_dH_dtheta = data.exp_dH
    # exp_theta = data.exp_theta
    # exp_p_po = data.exp_p_po
    # BET = data.BET
    rho = 1 / 18.016 * 100**3
    heat = 44000
    temperature = 298
    entropy = 44000 / 298

    def __init__(self, data):
        self.exp_dH_dtheta = data['exp_dH']
        self.exp_theta = data['exp_theta']
        self.exp_p_po = data['exp_p_po']
        self.BET = data['BET']

    def load_cust_fit(self, param=None):
        if param == None:
            return

        if c in param:
            self.c = param['c']

        if theta_physi in param:
            self.theta_physi = param['theta_physi']

        if theta_chemi in param:
            self.theta_chemi = param['theta_chemi']

        if b in param:
            self.b = param['b']

        if D in param:
            self.D = param['D']

        if d in param:
            self.d = param['d']

    def fit_theta(self):
        """ THETA """
        temp_p_po = np.take(self.exp_p_po, np.where(self.exp_p_po < 0.05))[0]

        if np.size(temp_p_po) <= 5:
            temp_p_po = np.where(self.exp_p_po < 0.1)

        if np.size(temp_p_po) <= 5:
            temp_p_po = np.where(self.exp_p_po < 0.2)

        if np.size(temp_p_po) <= 5:
            print("uhoh")
            # return "Data is not fittable"

        points = np.arange(0, np.size(temp_p_po))
        temp_exp_theta = np.array(np.take(self.exp_theta, points))

        theta_model = Model(theta_func)
        theta_params = theta_model.make_params()
        theta_params['b'].set(value=1e1, min=0e0, max=1e3)
        theta_params['theta_c'].set(
            value=(self.BET / 68654.22 / 2),
            min=(self.BET / 68654.22 / 2) / 10,
            max=(self.BET / 68654.22 / 2) * 10)

        theta_opt = theta_model.fit(
            temp_exp_theta, x=temp_p_po, params=theta_params)

        theta_chemi = theta_opt.best_values['theta_c']
        self.b = theta_opt.best_values['b']

        theta_two_model = Model(theta_two_func)
        theta_two_params = theta_two_model.make_params()
        theta_two_params['b'].set(value=self.b, min=self.b-1e-7, max=self.b+1e-7)
        theta_two_params['c'].set(value=5, min=0e0, max=1000)
        theta_two_params['theta_c'].set(value=theta_chemi, min=theta_chemi-1e-7, max=theta_chemi+1e-7)
        theta_two_params['theta_p'].set(value=(self.BET / 68654.22), min=(self.BET / 68654.22)-1e-7, max=(self.BET / 68654.22)+1e-7) 
        theta_two_opt = theta_two_model.fit(
            self.exp_theta,
            x=self.exp_p_po,
            params=theta_two_params)

        self.c = theta_two_opt.best_values['c']
        self.theta_physi = theta_two_opt.best_values['theta_p']
        self.theta_chemi = theta_two_opt.best_values['theta_c']
        self.b = theta_two_opt.best_values['b']

        plt.clf()
        plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
        plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))

        plt.plot(
            self.exp_p_po,
            theta_two_func(self.exp_p_po, self.b, self.theta_chemi, self.c,
                           self.theta_physi))
        plt.scatter(self.exp_p_po, self.exp_theta)

        image = BytesIO()
        plt.savefig(image, format='png')
        image_str = base64.b64encode(image.getvalue()).decode()
        html = '''data:image/png;base64,%s''' % (image_str)

        return html

    def fit_dH(self, dH_fittype):
        print("Sodijfsoidjfsodifjsodifjsoidfjsodijfsdoifj:")
        """ dH """
        index = np.where(self.exp_p_po < 0.015)
        dH_opt, dH_cov = curve_fit(
            dH_func, np.array(np.take(self.exp_theta, index)[0]),
            np.log(-np.array(np.take(self.exp_dH_dtheta, index)[0])))

        self.D = -np.exp(dH_opt[1]) + 44000

        if dH_fittype:
            dH_two_opt, dH_two_cov = curve_fit(
                dH_two_func,
                self.exp_theta,
                self.exp_dH_dtheta + self.heat, [self.D - 44000, 1E-3],
                bounds=([-np.inf, 0], [0, 1]))

            dH_two_rsq = utils.r_sq(
                self.exp_theta, self.exp_dH_dtheta + self.heat,
                dH_two_func(self.exp_theta, dH_two_opt[0], dH_two_opt[1]))

            self.D = dH_two_opt[0]
            self.d = dH_two_opt[1]
            self.E = 0
            self.e = 1
            self.f = 1
        else:
        #    dH_one_opt, dH_one_cov = curve_fit(
        #        dH_one_func,
        #        self.exp_theta,
        #        self.exp_dH_dtheta + self.heat,
        #        [self.D - 44000, 1e-3, 1e12, 1e-3, 1e-5],
        #        bounds=([-np.inf, 0, 0, 0, 0], [0, 1, np.inf, 1, 1])) 
            dH_one_model = Model(dH_one_func)
            dH_one_param = dH_one_model.make_params()
            dH_one_param['D'].set(value=self.D - 44000, min=-np.inf, max=0)
            dH_one_param['d'].set(value=1e-3, min=0, max=1)
            dH_one_param['E'].set(value=1e12, min=0, max=np.inf)
            dH_one_param['e'].set(value=1e-3, min=0, max=1)
            dH_one_param['f'].set(value=1e-5, min=0, max=1)
            dH_one_fit = dH_one_model.fit(self.exp_dH_dtheta + self.heat, x=self.exp_theta, params=dH_one_param)
            dH_one_opt = dH_one_fit.best_values
            self.D = dH_one_opt['D']
            self.d = dH_one_opt['d']
            self.E = dH_one_opt['E']
            self.e = dH_one_opt['e']
            self.f = dH_one_opt['f']

        self.BET = self.theta_physi * 68654.22

        radius = np.sqrt(self.BET / (4 * np.pi))
        new_radius = radius + 1E-6
        volume = 4 / 3 * np.pi * radius**3
        new_volume = 4 / 3 * np.pi * new_radius**3
        delta_mols = self.rho * (new_volume - volume)
        new_SA = 4 * np.pi * new_radius**2

        self.A = (new_SA - self.BET) / delta_mols

        self.regression_data = self.D * np.exp(
            -self.exp_theta / self.d) + self.E * (
                self.f * self.exp_theta - self.exp_theta**2) * np.exp(
                    -self.exp_theta / self.e) - 44000
        self.gamma_error = np.sum(
            np.abs((self.exp_dH_dtheta - self.regression_data) /
                   self.regression_data)) / np.size(self.exp_dH_dtheta)

        plt.clf()
        plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
        plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
        if dH_fittype:
            plt.plot(self.exp_theta,
                     dH_two_func(self.exp_theta, dH_two_opt[0], dH_two_opt[1]))
        else:
            plt.plot(self.exp_theta,
                     dH_one_func(self.exp_theta, self.D, self.d, self.E, self.e, self.f))
        plt.scatter(self.exp_theta, self.exp_dH_dtheta + self.heat)

        image = BytesIO()
        plt.savefig(image, format='png')
        image_str = base64.b64encode(image.getvalue()).decode()
        html = '''data:image/png;base64,%s''' % (image_str)

        return html

    def solve_differential(self):
        dx = 1e-17
        k = 1
        x = np.linspace(1 - 5.8e-17, 1 - 1e-12, k + 1)

        p_po = (np.linspace(x[0], x[1], int(abs(x[1] - x[0]) / dx))).T
        diff_ugas = 8.314 * 298 / p_po
        theta = self.theta_chemi * self.b * p_po**(0.5) / (
            1 + self.b * p_po**(0.5)) + self.theta_physi * self.c * p_po / (
                (1 - p_po) * (p_po * (self.c - 1) + 1))

        diff_theta = self.theta_chemi * self.b / (
            1 + self.b * p_po**
            (0.5))**2 / 2 / p_po**(0.5) + self.theta_physi * self.c * (
                p_po**2 * (self.c - 1) + 1) / (
                    (1 - p_po)**2 * (p_po * (self.c - 1) + 1)**2)

        diff_diff_theta = -self.theta_chemi * self.b * (
            3 * self.b * p_po**
            (0.5) + 1) / (4 * p_po**(3 / 2) * (1 + self.b * p_po**(0.5))**
                          (3)) + 2 * self.theta_physi * self.c * (
                              p_po**3 * (self.c - 1)**2 + 3 * p_po *
                              (self.c - 1) + 2 - self.c) / (
                                  (1 - p_po)**3 * (p_po * (self.c - 1) + 1)**3)

        diff_diff_H = -self.D / self.d * np.exp(-theta / self.d) + self.E * (
            self.f - 2 * theta) * np.exp(-theta / self.e) - self.E / self.e * (
                self.f * theta - theta**2) * np.exp(-theta / self.e)

        SA = self.BET + self.A * theta
        diff_SA = self.A
        diff_diff_gamma = np.zeros(np.size(p_po))
        diff_gamma = np.zeros(np.size(p_po))
        gamma = np.zeros(np.size(p_po))

        gamma[0] = 0.072
        diff_gamma[0] = 0.

        for i in range(1, np.size(p_po)):
            print(
                '\x1b[2K First Integral: {:2.1%} complete'.format(
                    ((i / np.size(p_po)))),
                end='\r')
            diff_gamma[i] = diff_gamma[i - 1] - dx * diff_diff_gamma[i - 1]
            gamma[i] = gamma[i - 1] - dx * diff_gamma[i - 1]
            diff_diff_gamma[i] = diff_theta[i] / SA[i] * (
                diff_diff_H[i] * diff_theta[i] + diff_ugas[i] - diff_gamma[i] *
                (2 * diff_SA - SA[i] / theta[i] -
                 SA[i] * diff_diff_theta[i] / diff_theta[i]**2))

        temp_gamma = gamma[np.size(gamma) - 2]
        temp_diff_gamma = diff_gamma[np.size(diff_gamma) - 2]

        dx = 5e-10
        k = 10000
        x = np.linspace(1 - 1e-12, 0, k + 1)

        res_p_po = (c_double * k)()
        res_theta = (c_double * k)()
        res_gamma = (c_double * k)()
        res_uads = (c_double * k)()

        return {
            'theta_chemi': self.theta_chemi,
            'theta_physi': self.theta_physi,
            'c': self.c,
            'D': self.D,
            'd': self.d,
            'E': self.E,
            'e': self.e,
            'A': self.A,
            'BET': self.BET,
            'entropy': self.entropy,
            'heat': self.heat,
            'temperature': self.temperature,
            'b': self.b,
            'temp_gamma': temp_gamma,
            'temp_diff_gamma': temp_diff_gamma,
            'f': self.f,
        }

        lib.numerical_integral(
            self.theta_chemi, self.theta_physi, self.c, self.D, self.d, self.E,
            self.e, self.A, self.BET, self.entropy, self.heat,
            self.temperature, self.b, temp_gamma, temp_diff_gamma, self.f,
            res_p_po, res_theta, res_gamma, res_uads)

        fin_p_po = np.zeros(k)
        fin_theta = np.zeros(k)
        fin_gamma = np.zeros(k)
        fin_uads = np.zeros(k)

        for i in range(0, 10000):
            fin_p_po[i] = res_p_po[i]
            fin_theta[i] = res_theta[i]
            fin_gamma[i] = res_gamma[i]
            fin_uads[i] = res_uads[i]

        # plt.scatter(self.exp_p_po, self.exp_theta * 6.0221413E5 / 2.0, c='r')
        # plt.plot(fin_p_po, fin_theta)
        # plt.ylim(
        #     ymax=np.amax(self.exp_theta * 6.0221413E5 / 2.0) * 1.2, ymin=0)
        # plt.xlim(xmax=np.amax(self.exp_p_po) * 1.2, xmin=0.0)

        gamma = fin_gamma[k - 1]

        # image = BytesIO()
        # plt.savefig(image, format='png')
        # image_str = base64.b64encode(image.getvalue()).decode()

        d = {
            'p_po': fin_p_po,
            'theta': fin_theta,
            'gamma': fin_gamma,
            'uads': fin_uads
        }
        df = pd.DataFrame(data=d)
        out = StringIO()
        df.to_csv(out,index=False)
        return out.getvalue(), gamma
