from flask import Flask, jsonify, render_template, request
from io import BytesIO
import base64

from surface_energy import Simulator
import numpy as np
import data

import pandas as pd

application = Flask(__name__)


@application.route('/_load_data')
def load_data():
    BET = request.args.get('BET', 0, type=float)

    str_theta = request.args.get('exp_theta', "0", type=str)
    str_p_po = request.args.get('exp_p_po', "0", type=str)
    str_dH = request.args.get('exp_dH', "0", type=str)

    exp_theta = np.array(np.fromstring(str_theta, sep='\n'))
    exp_p_po = np.array(np.fromstring(str_p_po, sep='\n'))
    exp_dH = np.array(np.fromstring(str_dH, sep='\n'))

    if (exp_theta.size == exp_p_po.size == exp_dH.size):
        return jsonify(result="VALID")
    else:
        return jsonify(result="ARRAY LENGTH MUST MATCH")


@application.route('/_fit_theta')
def fit_theta():
    BET = request.args.get('BET', 0, type=float)

    str_theta = request.args.get('exp_theta', "0", type=str)
    str_p_po = request.args.get('exp_p_po', "0", type=str)
    str_dH = request.args.get('exp_dH', "0", type=str)

    exp_theta = np.array(np.fromstring(str_theta, sep='\n'))
    exp_p_po = np.array(np.fromstring(str_p_po, sep='\n'))
    exp_dH = np.array(np.fromstring(str_dH, sep='\n'))

    d = {
        'exp_theta': exp_theta,
        'exp_p_po': exp_p_po,
        'exp_dH': exp_dH,
        'BET': BET
    }

    sim = Simulator(d)
    res = sim.fit_theta()

    ret = {
        "result": res,
        "fit_params": {
            "theta_chemi": sim.theta_chemi,
            "theta_physi": sim.theta_physi,
            "b": sim.b,
            "c": sim.c
        }
    }
    return jsonify(ret)


@application.route('/_fit_dh')
def fit_dh():
    BET = request.args.get('BET', 0, type=float)

    str_theta = request.args.get('exp_theta', "0", type=str)
    str_p_po = request.args.get('exp_p_po', "0", type=str)
    str_dH = request.args.get('exp_dH', "0", type=str)
    dh_method = request.args.get('dh_method', 1, type=int)
    dh_method = bool(dh_method)

    exp_theta = np.array(np.fromstring(str_theta, sep='\n'))
    exp_p_po = np.array(np.fromstring(str_p_po, sep='\n'))
    exp_dH = np.array(np.fromstring(str_dH, sep='\n'))

    d = {
        'exp_theta': exp_theta,
        'exp_p_po': exp_p_po,
        'exp_dH': exp_dH,
        'BET': BET
    }

    sim = Simulator(d)
    sim.fit_theta()
    res = sim.fit_dH(dh_method)
    ret = {
        "result": res,
        "fit_params": {
            "theta_chemi": sim.theta_chemi,
            "theta_physi": sim.theta_physi,
            "b": sim.b,
            "c": sim.c,
            "A": sim.A,
            "BET": sim.BET,
            "D": sim.D,
            "d": sim.d,
            "E": sim.E,
            "e": sim.e,
            "f": sim.f,
            "gamma_error": sim.gamma_error,
        }
    }
    return jsonify(ret)

@application.route('/_solve_diff')
def solve_diff():
    BET = request.args.get('BET', 0, type=float)

    str_theta = request.args.get('exp_theta', "0", type=str)
    str_p_po = request.args.get('exp_p_po', "0", type=str)
    str_dH = request.args.get('exp_dH', "0", type=str)
    dh_method = request.args.get('dh_method', 1, type=int)
    dh_method = bool(dh_method)

    exp_theta = np.array(np.fromstring(str_theta, sep='\n'))
    exp_p_po = np.array(np.fromstring(str_p_po, sep='\n'))
    exp_dH = np.array(np.fromstring(str_dH, sep='\n'))

    d = {
        'exp_theta': exp_theta,
        'exp_p_po': exp_p_po,
        'exp_dH': exp_dH,
        'BET': BET
    }

    sim = Simulator(d)
    sim.fit_theta()
    sim.fit_dH(dh_method)
    res = sim.solve_differential()
    return jsonify(result=res)


@application.route('/')
def index():
    return render_template('index.html')


if __name__ == "__main__":
    application.run(port=80)
