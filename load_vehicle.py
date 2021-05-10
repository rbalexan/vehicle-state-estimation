import numpy as np

def load_vehicle():
	veh = {
	"m": 1776,
	"g": 9.81,
	"rho": 1.225,
	"Iz": 2763.5,
	"a": 1.264,
	"b": 1.367,
	"L": 2.631,
	"Wf": 9052.3,
	"Wr": 8370.2,
	"hcg": .55,
	"cdA": .594,
	"frr": 0.015,
	"Fxmax_engine": 5500,
	"delta_limit": 0.5236,
	"brake_delay": 0.01,
	"steer_delay": 0.07,
	"brake_prop_front": 0.68,
	"brake_prop_rear": 0.32,
	"Re": 0.318}

	ftire = {
	"Ca_lin": 80000,
	"Cy": 120000,
	"mu_s": 0.9,
	"mu": 0.9}

	rtire = {
	"Ca_lin": 120000,
	"Cy": 200000,
	"mu_s": 0.91,
	"mu": 0.91}

	return veh, ftire, rtire


	