import numpy as np
import math
import random
import time
import datetime
from sp import Sp
from init import init
from save import save
from wmc_muscato import wmc
from annihilation import annihilation, annihilation_a

Nmax = 2_000_000

hbar = 6.62607015 / 1.602176634 / 2 / math.pi
m = 0.32
e = 1

Nx = 200
Np = 400
Nt = 400

xrange = np.array([-30, 30])
prange = np.array([-10, 10]) * hbar

initialParticles = 160_000

# fanni = 100

dx = (xrange[1] - xrange[0]) / Nx
dp = (prange[1] - prange[0]) / Np
dt = 0.05

x0 = -15
s0 = 2.852
p0 = 0.7 * hbar
xc = 0
a = 0.3
sigma = 1

def wFunction(x, p):

	return 1 / math.pi / hbar * np.exp(-0.5 * (((x - x0) / s0 )**2)) * np.exp(-2 * (s0 * (p - p0) / hbar)**2)

def potential(x):

	return e * a * np.exp(-0.5 * (x - xc)**2 / sigma**2)

def wKernel(x, p):

	return 2 * a * sigma * math.sqrt(2 * math.pi) / math.pi / hbar**2 * np.exp(-2 * (sigma * p / hbar)**2) * np.sin(2 * p * x / hbar)

cutoff = 8 * hbar
vwh = 2 * a * sigma * math.sqrt(2 * math.pi) / math.pi / hbar**2
gamma = vwh * cutoff

def main():

	random.seed(12345)
	
	sp = []

	x = np.arange(xrange[0] + dx / 2, xrange[1], dx)
	p = np.arange(prange[0] + dp / 2, prange[1], dp)
	t = np.arange(0, (Nt + 1) * dt, dt)

	print(f'Maximum number of particles allowed = {Nmax}\n')

	t0 = time.perf_counter()

	numberOfParticles = init(sp, Nmax, m, x, p, wFunction, initialParticles)

	save(sp, t, 0, x, p)

	t1 = t0
	t2 = time.perf_counter()
	print('Elapsed time is ' + str(t2 - t1) + ' seconds.')
	print('Total elapsed time is ' + str(datetime.timedelta(seconds=t2-t0)), flush=True)

	particleTime = 0

	for i in range(1, Nt + 1) :

		particleTime += dt

		print(f'\nStep {i} of {Nt} --- time = {particleTime}\n')

		numberOfParticles = wmc(sp, Nmax, x, p, dt, wKernel, cutoff, vwh, gamma)

		if numberOfParticles > 480_000 :
			
			sp = annihilation(sp, x, p)
		
		save(sp, t, i, x, p)

		t1 = t2
		t2 = time.perf_counter()
		print('Elapsed time is ' + str(t2 - t1) + ' seconds.')
		print('Total elapsed time is ' + str(datetime.timedelta(seconds=t2-t0)), flush=True)
	
	print('\nOutput files saved\n')

if __name__ == "__main__":
    main()
