import numpy as np
import math
import random
from sp import Sp

def propagate(sp, x0, dx, dt) :

	sp.x += sp.p / sp.m * dt
	sp.xcell = math.floor((sp.x - (x0 - dx / 2)) / dx)
	sp.time += dt

	return

def range_check(sp, Nx, Np) :

	if (sp.xcell >= 0 and sp.xcell < Nx) or (sp.pcell >= 0 and sp.pcell < Np) :
		
		return True

	else :

		return False

def create_particles(sp, n, currentParticles, newParticles, Nmax, dt, p0, dp, Vw, cutoff, vwh, Gamma) :

	if random.random() < (Gamma * dt) :
	
		ps = random.uniform(-cutoff, cutoff)

		if random.random() < math.fabs(Vw(sp[n].x, ps)) / vwh :
		
			currentParticles += 2
			newParticles += 2

			if currentParticles > Nmax :
			
				print('The number of particles has reached the limit.')
				exit(0)

			sp.append(None)
			sp[currentParticles-2] = Sp(currentParticles-2)
			sp[currentParticles-2].m = sp[n].m
			sp[currentParticles-2].x = sp[n].x
			sp[currentParticles-2].xcell = sp[n].xcell
			sp[currentParticles-2].p = sp[n].p + ps
			sp[currentParticles-2].pcell = math.floor((sp[currentParticles-2].p - (p0 - dp / 2)) / dp)
			sp[currentParticles-2].weight = sp[n].weight * np.sign(Vw(sp[n].x, ps))
			sp[currentParticles-2].time = sp[n].time
			sp[currentParticles-2].status = True

			sp.append(None)
			sp[currentParticles-1] = Sp(currentParticles-1)
			sp[currentParticles-1].m = sp[n].m
			sp[currentParticles-1].x = sp[n].x
			sp[currentParticles-1].xcell = sp[n].xcell
			sp[currentParticles-1].p = sp[n].p - ps
			sp[currentParticles-1].pcell = math.floor((sp[currentParticles-1].p - (p0 - dp / 2)) / dp)
			sp[currentParticles-1].weight = -sp[n].weight * np.sign(Vw(sp[n].x, ps))
			sp[currentParticles-1].time = sp[n].time
			sp[currentParticles-1].status = True

	return (currentParticles, newParticles)

def wmc(sp, Nmax, x, p, dt, Vw, cutoff, vwh, Gamma) :

	print('Evolving particles...')
	
	Nx, Np, x0, dx, p0, dp = (x.size, p.size, x[0], x[1] - x[0], p[0], p[1] - p[0])

	originalParticles = len(sp)
	
	currentParticles, newParticles, outsideParticles = originalParticles, 0, 0

	for n in range(originalParticles) :

		if sp[n].status == False :

			outsideParticles += 1
			continue
		
		propagate(sp[n], x0, dx, dt)
		
		if range_check(sp[n], Nx, Np) == False :

			sp[n].status = False
			outsideParticles += 1
			continue

		currentParticles, newParticles = create_particles(sp, n, currentParticles, newParticles, Nmax, dt, p0, dp, Vw, cutoff, vwh, Gamma)

	print(f'Number of created particles = {newParticles}')
	print(f'Total number of particles = {currentParticles}')

	for n in range(originalParticles, currentParticles) :

		if range_check(sp[n], Nx, Np) == False:
		
			sp[n].status = False
			outsideParticles += 1

	print(f'Number of outside particles = {outsideParticles}')

	return currentParticles
