import numpy as np
import random
from functools import cmp_to_key
from sp import Sp

def density(sp, Nx, Np, Fw) :

	for n in range(len(sp)) :

		i = sp[n].xcell
		j = sp[n].pcell

		if (i >= 0 and i < Nx) and (j >= 0 and j < Np) :

			Fw[i, j] += sp[n].weight

	return

def comp(sp1, sp2) :

	if sp1.status == True and sp2.status == False :
		return -1
	
	elif sp1.status == False and sp2.status == False :
		return 0
	
	elif sp1.status == True and sp2.status == True :

		if sp1.xcell < sp2.xcell :
			return -1
		
		elif sp1.xcell > sp2.xcell :
			return 1
		
		elif sp1.pcell < sp2.pcell :
			return -1
		
		elif sp1.pcell > sp2.pcell :
			return 1
		
		else :
			return 0
		
	else :
		return 1

def check_range(sp, Nx, Np) :

	for n in range(len(sp)) :

		if (sp[n].status == False) :
			continue

		else :
			if not (sp[n].xcell >= 0 and sp[n].xcell < Nx) or not (sp[n].pcell >= 0 and sp[n].pcell < Np) :

				sp[n].status = False

	return

def sparse_matrix(sp, sparse) :

	if sp[0].status != False :

		sparse.append([sp[0].xcell, sp[0].pcell, sp[0].weight])
	
	for i in range(1, len(sp)) :

		if sp[i].status == False :
			continue
		
		elif sp[i].xcell != sp[i - 1].xcell :

			if sparse[-1][2] == 0 :

				sparse.pop()

			sparse.append([sp[i].xcell, sp[i].pcell, sp[i].weight])
		
		elif sp[i].pcell != sp[i - 1].pcell :

			if sparse[-1][2] == 0 :

				sparse.pop()

			sparse.append([sp[i].xcell, sp[i].pcell, sp[i].weight])
		
		else :

			sparse[-1][2] += sp[i].weight

def annihilation(sp, x, p) :
	'''grid, O(n)'''

	print('\nAnnihilating particles...')
	
	Nx, Np, dx, dp = (x.size, p.size, x[1] - x[0], p[1] - p[0])

	print(f'# of particles before annihilation = {len(sp)}')
	
	Fw = np.zeros([x.size, p.size], dtype=int)
	
	density(sp, Nx, Np, Fw)

	count = 0

	for i in range(Nx) :
		for j in range(Np) :
			for k in range(np.absolute(Fw[i, j]).astype(int)) :

				sp[count].number = count

				sp[count].xcell = i
				sp[count].x = x[i] + (random.random() - 0.5) * dx
				sp[count].pcell = j
				sp[count].p = p[j] + (random.random() - 0.5) * dp
				sp[count].weight = np.sign(Fw[i, j])

				sp[count].status = True

				count += 1

	for i in range(count, len(sp)) :

		sp[i] = None

	while sp[-1] == None :

		sp.pop()
					
	print(f'# of particles after annihilation  = {count}\n')

	return sp

def annihilation_a(sp, x, p) :
	'''sort, O(n log n)'''

	print('\nAnnihilating particles...')

	Nx, Np, dx, dp = (x.size, p.size, x[1] - x[0], p[1] - p[0])

	print(f'# of particles before annihilation = {len(sp)}')

	check_range(sp, Nx, Np)

	sp.sort(key=cmp_to_key(comp))

	sparse = []

	sparse_matrix(sp, sparse)

	count = 0

	m = sp[0].m
	t = sp[0].time

	for i in range(len(sp)) :

		sp[i] = None

	for i in range(len(sparse)) :
		for n in range(np.absolute(sparse[i][2]).astype(int)) :

			sp[count] = Sp(count)

			sp[count].m = m
			sp[count].xcell = sparse[i][0]
			sp[count].x = x[sparse[i][0]] + (random.random() - 0.5) * dx
			sp[count].pcell = sparse[i][1]
			sp[count].p = p[sparse[i][1]] + (random.random() - 0.5) * dp
			sp[count].weight = np.sign(sparse[i][2])
			sp[count].time = t

			sp[count].status = True

			count += 1

	for i in range(count, len(sp)) :

		sp[i] = None

	while sp[-1] == None :

		sp.pop()

	print(f'# of particles after annihilation  = {len(sp)}\n')

	return sp
