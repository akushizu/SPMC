import numpy as np
import random
from sp import Sp

def init(sp, Nmax, m, x, p, fw, number) :
		
	print('Initializing fw...')

	Nx, Np, dx, dp = (x.size, p.size, x[1] - x[0], p[1] - p[0])
	
	actualNumber = 0
	count = 0

	for i in range(Nx) :
		for j in range(Np) :

			Fw = np.around(number * fw(x[i], p[j]) * dx * dp).astype(int)

			actualNumber += np.absolute(Fw).astype(int)

			if number > Nmax :

				print('The number of particles has reached the limit.')
				exit(0)

			for k in range(np.absolute(Fw)) :

				sp.append(None)
				sp[count] = Sp(count)

				sp[count].m = m
				sp[count].xcell = i
				sp[count].x = x[i] + (random.random() - 0.5) * dx
				sp[count].pcell = j
				sp[count].p = p[j] + (random.random() - 0.5) * dp
				sp[count].weight = np.sign(Fw)
				sp[count].time = 0

				sp[count].status = True

				count += 1

	if count == 0 :

		print('There are no particles in the system.')
		exit(0)

	print( f'Initial number of particles = {count}' )

	return count