import numpy as np
import multiprocessing

def saveRho(sp, rho, t, ti, x, Np) :

	Nx, Nt = (x.size, t.size)

	for n in range(len(sp)) :

		i = sp[n].xcell
		j = sp[n].pcell

		if (i >= 0 and i < Nx) and (j >= 0 and j < Np) :

			rho[i] += sp[n].weight

	if ti == 0 :
		with open('rho.dat', mode='w') as den :
			den.write( 'x\trho\n' )

	with open('rho.dat', mode='a') as den :
		for i in range(20) :
			den.write('-')

		den.write(f'\n# = {ti} of {Nt-1}\n')
		den.write('t = ' + str(t[ti]))
		den.write('\n\n')

		for i in range(Nx) :
			den.write(f'{x[i]:.3f}\t{rho[i]}\n')

		den.write('\n')

	return

def saveFw(sp, Fw, t, ti, x, p) :

	Nx, Np, Nt = (x.size, p.size, t.size)

	for n in range(len(sp)) :

		i = sp[n].xcell
		j = sp[n].pcell

		if (i >= 0 and i < Nx) and (j >= 0 and j < Np) :

			Fw[i, j] += sp[n].weight

	if ti == 0 :
		with open('Fw.dat', mode='w') as fw :
			fw.write('x\n')

			for i in range(Nx) :
				fw.write(f'{x[i]:.3f}\t')

			fw.write('\np\n')

			for i in range(Np) :
				fw.write(f'{p[i]:.3f}\t')

			fw.write('\nFw\n')

	with open('Fw.dat', mode='a') as fw :
		for i in range(20) :
			fw.write('-')

		fw.write(f'\n# = {ti} of {Nt-1}\n')
		fw.write('t = ' + str( t[ti] ))
		fw.write('\n\n')

		for j in range(Np - 1, -1, -1) :
			for i in range(Nx) :
				fw.write(f'{Fw[i, j]}\t')

			fw.write('\n')

		fw.write('\n\n')

	return

def savep(sp, t, ti, x, p) :

	Fw = np.zeros([x.size, p.size], dtype=int)
	rho = np.zeros(x.size, dtype=int)
	
	tRho = multiprocessing.Process(target=saveRho, args=(sp, rho, t, ti, x, p.size))
	tFw = multiprocessing.Process(target=saveFw, args=(sp, Fw, t, ti, x, p))
	
	tRho.start()
	tFw.start()

	tRho.join()
	tFw.join()

	return

def save(sp, t, ti, x, p) :

	Fw = np.zeros([x.size, p.size], dtype=int)
	rho = np.zeros(x.size, dtype=int)

	saveRho(sp, rho, t, ti, x, p.size)
	saveFw(sp, Fw, t, ti, x, p)
	
	return