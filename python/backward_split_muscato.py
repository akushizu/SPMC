# %%
import numpy as np
import math
import time
from scipy import sparse, fft
import matplotlib.pyplot as plt
import matplotlib.animation as animation

#%matplotlib widget

# %%

hb = 6.62607015 / 1.602176634 / 2 / math.pi
m = 0.32
e = 1
k0 = 0.7
s0 = 2.852
x0 = -15
xc = 0
a = 0.3
s = 1

#--------step size----------------------
dx = 0.005
dt = 0.001
l = (0.5 * 1j * hb / m) * (dt / dx**2)

x = np.arange(-40, 40 + dx, dx)
t = np.arange(0, 20 + dt, dt)

#------initial & boundaty condition------------------
w = np.sqrt(np.sqrt(1 / (2 * math.pi * s0**2))) * np.exp(-(x - x0)**2 / (4 * s0**2)) * np.exp(1j * k0 * (x - xc))

N = 1 #np.sqrt( dx * np.sum( np.abs(w)**2 ) )

w = w / N

#---------potential------------------------
Vw = e * a * np.exp(-(x - xc)**2 / (2 * s**2))

#---------backward difference------------------------
W = sparse.diags(( 1 + 2 * l ) * np.ones(len(x)) + (1j * dt / hb) * Vw.conj().T, 0, format='csc')
W += sparse.diags(-l * np.ones(len(x) - 1), 1, format='csc')
W += sparse.diags(-l * np.ones(len(x) - 1), -1, format='csc')

adata = [{} for i in range(int((len(t) - 1) / 50) + 1)]
adata[0] = {'t': t[0], 'x': x, 'Amp': w, 'time': 0}

t0 = time.perf_counter()

for i in range(1, len(t)):

    w1 = sparse.linalg.spsolve(W.T, w) # w1 * W = w
    N = 1 #np.sqrt(dx * np.sum(abs(w1)**2))
    w = w1 / N
    
    if i % 50 == 0:
        
        adata[int(i / 50)] = {'t': t[i], 'x': x, 'Amp': w, 'time': str(time.perf_counter() - t0)}
        

print('Elapsed time is ' + str(time.perf_counter() - t0) + ' seconds.')

#----------save----------------------------
np.savez('adata', adata=adata)

#------------------------------------------

# %%

x = adata[1]['x']
t = [adata[i]['t'] for i in range(len(adata))]

col = plt.rcParams['axes.prop_cycle'].by_key()['color']

px = 1 / plt.rcParams['figure.dpi']
fig = plt.figure(figsize=(1120*px, 840*px), dpi=150, layout='constrained')

v = animation.FFMpegWriter(fps=15, codec='ffv1')

with v.saving(fig, 'wob.avi', 150):
    for i in range(len(t)):
        fig.clf()

        ax = fig.add_subplot(2, 2, 1)
        ax.plot(x, 0.5*Vw, '--', color=col[1])
        ax.plot(x, np.abs(adata[i]['Amp'])**2, color=col[0])
        ax.set_xlim([x[0], x[-1]])
        ax.set_xlabel('x')
        ax.set_title('den, t = ' + str(t[i]) + 'fs')
        ax.xaxis.label.set_size(16)
        ax.yaxis.label.set_size(16)
        ax.title.set_size(16)
        ax.tick_params(axis='both', labelsize=16)
        
        ax = fig.add_subplot(2, 2, 2)
        ax.plot(x, Vw, '--', color=col[1])
        ax.plot(x, np.abs(adata[i]['Amp']), color=col[0])
        ax.set_xlim([x[0], x[-1]])
        ax.set_xlabel('x')
        ax.set_title('abs')
        ax.xaxis.label.set_size(16)
        ax.yaxis.label.set_size(16)
        ax.title.set_size(16)
        ax.tick_params(axis='both', labelsize=16)
        
        ax = fig.add_subplot(2, 2, 3)
        ax.plot(x, Vw, '--', color=col[1])
        ax.plot(x, np.abs(adata[i]['Amp']), 'k:', x , -np.abs(adata[i]['Amp']), 'k:')
        ax.plot(x, np.real(adata[i]['Amp']), color=col[0])
        ax.set_xlim([x[0], x[-1]])
        ax.set_xlabel('x')
        ax.set_title('real')
        ax.xaxis.label.set_size(16)
        ax.yaxis.label.set_size(16)
        ax.title.set_size(16)
        ax.tick_params(axis='both', labelsize=16)
        
        ax = fig.add_subplot(2, 2, 4)
        ax.plot(x, Vw , '--', color=col[1])
        ax.plot(x, np.abs(adata[i]['Amp']), 'k:', x, -np.abs(adata[i]['Amp']), 'k:')
        ax.plot(x, np.imag(adata[i]['Amp']), color=col[0])
        ax.set_xlim([x[0], x[-1]])
        ax.set_xlabel('x')
        ax.set_title('imag')
        ax.xaxis.label.set_size(16)
        ax.yaxis.label.set_size(16)
        ax.title.set_size(16)
        ax.tick_params(axis='both', labelsize=16)

        v.grab_frame()
        
# %%

hb = 6.62607015 / 1.602176634 / 2 / math.pi
m = 0.32
e = 1
k0 = 0.7
s0 = 2.852
x0 = -15
xc = 0
a = 0.3
s = 1

#--------step size----------------------
N = 16384
dt = 0.001

x = np.linspace(-40, 40, N)
dx = x[1] - x[0]
k = np.linspace(-math.pi / dx, math.pi / dx, N)
t = np.arange(0, 20 + dt, dt)

#------initial & boundaty condition------------------
psi = math.sqrt(math.sqrt(1 / (2 * math.pi * s0**2))) * np.exp(-(x - x0)**2 / (4 * s0**2)) * np.exp(1j * k0 * (x - xc))
V = e * a * np.exp(-(x - xc)**2 / (2 * s**2))

Ut = np.exp(-1j * (hb**2 * k**2 / m / 2) * dt / hb)
Uv = np.exp(-1j * V * dt / 2 / hb)

edata = [{} for i in range(int((len(t) - 1) / 50) + 1)]
edata[0] = {'t': t[0], 'x': x, 'Amp': psi, 'time': 0}

t0 = time.perf_counter()

for i in range(1, len(t)):
    
    psi = Uv * psi
    psih = fft.fftshift( fft.fft(psi) )
    psih = Ut * psih
    psi = fft.ifft( fft.ifftshift(psih) )
    psi = Uv * psi
    
    if i % 50 == 0:
        
        edata[int(i / 50)] = {'t': t[i], 'x': x, 'Amp': psi, 'time': str(time.perf_counter() - t0)}
        
print('Elapsed time is ' + str(time.perf_counter() - t0) + ' seconds.')

#----------save----------------------------
np.savez('edata', edata=edata)

#------------------------------------------

# %%

x = edata[1]['x']
t = [edata[i]['t'] for i in range(len(edata))]

col = plt.rcParams['axes.prop_cycle'].by_key()['color']

px = 1 / plt.rcParams['figure.dpi']
fig = plt.figure(figsize=(1120*px, 840*px), dpi=150, layout='tight')

v = animation.FFMpegWriter(fps=15, codec='ffv1')

with v.saving(fig, 'wos.avi', 150):
    for i in range(len(t)):
        fig.clf()

        ax = fig.add_subplot(2, 2, 1)
        ax.plot(x, 0.5*V, '--', color=col[1])
        ax.plot(x, np.abs(edata[i]['Amp'])**2, color=col[0])
        ax.set_xlim([x[0], x[-1]])
        ax.set_xlabel('x')
        ax.set_title('den, t = ' + str(t[i]) + 'fs')
        ax.xaxis.label.set_size(16)
        ax.yaxis.label.set_size(16)
        ax.title.set_size(16)
        ax.tick_params(axis='both', labelsize=16)
        
        ax = fig.add_subplot(2, 2, 2)
        ax.plot(x, V, '--', color=col[1])
        ax.plot(x, np.abs(edata[i]['Amp']), color=col[0])
        ax.set_xlim([x[0], x[-1]])
        ax.set_xlabel('x')
        ax.set_title('abs')
        ax.xaxis.label.set_size(16)
        ax.yaxis.label.set_size(16)
        ax.title.set_size(16)
        ax.tick_params(axis='both', labelsize=16)
        
        ax = fig.add_subplot(2, 2, 3)
        ax.plot(x, V, '--', color=col[1])
        ax.plot(x, np.abs(edata[i]['Amp']), 'k:', x, -np.abs(edata[i]['Amp']), 'k:')
        ax.plot(x, np.real(edata[i]['Amp']), color=col[0])
        ax.set_xlim([x[0], x[-1]])
        ax.set_xlabel('x')
        ax.set_title('real')
        ax.xaxis.label.set_size(16)
        ax.yaxis.label.set_size(16)
        ax.title.set_size(16)
        ax.tick_params(axis='both', labelsize=16)
        
        ax = fig.add_subplot(2, 2, 4)
        ax.plot(x, V, '--', color=col[1])
        ax.plot(x, np.abs(edata[i]['Amp']), 'k:', x, -np.abs(edata[i]['Amp']), 'k:')
        ax.plot(x, np.imag(edata[i]['Amp']), color=col[0])
        ax.set_xlim([x[0], x[-1]])
        ax.set_xlabel('x')
        ax.set_title( 'imag' )
        ax.xaxis.label.set_size(16)
        ax.yaxis.label.set_size(16)
        ax.title.set_size(16)
        ax.tick_params(axis='both', labelsize=16)
        
        v.grab_frame()

# %%

data = np.load('adata.npz', allow_pickle=True)
adata = data['adata']
data = np.load('edata.npz', allow_pickle=True)
edata = data['edata']

e = 1
xc = 0
a = 0.12
s = 1

x = np.arange(-40, 40 + 0.005, 0.005)
V = e * a * np.exp(-(x - xc)**2 / (2 * s**2))

px = 1 / plt.rcParams['figure.dpi']
fig = plt.figure(figsize=(560*px, 420*px), dpi=150)

col = plt.rcParams['axes.prop_cycle'].by_key()['color']

v = animation.FFMpegWriter(fps=15, codec='ffv1')

with v.saving(fig, 'woc.avi', 150):
    for i in range(len(adata)):
        fig.clf()
    
        ax = fig.add_subplot()
        ax.plot(x, 0.5*V, '--', color=col[8])
        p1, = ax.plot(adata[i]['x'], abs(adata[i]['Amp'])**2, color=col[0])
        p2, = ax.plot(edata[i]['x'], abs(edata[i]['Amp'])**2, color=col[1])
        ax.set_xlim([x[0], x[-1]])
        ax.set_xlabel('x')
        ax.set_title('den, t = ' + str(adata[i]['t']) + 'fs')
        ax.legend([p1, p2], ['backward difference', 'split operator'])

        v.grab_frame()

