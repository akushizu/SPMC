# %%
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib.cm import ScalarMappable

#%matplotlib widget

# %%

data = np.load('adata.npz', allow_pickle=True)
adata = data['adata']

v = animation.FFMpegWriter(fps=15, codec='ffv1')

px = 1 / plt.rcParams['figure.dpi']
fig = plt.figure(figsize=(560*px, 420*px), dpi=150)

with v.saving(fig, 'muscato_rho.avi', 100):
    with open('rho.dat', mode='r') as f:

        f.readline()

        N = 0

        st = f.readline()

        while st:
            _, _, i, _, n = f.readline().split()
            i = int(i)

            _, _, t = f.readline().split()
            t = float(t)

            f.readline()

            st = f.readline()
            C = np.array([], dtype=float)

            while st != '\n':
                C = np.vstack((C, np.asarray(st.split()).astype(float))) if C.size else np.asarray(st.split()).astype(float)

                st = f.readline()

            st = f.readline()

            x = C[:, 0]
            rho = C[:, 1]

            if i == 0:
                N = rho.sum()

            rho = rho / N / ( x[1] - x[0] )

            fig.clf()
            ax = fig.subplots()

            l1, = ax.plot(x, rho)
            ax.plot(adata[i]['x'], 0.05 * np.exp( -(adata[i]['x'] - 0)**2 / 2 / 1**2 ))
            l2, = ax.plot(adata[i]['x'], np.abs(adata[i]['Amp'])**2)
            ax.set_title(f'rho, t = {t} fs', fontname='DejaVu Sans', fontsize=12, fontweight='bold')
            ax.set_xlabel('x (nm)', fontname='DejaVu Sans', fontsize=12)
            ax.set_ylabel('density', fontname='DejaVu Sans', fontsize=12)
            ax.legend([l1, l2], ['Monte Carlo', 'backward difference'])
            ax.set_xlim(-30, 30)
            ax.tick_params(axis='both', labelfontfamily='DejaVu Sans', labelsize=12)

            v.grab_frame()

plt.savefig('rho.svg', format='svg')

# %%

v = animation.FFMpegWriter(fps=15, codec='ffv1')

px = 1 / plt.rcParams['figure.dpi']
fig = plt.figure(figsize=(560*px, 420*px), dpi=150)

with v.saving(fig, 'muscato_fw.avi', 100):
    with open('Fw.dat', mode='r') as f:

        f.readline()
        st = f.readline()
        x = np.asarray(st.split()).astype(float)

        f.readline()
        st = f.readline()
        p = np.flipud(np.asarray(st.split()).astype(float))

        f.readline()

        N = 0

        st = f.readline()

        while st:
            _, _, i, _, n = f.readline().split()
            i = int(i)

            _, _, t = f.readline().split()
            t = float(t)

            f.readline()

            st = f.readline()
            fw = np.array([], dtype=float)

            while st != '\n':
                fw = np.vstack((fw, np.asarray(st.split()).astype(float))) if fw.size else np.asarray(st.split()).astype(float)

                st = f.readline()

            f.readline()
            st = f.readline()

            if i == 0:
                N = fw.sum()

            fw = fw / N / ( x[1] - x[0] ) / ( p[0] - p[1] )

            [xx, pp] = np.meshgrid(x, p)

            fig.clf()
            ax = fig.subplots()

            cf = ax.contourf(xx, pp, fw, 9, vmin=-0.3, vmax=0.3, cmap='bwr')
            ax.contour(xx, pp, fw, 9, vmin=-0.3, vmax=0.3, colors='k', linewidths=0.5)

            cbar = fig.colorbar(ScalarMappable(norm=cf.norm, cmap=cf.cmap), ax=ax,
                                ticks=[-0.3, -0.2, -0.1, 0, 0.1, 0.2, 0.3], extend='both')
            ax.set_title(f'fw, t = {t} fs', fontname='DejaVu Sans', fontsize=12, fontweight='bold')
            ax.set_xlabel('x (nm)', fontname='DejaVu Sans', fontsize=12)
            ax.set_ylabel('p', fontname='DejaVu Sans', fontsize=12)
            ax.tick_params(axis='both', labelfontfamily='DejaVu Sans', labelsize=12)
            cbar.ax.tick_params(labelfontfamily='DejaVu Sans', labelsize=12)

            v.grab_frame()

plt.savefig('Fw.svg', format='svg')
