import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from matplotlib import rc
from mpl_toolkits import mplot3d


def α(s, z, A, a, r):
    return a * A * r - 2.0 * a * z - a**2 * s

def Sigm(v, vr):
    return 2.0 * r0 / (1 + np.exp(ρ * (vr - v)))

def LaNMM(t, u, p):
    p1, p2, p3, p4, p5 = p
    s = u[:5]
    z = u[5:]
    ds = z
    input_val = np.array([
        C[0]*s[1] + C[1]*s[2] + C[10]*s[3] + C[2]*A_AMPA/a_AMPA*(p1),
        C[3]*s[0],
        C[4]*s[0],
        C[5]*s[3] + C[6]*s[4] + C[11]*s[0] + A_AMPA/a_AMPA*(p4),
        C[8]*s[3] + C[9]*s[4] + C[12]*s[0]
    ])
    dz = α(s, z, B, b, Sigm(input_val, w0))
    return np.concatenate([ds, dz])

def fourier(sgnl,tipo):
	if tipo == 1:
		signal = pyr1(sgnl) - np.mean(pyr1(sgnl))
	elif tipo == 2:
		signal = pyr2(sgnl) - np.mean(pyr2(sgnl))
	else:
		print("Error, tipo must be 1 or 2.")
	t = [sgnl.t[-1],sgnl.t[0]]
	fft_result = np.fft.fft(signal)
	freqs = np.fft.fftfreq(len(signal), (t[1] - t[0])/len(signal))
	freqs = np.fft.fftshift(freqs)
	fft_result = np.fft.fftshift(fft_result)
	return freqs, np.abs(fft_result)

def pyr1(sol):
    return C[0]*sol.y[1] + C[1]*sol.y[2] + C[10]*sol.y[3]

def pyr2(sol):
    return C[5]*sol.y[3] + C[6]*sol.y[4] + C[11]*sol.y[0]

def pvs(sol):
	return C[8]*sol.y[3] + C[9]*sol.y[4] + C[12]*sol.y[0]

colores =["#5F0EFF","#FF470E"]

a_AMPA = 100.0 
a_GABAs = 50.0
a_GABAf = 220.0
r0 = 2.5 
va = 6.0
vb = 1.0
ρ = 0.56
A_AMPA = 3.25
A_GABAs = -22.0
A_GABAf = -30.0
B = np.array([A_AMPA, A_AMPA, A_GABAs, A_AMPA, A_GABAf])
b = np.array([a_AMPA, a_AMPA, a_GABAs, a_AMPA, a_GABAf])
C = np.array([108.0, 33.75, 1.0, 135.0, 33.75, 70.0, 300.0, 1.0, 200.0, 100.0, 80.0, 200.0, 30.0])
w0 = np.array([va, va, va, vb, va])



u0 = np.random.rand(10) 
t_span = (0.0, 20.0)

t_evl=np.linspace(10, 20, 10000)

p1 = 200
p2 = 0
p3 = 0
p4 = 90
p5 = 0
pp = [p1,p2,p3,p4,p5]


sol = solve_ivp(lambda t, y: LaNMM(t, y, pp), t_span, u0, method='RK23',t_eval = t_evl)#,rtol = 1e-10, atol = 1e-12)

fig,axs = plt.subplots(2,1,sharex=True)
axs[0].plot(sol.t, pyr1(sol), color=colores[0],label=r'$v_{P_1}$')
axs[1].plot(sol.t, pyr2(sol), color=colores[1],label=r'$v_{P_2}$')
axs[1].set_xlabel('Time')
axs[0].set_ylabel('Activity')
axs[1].set_ylabel('Activity')
axs[0].legend()
axs[1].legend()
axs[0].set_xlim(18,20)
axs[1].set_xlim(18,20)
plt.show()




lbls = ["P1","P2"]
fig, axx = plt.subplots(2, 1, figsize=(16,10),sharex=True)
for i in range(2):
	fft = fourier(sol,i+1)
	idxs = np.argwhere(fft[1] == np.max(fft[1])).T[0]
	axx[i].plot(fft[0], fft[1],color=colores[i],label=lbls[i])
	axx[i].set_xlabel('Frequency (Hz)')
	axx[i].set_ylabel('Magnitude')
	axx[i].set_xlim(0,100)
	axx[i].grid(True)



for ax in fig.get_axes():
	ax.label_outer()

plt.show()
