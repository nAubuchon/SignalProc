import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc
# rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
## for Palatino and other serif fonts use:
rc('font',**{'family':'serif','serif':['Palatino']})
rc('text', usetex=True)
params = {'text.usetex' : True,
          'font.size' : 16,
          'font.family' : 'lmodern',
          'text.latex.unicode': True,
          }
plt.rcParams.update(params)

def transferMagFunc(w, r, l, c):
    w0 = 1 / np.sqrt(l * c)
    re = 1 - (w / w0) ** 2
    im = (r * c * w) ** 2

    numer = np.sqrt(re**2)
    denom = np.sqrt(re**2 + im)

    return numer/denom

def transferFunc(w, r, l, c):
    w0 = 1 / np.sqrt(l * c)
    re = 1 - (w / w0) ** 2
    im = (r * c * 1j*w)

    numer = re
    denom = re + im

    return numer/denom

# Constants
N = 2**16
PI = np.pi
R = 921.0
L = 0.10497
C = 494e-9

# Measured Values From Scope
CH1 = np.genfromtxt('[PATH].csv', delimiter=',')
CH2 = np.genfromtxt('[PATH].csv', delimiter=',')
freq_M = np.genfromtxt('[PATH].csv', delimiter=',')
omega_M = freq_M * 2 * PI
phase_M = np.genfromtxt('[PATH].csv', delimiter=',')
phase_M[0:174] *= -1
H_M_bode = 20.0 * np.log10(CH2/CH1)
# phase_M_deg = phase_M * (180/PI)

# Theoretical Values
freq_T = np.linspace(20, 20000, N)
omega_T = freq_T * 2 * PI
H_T_mag = transferMagFunc(omega_T, R, L, C)
H_T_bode = 20*np.log10(H_T_mag)

#for phase
H_T = transferFunc(omega_T, R, L, C)
im = np.imag(H_T)
re = np.real(H_T)
phase_T = np.arctan2(im, re)
zero_line = np.zeros((len(freq_T),))

# Calculating quality factor and center/cutoff frequencies
Q = (1/R)*(np.sqrt(L/C))
f0 = (1/(2*PI*np.sqrt(L*C)))
a = 1/(4 * Q**2)
b = 1/(2 * Q)
f1 = f0 * (np.sqrt(1 + a) - b)
Wc1 = 2*PI*f1
f2 = f0 * (np.sqrt(1 + a) + b)
Wc2 = 2*PI*f2
BW = (f2 - f1)/2

print 'Q factor:', Q
print 'Center freq: ', 2*PI*f0
print 'Lower cutoff: %.1f Hz (%.1f rad)' % (f1, Wc1)
print 'Upper cutoff: %.1f Hz (%.1f rad)' % (f2, Wc2)
print 'Bandwidth: ', BW*2*PI
print '1/RC = %.1f' % (1/(R*C))
print 'R/L = %.1f' % (R/L)

print

# intersections #

# zeroline intersect f1
i1_x = 1/(R*C)
Hjw1 = transferMagFunc(i1_x, R, L, C)
i1_y = 20 * np.log10(Hjw1)
print 'i1 x,y: (', i1_x, ',', i1_y, ')'

# f2 intersect zeroline
i2_x = (R/L)
Hjw2 = transferMagFunc(i2_x, R, L, C)
i2_y = 20 * np.log10(Hjw2)
print 'i2 x,y: (', i2_x, ',', i2_y, ')',

x1 = -20*np.log10(omega_T) - 20*np.log10(R*C)
x2 = 20*np.log10(omega_T) - 20*np.log10(R/L)

# Plots #
ylims = (-30, 2)

# TRANSFER MAGNITUDE
fig, ax1 = plt.subplots()
ax2 = ax1.twinx()
p1, = ax1.plot(omega_T, H_T_mag, label=r'$\vert$ H($j\omega$) $\vert$', color='green', linewidth=2.0)
p2, = ax2.plot(omega_M, CH2, label=r'Scope Data', color='red', linewidth=2.0)
p3 = ax1.vlines(omega_T[np.argmin(H_T_bode)], 0, 1, color='green', linestyle='-', linewidth=0.75,
                label='Th. Center Freq\n(~%.2f Hz)' % freq_T[np.argmin(H_T_mag)])
p4 = ax1.vlines(omega_M[np.argmin(CH2)], 0, 1, color='red', linestyle='-', linewidth=1.00,
                label='M. Center Freq\n(~%.2f Hz)' % freq_M[np.argmin(CH2)])
plt.xscale('log')
plt.xlim(min(omega_T), max(omega_T))
ax1.set_ylim(ylims)
ax1.grid(b=True, which='major', color='black', linestyle='-', linewidth=0.2)
ax2.grid(b=True, which='major', color='red', linestyle='-', linewidth=0.2)
plt.legend(handles=[p1, p2, p3, p4], loc=3)
ax1.set_xlabel(r'Frequency (rad)')
ax1.set_ylabel(r'$\vert$ H($j\omega$) $\vert$')
ax2.set_ylabel(r'Circuit Output (V)')
plt.title(r'$Theoretical$ $Magnitude$   v   $Measured$ $V_{out}$')

# ASYMPTOTES
plt.figure('Bode Asymptotes')
plt.title(r'$Asymptotic$ $Approximations$ ($R=921$ $\Omega$)')
p1, = plt.plot(omega_M, H_M_bode, label='Measured', color='blue', linewidth=2.0)
p2, = plt.plot(omega_T, H_T_bode, label='Theoretical', color='green', linewidth=2.0)
plt.plot(omega_T, x1, color='black', linestyle='--', linewidth=1.5)
plt.plot(omega_T, x2, color='black', linestyle='--', linewidth=1.5)
plt.plot(omega_T, zero_line, label='Zero', color='black', linestyle='--', linewidth=1.5)
plt.plot(omega_T, zero_line, label='Zero', color='black', linestyle='--', linewidth=1.5)
p3 = plt.vlines(i1_x, -7.5, 2, color='red', linestyle='-', linewidth=2.0, label=r'$\frac{1}{RC}$ (%.1f rad)' % i1_x)
p4 = plt.vlines(i2_x, -7.5, 2, color='orange', linestyle='-', linewidth=2.0, label=r'$\frac{R}{L}$ (%.1f rad)' % i2_x)
plt.vlines(i2_x, -7.5, 2, color='orange', linestyle='-', linewidth=2.0, label=r'$\frac{R}{L}$ (%.1f rad)' % i2_x)
plt.xscale('log')
plt.xlim(min(omega_T), max(omega_T))
plt.ylim(ylims)
plt.grid(b=True, which='major', color='grey', linestyle=':', linewidth=1.0)
plt.grid(b=True, which='minor', color='red', linestyle=':', linewidth=1.0)
plt.legend(handles=[p1, p2, p3, p4], loc=3)
plt.xlabel('Frequency (rad)')
plt.ylabel('Magnitude (dB)')

# PHASE
fig, ax1 = plt.subplots()
ax2 = ax1.twinx()
p1, = ax1.plot(omega_M, phase_M, label='Measured', color='blue', linewidth=1.50)
p2, = ax1.plot(omega_T, phase_T, label='Theoretical', color='green', linewidth=1.50)
p3, = ax2.plot(omega_M, CH2, label=r'Scope Output', color='red', linewidth=2.0)
p4, = ax1.plot(omega_T, zero_line, color='black')
plt.xscale('log')
plt.xlim(min(omega_T), max(omega_T))
ax1.set_ylim(-1.6, 1.6)
plt.grid(b=True, which='minor', color='red', linestyle='--')
plt.legend(handles=[p1, p2, p3], loc=3)
ax1.set_xlabel(r'Frequency (rad)')
ax1.set_ylabel(r'Phase (rad)')
ax2.set_ylabel(r'Circuit Output (V)')
plt.title(r'$Phase$, $Theoretical$ v $Measured$')

# POLE-ZERO
plt.figure('Pole-Zero Plot')
plt.title(r'$Pole$-$Zero$ $Plot$', y=1.08)
plt.hlines(0, -10000, 10000, color='black', linestyle='-', linewidth=2.0)
plt.text(42, -20, r'$\sigma$', fontsize=35, horizontalalignment='center', verticalalignment='center')
plt.vlines(0, -6000, 6000, color='black', linestyle='-', linewidth=2.0)
plt.text(0, 5400, r'$j\omega$', fontsize=25, horizontalalignment='center', verticalalignment='center')
plt.plot(0, 4391, "o", markersize=20, markeredgecolor='red', markeredgewidth=2.0, markerfacecolor='#DCDCDC')
plt.annotate('(0, 4391j)', xy=(0, 4391), xytext=(-15, 4391),
             arrowprops=dict(arrowstyle="->", connectionstyle="arc3"))
plt.plot(0, -4391, "o", markersize=20, markeredgecolor='red', markeredgewidth=2.0, markerfacecolor='#DCDCDC')
plt.annotate('(0, -4391j)', xy=(0, -4391), xytext=(-15, -4391),
             arrowprops=dict(arrowstyle="->", connectionstyle="arc3"))
plt.plot(-43, 197, "x", markersize=10, markeredgecolor='red', markeredgewidth=2.0, markerfacecolor='white')
plt.annotate('(-43, 197j)', xy=(-43, 197), xytext=(-50, 1700),
             arrowprops=dict(arrowstyle="->", connectionstyle="arc3"))
plt.plot(-43, -197, "x", markersize=10, markeredgecolor='red', markeredgewidth=2.0, markerfacecolor='white')
plt.annotate('(-43, -197j)', xy=(-43, -197), xytext=(-50, -1700),
             arrowprops=dict(arrowstyle="->", connectionstyle="arc3"))
plt.axvspan(-42.8, 10000, color='#DCDCDC')
plt.xlim(-60, 40)
plt.ylim(-5000, 5000)

plt.show()
