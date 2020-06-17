try:
    import numpy as npy
except ImportError:
    print("Numpy module could not be found")

try:
	import matplotlib.pyplot as plt
	from matplotlib import colors
except ImportError:
    print("matplotlib pyplot and/or colors components could not be found")


# plot admittance matrix

cmap='cool'
fs1=18
fs2=16
fs3=12
ms1=6

d1=npy.genfromtxt('admmatR.dat')
d2=npy.genfromtxt('admmatI.dat')
a=d1**2+d2**2
a[npy.where(a<1.e-6)]=1.e-20

fig=plt.figure(num=None, figsize=(6, 6), dpi=150, facecolor='w', edgecolor='w')
ax=fig.add_axes([0.10,0.08,0.89,0.89])
plt.imshow(npy.log(a),cmap=plt.cm.BuPu)
for tick in ax.xaxis.get_major_ticks()+ax.yaxis.get_major_ticks():
    tick.label.set_fontsize(fs1) 

plt.savefig('admmat.pdf')

# plot voltage magniture and angles
d1 = npy.genfromtxt('statsol.dat')

fig=plt.figure(num=None, figsize=(6, 4), dpi=150, facecolor='w', edgecolor='w')
ax=fig.add_axes([0.13,0.14,0.85,0.84])
ax.plot(d1[:,1],'o',ms=ms1)
ax.set_ylim([0.9,1.1])
ax.set_yticks([0.9,1,1.1])
ax.set_xlabel('bus ID',fontsize=fs2)
ax.set_ylabel('voltage magnitude [p.u.]',fontsize=fs2)
for tick in ax.xaxis.get_major_ticks()+ax.yaxis.get_major_ticks():
    tick.label.set_fontsize(fs3) 

plt.savefig('vm.pdf')

fig=plt.figure(num=None, figsize=(6, 4), dpi=150, facecolor='w', edgecolor='w')
ax=fig.add_axes([0.13,0.14,0.85,0.84])
ax.plot(d1[:,2],'o',ms=ms1)
#ax.set_ylim([0.9,1.1])
ax.set_yticks([0,10,20,30,40])
ax.set_xlabel('bus ID',fontsize=fs2)
ax.set_ylabel('voltage angle [deg]',fontsize=fs2)
for tick in ax.xaxis.get_major_ticks()+ax.yaxis.get_major_ticks():
    tick.label.set_fontsize(fs3) 

plt.savefig('va.pdf')
