import numpy as np
import matplotlib.pyplot as plt
import scipy.io as sio
from miscellaneous import *

np.random.seed(2021)

subcatchments = sio.loadmat('../dataOCN/subcatchments.mat')


CTC = subcatchments['CTC'].astype('int')
X = subcatchments['X']
Y = subcatchments['Y']

NX = np.ceil(np.max(X)).astype('int')
NY = np.ceil(np.max(Y)).astype('int')

COLS = ["#ff595e", "#ff924c", "#ffca3a", "#8ac926", "#1982c4", "#6a4c93"]

CCR = np.random.randint(25, size=np.max(CTC))



CTD = np.zeros(len(CTC))
for i in range(len(CTD)):
    CTD[i] = CCR[CTC[i]-1]

CTM = np.zeros((NX, NY))
MSC = np.zeros((NX, NY))
XI = np.ceil(X)
YI = np.ceil(Y)

for i in range(len(X)):
    CTM[int(XI[i])-1, int(YI[i])-1] = CTD[i]
    MSC[int(XI[i])-1, int(YI[i])-1] = CTC[i]
    
FD = sio.loadmat('../dataOCN/FD.mat')

FD_A = FD['FD_A']
FD_X = FD['FD_X']
FD_Y = FD['FD_Y']
FD_downNode = FD['FD_downNode'].astype('int')



FD_nNodes = len(CTD)
FD_outlet = 1921-1
FD_A = FD_A.flatten()
FD_X = FD_X.flatten()
FD_Y = FD_Y.flatten()
FD_downNode = FD_downNode.flatten()

thrA = 120
cellsize = 1

CTA = np.zeros((NX, NY))

AG = sio.loadmat('../dataOCN/AG.mat')


AG_A = AG['A']

downNode = AG['downNode']

SC = sio.loadmat('../dataOCN/SC.mat')
SCX = SC['SCX']
SCY = SC['SCY']

for i in range(len(X)):
    CTA[int(XI[i])-1, int(YI[i])-1] = AG_A[CTC[i]-1]


ranks = np.random.permutation(len(AG_A))+1
P = zipf(ranks, 2, 250)

Nf = (AG_A * 1000) ** 0.3
hf = np.ones(len(CTC)) * 0.2

nNodes = len(SCX)




fig, ax = plt.subplots()
img = ax.imshow(CTM.T, cmap = 'Set1', aspect='auto', origin='lower')
plt.title('Watersheds')
for i in range(FD_nNodes):
    if i + 1 != FD_outlet:
        if FD_A[i] >= thrA and abs(FD_X[i] - FD_X[FD_downNode[i]-1]) <= cellsize and abs(FD_Y[i] - FD_Y[FD_downNode[i]-1]) <= cellsize:
            ax.plot([FD_X[i], FD_X[FD_downNode[i]-1]], [FD_Y[i], FD_Y[FD_downNode[i]-1]], linewidth=0.5+4.5*np.sqrt(FD_A[i]/(FD_nNodes*cellsize**2)), color='white')

for sc in range(len(P)):
    ax.plot(SCX[sc], SCY[sc], '.r', markersize=0.5+1.5*np.log(P[sc]))
    ax.text(SCX[sc], SCY[sc], str(sc+1), color='white')
    
drawborders(MSC)

ax.axis('off')
ax.set_xlim([0, NX])
ax.set_ylim([0, NY])
plt.show()


M = np.zeros((nNodes, nNodes))
DST = np.zeros((nNodes, nNodes))
LOC = 0.1 * np.ones((1, nNodes))
par_phi = 5

dist_matrix = np.sqrt((SCX[:, np.newaxis] - SCX) ** 2 + (SCY[:, np.newaxis] - SCY) ** 2)

M = np.where(np.eye(nNodes, dtype=bool), 0,  np.outer(P,P) / (dist_matrix ** par_phi))
M = M / np.sum(M, axis=0) * (1 - LOC)
np.fill_diagonal(M, LOC)

plt.figure()
plt.imshow(M)
plt.show()


HD = np.zeros((nNodes, nNodes))
HU = np.zeros((nNodes, nNodes))

nn_vals = np.arange(1, SC.nNodes+1)
downNode_nn = downNode[:, np.newaxis] == nn_vals

HD[nn_vals-1] = downNode_nn.T
HU[downNode_nn] = 1
