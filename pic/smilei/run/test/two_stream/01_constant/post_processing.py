import happi
import matplotlib.pyplot as plt

S=happi.Open()
uelm = S.Scalar('Uelm',data_log=True, ylabel = 'Log(Total electromagnetic)')
ne  = S.Field(0,'-Rho_eon1-Rho_eon2', xmin=0, xmax=1.05, ylabel = 'electron density')
ex  = S.Field(0,'Ex', xmin=0, xmax=1.05)
# phs = S.ParticleBinning(0,data_log=True)
phs = S.ParticleBinning(0,cmap="smilei_r")
happi.multiSlide(uelm,ne,ex,phs,shape=[2,2])

plt.show(block=True)