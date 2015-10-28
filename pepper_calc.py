import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import sys
sys.path.append('/home/jls/work/code/aa/lib')
try:
    reload(aa)
except NameError:
    import aa_py as aa
sys.path.append('/data/jls/galpy/')
try:
    reload(galpy.df_src.streamgapdf)
except NameError:
    import galpy.df_src.streamgapdf
from galpy.util import bovy_plot
import seaborn as sns
from matplotlib.ticker import MaxNLocator
labels = {'x':r'$x/{\rm kpc}$','y':r'$y/{\rm kpc}$','z':r'$z/{\rm kpc}$',
          'vx':r'$v_x/{\rm km\,s}^{-1}$','vy':r'$v_y/{\rm km\,s}^{-1}$','vz':r'$v_z/{\rm km\,s}^{-1}$',
          'dvx':r'$\Delta v_x/{\rm km\,s}^{-1}$','dvy':r'$\Delta v_y/{\rm km\,s}^{-1}$','dvz':r'$\Delta v_z/{\rm km\,s}^{-1}$'
          ,'R':r'$R/{\rm kpc}$','phi':r'$\phi/{\rm rad}$',
          'vR':r'$v_R/{\rm km\,s}^{-1}$','vphi':r'$v_\phi/{\rm km\,s}^{-1}$',
          'E':r'$E/({\rm km\,s}^{-1})^2$','Lz':r'$L_z/{\rm kpc\,km\,s}^{-1}$','E100':r'$E/({\rm 100\,km\,s}^{-1})^2$','dE100':r'$\Delta E/({\rm 100\,km\,s}^{-1})^2$',
          'JR':r'$J_R/{\rm kpc\,km\,s}^{-1}$','Jz':r'$J_z/{\rm kpc\,km\,s}^{-1}$',
          'OmR':r'$\Omega_R/{\rm kpc}^{-1}{\rm km\,s}^{-1}$','Omp':r'$\Omega_\phi/{\rm kpc}^{-1}{\rm km\,s}^{-1}$','Omz':r'$\Omega_z/{\rm kpc}^{-1}{\rm km\,s}^{-1}$',
          'angR':r'$\theta_R/{\rm rad}$','angp':r'$\theta_\phi/{\rm rad}$','angz':r'$\theta_z/{\rm rad}$',
          'dJR':r'$\Delta J_R/{\rm kpc\,km\,s}^{-1}$','dL':r'$\Delta L/{\rm kpc\,km\,s}^{-1}$','dLz':r'$\Delta L_z/{\rm kpc\,km\,s}^{-1}$','dJz':r'$\Delta J_z/{\rm kpc\,km\,s}^{-1}$',
          'dOmR':r'$\Delta \Omega_R/{\rm kpc}^{-1}{\rm km\,s}^{-1}$','dOmp':r'$\Delta \Omega_\phi/{\rm kpc}^{-1}{\rm km\,s}^{-1}$','dOmz':r'$\Delta \Omega_z/{\rm kpc}^{-1}{\rm km\,s}^{-1}$','dOm':r'$\Delta \Omega/{\rm kpc}^{-1}{\rm km\,s}^{-1}$',
          'dOmRG':r'$\Delta \Omega_R/{\rm Gyr}^{-1}$','dOmpG':r'$\Delta \Omega_\phi/{\rm Gyr}^{-1}$','dOmzG':r'$\Delta \Omega_z/{\rm Gyr}^{-1}$','dOmG':r'$\Delta \Omega/{\rm Gyr}^{-1}$',
            'OmRG':r'$\Omega_R/{\rm Gyr}^{-1}$','OmpG':r'$\Omega_\phi/{\rm Gyr}^{-1}$','OmzG':r'$\Omega_z/{\rm Gyr}^{-1}$','OmG':r'$\Omega/{\rm Gyr}^{-1}$',
            'dangR':r'$\Delta \theta_R/{\rm rad}$','dangp':r'$\Delta \theta_\phi/{\rm rad}$','dangz':r'$\Delta \theta_z/{\rm rad}$','dang':r'$\Delta \theta/{\rm rad}$','dangpar':r'$\Delta \theta_{||}/{\rm rad}$','dangperp':r'$\Delta \theta_{\perp}/{\rm rad}$','angpar':r'$\theta_{||}/{\rm rad}$'
          ,'dOmg':r'$\delta \Omega^g/{\rm kpc}^{-1}{\rm km\,s}^{-1}$','dOmgG':r'$\delta \Omega^g/{\rm Gyr}^{-1}$'
          ,'dangg':r'$\delta \theta^g/{\rm rad}$','dvg':r'$\delta v^g/{\rm km\,s}^{-1}$',
          'dJRg':r'$\delta J^g_R/{\rm kpc\,km\,s}^{-1}$','dLg':r'$\delta L^g/{\rm kpc\,km\,s}^{-1}$','dLzg':r'$\delta L^g_z/{\rm kpc\,km\,s}^{-1}$','dJzg':r'$\delta J^g_z/{\rm kpc\,km\,s}^{-1}$',
          'dOmparG':r'$\Delta \Omega_{||}/{\rm Gyr}^{-1}$','dOmperpG':r'$\Delta \Omega_{\perp}/{\rm Gyr}^{-1}$'}
a = sns.plotting_context('ipython')
data_folder = '/data/jls/stream_gaps/tilted/'
plots_folder = '../paper/plots/'
kms2kpcMyr=1./977.775320024919

def find_scale_radius(mass):
    return (mass/10.**8.5)**(1./2.5)

def bombard(data,N,gc,pot,sigmaw,M=10**5,sigmab=0.01,tevolve=1000.,output_plot=False,deltav=False):
    w=np.array([6.822,132.77,149.417])
    rs=0.625
    b=0
    sb = sub_halo(M,b,rs,w,gc)
    afk_r=angle_freq_kicks(data,sb,pot,output_plots=output_plot)
    data2 = data.copy()
    data2['angR_now']=data['angR_y']
    data2['angp_now']=data['angp_y']
    data2['angz_now']=data['angz_y']
    data2['OmR_now']=data['OmR_y']
    data2['Omp_now']=data['Omp_y']
    data2['Omz_now']=data['Omz_y']
    vdotv = np.zeros(N)
    xdist = np.zeros((N,3))
    for n in range(N):
        indxf = np.random.uniform()*(len(afk_r.stream_track)-1)
        indx,indxu = np.floor(indxf),np.ceil(indxf)
        centre_hit = afk_r.stream_track[indx]+(indxf-indx)*(afk_r.stream_track[indxu]-afk_r.stream_track[indx])
        w = np.random.normal(loc=0.,scale=sigmaw)
        if(deltav):
        	indx=deltav
        	centre_hit = afk_r.stream_track[indx]
	        w = sigmaw
        b = np.random.normal(loc=0.,scale=sigmab)
        rs = find_scale_radius(M)
        sb = sub_halo(M,b,rs,w,centre_hit)
        vecb = np.cross(w,centre_hit[3:])
        vecb = b*vecb/np.sqrt(np.sum(vecb**2))
        xdist[n]=centre_hit[:3]+vecb
        afk_r.new_kick(sb)
        data2 = afk_r.kick_unperturbed(data2)
        vdotv[n]=np.max(afk_r.vdotdv)
        data2['angR_now']=np.fmod(data2['angR_now']+data2['angR_k'],2.*np.pi)
        data2['angp_now']=np.fmod(data2['angp_now']+data2['angp_k'],2.*np.pi)
        data2['angz_now']=np.fmod(data2['angz_now']+data2['angz_k'],2.*np.pi)
        data2['OmR_now']=(data2['OmR_now']+data2['OmR_k'])
        data2['Omp_now']=(data2['Omp_now']+data2['Omp_k'])
        data2['Omz_now']=(data2['Omz_now']+data2['Omz_k'])
    data2['tRnow']=data2['angR_now']+tevolve*kms2kpcMyr*data2['OmR_now']
    data2['tpnow']=data2['angp_now']+tevolve*kms2kpcMyr*data2['Omp_now']
    data2['tznow']=data2['angz_now']+tevolve*kms2kpcMyr*data2['Omz_now']
    data2['tRnow_u']=data2['angR_y']+tevolve*kms2kpcMyr*data2['OmR_y']
    data2['tpnow_u']=data2['angp_y']+tevolve*kms2kpcMyr*data2['Omp_y']
    data2['tznow_u']=data2['angz_y']+tevolve*kms2kpcMyr*data2['Omz_y']
    data2['tparnow_u']=afk_r.n[0]*data2['tRnow_u']+afk_r.n[1]*data2['tpnow_u']+afk_r.n[2]*data2['tznow_u']
    data2['tparnow']=afk_r.n[0]*data2['tRnow']+afk_r.n[1]*data2['tpnow']+afk_r.n[2]*data2['tznow']
    data2['Omparnow_u']=afk_r.n[0]*data2['OmR_y']+afk_r.n[1]*data2['Omp_y']+afk_r.n[2]*data2['Omz_y']
    data2['Omparnow']=afk_r.n[0]*data2['OmR_now']+afk_r.n[1]*data2['Omp_now']+afk_r.n[2]*data2['Omz_now']
    return data2,vdotv,xdist,tevolve


def bombard_plot(bombarded_peri,vv,xd,time):
    f,a = plt.subplots(2,4,figsize=(13,7))
    plt.subplots_adjust(wspace=0.,hspace=0.35)
    n,b,p = a[0][0].hist(bombarded_peri.tRnow,histtype='step',bins=150,lw=1)
    n2,b2,p2 = a[0][0].hist(bombarded_peri.tRnow_u,histtype='step',bins=b,lw=1)
    a[0][0].plot(.5*(b[1:]+b[:-1]),n-n2+1500.)
    n,b,p = a[0][1].hist(bombarded_peri.tpnow,histtype='step',bins=150,lw=1)
    n2,b2,p2 = a[0][1].hist(bombarded_peri.tpnow_u,histtype='step',bins=b,lw=1)
    a[0][1].plot(.5*(b[1:]+b[:-1]),n-n2+1500.)
    n,b,p = a[0][2].hist(bombarded_peri.tznow,histtype='step',bins=150,lw=1)
    n2,b2,p2 = a[0][2].hist(bombarded_peri.tznow_u,histtype='step',bins=b,lw=1)
    a[0][2].plot(.5*(b[1:]+b[:-1]),n-n2+1500.)
    n,b,p = a[0][3].hist(bombarded_peri.tparnow,histtype='step',bins=150,lw=1)
    n2,b2,p2 = a[0][3].hist(bombarded_peri.tparnow_u,histtype='step',bins=b,lw=1)
    a[0][3].plot(.5*(b[1:]+b[:-1]),n-n2+1500.)
    n,b,p = a[1][0].hist(bombarded_peri.OmR_now,histtype='step',bins=150,lw=1)
    n2,b2,p2 = a[1][0].hist(bombarded_peri.OmR_y,histtype='step',bins=b,lw=1)
    a[1][0].plot(.5*(b[1:]+b[:-1]),n-n2+1500.)
    n,b,p = a[1][1].hist(bombarded_peri.Omp_now,histtype='step',bins=150,lw=1)
    n2,b2,p2 = a[1][1].hist(bombarded_peri.Omp_y,histtype='step',bins=b,lw=1)
    a[1][1].plot(.5*(b[1:]+b[:-1]),n-n2+1500.)
    n,b,p = a[1][2].hist(bombarded_peri.Omz_now,histtype='step',bins=150,lw=1)
    n2,b2,p2 = a[1][2].hist(bombarded_peri.Omz_y,histtype='step',bins=b,lw=1)
    a[1][2].plot(.5*(b[1:]+b[:-1]),n-n2+1500.)
    n,b,p = a[1][3].hist(bombarded_peri.Omparnow,histtype='step',bins=150,lw=1)
    n2,b2,p2 = a[1][3].hist(bombarded_peri.Omparnow_u,histtype='step',bins=b,lw=1)
    a[1][3].plot(.5*(b[1:]+b[:-1]),n-n2+1500.)
    for i in range(1,4):
        a[1][i].set_yticklabels([])
        a[0][i].set_yticklabels([])
        a[0][i-1].xaxis.set_major_locator(MaxNLocator(5,prune='upper'))
        a[1][i-1].xaxis.set_major_locator(MaxNLocator(5,prune='upper'))
    a[0][3].xaxis.set_major_locator(MaxNLocator(5))
    a[1][3].xaxis.set_major_locator(MaxNLocator(5))
    a[0][0].set_xlabel(labels['angR'])
    a[0][1].set_xlabel(labels['angp'])
    a[0][2].set_xlabel(labels['angz'])
    a[1][0].set_xlabel(labels['OmR'])
    a[1][1].set_xlabel(labels['Omp'])
    a[1][2].set_xlabel(labels['Omz'])
    a[0][3].set_xlabel(r'$\theta_{||}/\mathrm{rad}$')
    a[1][3].set_xlabel(r'$\Omega_{||}/\mathrm{kpc}^{-1}\,\mathrm{km\,s}^{-1}$')
    a[0][0].set_ylabel(r'$\mathrm{d}N/\mathrm{d}\theta$')
    a[1][0].set_ylabel(r'$\mathrm{d}N/\mathrm{d}\Omega$')
    a[0][0].text(0.,1.,r'$t='+str(time/1000.)+r'\,\mathrm{Gyr}$',horizontalalignment='left',verticalalignment='bottom',transform=a[0][0].transAxes,fontsize=22)

    f,a = plt.subplots(1,3,figsize=(13,4))
    plt.subplots_adjust(wspace=0.3)
    a[0].plot(bombarded_peri.tparnow,bombarded_peri.Omparnow,'.',ms=1)
    a[0].set_xlabel(r'$\theta_{||}/\mathrm{rad}$')
    a[0].set_ylabel(r'$\Omega_{||}/\mathrm{kpc}^{-1}\,\mathrm{km\,s}^{-1}$')
    q = a[1].hist(vv,histtype='step',lw=3)
    a[1].set_xlabel(r'${\bf v}\cdot\delta{\bf v}$')
    a[1].set_ylabel(r'$\mathrm{d}N/\mathrm{d}({\bf v}\cdot\delta{\bf v})$')
    a[2].plot(bombarded_peri.x_y,bombarded_peri.y_y,'k.',ms=1)
    q = a[2].scatter(xd.T[0],xd.T[1],c=vv,lw=0,s=50,cmap=plt.cm.jet)
    a[2].set_aspect('equal')
    a[2].set_xlabel(labels['x'])
    a[2].set_ylabel(labels['y'])

def plot_power(data,ax,axcorr=None):
	angpar=np.histogram(data,bins=400)
	FFT = np.fft.fft(angpar[0])
	deltaX = (angpar[1][-1]-angpar[1][0])/len(angpar)
	Xrange = np.fft.fftfreq(angpar[0].size, d=deltaX)
	ax.plot(Xrange[:len(Xrange)/2],(FFT*np.conj(FFT))[:len(Xrange)/2])
	ax.semilogy()
	ax.set_xlim(0.,np.max(Xrange))
	if(axcorr):
		IFFT=np.fft.fftshift(np.fft.ifft(FFT*np.conj(FFT)))*deltaX
		axcorr.plot((angpar[1][1:]+angpar[1][:-1])/2.,IFFT)

def plot_diff_power(data,data2,ax,axcorr=None):
	angpar=np.histogram(data,bins=400)
	angparu=np.histogram(data2,bins=angpar[1])
	FFT = np.fft.fft(angpar[0]-angparu[0])
	deltaX = (angpar[1][-1]-angpar[1][0])/len(angpar)
	Xrange = np.fft.fftfreq(angpar[0].size, d=deltaX)
	ax.plot(Xrange[:len(Xrange)/2],(FFT*np.conj(FFT))[:len(Xrange)/2])
	ax.semilogy()
	ax.set_xlim(0.,np.max(Xrange))
	if(axcorr):
		IFFT=np.fft.fftshift(np.fft.ifft(FFT*np.conj(FFT)))*deltaX
		bins = (angpar[1][1:]+angpar[1][:-1])/2.
		bins = bins-bins[len(bins)/2]
		axcorr.plot(bins,IFFT)

def bombard_ft_plot(bombarded_peri):
	f,a = plt.subplots(4,2,figsize=(13,15))
	plt.subplots_adjust(wspace=0.2,hspace=0.35)
	plot_power(bombarded_peri['tparnow_u'],ax=a[0][0],axcorr=a[0][1])
	plot_power(bombarded_peri['tparnow'],ax=a[0][0],axcorr=a[0][1])
	a[0][0].set_xlabel(r'$k_\theta/2\pi/\mathrm{rad}^{-1}$')
	a[0][0].set_ylabel(r'Power')
	a[0][1].set_xlabel(r'$\theta/\mathrm{rad}$')
	a[0][1].set_ylabel(r'Corr.')
	plot_power(bombarded_peri['Omparnow_u'],ax=a[1][0],axcorr=a[1][1])
	plot_power(bombarded_peri['Omparnow'],ax=a[1][0],axcorr=a[1][1])
	a[1][0].set_xlabel(r'$k_\Omega/2\pi/\mathrm{Gyr}$')
	a[1][1].set_xlabel(r'$\Omega/\mathrm{Gyr}^{-1}$')
	a[1][0].set_ylabel(r'Power')
	a[1][1].set_ylabel(r'Corr.')
	plot_diff_power(bombarded_peri['tparnow'],bombarded_peri['tparnow_u'],ax=a[2][0],axcorr=a[2][1])
	a[2][0].set_xlabel(r'$k_\theta/2\pi/\mathrm{rad}^{-1}$')
	a[2][0].set_ylabel(r'Power')
	a[2][1].set_xlabel(r'$\theta/\mathrm{rad}$')
	a[2][1].set_ylabel(r'Corr.')
	plot_diff_power(bombarded_peri['Omparnow'],bombarded_peri['Omparnow_u'],ax=a[3][0],axcorr=a[3][1])
	a[3][0].set_xlabel(r'$k_\Omega/2\pi/\mathrm{Gyr}$')
	a[3][1].set_xlabel(r'$\Omega/\mathrm{Gyr}^{-1}$')
	a[3][0].set_ylabel(r'Power')
	a[3][1].set_ylabel(r'Corr.')

def load_denis_simulation_from_file():
	merged_peri=pd.read_csv(data_folder+'merged_peri.aa',sep=' ')
	merged_peri['angR_y']=merged_peri['angR_y']-2.*np.pi*(merged_peri['angR_y']>np.pi)
	merged_peri['angz_y']=merged_peri['angz_y']-2.*np.pi*(merged_peri['angz_y']>np.pi)
	merged_peri['dangz_y']=merged_peri['dangz_y']+2.*np.pi*(merged_peri['dangz_y']<-np.pi)
	merged_peri['angR0']=merged_peri['angR0']+2.*np.pi*(merged_peri['angR0']<-np.pi)
	merged_peri['angp0']=merged_peri['angp0']+2.*np.pi*(merged_peri['angp0']<0.)
	merged_peri['angz0']=merged_peri['angz0']+2.*np.pi*(merged_peri['angz0']<-np.pi)
	merged_peri['angz0']=merged_peri['angz0']+2.*np.pi*(merged_peri['angz0']<np.pi)
	return merged_peri

from scipy.interpolate import LSQUnivariateSpline, interp1d
from scipy.stats import linregress

class sub_halo:
    def __init__(self,M,b,rs,w,gc):
        G=4.300918e-6
        self.M=M
        self.GM = G*M
        self.b = b
        self.rs = rs
        self.w = w
        self.magw = np.sqrt(np.sum(w**2))
        self.gap_centre_impact = gc

class angle_freq_kicks:

    def dOmAngdv(self,x,dv):
        A=self.AA.angles(x)
        dOdv=np.zeros((3,3))
        dAdv=np.zeros((3,3))
        AOu=np.zeros((3,6))
        AOd=np.zeros((3,6))
        xu,xd=np.copy(x),np.copy(x)
        xu[3]+=dv
        xd[3]-=dv
        AOu[0]=self.AA.angles(xu)
        AOd[0]=self.AA.angles(xd)
        xu,xd=np.copy(x),np.copy(x)
        xu[4]+=dv
        xd[4]-=dv
        AOu[1]=self.AA.angles(xu)
        AOd[1]=self.AA.angles(xd)
        xu,xd=np.copy(x),np.copy(x)
        xu[5]+=dv
        xd[5]-=dv
        AOu[2]=self.AA.angles(xu)
        AOd[2]=self.AA.angles(xd)
        dAdv=AOu.T[:3].T-AOd.T[:3].T
        dOdv=AOu.T[3:6].T-AOd.T[3:6].T
        l,e = np.linalg.eigh(dAdv.T/2./dv)
        return A[:3],A[3:],dAdv.T/2./dv,dOdv.T/2./dv,e[-2]


    def plot_kicks(self):
        f,a=plt.subplots(1,3,figsize=(15,5));plt.subplots_adjust(wspace=0.4);
        a[0].plot(self.deltaang_par,self.dv.T[0],label=r'$\delta v^g_x$');
        a[0].plot(self.deltaang_par,self.dv.T[1],label=r'$\delta v^g_y$');
        a[0].plot(self.deltaang_par,self.dv.T[2],label=r'$\delta v^g_z$');
        a[0].set_xlabel(labels['dangpar']);a[0].set_ylabel(labels['dvg']);a[0].legend(loc=3)

        a[1].plot(self.deltaang_par,self.dA.T[0],label=r'$\delta \theta^g_R$');
        # a[1].plot(afk.deltaang_par,afk.angR_k(afk.deltaang_par));
        a[1].plot(self.deltaang_par,self.dA.T[1],label=r'$\delta \theta^g_\phi$');
        # a[1].plot(afk.deltaang_par,afk.angp_k(afk.deltaang_par));
        a[1].plot(self.deltaang_par,self.dA.T[2],label=r'$\delta \theta^g_z$')
        # a[1].plot(afk.deltaang_par,afk.angz_k(afk.deltaang_par));
        a[1].set_xlabel(labels['dangpar']);a[1].set_ylabel(labels['dangg']);a[1].legend(loc=3)

        a[2].plot(self.deltaang_par,self.dO.T[0],label=r'$\delta \Omega^g_R$');
        # a[2].plot(afk.deltaang_par,afk.OmR_k(afk.deltaang_par));
        a[2].plot(self.deltaang_par,self.dO.T[1],label=r'$\delta \Omega^g_\phi$');
        # a[2].plot(afk.deltaang_par,afk.Omp_k(afk.deltaang_par));
        a[2].plot(self.deltaang_par,self.dO.T[2],label=r'$\delta \Omega^g_z$')
        # a[2].plot(afk.deltaang_par,afk.Omz_k(afk.deltaang_par));
        a[2].set_xlabel(labels['dangpar']);a[2].set_ylabel(labels['dOmg']);a[2].legend(loc=2)


    def __init__(self,merged_peri,sub_halo_props,pot,num=230,output_plots=False):
        self.AA = aa.Actions_Genfunc(pot,"axisymmetric")

        merged_peri['phi_int']=merged_peri['phi_y']-2.*np.pi*(merged_peri['phi_y']>0.)

        ## Build spline of track
        xs = np.linspace(np.min(merged_peri['phi_int']),np.max(merged_peri['phi_int']),num)
        knots = np.linspace(np.min(merged_peri['phi_int']),np.max(merged_peri['phi_int']),12)[1:-1]

        sx=LSQUnivariateSpline(merged_peri['phi_int'][np.argsort(merged_peri['phi_int'])].values,merged_peri['x_y'][np.argsort(merged_peri['phi_int'])].values,knots)

        s=LSQUnivariateSpline(merged_peri['phi_int'][np.argsort(merged_peri['phi_int'])].values,merged_peri['y_y'][np.argsort(merged_peri['phi_int'])].values,knots)

        sz=LSQUnivariateSpline(merged_peri['phi_int'][np.argsort(merged_peri['phi_int'])].values,merged_peri['z_y'][np.argsort(merged_peri['phi_int'])].values,knots)

        s1=LSQUnivariateSpline(merged_peri['phi_int'][np.argsort(merged_peri['phi_int'])].values,merged_peri['vx_y'][np.argsort(merged_peri['phi_int'])].values,knots)

        s2=LSQUnivariateSpline(merged_peri['phi_int'][np.argsort(merged_peri['phi_int'])].values,merged_peri['vy_y'][np.argsort(merged_peri['phi_int'])].values,knots)

        s2z=LSQUnivariateSpline(merged_peri['phi_int'][np.argsort(merged_peri['phi_int'])].values,merged_peri['vz_y'][np.argsort(merged_peri['phi_int'])].values,knots)

        stream_track = np.ones((len(xs),6))*1e-5
        for i,j in enumerate(xs):
            stream_track[i][0]=sx(j)
            stream_track[i][1]=s(j)
            stream_track[i][2]=sz(j)
            stream_track[i][3]=s1(j)
            stream_track[i][4]=s2(j)
            stream_track[i][5]=s2z(j)

        dOdAdv=map(lambda i:self.dOmAngdv(i,0.5),stream_track)
        self.dOdAdv = dOdAdv
        self.subhalo=sub_halo_props
        b = sub_halo_props.b
        w = sub_halo_props.w
        GM = sub_halo_props.GM
        rs = sub_halo_props.rs
        gci = sub_halo_props.gap_centre_impact
        self.dv = galpy.df_src.streamgapdf.impulse_deltav_plummer_curvedstream(stream_track.T[3:].T,stream_track.T[:3].T,b,w,gci[:3],gci[3:],GM,rs)
        self.dA = np.array([map(lambda i,j:np.dot(i[2],j),dOdAdv,self.dv)])[0]
        self.dO = np.array([map(lambda i,j:np.dot(i[3],j),dOdAdv,self.dv)])[0]
        self.O =  np.array([map(lambda i:i[1],dOdAdv)])[0]
        self.A = np.array([map(lambda i:i[0],dOdAdv)])[0]
        self.e = np.array([map(lambda i:i[4],dOdAdv)])[0]
        for i in [0,2]:
            self.A.T[i]=self.A.T[i]-2.*np.pi*(self.A.T[i]>np.pi)

        slope, intercept, r_value, p_value, std_err = linregress(self.O.T[0],self.O.T[1])
        slope2, intercept2, r_value, p_value, std_err = linregress(self.O.T[0],self.O.T[2])
        self.n = np.array([1.,slope,slope2])
        self.n = self.n/np.sqrt(np.sum(self.n**2))
        self.acts_centre = self.AA.actions(gci)
        self.angs_centre = self.AA.angles(gci)
        for i in [0,2]:
            self.angs_centre[i]=self.angs_centre.T[i]-2.*np.pi*(self.angs_centre.T[i]>np.pi)
        self.ndotac = np.dot(self.n,self.angs_centre[:3])
        self.ndotfc = np.dot(self.n,self.angs_centre[3:])
        self.deltaang_par=np.dot(self.A-self.angs_centre[:3],self.n)
#         self.deltaang_par=self.deltaang_par-2.*np.pi*(self.deltaang_par>np.pi)
        self.vdotdv = np.array(map(lambda i:np.dot(i[0],i[1]),zip(stream_track.T[3:].T,self.dv)))
        self.stream_track=stream_track
        self.angR_k=interp1d(self.deltaang_par,self.dA.T[0])
        self.angp_k=interp1d(self.deltaang_par,self.dA.T[1])
        self.angz_k=interp1d(self.deltaang_par,self.dA.T[2])
        self.OmR_k=interp1d(self.deltaang_par,self.dO.T[0])
        self.Omp_k=interp1d(self.deltaang_par,self.dO.T[1])
        self.Omz_k=interp1d(self.deltaang_par,self.dO.T[2])
        if(output_plots):
			f,a = plt.subplots(2,6,figsize=(20,7));plt.subplots_adjust(wspace=0.7,hspace=0.3)
			hi=a[0][0].plot(np.sort(merged_peri['phi_int']),merged_peri['x_y'][np.argsort(merged_peri['phi_int'])],'.',ms=3)
			a[0][0].plot(xs,sx(xs));a[0][0].set_xlabel(labels['x']);a[0][0].set_ylabel(labels['x']);
			hi=a[0][1].plot(np.sort(merged_peri['phi_int']),merged_peri['y_y'][np.argsort(merged_peri['phi_int'])],'.',ms=3)
			a[0][1].plot(xs,s(xs));a[0][1].set_xlabel(labels['x']);a[0][1].set_ylabel(labels['y']);
			hi=a[0][2].plot(np.sort(merged_peri['phi_int']),merged_peri['z_y'][np.argsort(merged_peri['phi_int'])],'.',ms=3)
			a[0][2].plot(xs,sz(xs));a[0][2].set_xlabel(labels['x']);a[0][2].set_ylabel(labels['z']);
			hi=a[0][3].plot(np.sort(merged_peri['phi_int']),merged_peri['vx_y'][np.argsort(merged_peri['phi_int'])],'.',ms=3)
			a[0][3].plot(xs,s1(xs));a[0][3].set_xlabel(labels['x']);a[0][3].set_ylabel(labels['vx']);
			hi=a[0][4].plot(np.sort(merged_peri['phi_int']),merged_peri['vy_y'][np.argsort(merged_peri['phi_int'])],'.',ms=3)
			a[0][4].plot(xs,s2(xs));a[0][4].set_xlabel(labels['x']);a[0][4].set_ylabel(labels['vy']);
			hi=a[0][5].plot(np.sort(merged_peri['phi_int']),merged_peri['vz_y'][np.argsort(merged_peri['phi_int'])],'.',ms=3)
			a[0][5].plot(xs,s2z(xs));a[0][5].set_xlabel(labels['x']);a[0][5].set_ylabel(labels['vz']);
			a[0][0].plot(xs,stream_track.T[0])
			a[0][1].plot(xs,stream_track.T[1])
			a[0][2].plot(xs,stream_track.T[2])
			a[0][3].plot(xs,stream_track.T[3])
			a[0][4].plot(xs,stream_track.T[4])
			a[0][5].plot(xs,stream_track.T[5])
			a[1][0].plot(self.A.T[0],self.A.T[1])
			a[1][1].plot(self.A.T[0],self.A.T[2])
			a[1][2].plot(self.O.T[0],self.O.T[1])
			a[1][3].plot(self.O.T[1],self.O.T[2])
			a[1][2].plot(self.O.T[0],self.n[1]*self.O.T[0]+intercept)
			a[1][3].plot(self.n[1]*self.O.T[0]+intercept,self.n[2]*self.O.T[0]+intercept2)
			a[1][4].plot(self.A.T[0])
			self.plot_kicks()

    def kick(self,ang):
        return self.angR_k(ang),self.angp_k(ang),self.angz_k(ang),self.OmR_k(ang),self.Omp_k(ang),self.Omz_k(ang)

    def kick_unperturbed(self,data2):
        data=data2.copy()
        data['dangpar']=np.dot(self.n,np.array([data['angR_y'].T,data['angp_y'].T,data['angz_y'].T]))-self.ndotac
        mm = np.maximum(np.min(self.deltaang_par),np.minimum(np.max(self.deltaang_par),data['dangpar']))
        data['angR_k']=self.angR_k(mm)
        data['angp_k']=self.angp_k(mm)
        data['angz_k']=self.angz_k(mm)
        data['OmR_k']=self.OmR_k(mm)
        data['Omp_k']=self.Omp_k(mm)
        data['Omz_k']=self.Omz_k(mm)
        return data

    def new_kick(self,sub_halo_props):
        self.subhalo=sub_halo_props
        b = sub_halo_props.b
        w = sub_halo_props.w
        GM = sub_halo_props.GM
        rs = sub_halo_props.rs
        gci = sub_halo_props.gap_centre_impact
        self.dv = galpy.df_src.streamgapdf.impulse_deltav_plummer_curvedstream(self.stream_track.T[3:].T,self.stream_track.T[:3].T,b,w,gci[:3],gci[3:],GM,rs)
        self.dA = np.array([map(lambda i,j:np.dot(i[2],j),self.dOdAdv,self.dv)])[0]
        self.dO = np.array([map(lambda i,j:np.dot(i[3],j),self.dOdAdv,self.dv)])[0]
        self.angs_centre = self.AA.angles(gci)
        for i in [0,2]:
            self.angs_centre[i]=self.angs_centre.T[i]-2.*np.pi*(self.angs_centre.T[i]>np.pi)
        self.ndotac = np.dot(self.n,self.angs_centre[:3])
        self.ndotfc = np.dot(self.n,self.angs_centre[3:])
        self.deltaang_par=np.dot(self.A-self.angs_centre[:3],self.n)
        self.vdotdv = np.array(map(lambda i:np.dot(i[0],i[1]),zip(self.stream_track.T[3:].T,self.dv)))
        self.angR_k=interp1d(self.deltaang_par,self.dA.T[0])
        self.angp_k=interp1d(self.deltaang_par,self.dA.T[1])
        self.angz_k=interp1d(self.deltaang_par,self.dA.T[2])
        self.OmR_k=interp1d(self.deltaang_par,self.dO.T[0])
        self.Omp_k=interp1d(self.deltaang_par,self.dO.T[1])
        self.Omz_k=interp1d(self.deltaang_par,self.dO.T[2])
