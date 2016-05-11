import os, os.path
import copy
import pickle
import numpy
import tqdm
import subprocess
import matplotlib
matplotlib.use('Agg')
from galpy.util import bovy_plot, bovy_conversion
from matplotlib import pyplot
import seaborn as sns
import simulate_streampepper
from streampepperdf import streampepperdf
from optparse import OptionParser
from gd1_util import R0, V0
def get_options():
    usage = "usage: %prog [options]"
    parser = OptionParser(usage=usage)
    # make frames or combine into movie?
    parser.add_option("--movie",action="store_true", 
                      dest="movie",default=False,
                      help="If set, combine frames into a movie")
    # Output
    parser.add_option("-o",dest='outputfilename',default=None,
                      help="Name of the file that will hold the movie")
    parser.add_option("--base",dest='basefilename',default=None,
                      help="Basename of the frame files (incl. directory)")
    parser.add_option("--skip",action="store_true", 
                      dest="skip",default=False,
                      help="If set, skip existing frames")
    # Simulation options
    parser.add_option("-n",dest='nparticles',default=10000,type='int',
                      help="Number of particles to sample for each tail")
    parser.add_option("--nsnap",dest='nsnap',default=1024,type='int',
                      help="Number of snapshots to produce (must match a gd1pappernsnap.pkl file and gd1peppernsnap_trailing file)")
    # Single impact?
    parser.add_option("--single",action="store_true", 
                      dest="single",default=False,
                      help="If set, perform a single, large impact")
    parser.add_option("--singletimpact",
                      dest='singletimpact',default=4.,type='float',
                      help="Time of impact (in Gyr) for the single impact case")
    parser.add_option("--singlemimpact",
                      dest='singlemimpact',default=10.,type='float',
                      help="Mass of impact (in 10^7 msun) for the single impact case")
    return parser

def create_frames(options,args):
    # First reload the model
    with open('gd1pepper%isampling.pkl' % options.nsnap,'rb') as savefile:
        sdf_pepper_leading= pickle.load(savefile)
    with open('gd1pepper%isampling_trailing.pkl' % options.nsnap,'rb') as savefile:
        sdf_pepper_trailing= pickle.load(savefile)
    # Output times
    timpacts= sdf_pepper_leading._uniq_timpact
    # Sample unperturbed aAt
    numpy.random.seed(1)
    Oml,anglel,dtl= super(streampepperdf,sdf_pepper_leading)._sample_aAt(\
        options.nparticles)
    Omt,anglet,dtt= super(streampepperdf,sdf_pepper_trailing)._sample_aAt(\
        options.nparticles)
    # Setup progenitor
    prog= sdf_pepper_leading._progenitor().flip()
    prog.integrate(numpy.linspace(0.,9./bovy_conversion.time_in_Gyr(V0,R0),
                                  10001),sdf_pepper_leading._pot)
    prog.flip()
    # Setup impacts
    if options.single:
        # Hit the leading arm and the trailing arm 1 Gyr later
        m= options.singlemimpact/bovy_conversion.mass_in_1010msol(V0,R0)/1000.
        t= timpacts[\
            numpy.argmin(\
                numpy.fabs(\
                    numpy.array(timpacts)\
                        -options.singletimpact\
                        /bovy_conversion.time_in_Gyr(V0,R0)))]
        sdf_pepper_leading.set_impacts(\
            impactb=[0.5*simulate_streampepper.rs(options.singlemimpact*10.**7.)],
            subhalovel=numpy.array([[-25.,155.,30.]])/V0,
            impact_angle=[0.2],
            timpact=[t],
            GM=[m],rs=[simulate_streampepper.rs(options.singlemimpact*10.**7.)])
        # Trailing
        m= options.singlemimpact/bovy_conversion.mass_in_1010msol(V0,R0)/1000.
        t= timpacts[\
            numpy.argmin(\
                numpy.fabs(\
                    numpy.array(timpacts)\
                        -(options.singletimpact+1.)\
                        /bovy_conversion.time_in_Gyr(V0,R0)))]
        sdf_pepper_trailing.set_impacts(\
            impactb=[1.*simulate_streampepper.rs(options.singlemimpact*10.**7.)],
            subhalovel=numpy.array([[-25.,155.,30.]])/V0,
            impact_angle=[-0.3],
            timpact=[t],
            GM=[m],rs=[simulate_streampepper.rs(options.singlemimpact*10.**7.)])
    else:
        # Hit both with zero
        sdf_pepper_leading.set_impacts(\
            impactb=[0.],
            subhalovel=numpy.array([[-25.,155.,30.]])/V0,
            impact_angle=[0.2],
            timpact=[timpacts[0]],
            GM=[0.],rs=[1.])
        sdf_pepper_trailing.set_impacts(\
            impactb=[0.],
            subhalovel=numpy.array([[-25.,155.,30.]])/V0,
            impact_angle=[-0.2],
            timpact=[timpacts[0]],
            GM=[0.],rs=[1.])
    # Now make all frames
    bovy_plot.bovy_print(fig_height=3.,fig_width=7.)
    xlabel= r'$X_{\mathrm{orb}}\,(\mathrm{kpc})$'
    ylabel= r'$Y_{\mathrm{orb}}\,(\mathrm{kpc})$'
    xrange= [-12.,12.]
    yrange= [-1.5,.3]
    TL= _projection_orbplane(prog)
    for ii in tqdm.tqdm(range(len(timpacts))):
        if options.skip and os.path.exists(options.basefilename+'_%s.png'\
                                    % str(ii).zfill(5)):
            continue
        timpact= timpacts[-ii-1]
        overplot= False
        for sdf_pepper,Om,angle,dt \
                in zip([sdf_pepper_leading,sdf_pepper_trailing],
                       [Oml,Omt],[anglel,anglet],[dtl,dtt]):
            tOm= copy.deepcopy(Om)
            tangle= copy.deepcopy(angle)
            tdt= copy.deepcopy(dt)
            if timpact > sdf_pepper._timpact[-1]:
                # No impact yet
                tangle-= tOm*timpact
            else:
                tangle-= tOm*timpact
                # Apply all kicks relevant for this impact 
                # (copied from streampepperdf)
                dangle_at_impact= angle\
                    -numpy.tile(sdf_pepper._progenitor_angle.T,
                                (options.nparticles,1)).T\
                                -(tOm-numpy.tile(sdf_pepper._progenitor_Omega.T,
                                                (options.nparticles,1)).T)\
                                                *sdf_pepper._timpact[-1]
                dangle_par_at_impact=\
                    numpy.dot(dangle_at_impact.T,
                              sdf_pepper._dsigomeanProgDirection)\
                              *sdf_pepper._sgapdfs[-1]._gap_sigMeanSign
                dOpar= numpy.dot((tOm-numpy.tile(sdf_pepper._progenitor_Omega.T,
                                                (options.nparticles,1)).T).T,
                                 sdf_pepper._dsigomeanProgDirection)\
                                 *sdf_pepper._sgapdfs[-1]._gap_sigMeanSign
                relevant_timpact= sdf_pepper._timpact[\
                    sdf_pepper._timpact > timpact]
                for kk,ti in enumerate(relevant_timpact[::-1]):
                    # Calculate and apply kicks (points not yet released have 
                    # zero kick)
                    dOr= sdf_pepper._sgapdfs[-kk-1]._kick_interpdOr(dangle_par_at_impact)
                    dOp= sdf_pepper._sgapdfs[-kk-1]._kick_interpdOp(dangle_par_at_impact)
                    dOz= sdf_pepper._sgapdfs[-kk-1]._kick_interpdOz(dangle_par_at_impact)
                    tOm[0,:]+= dOr
                    tOm[1,:]+= dOp
                    tOm[2,:]+= dOz
                    if kk < len(relevant_timpact)-1:
                        run_to_timpact= relevant_timpact[::-1][kk+1]
                    else:
                        run_to_timpact= timpact
                    tangle[0,:]+=\
                        sdf_pepper._sgapdfs[-kk-1]._kick_interpdar(dangle_par_at_impact)\
                        +dOr*(ti-timpact)
                    tangle[1,:]+=\
                        sdf_pepper._sgapdfs[-kk-1]._kick_interpdap(dangle_par_at_impact)\
                        +dOp*(ti-timpact)
                    tangle[2,:]+=\
                        sdf_pepper._sgapdfs[-kk-1]._kick_interpdaz(dangle_par_at_impact)\
                        +dOz*(ti-timpact)
                    # Update parallel evolution
                    dOpar+=\
                        sdf_pepper._sgapdfs[-kk-1]._kick_interpdOpar(dangle_par_at_impact)
                    dangle_par_at_impact+= dOpar*(ti-run_to_timpact)
            # Convert to RvR coordinates for this time
            coorddf= copy.deepcopy(sdf_pepper._sgapdfs_coordtransform[timpact])
            coorddf._interpolate_stream_track_kick_aA()
            coorddf._interpolatedObsTrack= coorddf._kick_interpolatedObsTrack
            coorddf._ObsTrack= coorddf._gap_ObsTrack
            coorddf._interpolatedObsTrackXY= coorddf._kick_interpolatedObsTrackXY
            coorddf._ObsTrackXY= coorddf._gap_ObsTrackXY
            coorddf._allinvjacsTrack= coorddf._gap_allinvjacsTrack
            coorddf._interpolatedObsTrackAA= coorddf._kick_interpolatedObsTrackAA
            coorddf._ObsTrackAA= coorddf._gap_ObsTrackAA
            coorddf._nTrackChunks= coorddf._nTrackChunksImpact
            coorddf._thetasTrack= coorddf._gap_thetasTrack
            coorddf._interpolatedThetasTrack= coorddf._kick_interpolatedThetasTrack
            coorddf._progenitor_angle-= coorddf._progenitor_Omega*timpact
            coorddf._progenitor_angle= coorddf._progenitor_angle % (2.*numpy.pi)
            tangle= tangle % (2.*numpy.pi)
            RvR= coorddf._approxaAInv(tOm[0,:],tOm[1,:],tOm[2,:],
                                      tangle[0,:],tangle[1,:],tangle[2,:],
                                      interp=True)
            cindx= numpy.array([coorddf._find_closest_trackpointaA(tOm[0,ll],tOm[1,ll],tOm[2,ll],
                                      tangle[0,ll],tangle[1,ll],tangle[2,ll],
                                      interp=True)
                          for ll in range(len(tOm[0]))],dtype='int')
            # Progenitor and its orbit at the current time
            cprog= prog(timpact)
            cprog.integrate(numpy.linspace(0.,3.,101),sdf_pepper._pot)
            cprogf= cprog.flip()
            cprogf.integrate(numpy.linspace(0.,3.,101),sdf_pepper._pot)
            # compute the orbit and rotate everything such that the derivative
            # of the orbit points along X
            tvec= numpy.empty((3,2))
            tvec[0,0]= cprog.x(numpy.linspace(0.,3.,101)[1])
            tvec[1,0]= cprog.y(numpy.linspace(0.,3.,101)[1])
            tvec[2,0]= cprog.z(numpy.linspace(0.,3.,101)[1])
            tvec[0,1]= cprogf.x(numpy.linspace(0.,3.,101)[1])
            tvec[1,1]= cprogf.y(numpy.linspace(0.,3.,101)[1])
            tvec[2,1]= cprogf.z(numpy.linspace(0.,3.,101)[1])
            tx= numpy.dot(TL,tvec)[0]
            ty= numpy.dot(TL,tvec)[1]
            dx= tx[1]-tx[0]
            dy= ty[1]-ty[0]
            mag= numpy.sqrt(dx**2.+dy**2.)
            dx/= mag
            dy/= mag
            rot= numpy.array([[dx,dy],[-dy,dx]])
            # Plot
            indx= tdt > timpact
            # Rotate to 'orbital plane'
            tvec= numpy.empty((3,options.nparticles))
            tvec[0,:]= RvR[0]*numpy.cos(RvR[5])
            tvec[1,:]= RvR[0]*numpy.sin(RvR[5])
            tvec[2,:]= RvR[3]
            tx= numpy.dot(TL,tvec)[0]
            ty= numpy.dot(TL,tvec)[1]
            tpx= numpy.dot(TL,[cprog.x(),cprog.y(),cprog.z()])[0]
            tpy= numpy.dot(TL,[cprog.x(),cprog.y(),cprog.z()])[1]
            plotx= numpy.dot(rot,numpy.array([(tx[indx]-tpx)*R0,
                                               (ty[indx]-tpy)*R0]))[0]
            ploty= numpy.dot(rot,numpy.array([(tx[indx]-tpx)*R0,
                                               (ty[indx]-tpy)*R0]))[1]
            txrange=\
                [xrange[0]\
                     *(9.-timpact*bovy_conversion.time_in_Gyr(V0,R0))/9.-1.,
                 xrange[1]\
                     *(9.-timpact*bovy_conversion.time_in_Gyr(V0,R0))/9.+1.]
            bovy_plot.bovy_plot(plotx,
                                ploty,
                                '.',ms=2.,
                                color=sns.color_palette("colorblind")[0],
                                alpha=0.2,
                                xlabel=xlabel,
                                ylabel=ylabel,
                                xrange=txrange,
                                yrange=yrange,
                                zorder=4,
                                overplot=overplot)
            overplot= True
        # Plot progenitor orbit
        for tp in [cprog,cprogf]:
            tvec= numpy.empty((3,101))
            tvec[0]= tp.x(numpy.linspace(0.,3.,101))
            tvec[1]= tp.y(numpy.linspace(0.,3.,101))
            tvec[2]= tp.z(numpy.linspace(0.,3.,101))
            tx= numpy.dot(TL,tvec)[0]
            ty= numpy.dot(TL,tvec)[1]
            plotx= numpy.dot(rot,numpy.array([(tx-tpx)*R0,
                                              (ty-tpy)*R0]))[0]
            ploty= numpy.dot(rot,numpy.array([(tx-tpx)*R0,
                                              (ty-tpy)*R0]))[1]
            bovy_plot.bovy_plot(plotx,ploty,
                                color=sns.color_palette('colorblind')[2],
                                lw=2.,zorder=0,
                                overplot=True)
        pyplot.subplots_adjust(bottom=0.175,left=0.11,right=0.965,top=0.95)
        bovy_plot.bovy_end_print(options.basefilename+'_%s.png'\
                                     % str(ii).zfill(5))
    return None

def _projection_orbplane(prog):
    L= prog.L()[0]
    Lx= L[0]
    Ly= L[1]
    Lz= L[2]
    L= numpy.sqrt(Lx**2.+Ly**2.+Lz**2.)
    Lx/= L
    Ly/= L
    Lz/= L
    Txz= numpy.zeros((3,3))
    Tz= numpy.zeros((3,3))
    Txz[0,0]= Lx/numpy.sqrt(Lx**2.+Ly**2.)
    Txz[1,1]= Lx/numpy.sqrt(Lx**2.+Ly**2.)
    Txz[1,0]= Ly/numpy.sqrt(Lx**2.+Ly**2.)
    Txz[0,1]= -Ly/numpy.sqrt(Lx**2.+Ly**2.)
    Txz[2,2]= 1.
    Tz[0,0]= Lz
    Tz[1,1]= 1.
    Tz[2,2]= Lz
    Tz[2,0]= -numpy.sqrt(Lx**2.+Ly**2.)
    Tz[0,2]= numpy.sqrt(Lx**2.+Ly**2.)
    TL= numpy.dot(Tz,Txz)
    return TL

def create_movie(options,args):
    framerate= 25
    bitrate= 1000000
    try:
        subprocess.check_call(['ffmpeg',
                               '-i',
                               options.basefilename+'_%05d.png',
                               '-y',
                               '-framerate',str(framerate),
                               '-r',str(framerate),
                               '-b', str(bitrate),
                               options.outputfilename])
    except subprocess.CalledProcessError:
        print "'ffmpeg' failed"
    return None  

if __name__ == '__main__':
    parser= get_options()
    options,args= parser.parse_args()
    if not options.movie:
        create_frames(options,args)
    else:
        create_movie(options,args)
