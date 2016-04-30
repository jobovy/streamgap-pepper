# run_pal5_abc.py: simple ABC method for constraining Nsubhalo from Pal 5 data
import os, os.path
import glob
import csv
import time
import pickle
from optparse import OptionParser
import numpy
from numpy.polynomial import Polynomial
from scipy import interpolate, signal
from galpy.util import save_pickles, bovy_conversion, bovy_coords
import simulate_streampepper
import pal5_util
from gd1_util import R0,V0
_DATADIR= os.getenv('DATADIR')
def get_options():
    usage = "usage: %prog [options]"
    parser = OptionParser(usage=usage)
    # stream
    parser.add_option("-s",dest='streamsavefilename',
                      default=None,
                      help="Filename to save the streampepperdf object in")
    # savefilenames
    parser.add_option("--outdens",dest='outdens',default=None,
                      help="Name of the output file for the density")
    parser.add_option("--outomega",dest='outomega',default=None,
                      help="Name of the output file for the mean Omega")
    parser.add_option("-o","--abcfile",dest='abcfile',default=None,
                      help="Name of the output file for the ABC")
    parser.add_option("-b","--batch",dest='batch',default=None,
                      type='int',
                      help="If running batches of ABC simulations, batch number")
    # Parameters of the subhalos simulation
    parser.add_option("-t","--timpacts",dest='timpacts',default='256sampling',
                      help="Impact times in Gyr to consider; should be a comma separated list")
    parser.add_option("-X",dest='Xrs',default=5.,
                      type='float',
                      help="Number of times rs to consider for the impact parameter")
    parser.add_option("-l",dest='length_factor',default=1.,
                      type='float',
                      help="length_factor input to streampepperdf (consider impacts to length_factor x length)")
    parser.add_option("-M",dest='mass',default='5,9',
                      help="Mass or mass range to consider; given as log10(mass)")
    parser.add_option("--rsfac",dest='rsfac',default=1.,type='float',
                      help="Use a r_s(M) relation that is a factor of rsfac different from the fiducial one")
    parser.add_option("--plummer",action="store_true", 
                      dest="plummer",default=False,
                      help="If set, use a Plummer DM profile rather than Hernquist")
    parser.add_option("--age",dest='age',default=5.,type='float',
                      help="Age of the stream in Gyr")
    # Parallel angles at which to compute stuff
    parser.add_option("--amin",dest='amin',default=0.,
                      type='float',
                      help="Minimum parallel angle to consider")
    parser.add_option("--amax",dest='amax',default=1.75,
                      type='float',
                      help="Maximum parallel angle to consider (default: 2*meandO*mintimpact)")
    parser.add_option("--da",dest='dapar',default=0.01,
                      type='float',
                      help="Step in apar to use")
    # Data handling and continuum normalization
    parser.add_option("--polydeg",dest='polydeg',default=0,
                      type='int',
                      help="Polynomial order to fit to smooth stream density")
    parser.add_option("--minxi",dest='minxi',default=4.05,
                      type='float',
                      help="Minimum xi to consider")   
    parser.add_option("--maxxi",dest='maxxi',default=14.35,
                      type='float',
                      help="Maximum xi to consider")   
    parser.add_option("--nerrsim",dest='nerrsim',default=100,
                      type='int',
                      help="Simulate this many realizations of the errors per rate simulation")   
    parser.add_option("-m",dest='mockfilename',
                      default=None,
                      help="If set, filename of a mock Pal 5 simulation to use instead of real data")
    # Parameters of the ABC simulation
    parser.add_option("--ratemin",dest='ratemin',default=-1.,
                      type='float',
                      help="Minimum rate compared to CDM expectation; in log10")
    parser.add_option("--ratemax",dest='ratemax',default=1.,
                      type='float',
                      help="Maximum rate compared to CDM expectation; in log10")
    parser.add_option("-n","--nsamples",dest='nsamples',default=100,
                      type='int',
                      help="Number of simulations to run")
    parser.add_option("-r","--recompute",action="store_true", 
                      dest="recompute",default=False,
                      help="If set, do not run simulations, but recompute the tatistics for existing densities")
    return parser

def load_abc(filename):
    """
    NAME:
       load_abc
    PURPOSE:
       Load all ABC runs for a given filename (all batches)
    INPUT:
       filename - filename w/o batch
    OUTPUT:
       array with ABC outputs
    HISTORY:
       2016-04-10 - Written - Bovy (UofT)
    """
    allfilenames= glob.glob(filename.replace('.dat','.*.dat'))
    out= numpy.loadtxt(filename,delimiter=',')
    for fname in allfilenames:
        out= numpy.vstack((out,numpy.loadtxt(fname,delimiter=',')))
    return out

# Convert track to xi, eta
def convert_dens_to_obs(sdf_pepper,apars,
                        dens,mO,dens_smooth,minxi=0.25,maxxi=14.35):
    """
    NAME:
        convert_dens_to_obs
    PURPOSE:
        Convert track to observed coordinates
    INPUT:
        sdf_pepper - streampepperdf object
        apars - parallel angles
        dens - density(apars)
        dens_smooth - smooth density(apars)
        mO= (None) mean parallel frequency (1D) 
            [needs to be set to get density on same grid as track]
        minxi= (0.25) minimum xi to consider
    OUTPUT:
        (xi,dens/smooth)
    """
    mT= sdf_pepper.meanTrack(apars,_mO=mO,coord='lb')
    mradec= bovy_coords.lb_to_radec(mT[0],mT[1],degree=True)
    mxieta= pal5_util.radec_to_pal5xieta(mradec[:,0],mradec[:,1],degree=True)
    outll= numpy.arange(minxi,maxxi,0.1)
    # Interpolate density
    ipll= interpolate.InterpolatedUnivariateSpline(mxieta[:,0],apars)
    ipdens= interpolate.InterpolatedUnivariateSpline(apars,dens/dens_smooth)
    return (outll,ipdens(ipll(outll)))

def setup_densOmegaWriter(apar,options):
    outdens= options.outdens
    outomega= options.outomega
    if not options.batch is None:
        outdens= outdens.replace('.dat','.%i.dat' % options.batch)
    if not options.batch is None:
        outomega= outomega.replace('.dat','.%i.dat' % options.batch)
    if os.path.exists(outdens):
        # First read the file to check apar
        apar_file= numpy.genfromtxt(outdens,delimiter=',',max_rows=1)
        assert numpy.amax(numpy.fabs(apar_file-apar)) < 10.**-5., 'apar according to options does not correspond to apar already in outdens'
        apar_file= numpy.genfromtxt(outomega,delimiter=',',max_rows=1)
        assert numpy.amax(numpy.fabs(apar_file-apar)) < 10.**-5., 'apar according to options does not correspond to apar already in outomega'
        csvdens= open(outdens,'a')
        csvomega= open(outomega,'a')       
        denswriter= csv.writer(csvdens,delimiter=',')
        omegawriter= csv.writer(csvomega,delimiter=',')
    else:
        csvdens= open(outdens,'w')
        csvomega= open(outomega,'w')
        denswriter= csv.writer(csvdens,delimiter=',')
        omegawriter= csv.writer(csvomega,delimiter=',')
        # First write apar
        denswriter.writerow([a for a in apar])
        omegawriter.writerow([a for a in apar])
        csvdens.flush()
        csvomega.flush()
    return (denswriter,omegawriter,csvdens,csvomega)

def process_pal5_densdata(options):
    # Read and prep data
    backg= 400.
    data= numpy.loadtxt('data/ibata_fig7b_raw.dat',delimiter=',')
    sindx= numpy.argsort(data[:,0])
    data= data[sindx]
    data_lowerr= numpy.loadtxt('data/ibata_fig7b_rawlowerr.dat',delimiter=',')
    sindx= numpy.argsort(data_lowerr[:,0])
    data_lowerr= data_lowerr[sindx]
    data_uperr= numpy.loadtxt('data/ibata_fig7b_rawuperr.dat',delimiter=',')
    sindx= numpy.argsort(data_uperr[:,0])
    data_uperr= data_uperr[sindx]
    data_err= 0.5*(data_uperr-data_lowerr)
    # CUTS
    indx= (data[:,0] > options.minxi-0.05)*(data[:,0] < options.maxxi)
    data= data[indx]
    data_lowerr= data_lowerr[indx]
    data_uperr= data_uperr[indx]
    data_err= data_err[indx]
    # Compute power spectrum
    tdata= data[:,1]-backg
    pp= Polynomial.fit(data[:,0],tdata,deg=options.polydeg,w=1./data_err[:,1])
    tdata/= pp(data[:,0])
    ll= data[:,0]
    py= signal.csd(tdata,tdata,fs=1./(ll[1]-ll[0]),scaling='spectrum',
                   nperseg=len(ll))[1]
    py= py.real
    return (numpy.sqrt(py*(ll[-1]-ll[0])),data_err[:,1]/pp(data[:,0]))

def process_mock_densdata(options):
    # Read and prep data for mocks
    xvid= numpy.loadtxt(options.mockfilename)
    xv= xvid[:,:6]
    xv= xv[numpy.argsort(xvid[:,6])]
    XYZ= bovy_coords.galcenrect_to_XYZ(xv[:,0],xv[:,1],xv[:,2],
                                       Xsun=R0,Zsun=0.025)
    lbd= bovy_coords.XYZ_to_lbd(XYZ[0],XYZ[1],XYZ[2],degree=True)
    radec= bovy_coords.lb_to_radec(lbd[:,0],lbd[:,1],degree=True)
    xieta= pal5_util.radec_to_pal5xieta(radec[:,0],radec[:,1],degree=True)
    # make sure the progenitor is at (0,0)
    xieta[:,0]-= numpy.median(xieta[:,0])
    xieta[:,1]-= numpy.median(xieta[:,1])
    h,e= numpy.histogram(xieta[:,0],range=[0.2,14.3],bins=141)
    xdata= numpy.arange(0.25,14.35,0.1)
    # Compute power spectrum
    tdata= h-0.
    pp= Polynomial.fit(xdata,tdata,deg=options.polydeg,w=1./numpy.sqrt(h+1.))
    tdata/= pp(xdata)
    ll= xdata
    px, py= signal.csd(tdata,tdata,fs=1./(ll[1]-ll[0]),scaling='spectrum',
                       nperseg=len(ll))
    py= py.real
    px= 1./px
    py= numpy.sqrt(py*(ll[-1]-ll[0]))
    return (numpy.sqrt(py*(ll[-1]-ll[0])),numpy.sqrt(h+1.)/pp(xdata))

def pal5_abc(sdf_pepper,sdf_smooth,options):
    """
    """
    # Setup apar grid
    apar= numpy.arange(options.amin,options.amax,options.dapar)
    dens_unp= numpy.array([sdf_smooth._density_par(a) for a in apar])
    if options.recompute:
        # Load density and omega from file
        outdens= options.outdens
        outomega= options.outomega
        if not options.batch is None:
            outdens= outdens.replace('.dat','.%i.dat' % options.batch)
        if not options.batch is None:
            outomega= outomega.replace('.dat','.%i.dat' % options.batch)
        densdata= numpy.genfromtxt(outdens,delimiter=',',skip_header=1)
        omegadata= numpy.genfromtxt(outomega,delimiter=',',skip_header=1)
        nd= 0
    else:
        # Setup saving of the densities and mean Omegas
        denswriter, omegawriter, csvdens, csvomega=\
            setup_densOmegaWriter(apar,options)
        # Setup sampling
        massrange= simulate_streampepper.parse_mass(options.mass)
        rs= simulate_streampepper.rs
        sample_GM= lambda: (10.**((-0.5)*massrange[0])\
                            +(10.**((-0.5)*massrange[1])\
                              -10.**((-0.5)*massrange[0]))\
                            *numpy.random.uniform())**(1./(-0.5))\
            /bovy_conversion.mass_in_msol(V0,R0)
        sample_rs= lambda x: rs(x*bovy_conversion.mass_in_1010msol(V0,R0)*10.**10.,
                                plummer=options.plummer)
        rate_range= numpy.arange(massrange[0]+0.5,massrange[1]+0.5,1)
        cdmrate= numpy.sum([simulate_streampepper.\
                            dNencdm(sdf_pepper,10.**r,Xrs=options.Xrs,
                                    plummer=options.plummer,
                                    rsfac=options.rsfac)
                            for r in rate_range])
        print "Using an overall CDM rate of %f" % cdmrate
    # Load Pal 5 data to compare to
    if options.mockfilename is None:
        power_data, data_err= process_pal5_densdata(options)
    else:
        power_data, data_err= process_mock_densdata(options)
    # Run ABC
    while True:
        if not options.recompute:
            # Simulate a rate
            l10rate= (numpy.random.uniform()*(options.ratemax-options.ratemin)
                      +options.ratemin)
            rate= 10.**l10rate*cdmrate
            print l10rate, rate
            # Simulate
            sdf_pepper.simulate(rate=rate,sample_GM=sample_GM,sample_rs=sample_rs,
                                Xrs=options.Xrs)
            # Compute density and meanOmega and save
            try:
                densOmega= numpy.array([\
                    sdf_pepper._densityAndOmega_par_approx(a) for a in apar]).T
            except IndexError: # no hit
                dens= numpy.array([sdf_smooth._density_par(a) for a in apar])
                omega= numpy.array([sdf_smooth.meanOmega(a,oned=True) for a in apar])
            else:
                dens= densOmega[0]
                omega= densOmega[1]
            write_dens= [l10rate]
            write_omega= [l10rate]
            write_dens.extend(list(dens))
            write_omega.extend(list(omega))
            denswriter.writerow(write_dens)
            omegawriter.writerow(write_omega)
            csvdens.flush()
            csvomega.flush()
        else:
            if nd >= len(densdata): break
            l10rate= densdata[nd,0]
            dens = densdata[nd,1:]
            omega= omegadata[nd,1:]
            nd+= 1
        # Convert density to observed density
        xixi,dens= convert_dens_to_obs(sdf_pepper,apar,
                                        dens,omega,dens_unp,
                                        minxi=options.minxi,
                                        maxxi=options.maxxi)
        # Add errors (Rao-Blackwellize...)
        for ee in range(options.nerrsim):
            tdens= dens+numpy.random.normal(size=len(xixi))*data_err
            # Compute power spectrum
            tcsd= signal.csd(tdens,tdens,fs=1./(xixi[1]-xixi[0]),
                             scaling='spectrum',nperseg=len(xixi))[1].real
            power= numpy.sqrt(tcsd*(xixi[-1]-xixi[0]))
            yield (l10rate,
                   numpy.fabs(power[1]-power_data[1]),
                   numpy.fabs(power[2]-power_data[2]),
                   numpy.fabs(power[3]-power_data[3]),
                   numpy.fabs(numpy.log(numpy.mean(tdens[7:17])\
                                            /numpy.mean(tdens[107:117]))),
                   ee)

def abcsims(sdf_pepper,sdf_smooth,options):
    """
    NAME:
       abcsims
    PURPOSE:
       Run a bunch of ABC simulations
    INPUT:
       sdf_pepper - streampepperdf object to compute peppering
       sdf_smooth - streamdf object for smooth stream
       options - the options dictionary
    OUTPUT:
       (none; just saves the simulations to a file)
    HISTORY:
       2016-04-08 - Written - Bovy (UofT)
    """
    print("Running ABC sims ...")
    abcfile= options.abcfile
    if not options.batch is None:
        abcfile= abcfile.replace('.dat','.%i.dat' % options.batch)
    if os.path.exists(abcfile):
        # First read the file to check apar
        csvabc= open(abcfile,'a')
        abcwriter= csv.writer(csvabc,delimiter=',')
    else:
        csvabc= open(abcfile,'w')
        abcwriter= csv.writer(csvabc,delimiter=',')
    nit= 0
    for sim in pal5_abc(sdf_pepper,sdf_smooth,options):
        abcwriter.writerow(list(sim)[:-1])
        csvabc.flush()
        nit+= 1
        if nit >= options.nerrsim*options.nsamples: break
    return None

def recompute(sdf_pepper,sdf_smooth,options):
    """
    NAME:
       recompute
    PURPOSE:
       Recompute the ABC summaries for existing simulations
    INPUT:
       sdf_pepper - streampepperdf object to compute peppering
       sdf_smooth - streamdf object for smooth stream
       options - the options dictionary
    OUTPUT:
       (none; just saves the simulations to a file)
    HISTORY:
       2016-04-14 - Written - Bovy (UofT)
    """
    print("Recomputing ABC sims ...")
    abcfile= options.abcfile
    if not options.batch is None:
        abcfile= abcfile.replace('.dat','.%i.dat' % options.batch)
    if os.path.exists(abcfile):
        raise IOError("ERROR: abcfile already exists, would be overridden...")
    else:
        csvabc= open(abcfile,'w')
        abcwriter= csv.writer(csvabc,delimiter=',')
    for sim in pal5_abc(sdf_pepper,sdf_smooth,options):
        abcwriter.writerow(list(sim)[:-1])
        csvabc.flush()
    return None

if __name__ == '__main__':
    parser= get_options()
    options,args= parser.parse_args()
    # Setup the streampepperdf object
    if not os.path.exists(options.streamsavefilename):
        timpacts= simulate_streampepper.parse_times(\
            options.timpacts,options.age)
        sdf_smooth= pal5_util.setup_pal5model(age=options.age)
        sdf_pepper= pal5_util.setup_pal5model(timpact=timpacts,
                                              hernquist=not options.plummer,
                                              age=options.age,
                                              length_factor=options.length_factor)
        save_pickles(options.streamsavefilename,sdf_smooth,sdf_pepper)
    else:
        with open(options.streamsavefilename,'rb') as savefile:
            sdf_smooth= pickle.load(savefile)
            sdf_pepper= pickle.load(savefile)
    if options.recompute:
        recompute(sdf_pepper,sdf_smooth,options)
    else:
        abcsims(sdf_pepper,sdf_smooth,options)
