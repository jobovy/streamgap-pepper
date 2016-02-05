# script to run simulations of stream peppering
import os, os.path
import csv
import time
import numpy
from optparse import OptionParser
from galpy.util import bovy_conversion
import gd1_util
from gd1_util import R0,V0
_DATADIR= os.getenv('DATADIR')
def get_options():
    usage = "usage: %prog [options]"
    parser = OptionParser(usage=usage)
    # stream
    parser.add_option("-s","--stream",dest='stream',default='gd1like',
                      help="Stream to consider")
    # savefilenames
    parser.add_option("--outdens",dest='outdens',default=None,
                      help="Name of the output file for the density")
    parser.add_option("--outomega",dest='outomega',default=None,
                      help="Name of the output file for the mean Omega")
    # Parameters of this simulation
    parser.add_option("-t","--timpacts",dest='timpacts',default=None,
                      help="Impact times in Gyr to consider; should be a comma separated list")
    parser.add_option("-X",dest='Xrs',default=3.,
                      type='float',
                      help="Number of times rs to consider for the impact parameter")
    parser.add_option("-M",dest='mass',default='6.5',
                      help="Mass or mass range to consider; given as log10(mass)")
    parser.add_option("--timescdm",dest='timescdm',default=1.,
                      help="Use a rate that is timescdm times the CDM prediction")
    parser.add_option("--rsfac",dest='rsfac',default=1.,type='float',
                      help="Use a r_s(M) relation that is a factor of rsfac different from the fiducial one")
    parser.add_option("--plummer",action="store_true", 
                      dest="plummer",default=False,
                      help="If set, use a Plummer DM profile rather than Hernquist")
    # Parallel angles at which to compute stuff
    parser.add_option("--amin",dest='amin',default=None,
                      type='float',
                      help="Minimum parallel angle to consider")
    parser.add_option("--amax",dest='amax',default=None,
                      type='float',
                      help="Maximum parallel angle to consider (default: 2*meandO*mintimpact)")
    parser.add_option("--da",dest='dapar',default=0.003,
                      type='float',
                      help="Step in apar to use")
    parser.add_option("--dt",dest='dt',default=60.,
                      type='float',
                      help="Number of minutes to run simulations for")
    return parser

def parse_times(times):
    if 'sampling' in times:
        nsam= int(times.split('sampling')[0])
        return [float(ti)/bovy_conversion.time_in_Gyr(V0,R0)
                for ti in numpy.arange(1,nsam+1)/(nsam+1.)*9.]
    return [float(ti)/bovy_conversion.time_in_Gyr(V0,R0)
            for ti in times.split(',')]
def parse_mass(mass):   
    return [float(m) for m in mass.split(',')]

# Functions to sample
def nsubhalo(m):
    return 0.3*(10.**6.5/m)
def rs(m,plummer=False,rsfac=1.):
    if plummer:
        return 1.62*rsfac/R0*(m/10.**8.)**0.5
    else:
        return 1.05*rsfac/R0*(m/10.**8.)**0.5
def dNencdm(sdf_pepper,m,Xrs=3.,plummer=False,rsfac=1.):
    return sdf_pepper.subhalo_encounters(\
        sigma=120./220.,nsubhalo=nsubhalo(m),
        bmax=Xrs*rs(m,plummer=plummer,rsfac=rsfac))

# Function to run the simulations
def run_simulations(sdf_pepper,sdf_smooth,options):
    # Setup apar grid
    apar= numpy.arange(options.amin,options.amax,options.dapar)
    # Check whether the output files already exist and if so, get the amin, amax, da from them
    if os.path.exists(options.outdens):
        # First read the file to check apar
        apar_file= numpy.genfromtxt(options.outdens,delimiter=',',max_rows=1)
        print numpy.amax(numpy.fabs(apar_file-apar))
        assert numpy.amax(numpy.fabs(apar_file-apar)) < 10.**-5., 'apar according to options does not correspond to apar already in outdens'
        apar_file= numpy.genfromtxt(options.outomega,delimiter=',',max_rows=1)
        print numpy.amax(numpy.fabs(apar_file-apar))
        assert numpy.amax(numpy.fabs(apar_file-apar)) < 10.**-5., 'apar according to options does not correspond to apar already in outomega'
        csvdens= open(options.outdens,'a')
        csvomega= open(options.outomega,'a')       
        denswriter= csv.writer(csvdens,delimiter=',')
        omegawriter= csv.writer(csvomega,delimiter=',')
    else:
        csvdens= open(options.outdens,'w')
        csvomega= open(options.outomega,'w')
        denswriter= csv.writer(csvdens,delimiter=',')
        omegawriter= csv.writer(csvomega,delimiter=',')
        # First write apar and the smooth calculations
        denswriter.writerow([a for a in apar])
        omegawriter.writerow([a for a in apar])
        if sdf_smooth is None and options.stream.lower() == 'gd1like':
            sdf_smooth= gd1_util.setup_gd1model()
        dens_unp= [sdf_smooth._density_par(a) for a in apar]
        denswriter.writerow(dens_unp)
        omega_unp= [sdf_smooth.meanOmega(a,oned=True) for a in apar]
        omegawriter.writerow(omega_unp)
        csvdens.flush()
        csvomega.flush()
    # Parse mass
    massrange= parse_mass(options.mass)
    if len(massrange) == 1:
        sample_GM= lambda: 10.**(massrange[0]-10.)\
            /bovy_conversion.mass_in_1010msol(V0,R0)
        rate= options.timescdm*dNencdm(sdf_pepper,
                                       10.**massrange[0],Xrs=options.Xrs,
                                       plummer=options.plummer,
                                       rsfac=options.rsfac)
    elif len(massrange) == 2:
        # Sample from power-law
        sample_GM= lambda: (10.**-(massrange[0]/2.)\
                                   +(10.**-(massrange[1]/2.)\
                                         -10.**(-massrange[0]/2.))\
                                   *numpy.random.uniform())**-2.\
                                   /bovy_conversion.mass_in_msol(V0,R0)
        rate_range= numpy.arange(massrange[0]+0.5,massrange[1]+0.5,1)
        rate= options.timescdm\
            *numpy.sum([dNencdm(sdf_pepper,10.**r,Xrs=options.Xrs,
                                plummer=options.plummer,rsfac=options.rsfac)
                        for r in rate_range])
    print "Using an overall rate of %f" % rate
    sample_rs= lambda x: rs(x*bovy_conversion.mass_in_1010msol(V0,R0)*10.**10.,
                            plummer=options.plummer,rsfac=options.rsfac)
    # Simulate
    start= time.time()
    while time.time() < (start+options.dt*60.):
        sdf_pepper.simulate(rate=rate,sample_GM=sample_GM,sample_rs=sample_rs,
                            Xrs=options.Xrs)
        # Compute density and meanOmega and save
        densOmega= numpy.array([sdf_pepper._densityAndOmega_par_approx(a)
                                for a in apar]).T
        denswriter.writerow(list(densOmega[0]))
        omegawriter.writerow(list(densOmega[1]))
        csvdens.flush()
        csvomega.flush()
    csvdens.close()
    csvomega.close()
    return None
        
if __name__ == '__main__':
    parser= get_options()
    options,args= parser.parse_args()
    # Setup 
    if options.stream.lower() == 'gd1like':
        timpacts= parse_times(options.timpacts)
        sdf_pepper= gd1_util.setup_gd1model(timpact=timpacts,
                                            hernquist=not options.plummer)
    # Need smooth?
    if options.amax is None or options.amin is None:
        sdf_smooth= gd1_util.setup_gd1model()
    if options.amax is None:
        options.amax= sdf_smooth.length()+options.dapar
    else:
        sdf_smooth= None
    if options.amin is None:
        options.amin= 2.*sdf_smooth.meanOmega(0.1,oned=True)\
            *numpy.amin(numpy.array(timpacts))
    run_simulations(sdf_pepper,sdf_smooth,options)
