# custom_stripping_df.py: streamgapdf for Jason's custom stripping time distribution
import numpy
import galpy.df_src.streamgapdf
import streampepperdf
import galpy.df_src.streamdf
from galpy.util import bovy_conversion
class streamdf_jason(galpy.df_src.streamdf.streamdf):
    def _sample_aAt(self,n):
        return custom_sample_aAt(self,n)
class streamgapdf_jason(galpy.df_src.streamgapdf.streamgapdf):
    def _sample_aAt(self,n):
        return custom_sample_aAt_gap(self,n)
class streamgapdf_jason_onlyopar(galpy.df_src.streamgapdf.streamgapdf):
    """Same as above, but only apply dOmega kick"""
    def _sample_aAt(self,n):
        return custom_sample_aAt_gap_onlyopar(self,n)
class streampepperdf_jason(streampepperdf.streampepperdf):
    def _sample_aAt(self,n):
        return custom_sample_aAt_pepper(self,n)
class streampepperdf_jason_onlyopar(streampepperdf.streampepperdf):
    """Same as above, but only apply dOmega kick"""
    def _sample_aAt(self,n):
        return custom_sample_aAt_pepper_onlyopar(self,n)
# Pal 5
class streamdf_pal5_jason(galpy.df_src.streamdf.streamdf):
    def _sample_aAt(self,n):
        return custom_sample_aAt_pal5(self,n)

def custom_sample_aAt(sdf,n):
    """Custom time and frequency sampling based on Jason's investigations of Denis' simulations"""
    # Sample pericenter vs. continuous
    w= 0.28
    atperi= numpy.random.binomial(1,1.-w,size=n)
    # Setup output
    dO1s= numpy.zeros(n)
    dt= numpy.zeros(n)
    # First go through the pericenter ones: 
    # find all of the pericenter passages up to tdisrupt
    nperi= sdf._tdisrupt/sdf._progenitor_Omega[0]+2
    tps= sdf._progenitor_angle[0]/sdf._progenitor_Omega[0]\
        +numpy.arange(nperi)*2.*numpy.pi/sdf._progenitor_Omega[0]
    tps= tps[tps < sdf._tdisrupt]
    # Sample these with a linearly-increasing probability
    dt[atperi==1]= tps[numpy.random.choice(len(tps),
                                           size=numpy.sum(atperi),
                                           p=tps/numpy.sum(tps))]
    # Generate frequencies from the broad Gaussian
    dO1s[atperi==1]= numpy.random.normal(size=numpy.sum(atperi))\
        *0.049+0.306
    # Next do the continuous part
    dt[atperi==0]= (numpy.random.uniform(size=n-numpy.sum(atperi)))**(1./2.)\
        *sdf._tdisrupt
    dO1s[atperi==0]= numpy.random.normal(size=n-numpy.sum(atperi))\
        *0.023+0.243
    # Put everything together and rotate etc.
    dO1s/= bovy_conversion.freq_in_kmskpc(sdf._vo,sdf._ro)
    dO1s= numpy.array(dO1s)*sdf._sigMeanSign
    dO2s= numpy.random.normal(size=n)*numpy.sqrt(sdf._sortedSigOEig[1])
    dO3s= numpy.random.normal(size=n)*numpy.sqrt(sdf._sortedSigOEig[0])
    #Rotate into dOs in R,phi,z coordinates
    dO= numpy.vstack((dO3s,dO2s,dO1s))
    dO= numpy.dot(sdf._sigomatrixEig[1][:,sdf._sigomatrixEigsortIndx],
                  dO)
    Om= dO+numpy.tile(sdf._progenitor_Omega.T,(n,1)).T
    #Also generate angles
    da= numpy.random.normal(size=(3,n))*sdf._sigangle
    #Integrate the orbits relative to the progenitor
    da+= dO*numpy.tile(dt,(3,1))
    angle= da+numpy.tile(sdf._progenitor_angle.T,(n,1)).T
    return (Om,angle,dt)

def custom_sample_aAt_pal5(sdf,n):
    """Custom time and frequency sampling based on Jason's investigations of Denis' simulation of Pal 5's disruption"""
    # Sample pericenter vs. continuous
    w= 0.267
    atperi= numpy.random.binomial(1,1.-w,size=n)
    # Setup output
    dO1s= numpy.zeros(n)
    dt= numpy.zeros(n)
    # First go through the pericenter ones: 
    # find all of the pericenter passages up to tdisrupt
    nperi= sdf._tdisrupt/sdf._progenitor_Omega[0]+2
    tps= sdf._progenitor_angle[0]/sdf._progenitor_Omega[0]\
        +numpy.arange(nperi)*2.*numpy.pi/sdf._progenitor_Omega[0]
    tps= tps[tps < sdf._tdisrupt]
    # Sample these with a uniform probability
    dt[atperi==1]= tps[numpy.random.choice(len(tps),
                                           size=numpy.sum(atperi))]
    # Generate frequencies from the broad Gaussian
    dO1s[atperi==1]= numpy.random.normal(size=numpy.sum(atperi))\
        *0.072+0.370
    # Next do the continuous part
    dt[atperi==0]= (numpy.random.uniform(size=n-numpy.sum(atperi)))\
        *sdf._tdisrupt
    dO1s[atperi==0]= numpy.random.normal(size=n-numpy.sum(atperi))\
        *0.033+0.277
    # Put everything together and rotate etc.
    dO1s/= bovy_conversion.freq_in_kmskpc(sdf._vo,sdf._ro)
    dO1s= numpy.array(dO1s)*sdf._sigMeanSign
    dO2s= numpy.random.normal(size=n)*numpy.sqrt(sdf._sortedSigOEig[1])
    dO3s= numpy.random.normal(size=n)*numpy.sqrt(sdf._sortedSigOEig[0])
    #Rotate into dOs in R,phi,z coordinates
    dO= numpy.vstack((dO3s,dO2s,dO1s))
    dO= numpy.dot(sdf._sigomatrixEig[1][:,sdf._sigomatrixEigsortIndx],
                  dO)
    Om= dO+numpy.tile(sdf._progenitor_Omega.T,(n,1)).T
    #Also generate angles
    da= numpy.random.normal(size=(3,n))*sdf._sigangle
    #Integrate the orbits relative to the progenitor
    da+= dO*numpy.tile(dt,(3,1))
    angle= da+numpy.tile(sdf._progenitor_angle.T,(n,1)).T
    return (Om,angle,dt)

def custom_sample_aAt_gap(sdf,n):
    Om,angle,dt= custom_sample_aAt(sdf,n)
    # Copied from streamgapdf
    # Now rewind angles by timpact, apply the kicks, and run forward again
    dangle_at_impact= angle-numpy.tile(sdf._progenitor_angle.T,(n,1)).T\
        -(Om-numpy.tile(sdf._progenitor_Omega.T,(n,1)).T)*sdf._timpact
    dangle_par_at_impact= numpy.dot(dangle_at_impact.T,
                                    sdf._dsigomeanProgDirection)\
                                    *sdf._gap_sigMeanSign
    # Calculate and apply kicks (points not yet released have zero kick)
    dOr= sdf._kick_interpdOr(dangle_par_at_impact)
    dOp= sdf._kick_interpdOp(dangle_par_at_impact)
    dOz= sdf._kick_interpdOz(dangle_par_at_impact)
    Om[0,:]+= dOr
    Om[1,:]+= dOp
    Om[2,:]+= dOz
    angle[0,:]+=\
        sdf._kick_interpdar(dangle_par_at_impact)+dOr*sdf._timpact
    angle[1,:]+=\
        sdf._kick_interpdap(dangle_par_at_impact)+dOp*sdf._timpact
    angle[2,:]+=\
        sdf._kick_interpdaz(dangle_par_at_impact)+dOz*sdf._timpact
    return (Om,angle,dt)
        
def custom_sample_aAt_gap_onlyopar(sdf,n):
    Om,angle,dt= custom_sample_aAt(sdf,n)
    # Copied from streamgapdf
    # Now rewind angles by timpact, apply the kicks, and run forward again
    dangle_at_impact= angle-numpy.tile(sdf._progenitor_angle.T,(n,1)).T\
        -(Om-numpy.tile(sdf._progenitor_Omega.T,(n,1)).T)*sdf._timpact
    dangle_par_at_impact= numpy.dot(dangle_at_impact.T,
                                    sdf._dsigomeanProgDirection)\
                                    *sdf._gap_sigMeanSign
    # Calculate and apply kicks (points not yet released have zero kick)
    dO= (numpy.tile(sdf._kick_interpdOpar(dangle_par_at_impact),
                    (3,1))
         *numpy.tile(sdf._dsigomeanProgDirection,(n,1)).T
         *sdf._gap_sigMeanSign)
    Om+= dO 
    angle+= dO*sdf._timpact
    return (Om,angle,dt)
        
def custom_sample_aAt_pepper(sdf,n):
    Om,angle,dt= custom_sample_aAt(sdf,n)
    # Copied from streampepperdf
    # Now rewind angles to the first impact, then apply all kicks,          
    # and run forward again                                                 
    dangle_at_impact= angle-numpy.tile(sdf._progenitor_angle.T,(n,1)).T\
        -(Om-numpy.tile(sdf._progenitor_Omega.T,(n,1)).T)\
        *sdf._timpact[-1]
    dangle_par_at_impact=\
        numpy.dot(dangle_at_impact.T,
                  sdf._dsigomeanProgDirection)\
                  *sdf._sgapdfs[-1]._gap_sigMeanSign
    dOpar= numpy.dot((Om-numpy.tile(sdf._progenitor_Omega.T,(n,1)).T).T,
                     sdf._dsigomeanProgDirection)\
                     *sdf._sgapdfs[-1]._gap_sigMeanSign
    for kk,timpact in enumerate(sdf._timpact[::-1]):
        # Calculate and apply kicks (points not yet released have           
        # zero kick)                                                        
        dOr= sdf._sgapdfs[-kk-1]._kick_interpdOr(dangle_par_at_impact)
        dOp= sdf._sgapdfs[-kk-1]._kick_interpdOp(dangle_par_at_impact)
        dOz= sdf._sgapdfs[-kk-1]._kick_interpdOz(dangle_par_at_impact)
        Om[0,:]+= dOr
        Om[1,:]+= dOp
        Om[2,:]+= dOz
        if kk < len(sdf._timpact)-1:
            run_to_timpact= sdf._timpact[::-1][kk+1]
        else:
            run_to_timpact= 0.
        angle[0,:]+=\
            sdf._sgapdfs[-kk-1]._kick_interpdar(dangle_par_at_impact)\
            +dOr*timpact
        angle[1,:]+=\
            sdf._sgapdfs[-kk-1]._kick_interpdap(dangle_par_at_impact)\
            +dOp*timpact
        angle[2,:]+=\
            sdf._sgapdfs[-kk-1]._kick_interpdaz(dangle_par_at_impact)\
            +dOz*timpact
        # Update parallel evolution                                         
        dOpar+=\
            sdf._sgapdfs[-kk-1]._kick_interpdOpar(dangle_par_at_impact)
        dangle_par_at_impact+= dOpar*(timpact-run_to_timpact)
    return (Om,angle,dt)

def custom_sample_aAt_pepper_onlyopar(sdf,n):
    Om,angle,dt= custom_sample_aAt(sdf,n)
    # Copied from streampepperdf
    # Now rewind angles to the first impact, then apply all kicks,          
    # and run forward again                                                 
    dangle_at_impact= angle-numpy.tile(sdf._progenitor_angle.T,(n,1)).T\
        -(Om-numpy.tile(sdf._progenitor_Omega.T,(n,1)).T)\
        *sdf._timpact[-1]
    dangle_par_at_impact=\
        numpy.dot(dangle_at_impact.T,
                  sdf._dsigomeanProgDirection)\
                  *sdf._sgapdfs[-1]._gap_sigMeanSign
    dOpar= numpy.dot((Om-numpy.tile(sdf._progenitor_Omega.T,(n,1)).T).T,
                     sdf._dsigomeanProgDirection)\
                     *sdf._sgapdfs[-1]._gap_sigMeanSign
    for kk,timpact in enumerate(sdf._timpact[::-1]):
        # Calculate and apply kicks (points not yet released have           
        # zero kick)                                                        
        dO= (numpy.tile(\
                sdf._sgapdfs[-kk-1]._kick_interpdOpar(dangle_par_at_impact),
                (3,1))
             *numpy.tile(sdf._dsigomeanProgDirection,(n,1)).T
             *sdf._sgapdfs[-kk-1]._gap_sigMeanSign)
        Om+= dO 
        if kk < len(sdf._timpact)-1:
            run_to_timpact= sdf._timpact[::-1][kk+1]
        else:
            run_to_timpact= 0.
        angle+= dO*timpact
        # Update parallel evolution                                         
        dOpar+=\
            sdf._sgapdfs[-kk-1]._kick_interpdOpar(dangle_par_at_impact)
        dangle_par_at_impact+= dOpar*(timpact-run_to_timpact)
    return (Om,angle,dt)

