from galpy.df import streamdf
from galpy.df_src.streampepperdf import streampepperdf
from galpy.orbit import Orbit
from galpy.potential import LogarithmicHaloPotential
from galpy.actionAngle import actionAngleIsochroneApprox
from galpy.util import bovy_conversion #for unit conversions
R0, V0= 8., 220.
def setup_gd1model(leading=True,
                   timpact=None,
                   hernquist=True):
    lp= LogarithmicHaloPotential(normalize=1.,q=0.9)
    aAI= actionAngleIsochroneApprox(pot=lp,b=0.8)
    obs= Orbit([1.56148083,0.35081535,-1.15481504,0.88719443,
                -0.47713334,0.12019596])
    sigv= 0.365/2. #km/s, /2 bc tdis x2
    if timpact is None:
        sdf= streamdf(sigv/220.,progenitor=obs,pot=lp,aA=aAI,leading=leading,
                      nTrackChunks=11,
                      tdisrupt=9./bovy_conversion.time_in_Gyr(V0,R0),
                      Vnorm=V0,Rnorm=R0)
    else:
        sdf= streampepperdf(sigv/220.,progenitor=obs,pot=lp,aA=aAI,
                            leading=leading,
                            nTrackChunks=11,
                            tdisrupt=9./bovy_conversion.time_in_Gyr(V0,R0),
                            Vnorm=V0,Rnorm=R0,
                            timpact=timpact,
                            spline_order=1,
                            hernquist=hernquist)
    return sdf
