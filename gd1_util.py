from galpy.df import streamdf, streamgapdf
from streampepperdf import streampepperdf
from galpy.orbit import Orbit
from galpy.potential import LogarithmicHaloPotential
from galpy.actionAngle import actionAngleIsochroneApprox
from galpy.util import bovy_conversion #for unit conversions
R0, V0= 8., 220.
def setup_gd1model(leading=True,
                   timpact=None,
                   hernquist=True,
                   age=9.,
                   singleImpact=False,
                   length_factor=1.,
                   **kwargs):
    lp= LogarithmicHaloPotential(normalize=1.,q=0.9)
    aAI= actionAngleIsochroneApprox(pot=lp,b=0.8)
    obs= Orbit([1.56148083,0.35081535,-1.15481504,0.88719443,
                -0.47713334,0.12019596])
    sigv= 0.365/2.*(9./age) #km/s, /2 bc tdis x2, adjust for diff. age
    if timpact is None:
        sdf= streamdf(sigv/220.,progenitor=obs,pot=lp,aA=aAI,leading=leading,
                      nTrackChunks=11,
                      tdisrupt=age/bovy_conversion.time_in_Gyr(V0,R0),
                      Vnorm=V0,Rnorm=R0)
    elif singleImpact:
        sdf= streamgapdf(sigv/220.,progenitor=obs,pot=lp,aA=aAI,
                         leading=leading,
                         nTrackChunks=11,
                         tdisrupt=age/bovy_conversion.time_in_Gyr(V0,R0),
                         Vnorm=V0,Rnorm=R0,
                         timpact=timpact,
                         spline_order=3,
                         hernquist=hernquist,**kwargs)
    else:
        sdf= streampepperdf(sigv/220.,progenitor=obs,pot=lp,aA=aAI,
                            leading=leading,
                            nTrackChunks=101,
                            tdisrupt=age/bovy_conversion.time_in_Gyr(V0,R0),
                            Vnorm=V0,Rnorm=R0,
                            timpact=timpact,
                            spline_order=1,
                            hernquist=hernquist,
                            length_factor=length_factor)
    sdf.turn_physical_off()
    return sdf
