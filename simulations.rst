Commands to run stream-pepper simulation
==========================================

Single times:
-------------

Commands like::

	 python simulate_streampepper.py --outdens=/Users/bovy/data/streamgap-pepper/gd1_onetime/gd1_t1.0_6.5_dens.dat --outomega=/Users/bovy/data/streamgap-pepper/gd1_onetime/gd1_t1.0_6.5_omega.dat -t 1. -M 6.5 --dt=10.

and::

	 python simulate_streampepper.py --outdens=/Users/bovy/data/streamgap-pepper/gd1_onetime/gd1_t0.81peri_5.5_dens.dat --outomega=/Users/bovy/data/streamgap-pepper/gd1_onetime/gd1_t0.81peri_5.5_omega.dat -t 0.81 -M 5.5 --dt=10.

Multiple times:
----------------

Commands like::

	 python simulate_streampepper.py --outdens=$DATADIR/bovy/streamgap-pepper/gd1_multtime/gd1_t2sampling_6.5_dens.dat --outomega=$DATADIR/bovy/streamgap-pepper/gd1_multtime/gd1_t2sampling_6.5_omega.dat -t 2sampling -M 6.5 --dt=10.


Determine X times rs factor:
----------------------------

Commands like::

	 python simulate_streampepper.py --outdens=$DATADIR/bovy/streamgap-pepper/gd1_multtime/gd1_t64sampling_X5_6.5_dens.dat --outomega=$DATADIR/bovy/streamgap-pepper/gd1_multtime/gd1_t64sampling_X5_6.5_omega.dat -t 64sampling -M 6.5 --dt=10. -X 5.

Full mass range
---------------

Do::

	python simulate_streampepper.py --outdens=$DATADIR/bovy/streamgap-pepper/gd1_multtime/gd1_t64sampling_X5_5-9_dens.dat --outomega=$DATADIR/bovy/streamgap-pepper/gd1_multtime/gd1_t64sampling_X5_5-9_omega.dat -t 64sampling -M 5,9 --dt=350. -X 5.

Plummer
-------

Do::

	python simulate_streampepper.py --outdens=$DATADIR/bovy/streamgap-pepper/gd1_multtime/gd1_t64sampling_X5_plum_5-9_dens.dat --outomega=$DATADIR/bovy/streamgap-pepper/gd1_multtime/gd1_t64sampling_X5_plum_5-9_omega.dat -t 64sampling -M 5,9 --dt=350. -X 5. --plummer
