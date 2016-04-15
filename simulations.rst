Commands to run stream-pepper simulation
==========================================

Single times:
-------------

Commands like::

	 python simulate_streampepper.py --outdens=/data/bovy/streamgap-pepper/gd1_onetime/gd1_t1.0_6.5_dens.dat --outomega=/data/bovy/streamgap-pepper/gd1_onetime/gd1_t1.0_6.5_omega.dat -t 1. -M 6.5 -n 1002 -X 5.

and::

	 python simulate_streampepper.py --outdens=/data/bovy/streamgap-pepper/gd1_onetime/gd1_t0.81peri_6.5_dens.dat --outomega=/data/bovy/streamgap-pepper/gd1_onetime/gd1_t0.81peri_6.5_omega.dat -t 0.81 -M 6.5 -n 1002 -X 5.

Multiple times:
----------------

Commands like::

	 python simulate_streampepper.py --outdens=$DATADIR/bovy/streamgap-pepper/gd1_multtime/gd1_t2sampling_X5_6.5_dens.dat --outomega=$DATADIR/bovy/streamgap-pepper/gd1_multtime/gd1_t2sampling_X5_6.5_omega.dat -t 2sampling -M 6.5 -n 1002 -X 5.


Determine X times rs factor:
----------------------------

Commands like::

	 python simulate_streampepper.py --outdens=$DATADIR/bovy/streamgap-pepper/gd1_multtime/gd1_t64sampling_X5_6.5_dens.dat --outomega=$DATADIR/bovy/streamgap-pepper/gd1_multtime/gd1_t64sampling_X5_6.5_omega.dat -t 64sampling -M 6.5 -n 1002 -X 5.

Length factor
--------------

Commands like::

	 python simulate_streampepper.py --outdens=$DATADIR/bovy/streamgap-pepper/gd1_multtime/gd1_t64sampling_lf1p25_6.5_dens.dat --outomega=$DATADIR/bovy/streamgap-pepper/gd1_multtime/gd1_t64sampling_lf1p25_6.5_omega.dat -t 64sampling -M 6.5 -n 1002 -X 5. -l 1.25 --timescdm=1.25

Full mass range
---------------

Do::

	python simulate_streampepper.py --outdens=$DATADIR/bovy/streamgap-pepper/gd1_multtime/gd1_t64sampling_X5_5-9_dens.dat --outomega=$DATADIR/bovy/streamgap-pepper/gd1_multtime/gd1_t64sampling_X5_5-9_omega.dat -t 64sampling -M 5,9 -n 1002 -X 5.

Plummer and different rs for Hernquist
---------------------------------------

Commands like::

	python simulate_streampepper.py --outdens=$DATADIR/bovy/streamgap-pepper/gd1_multtime/gd1_t64sampling_X3p24_plum_5-9_dens.dat --outomega=$DATADIR/bovy/streamgap-pepper/gd1_multtime/gd1_t64sampling_X3p24_plum_5-9_omega.dat -t 64sampling -M 5,9 -n 1002 -X 3.24 --plummer
	python simulate_streampepper.py --outdens=$DATADIR/bovy/streamgap-pepper/gd1_multtime/gd1_t64sampling_X2_rsfac2p5_5-9_dens.dat --outomega=$DATADIR/bovy/streamgap-pepper/gd1_multtime/gd1_t64sampling_X2_rsfac2p5_5-9_omega.dat -t 64sampling -M 5,9 -n 1002 -X 2. --rsfac=2.5
	python simulate_streampepper.py --outdens=$DATADIR/bovy/streamgap-pepper/gd1_multtime/gd1_t64sampling_X12p5_rsfacp4_5-9_dens.dat --outomega=$DATADIR/bovy/streamgap-pepper/gd1_multtime/gd1_t64sampling_X12p5_rsfacp4_5-9_omega.dat -t 64sampling -M 5,9 -n 1002 -X 12.5 --rsfac=0.4

Higher or lower rate than CDM
------------------------------

Commands like::

	 python simulate_streampepper.py --outdens=$DATADIR/bovy/streamgap-pepper/gd1_multtime/gd1_t64sampling_X5_cdm3_5-9_dens.dat --outomega=$DATADIR/bovy/streamgap-pepper/gd1_multtime/gd1_t64sampling_X5_cdm3_5-9_omega.dat -t 64sampling -M 5,9 -n 1002 -X 5. --timescdm=3.
	 python simulate_streampepper.py --outdens=$DATADIR/bovy/streamgap-pepper/gd1_multtime/gd1_t64sampling_X5_cdmp33_5-9_dens.dat --outomega=$DATADIR/bovy/streamgap-pepper/gd1_multtime/gd1_t64sampling_X5_cdmp33_5-9_omega.dat -t 64sampling -M 5,9 -n 1002 -X 5. --timescdm=0.33333333333
	 python simulate_streampepper.py --outdens=$DATADIR/bovy/streamgap-pepper/gd1_multtime/gd1_t64sampling_X5_cdm10_5-9_dens.dat --outomega=$DATADIR/bovy/streamgap-pepper/gd1_multtime/gd1_t64sampling_X5_cdm10_5-9_omega.dat -t 64sampling -M 5,9 -n 1002 -X 5. --timescdm=10.
	 python simulate_streampepper.py --outdens=$DATADIR/bovy/streamgap-pepper/gd1_multtime/gd1_t64sampling_X5_cdmp1_5-9_dens.dat --outomega=$DATADIR/bovy/streamgap-pepper/gd1_multtime/gd1_t64sampling_X5_cdmp1_5-9_omega.dat -t 64sampling -M 5,9 -n 1002 -X 5. --timescdm=0.1

Different age
--------------

Commands like::

	 python simulate_streampepper.py --outdens=$DATADIR/bovy/streamgap-pepper/gd1_multtime/gd1_t64sampling_X5_age4p5_5-9_dens.dat --outomega=$DATADIR/bovy/streamgap-pepper/gd1_multtime/gd1_t64sampling_X5_age4p5_5-9_omega.dat -t 64sampling -M 5,9 -n 1002 -X 5. --age=4.5
	 python simulate_streampepper.py --outdens=$DATADIR/bovy/streamgap-pepper/gd1_multtime/gd1_t64sampling_X5_age4p5_cdm2_5-9_dens.dat --outomega=$DATADIR/bovy/streamgap-pepper/gd1_multtime/gd1_t64sampling_X5_age4p5_cdm2_5-9_omega.dat -t 64sampling -M 5,9 -n 1002 -X 5. --age=4.5 --timescdm=2.
	 python simulate_streampepper.py --outdens=$DATADIR/bovy/streamgap-pepper/gd1_multtime/gd1_t64sampling_X5_age13p5_5-9_dens.dat --outomega=$DATADIR/bovy/streamgap-pepper/gd1_multtime/gd1_t64sampling_X5_age13p5_5-9_omega.dat -t 64sampling -M 5,9 -n 1002 -X 5. --age=13.5
	 python simulate_streampepper.py --outdens=$DATADIR/bovy/streamgap-pepper/gd1_multtime/gd1_t64sampling_X5_age13p5_cdmp66_5-9_dens.dat --outomega=$DATADIR/bovy/streamgap-pepper/gd1_multtime/gd1_t64sampling_X5_age13p5_cdmp66_5-9_omega.dat -t 64sampling -M 5,9 -n 1002 -X 5. --age=13.5 --timescdm=.666666666666

CDM spectrum w/ cut-off
-----------------------

Commands like::

	 python simulate_streampepper.py --outdens=$DATADIR/bovy/streamgap-pepper/gd1_multtime/gd1_t64sampling_X5_cutoff5p5_5-9_dens.dat --outomega=$DATADIR/bovy/streamgap-pepper/gd1_multtime/gd1_t64sampling_X5_cutoff5p5_5-9_omega.dat -t 64sampling -M 5,9 -n 1002 -X 5. --cutoff=5.5
	 python simulate_streampepper.py --outdens=$DATADIR/bovy/streamgap-pepper/gd1_multtime/gd1_t64sampling_X5_cutoff6p5_5-9_dens.dat --outomega=$DATADIR/bovy/streamgap-pepper/gd1_multtime/gd1_t64sampling_X5_cutoff6p5_5-9_omega.dat -t 64sampling -M 5,9 -n 1002 -X 5. --cutoff=6.5
	 python simulate_streampepper.py --outdens=$DATADIR/bovy/streamgap-pepper/gd1_multtime/gd1_t64sampling_X5_cutoff7p5_5-9_dens.dat --outomega=$DATADIR/bovy/streamgap-pepper/gd1_multtime/gd1_t64sampling_X5_cutoff7p5_5-9_omega.dat -t 64sampling -M 5,9 -n 1002 -X 5. --cutoff=7.5

Pal 5
------

Do::

	python simulate_streampepper.py --outdens=$DATADIR/bovy/streamgap-pepper/pal5_multtime/pal5_t64sampling_X5_5-9_dens.dat --outomega=$DATADIR/bovy/streamgap-pepper/pal5_multtime/pal5_t64sampling_X5_5-9_omega.dat -t 64sampling -M 5,9 -n 1002 -X 5. -s pal5like --age=5. --amax=1.75 --da=0.01

