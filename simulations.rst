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

Length factor
--------------

Commands like::

	 python simulate_streampepper.py --outdens=$DATADIR/bovy/streamgap-pepper/gd1_multtime/gd1_t64sampling_lf1p25_6.5_dens.dat --outomega=$DATADIR/bovy/streamgap-pepper/gd1_multtime/gd1_t64sampling_lf1p25_6.5_omega.dat -t 64sampling -M 6.5 --dt=10. -X 5. -l 1.25 --timescdm=1.25

Full mass range
---------------

Do::

	python simulate_streampepper.py --outdens=$DATADIR/bovy/streamgap-pepper/gd1_multtime/gd1_t64sampling_X5_5-9_dens.dat --outomega=$DATADIR/bovy/streamgap-pepper/gd1_multtime/gd1_t64sampling_X5_5-9_omega.dat -t 64sampling -M 5,9 --dt=350. -X 5.

Plummer and different rs for Hernquist
---------------------------------------

Commands like::

	python simulate_streampepper.py --outdens=$DATADIR/bovy/streamgap-pepper/gd1_multtime/gd1_t64sampling_X3p24_plum_5-9_dens.dat --outomega=$DATADIR/bovy/streamgap-pepper/gd1_multtime/gd1_t64sampling_X3p24_plum_5-9_omega.dat -t 64sampling -M 5,9 --dt=1000. -X 3.24 --plummer
	python simulate_streampepper.py --outdens=$DATADIR/bovy/streamgap-pepper/gd1_multtime/gd1_t64sampling_X2_rsfac2p5_5-9_dens.dat --outomega=$DATADIR/bovy/streamgap-pepper/gd1_multtime/gd1_t64sampling_X2_rsfac2p5_5-9_omega.dat -t 64sampling -M 5,9 --dt=1000. -X 2. --rsfac=2.5
	python simulate_streampepper.py --outdens=$DATADIR/bovy/streamgap-pepper/gd1_multtime/gd1_t64sampling_X12p5_rsfacp4_5-9_dens.dat --outomega=$DATADIR/bovy/streamgap-pepper/gd1_multtime/gd1_t64sampling_X12p5_rsfacp4_5-9_omega.dat -t 64sampling -M 5,9 --dt=1000. -X 12.5 --rsfac=0.4

Higher or lower rate than CDM
------------------------------

Commands like::

	 python simulate_streampepper.py --outdens=$DATADIR/bovy/streamgap-pepper/gd1_multtime/gd1_t64sampling_X5_cdm3_5-9_dens.dat --outomega=$DATADIR/bovy/streamgap-pepper/gd1_multtime/gd1_t64sampling_X5_cdm3_5-9_omega.dat -t 64sampling -M 5,9 --dt=350. -X 5. --timescdm=3.
	 python simulate_streampepper.py --outdens=$DATADIR/bovy/streamgap-pepper/gd1_multtime/gd1_t64sampling_X5_cdmp33_5-9_dens.dat --outomega=$DATADIR/bovy/streamgap-pepper/gd1_multtime/gd1_t64sampling_X5_cdmp33_5-9_omega.dat -t 64sampling -M 5,9 --dt=350. -X 5. --timescdm=0.33333333333
	 python simulate_streampepper.py --outdens=$DATADIR/bovy/streamgap-pepper/gd1_multtime/gd1_t64sampling_X5_cdm10_5-9_dens.dat --outomega=$DATADIR/bovy/streamgap-pepper/gd1_multtime/gd1_t64sampling_X5_cdm10_5-9_omega.dat -t 64sampling -M 5,9 --dt=1000. -X 5. --timescdm=10.
	 python simulate_streampepper.py --outdens=$DATADIR/bovy/streamgap-pepper/gd1_multtime/gd1_t64sampling_X5_cdmp1_5-9_dens.dat --outomega=$DATADIR/bovy/streamgap-pepper/gd1_multtime/gd1_t64sampling_X5_cdmp1_5-9_omega.dat -t 64sampling -M 5,9 --dt=1000. -X 5. --timescdm=0.1

Different age
--------------

Commands like::

	 python simulate_streampepper.py --outdens=$DATADIR/bovy/streamgap-pepper/gd1_multtime/gd1_t64sampling_X5_age4p5_5-9_dens.dat --outomega=$DATADIR/bovy/streamgap-pepper/gd1_multtime/gd1_t64sampling_X5_age4p5_5-9_omega.dat -t 64sampling -M 5,9 --dt=900. -X 5. --age=4.5
	 python simulate_streampepper.py --outdens=$DATADIR/bovy/streamgap-pepper/gd1_multtime/gd1_t64sampling_X5_age4p5_cdm2_5-9_dens.dat --outomega=$DATADIR/bovy/streamgap-pepper/gd1_multtime/gd1_t64sampling_X5_age4p5_cdm2_5-9_omega.dat -t 64sampling -M 5,9 --dt=900. -X 5. --age=4.5 --timescdm=2.

CDM spectrum w/ cut-off
-----------------------

Commands like::

	 python simulate_streampepper.py --outdens=$DATADIR/bovy/streamgap-pepper/gd1_multtime/gd1_t64sampling_X5_cutoff5p5_5-9_dens.dat --outomega=$DATADIR/bovy/streamgap-pepper/gd1_multtime/gd1_t64sampling_X5_cutoff5p5_5-9_omega.dat -t 64sampling -M 5,9 --dt=800. -X 5. --cutoff=5.5
	 python simulate_streampepper.py --outdens=$DATADIR/bovy/streamgap-pepper/gd1_multtime/gd1_t64sampling_X5_cutoff6p5_5-9_dens.dat --outomega=$DATADIR/bovy/streamgap-pepper/gd1_multtime/gd1_t64sampling_X5_cutoff6p5_5-9_omega.dat -t 64sampling -M 5,9 --dt=800. -X 5. --cutoff=6.5
	 python simulate_streampepper.py --outdens=$DATADIR/bovy/streamgap-pepper/gd1_multtime/gd1_t64sampling_X5_cutoff7p5_5-9_dens.dat --outomega=$DATADIR/bovy/streamgap-pepper/gd1_multtime/gd1_t64sampling_X5_cutoff7p5_5-9_omega.dat -t 64sampling -M 5,9 --dt=800. -X 5. --cutoff=7.5

CDM spectrum w/ different exponent
-----------------------------------

Commands like::

	 python simulate_streampepper.py --outdens=$DATADIR/bovy/streamgap-pepper/gd1_multtime/gd1_t64sampling_X5_massexpm1p5_5-9_dens.dat --outomega=$DATADIR/bovy/streamgap-pepper/gd1_multtime/gd1_t64sampling_X5_massexpm1p5_5-9_omega.dat -t 64sampling -M 5,9 --dt=800. -X 5. --massexp=-1.5
	 python simulate_streampepper.py --outdens=$DATADIR/bovy/streamgap-pepper/gd1_multtime/gd1_t64sampling_X5_massexpm2p5_5-9_dens.dat --outomega=$DATADIR/bovy/streamgap-pepper/gd1_multtime/gd1_t64sampling_X5_massexpm2p5_5-9_omega.dat -t 64sampling -M 5,9 --dt=800. -X 5. --massexp=-2.5

Pal 5
------

Do::

	python simulate_streampepper.py --outdens=$DATADIR/bovy/streamgap-pepper/pal5_multtime/pal5_t64sampling_X5_5-9_dens.dat --outomega=$DATADIR/bovy/streamgap-pepper/pal5_multtime/pal5_t64sampling_X5_5-9_omega.dat -t 64sampling -M 5,9 --dt=350. -X 5. -s pal5like --age=5.

