Building the SpecFlat Product
=============================

The specflat product is not direcly part of the `idlspec2d <https://github.com/sdss/idlspec2d>`_) product, it is a dependency of it. It contains two different types of products:

Monthly Calibration Products
""""""""""""""""""""""""""""
The Monthly Calibration Products can be built with the monthly calibration sets of bias, darks, and flats taken by the observers. The set of required calibrations differ slightly between the observatories, but in general, it contails ~25 bias frames, ~3 x 900s dark frames, and 5 flat frames (taken with each flat field screen available). The Flats are primiarly used to help build trace flats and the flatlib. The bias and darks are used to build pixel bias frames, and bad pixel masks. These can all be build using the `MontlyCheck.py` command in the specflat product, with the exposure ranges of the darks and flats given as below:

::
    MonthlyCheck.py --mjd 60584 --lco --biasrange 00040934 00040958 --darkrange 00040959 00040961
    MonthlyCheck.py --mjd 60584 --apo --biasrange 00040934 00040958 --darkrange 00040959 00040961

These can be compaired visually and graphically with the previous sets using the `Analysis.py` command:

::
    Analysis.py --bias --tagged_loaded
    Analysis.py --bias --lco --tagged_loaded
    Analysis.py --bpm --tagged_loaded
    Analysis.py --bpm --lco --tagged_loaded



Fiber Flat Products
"""""""""""""""""""
The Fiber Flat Products can be build with leaky Lossy fiber pixel flats taken with the custom lossy fiber light source that mounts in place of the BOSS Fiber slithead. The exact procedures differ by telescope, but in general, 2 different sets are taken, optimizing the signal for the 2 cameras. A set of shorter frames are taken, optimized for the red, followed by a longer set optimized for the blue (and staturating the red).  These can all be build using the `useLossyFlats.py` command in the specflat product, with the exposure ranges of the pixel flats and an unsturated exposure for measuring gain as given below:


::
    useLossyFlats.py --mjd 59768 --expid1 344228   --expid2 344279   --APO --gain_expid 344228 --blue_exptime 500  --red_exptime 150
    useLossyFlats.py --mjd 60516 --expid1 00037103 --expid2 00037187 --LCO --gain_expid 00037103 --blue_exptime 300 --red_exptime 45


These can be compaired visually and graphically with the previous sets using the `Analysis.py` command:

::
    Analysis.py --gain
    Analysis.py --pixflat --tagged_loaded
    Analysis.py --pixflat --lco --tagged_loaded


