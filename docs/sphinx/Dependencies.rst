.. title:: IDLspec2D: Dependencies

Pipeline Dependencies
=====================

Setup Requirements
-------------------

A number of system environmental variables and paths are required to run the BOSS DRP.

Running at the University of Utah CHPC
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
At Utah, modules (`github.com/sdss/sdss_modules <https://github.com/sdss/sdss_modules>`_) are used to maintain the the required packages and paths. For the BOSS pipeline, the bhm modules are used to manage these for the BOSS pipeline.

Environmental Variables
^^^^^^^^^^^^^^^^^^^^^^^
::

    export DATABASE_PROFILE="pipelines"
    export RUN2D="|idlspec2d_version|"
    export RUN1D="|idlspec2d_version|"
    export BOSS_SPECTRO_REDUX="/uufs/chpc.utah.edu/common/home/sdss50/sdsswork/bhm/boss/spectro/redux"
    export BOSS_SPECTRO_DATA_N="/uufs/chpc.utah.edu/common/home/sdss50/sdsswork/bhm/boss/spectro/apo"
    export BOSS_SPECTRO_DATA_S="/uufs/chpc.utah.edu/common/home/sdss50/sdsswork/bhm/boss/spectro/lco"
    export SDHDRFIX_DIR="/uufs/chpc.utah.edu/common/home/sdss50/software/git/sdss/sdsscore_test/main"
    export IDLSPEC2D_DIR="/uufs/chpc.utah.edu/common/home/sdss50/software/git/sdss/idlspec2d/$RUN2D"
    export SDSSCORE_DIR="/uufs/chpc.utah.edu/common/home/sdss50/software/git/sdss/sdsscore_test/main"

.. _Paths:

Paths
^^^^^
::

    export IDL_PATH="$IDL_PATH:$IDLSPEC2D_DIR"
    export PATH="$PATH:$IDLSPEC2D_DIR/bin"
    export PYTHONPATH="$PYTHONPATH:$IDLSPEC2D_DIR/python"

Dependencies
^^^^^^^^^^^^

* idl
* python(3.7-3.11)
* SDSS Collaboration Package Dependencies
    * `idlutils <https://github.com/sdss/idlutils>`_: idlutils is a collection of IDL functions and routines used by a variety of SDSS software.
    * `sdssdb <https://github.com/sdss/sdssdb/>`_: sdssdb contains the source catalogs, targeting catalogs, and operational databases.
    * `semaphore <https://github.com/sdss/semaphore>`_: provides codes and reference files to decode and understand the sdss_target_flags
    * `sdss-tree <https://github.com/sdss/tree>`_: provides environment variables to manage paths for SDSS data as they organized on the Science Archive Server (SAS)
    * `sdss_access <https://github.com/sdss/sdss_access>`_: provides a convenient way of navigating local and remote file system paths from the Science Archive Server (SAS)
    * `pyVista <https://github.com/holtzmanjon/pyvista>`_:
    * `SDSS Slurm <https://github.com/sdss/slurm>`_: Required for use of uurundaily at Utah, all other command can be manually (at Utah or elsewhere) run without access to the slurm manager provided by this package
    * `sdsstools <https://github.com/sdss/sdsstools>`_:
* SDSS Product Dependencies
    * `elodie <https://svn.sdss.org/public/data/eboss/elodie/>`_: A database of high and medium-resolution stellar spectra (Prugniel+, 2001) used by spec1d to classify spectra and determine stellar parameters.
    * `dust <https://svn.sdss.org/public/data/sdss/catalogs/dust/>`_: A catalog of dust extinction models, including the SFD model.
    * `speclog <https://svn.sdss.org/public/data/sdss/speclog/trunk/>`_: speclog is an SDSS product that contains information about SDSS BOSS plate operations including seeing measured by the guides (guiderMon-{MJD}.par, plate plug maps (plPlugMapM-{plateid}-{mjd}-{plugid}.par, and plate header correction files to change the header exposure values (sdHdrFix-{mjd}.par)
    * `platelist <https://svn.sdss.org/public/data/sdss/platelist/trunk/>`_: platelist is an SDSS product that contains information on the plate designs and plugging. The plateHoles files include additional metadata associated with the targets on a plate
    * `specflat <https://svn.sdss.org/public/data/sdss/specflat/>`_: specflat is an SDSS product that contains master calibration frames and bad pixel masks for use in the idlspec2d pipeline.
    * `gaia/dr2 <https://cdn.gea.esac.esa.int/Gaia/gdr2/>`_: idlspec2d utilizes gaia_source/csv to calculate the distance to standard stars from GAIA DR2 proper motion.
* External Dependencies
    * `pyDL <https://pydl.readthedocs.io/en/latest/index.html>`_: a package that consists of python replacements for IDL function, both built-in and from external astronomical libraries
    * `dustmaps <https://github.com/gregreen/dustmaps>`_: provides a unified interface for several 2D and 3D maps of interstellar dust reddening and extinction. idlspec2d makes use of the Bayestar 2015 dustmaps (`Green, Schlafly, Finkbeiner et al. 2015 <https://ui.adsabs.harvard.edu/abs/2015ApJ...810...25G>`_)
    * `PyXCSAO <https://github.com/mkounkel/pyxcsao>`_: a python package designed to replicate the functionality of `IRAF XCSAO <http://tdc-www.harvard.edu/iraf/rvsao/xcsao/xcsao.html>`_.
    * `numpy <https://numpy.org/>`_: a standard Python package for arrays and high-level mathematical functions
    * `astropy <https://www.astropy.org/>`_: a collection of astronomy packages written in Python
    * `matplotlib <https://matplotlib.org/>`_: a python plotting library
    * `healpy <https://healpy.readthedocs.io/en/latest/>`_: a Python package based on the Hierarchical Equal Area isoLatitude Pixelization (HEALPix) scheme
    * `tqdm <https://tqdm.github.io/>`_: a progress bar for Python
    * `pandas <https://pandas.pydata.org/>`_: a python package designed for data manipulation and analysis
    * `h5py <https://www.h5py.org/>`_: a python interface between numpy and HDF5 data
    * `scipy <https://scipy.org/>`_: a python package for scientific and technical computing
    * `pillow <https://pillow.readthedocs.io/en/stable/index.html>`_: a python image file processing library
    * `termcolor <https://pypi.org/project/termcolor/>`_: a python package for color formatting of terminal outputs (not required but recommended)

