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
The required environmental variables for the `idlspec2d <https://github.com/sdss/idlspec2d>`_ package and data products are given below. This does not include the paths required for `idlutils <https://github.com/sdss/idlutils>`_, as these are described in its install instructions on `https://www.sdss.org/dr18/software/packages/idlutils <https://www.sdss.org/dr18/software/packages/idlutils>`_. The paths are supplied as they currently exist at the University of Utah CHPC (where the pipeline is routinely run), but can be used as a template for your enviroment.

.. code-block:: shell

    export DATABASE_PROFILE="pipelines"
    export RUN2D="|idlspec2d_version|"
    export RUN1D="|idlspec2d_version|"
    export BOSS_SPECTRO_REDUX="/uufs/chpc.utah.edu/common/home/sdss50/sdsswork/bhm/boss/spectro/redux"
    export BOSS_SPECTRO_DATA_N="/uufs/chpc.utah.edu/common/home/sdss50/sdsswork/bhm/boss/spectro/apo"
    export BOSS_SPECTRO_DATA_S="/uufs/chpc.utah.edu/common/home/sdss50/sdsswork/bhm/boss/spectro/lco"
    export SDHDRFIX_DIR="/uufs/chpc.utah.edu/common/home/sdss50/software/git/sdss/sdsscore_test/main"
    export IDLSPEC2D_DIR="/uufs/chpc.utah.edu/common/home/sdss50/software/git/sdss/idlspec2d/$RUN2D"
    export SDSSCORE_DIR="/uufs/chpc.utah.edu/common/home/sdss50/software/git/sdss/sdsscore_test/main"
    export BOSS_DRP_DAILY_DIR="$HOME/daily/"
    export BOSS_QA_DIR=$BOSS_SPECTRO_REDUX/
    export BOSS_DRP_EMAIL_DOMAIN='chpc.utah.edu'

    export ELODIE_DIR="/uufs/chpc.utah.edu/common/home/sdss09/software/eboss/NULL/elodie/v1_3"
    export SPECLOG_DIR="/uufs/chpc.utah.edu/common/home/sdss09/software/svn.sdss.org/data/sdss/speclog/trunk"
    export PLATELIST_DIR="/uufs/chpc.utah.edu/common/home/sdss09/software/svn.sdss.org/data/sdss/platelist/trunk"
    export SPECFLAT_DIR="/uufs/chpc.utah.edu/common/home/sdss50/software/git/sdss/specflat/master"

.. _Paths:
Paths
^^^^^
.. code-block:: shell

    export IDL_PATH="$IDL_PATH:$IDLSPEC2D_DIR"
    export PATH="$PATH:$IDLSPEC2D_DIR/bin"
    export PYTHONPATH="$PYTHONPATH:$IDLSPEC2D_DIR/python"

Compiling Internal Functions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^
The idlspec2d package contains a shared object library and C scripts. If you have issues with these, you can run ``evilmake all`` within the main idlspec2d directory ($IDLSPEC2D_DIR). ``evilmake`` is included within idlutils, so follow the directions on (`sdss.org/dr18/software/packages/idlutils/ <https://www.sdss.org/dr18/software/packages/idlutils/>`_) to install idlutils first.

Compiling Internal Documentation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
This documentation can be compiled locally by navigating to ``$IDLSPEC2D_DIR/docs/sphinx`` and run ``make all`` (This can also be run with ``evilmake all`` in the same directory or ``evilmake doc`` within ``$IDLSPEC2D_DIR``). This will compile the html, singlehtml, and pdf forms of the documentation within the associated subdirectoires of ``$IDLSPEC2D_DIR/docs/sphinx/_build``.


Downloading File Dependencies
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Some of the dependencies listed below contain files rather then code and can be obtained by the following commands (assuming the Environmental Variables are set as described above):

.. code-block:: shell
    
    cd "$(dirname "$ELODIE_DIR")" && svn co https://svn.sdss.org/public/data/eboss/elodie/v1_3 v1_3
    cd "$(dirname "$SPECLOG_DIR")" && svn co https://svn.sdss.org/public/data/sdss/speclog/trunk trunk
    cd "$(dirname "$PLATELIST_DIR")" && svn co https://svn.sdss.org/public/data/sdss/platelist/trunk trunk
    cd "$(dirname "$SPECFLAT_DIR")" && git clone https://github.com/sdss/specflat.git
    
    if [ -n "$SDSSCORE_DIR" ]; then
        cd "$(dirname "$SDSSCORE_DIR")" && git clone git@github.com:sdsscore.git
        cd sdsscore
        git submodule init
        git submodule update --recursive
    fi


.. note::
    svn can be obtained by following the directions on `https://subversion.apache.org/packages.html <https://subversion.apache.org/packages.html>`_

Raw, intermediate, or Reduced Products
""""""""""""""""""""""""""""""""""""""
The raw, intermediate, and reduced products can be downloaded from the SAS to use in the pipeline. SDSS supplied a number of ways to do this depending on volumes. Bulk downloading is described on `https://www.sdss.org/dr18/data_access/bulk/ <https://www.sdss.org/dr18/data_access/bulk>`_, while smaller downloads can be done with `sdss_access <https://sdss-access.readthedocs.io/en/latest/intro.html`_

Dustmaps for readfibermaps
""""""""""""""""""""""""""
If you plan on running :ref:`readfibermaps<readfibermaps>` yourself and not download the previously run files from the SAS (Note: downloading them is the suggested path since it requires access to internal databases to completely match what is produced by the pipeline), pre-caching the dust maps is highly recommended. If you would like to manually control the location of these cache, you can follow the directions in the `dustmaps <https://dustmaps.readthedocs.io/en/latest/installation.html#custom-configuration-file-location-optional>`_ documentation. Running the following will cache the maps used by the pipeline.

.. code-block:: python
    
    import dustmaps.sfd
    dustmaps.sfd.fetch()

    import dustmaps.bayestar
    dustmaps.bayestar.fetch(version='bayestar2015')
    
    import dustmaps.edenhofer2023
    dustmaps.edenhofer2023.fetch(fetch_2kpc=True)
    dustmaps.edenhofer2023.fetch(fetch_2kpc=False)

The default dust map used by the pipeline is a merge of SFD, Bayestar2015, and simple_dust_2023. Simple_dust_2023 is an unpublished proprietary dust map. While `idlspec2d <https://github.com/sdss/idlspec2d>`_ contains the code to use it, it will default back to edenhofer2023 in its place (though edenhofer2023 is a much heavier dustmap).

Dependencies
-------------------

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
    * `sdsscore <https://github.com/sdss/sdsscore/>`_: sdsscore is an SDSS-V product that contains the FPS fiber configuration files and header correction files to change the header exposure values (sdHdrFix-{mjd}.par)
    * `elodie <https://svn.sdss.org/public/data/eboss/elodie/>`_: A database of high and medium-resolution stellar spectra (Prugniel+, 2001) used by spec1d to classify spectra and determine stellar parameters.
    * `speclog <https://svn.sdss.org/public/data/sdss/speclog/trunk/>`_: speclog is an SDSS product that contains information about SDSS BOSS plate operations including seeing measured by the guides (guiderMon-{MJD}.par, plate plug maps (plPlugMapM-{plateid}-{mjd}-{plugid}.par, and plate header correction files to change the header exposure values (sdHdrFix-{mjd}.par)
    * `platelist <https://svn.sdss.org/public/data/sdss/platelist/trunk/>`_: platelist is an SDSS product that contains information on the plate designs and plugging. The plateHoles files include additional metadata associated with the targets on a plate
    * `specflat <https://svn.sdss.org/public/data/sdss/specflat/>`_: specflat is an SDSS product that contains master calibration frames and bad pixel masks for use in the idlspec2d pipeline.
* Deprecated SDSS Product Dependencies
    * `dust <https://svn.sdss.org/public/data/sdss/catalogs/dust/>`_: A catalog of dust extinction models, including the SFD model (deprecated).
    * `gaia/dr2 <https://cdn.gea.esac.esa.int/Gaia/gdr2/>`_: idlspec2d utilizes gaia_source/csv to calculate the distance to standard stars from GAIA DR2 proper motion.
* External Dependencies
    * `pyDL <https://pydl.readthedocs.io/en/latest/index.html>`_: a package that consists of python replacements for IDL function, both built-in and from external astronomical libraries
    * `dustmaps <https://github.com/gregreen/dustmaps>`_: provides a unified interface for several 2D and 3D maps of interstellar dust reddening and extinction. idlspec2d makes use of the Bayestar 2015 dustmaps (`Green, Schlafly, Finkbeiner et al. 2015 <https://ui.adsabs.harvard.edu/abs/2015ApJ...810...25G>`_)
    * `PyXCSAO <https://github.com/mkounkel/pyxcsao>`_: a python package designed to replicate the functionality of `IRAF XCSAO <http://tdc-www.harvard.edu/iraf/rvsao/xcsao/xcsao.html>`_.
    * `numpy <https://numpy.org/>`_: a standard Python package for arrays and high-level mathematical functions
    * `astropy(<7.0) <https://www.astropy.org/>`_: a collection of astronomy packages written in Python
    * `matplotlib <https://matplotlib.org/>`_: a python plotting library
    * `healpy <https://healpy.readthedocs.io/en/latest/>`_: a Python package based on the Hierarchical Equal Area isoLatitude Pixelization (HEALPix) scheme
    * `tqdm <https://tqdm.github.io/>`_: a progress bar for Python
    * `pandas <https://pandas.pydata.org/>`_: a python package designed for data manipulation and analysis
    * `h5py <https://www.h5py.org/>`_: a python interface between numpy and HDF5 data
    * `scipy <https://scipy.org/>`_: a python package for scientific and technical computing
    * `pillow <https://pillow.readthedocs.io/en/stable/index.html>`_: a python image file processing library
    * `jinja2 <https://jinja.palletsprojects.com/en/3.1.x/>`_: a python templating engine used to build the HTMLs produced by the pipeline
    * `termcolor <https://pypi.org/project/termcolor/>`_: a python package for color formatting of terminal outputs (not required but recommended)
    * `plotly <https://plotly.com/python/>`_: a python package for interactive plots (not required for the core pipeline but used by some of the supplementary tools)
    * `psutil <https://psutil.readthedocs.io/en/latest/>`_: A cross-platform libary for system monitoring and process running via Python (only used by SOS)
    * `GitPython <https://gitpython.readthedocs.io/en/stable/>`_: A python libary used for git iteractions (only used by SOS - option)
