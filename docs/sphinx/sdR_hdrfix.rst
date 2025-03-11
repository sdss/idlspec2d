
Fixing Raw Fits Headers & Specifying Alternative Calibrations
=============================================================


Fixing Raw Fits Headers
-----------------------
At times the raw fits headers have to be updated after the data is taken. As the raw frames are immediatly frozen in their initial state, this is is handled via alternative yanny sdHdrfix-<MJD>.par files. In the era of plate operations, these files were stored in the `speclog <https://svn.sdss.org/public/data/sdss/speclog/trunk/>`_ data product, while in FPS operations, these are stored in the `sdsscore <https://github.com/sdss/sdsscore/>`_ data product. The files are read in as part of the pipeline and can be used to correct individual values, or to flag exposures as bad or test, which are then excluded by the pipeline from proceessing.

sdR_hdrfix
^^^^^^^^^^
The :ref:`sdR_hdrfix<sdR_hdrfix>` command is used by the observers and pipeline team to create the sdHdrfix yanny files. After the script is run, the created (or appended) sdHdrfix-<MJD>.par files should be manually added the `sdsscore <https://github.com/sdss/sdsscore/>` repo. The can be done by navigating to $SDHDRFIX_DIR/<obs>/sdHdrfix and `git add .`.


sphdrfix.pro
^^^^^^^^^^^^
The `sphdrfix.pro` command is used internally by the pipeline and contains internal documentation if one should desire to run it manually, but it can be simply run with `sphdrfix, filename, hdr` within idl, supplying the filename and raw fits header, and it returns the modified header.


Sphdrfix.py
^^^^^^^^^^^
`Sphdrfix.py` contains the `Sphdrfix` Python class that is used by the Python portions of pipeline, where the class is initialized with the MJD (and observatory) of the interest (and other other optional parameters). The class instance function of `sphdrfix.fix(filename, hdr)`, supplying the filename and raw fits header, can be used to which returns the modified header.


Specifying Alternative Calibrations
-----------------------------------
In the plate era, there is a hard requirement of a flat and arc with each field, however in the FPS era this requirement was slightly relaxed. Flats are now only taken several times a night, while arcs are still nominally required for each field. However, on a case-by-case basis these FPS requirements can be relaxed. This is handled by the :ref:`flag_manual_cal<flag_manual_cal>` command, where either the associated arc requirement is turned off, or a particalar arc or flat can be flagged to be used for a field-mjd. These are stored as spManCal.par next to sdHdrfix-<MJD>.par files in in the `sdsscore <https://github.com/sdss/sdsscore/>`_ data product and are read in the field-mjd plan creation step.

