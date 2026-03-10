.. title:: BOSS DRP Output File Tree 

BOSS DRP Output File Tree 
=========================
Historically, the BOSS DRP output files have been stored in a fairly flat directory tree (shown below) with all non-image files stored in a single RUN2D directory, and images stored in a seperate directory (images/RUN2D) next to the RUN2D directory. 
Within the main RUN2D directory, there was seperate directories for each plate/field, and within those directories there were directories for the RUN1D pipeline outputs, and spSpec intermediate files. Next to these plate/field directories,
there where all of the summary files and a spectra directory that contained all of the specLite and specFull outputs. In SDSS-V, this structure became a bit more complicated with the addtion of the different coadding schemes and their assocaited files. 

.. figure:: _static/tree/legacy_tree.png
    :alt: Legacy BOSS DRP Output File Tree
    :width: 90%
    :align: center

    Organization of the Pre-/Early SDSS-V BOSS DRP output directory structure.

However, this file structure started to become unwieldy with the increased data volume of SDSS-V, and the shere number of field/plates observed in the FPS era. In order to make the file structure more managable, 
the output files are now stored in a more hierarchical directory tree (shown below), with 5 top level folders (with in the RUN2D directory) for the different types of files (summary, spectra, images, fields, and traces). 
Within the fields directory there is a an additional layer of orgranzation, with FIELDGROUP directories that contain field/plate IDS with the same thousand digit. The custom coadds are stored in a seperate FIELDGROUP level directory 
(`allepoch` for SDSS-V) within the fields directory with each each observatory (and in theory allsky) having its own field level directory with in the `allepoch` directory. Within each field level directory, the schema remains roughly the same as before,
but with the addition of an epoch subdirectory for the field epoch coadds, which has its own subdirectories for the spSpec and RUN1D outputs. The summary files are stored in the summary directory, with seperate directories associated with each 
coadding schema (daily, epoch, allepoch). The spectra directory is organized by coadding schema, with the old spectra file tree within each (Full/Lite). However, as with the field folder, there is th additional layer for FIELDGROUPs. The image directory
was also organized by Coadding schema and FIELDGROUPs, as well as being moved within the RUN2D directory next to the other products. The latest versions of `sdss_tree <https://github.com/sdss/tree>`_ and `sdss_access <https://github.com/sdss/sdss_access>`_ have be 
updated to dynamically build the tree paths dependent on the RUN2D of the BOSS DRP used to produce the files.


.. figure:: _static/tree/tree.png
    :alt: BOSS DRP Output File Tree
    :width: 800px
    :align: center

    Updated Organization of the SDSS-V BOSS DRP output directory structure. (`Simplified version of the updated structure <tree_simple.html>`_).

