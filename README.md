# idlspec2d

The SDSS BOSS Spectrograph data reduction pipeline and associated tools.

## Should I be using the SVN version instead?

Until the end of SDSS-IV `idlspec2d` was maintained as a [SVN repository](https://svn.sdss.org/public/repo/eboss/idlspec2d/). SVN development happened as `v5_x` (and earlier), while Git development uses the `6.x` (and later) series. The SVN versions were not migrated to Github.

So, what version should you be using? If you're reducing SDSS-IV data, you probably want to keep using the `v5` branch. If you're reducing SDSS-V (or later) you *should* be using the GitHub repository.

## Using IDLspec2D

Starting with the v6_2_x tag series, the documentation for running the pipeline was added to [ReadtheDocs](https://sdss-idlspec2d.readthedocs.io/en/latest/). Prior to this, some (incomplete) documentation existed within the package itself. 

## Installing IDLspec2d

Starting with the v6_3_x tag series, the installation of the required Python dependencies has been simplified with the use of uv. However, some manual setup is still required due to the need of environmental variables and data. The [Dependencies](https://sdss-idlspec2d.readthedocs.io/en/latest/Dependencies.html) page on the [ReadtheDocs](https://sdss-idlspec2d.readthedocs.io/en/latest/Dependencies.html) describes the other requirements. 

Within the UV install, there are a number of optional extra installs:
 - dev: Install the (editable) development versions (from github) of the sdss packages **[Internal]**
 - dev_db: Install the (editable) development version (from github) sdssdb **[Internal]**
 - utah: Adds additional tools (such as sdssdb and juptyter used internally at Utah) **[Internal]**
 - sdss_slurm: Adds proprietary slurm tools used internally at Utah **[Internal]**
 - sos: Adds additional tools (such as sdssdb and ipython used internally at the mountains) **[Internal]**
 - sos_pin: Similar to sos, but also fixes versions of dependencies due to OS version of mountain **[Internal]**
 - docs: Adds the dependencies required to build the Sphinx documenation (either locally or ReadtheDocs)