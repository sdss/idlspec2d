.. title:: IDLspec2D: The SDSS BOSS Data Reduction Pipeline

IDLspec2D: The SDSS BOSS Data Reduction Pipeline
================================================

The BOSS DRP (officially known as `idlspec2d <https://github.com/sdss/idlspec2d>`_) is really a series of steps run in sequence. There are some initial steps to produce the plans and supplementary information required. After which, uubatchpbs is used to produce the redux scripts. These scripts are really a wrapper script designed to run each of the individual steps on a field-mjd basis by the slurm manager at CHPC (Utah) or manually on any computer. Starting with the v6_2_x version of the pipeline, the internal python commands have been organized into a boss_drp package within the IDLspec2D GitHub repo.

The Pipeline can be run in various ways, where the catchup method is designed to run a large number of MJDs in a short time and the daily run method is designed to run for 1 (or a few) MJDs.


Major changes since version v6_1_X
----------------------------------

The idlspec2d package received a cleanup of old unused scripts and reorganized the internal python functions into an internal boss_drp python package. The scripts and files removed from the newer versions can be found in earlier tags on the `github repo <https://github.com/sdss/idlspec2d>` or `svn repo <https://svn.sdss.org/public/repo/sdss/idlspec2d>`
