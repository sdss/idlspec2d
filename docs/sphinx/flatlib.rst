Building the Flatlib QA Library
===============================

The Flatlib QA libary is a tool to monitor fiber throughput, as measured by the trace flats as a function of different parameters. These parameters include rotator position, telescope altitiude, FPS beta arm angle, and time. The time series analysis is especially important for monitoring low throughput or broken fibers to me replaced or repaired in upcoming engineering runs.

Building the Libary
"""""""""""""""""""
The flatlib analysis commands are all included as sub commands in the :ref:`flatlib<flatlib>` script, with the step by step process outlined below, :

::

    flatlib reduce --fps --link_all --link_traceflat --no_run
    flatlib build
    flatlib analyze
    
These steps assume that you have viable spFlat and spTraceFlat files produced by the normal pipeline. If you drop the `--no_run` flag, then it use the cluster to build the missing spFlat files. If you want to do a detailed analysis of any of the parameters, this is advisable, but for monitoring of fiber throughput, this is unnecessary.


.. note::
    if you want to build for lco add `-l` flag to every command

Quick Run for TimeSeries Analysis can be run via
::

    flatlib reduce --fps --link_all --link_traceflat --no_run
    flatlib build
    flatlib csv
    flatlib timeSeries


