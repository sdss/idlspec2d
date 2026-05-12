Using Configuration Files
=========================

Starting with v6_2_2, the use of configuration files have been added to the BOSS DRP. Two different configuration files are used at all times:

- queue.yml
- boss_drp.yml

The queue.yml file contains the necessary information for submitting the pipeline jobs to a slurm or pbs queue. 
The boss_drp.yml contains the paramters and flags for the various steps of the pipeline.

The default versions of these are stored in the idlspec2d repository in boss_drp/etc. To use these, the user can copy them to their own 
~/.config/sdss/boss_drp/boss_drp.yml, or set the environment variable BOSS_DRP_CONFIG_PATH to point to the location of the config file. 
**BOSS_DRP_PIPE_CONFIG_PATH** and **BOSS_DRP_QUEUE_CONFIG_PATH** can also be set to a directory, in which case the pipeline will look for boss_drp.yml and queue.yml in that directory. 
Additionally, these config files can be specified on the command line using the following flags:

- `--pipe_config` - The name of the pipeline config file (e.g., 'boss_drp') - this overrides the default 'boss_drp/boss_drp.yml' file in your user directory. This is used to determine the parameters for the various steps of the pipeline.
- `--queue_config` - Sets the config name (e.g., 'tagged_daily') within the queue.yml file, which is used to determine the parameters for submitting jobs to the queue. 
- `--pipe_config_file` - This is the path to a specific pipeline config file, which will override any **BOSS_DRP_PIPE_CONFIG_PATH** settings.
- `--queue_config_file` - This is the path to a specific queue config file, which will override any **BOSS_DRP_QUEUE_CONFIG_PATH** settings.

If pipe_config_file or queue_config_file are specified, the associated environment variables **BOSS_DRP_PIPE_CONFIG_PATH** or **BOSS_DRP_QUEUE_CONFIG_PATH** will be cleared 
to avoid confusion, but if they are not specified, then the pipeline will look for config files in the directories specified by those environment variables 
(or the default ~/.config/sdss/boss_drp/ if the environment variables are not set).


Viewing and Editing
-------------------

The `boss_drp config` subcommand provides an interface to view and edit the configuration files, however it can also be done manually. 