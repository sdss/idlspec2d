import os
import warnings
from datetime import datetime
from jinja2 import Template
import shutil

import boss_drp
from boss_drp.utils.splog import splog
from boss_drp.run.monitor_job import monitor_job

os.environ['BOSS_DPR_QUEUE_TYPE'] = 'Slurm'


# -------------------------------
# Warning for SLURM fallback
# -------------------------------
class SlurmWarning(Warning):
    def __init__(self, message):
        self.message = message
    def __str__(self):
        return repr(self.message)

# Try importing SLURM queue if needed
if os.getenv('BOSS_DPR_QUEUE_TYPE', 'SDSS_CHPC') == 'SDSS_CHPC':
    try:
        from slurm import queue
        os.environ['BOSS_DPR_QUEUE_TYPE'] = 'SDSS_CHPC'

    except Exception:
        warnings.warn(
            'No slurm package installed: printing command to STDOUT/logfile for manual run',
            SlurmWarning
        )
        os.environ['BOSS_DPR_QUEUE_TYPE'] = 'NoCluster'

# ---------------------------------
# Check if GNU Parallel is avaiable
# ---------------------------------
parallel_path = shutil.which("parallel")


# -------------------------------
# Base Queue class with registry
# -------------------------------
class Queue:
    registry = {}

    @classmethod
    def register(cls, name):
        """Decorator to register subclasses by name."""
        def decorator(subcls):
            cls.registry[name] = subcls
            return subcls
        return decorator

    def __new__(cls, config, *args, **kwargs):
        # Dispatch only if constructing base class
        if cls is Queue:
            key = os.getenv("BOSS_DPR_QUEUE_TYPE", "NoCluster")

            # If using Scheduler and type is None, fall back to NoCluster
            if key == "Scheduler" and config.get('scheduler') is None:
                splog.info('Queue is set to use Scheduler mode, but no scheduler defined in config... reverting to NoCluster Queue')
                key = "NoCluster"

            key_lower = key.lower()
            registry_lower = {k.lower(): v for k, v in cls.registry.items()}

            if key_lower not in registry_lower:
                raise RuntimeError(f"Unknown queue type {key}. Available: {list(cls.registry)}")
            impl = registry_lower[key_lower]

#            # --- LOG which subclass was selected ---
#            splog.info(f"Queue dispatch: Using implementation '{key}' -> {impl.__name__}")
            
            return super().__new__(impl)  # construct instance of the selected subclass
        else:
            return super().__new__(cls)

    def __init__(self, config, *args, **kwargs):
        """Base __init__, called by subclasses."""
        self.type = "base"
        self.config = config
        self._queue = None
        self._ntask = 0
        self.key = ""

    # -------------------------------
    # Methods to be overridden
    # -------------------------------
    def create(self, *args, **kwargs):
        pass

    def append(self, cmd, outfile=None, errfile=None):
        pass

    def commit(self, *args, submit=True, **kwargs):
        pass

    def monitor_job(self, monitor=True, *args, **kwargs):
        return False

# -------------------------------
# SLURM Implementation
# -------------------------------
@Queue.register("SDSS_CHPC")
class SDSS_CHPC(Queue):
    def __init__(self, config, *args, **kwargs):
        super().__init__(config, *args, **kwargs)
        self.type = "slurm"
        # Only initialize _queue if SLURM package is available
        if not self.config.get("no_write"):
            if "queue" in globals():
                self._queue = queue(*args, **kwargs)
            else:
                self._queue = None
        else:
            self._queue = None

    def create(self, *args, **kwargs):
        
        if self.config.get("no_write"):
            return

        if self.config.get("queue_sub_dir"):
            os.environ["SLURM_SCRATCH_DIR"] = self.config["queue_sub_dir"]

        if os.getenv("SLURM_SCRATCH_DIR") is None:
            os.environ["SLURM_SCRATCH_DIR"] = os.getcwd()

        if self._queue:
            self._queue.create(*args, **kwargs)
            self.key = self._queue.key

    def append(self, cmd, outfile=None, errfile=None):
        if self._queue:
            self._ntask += 1
            self._queue.append(cmd, outfile=outfile, errfile=errfile)

    def commit(self, *args, submit=None, **kwargs):
        if not self._queue:
            return
        if submit is None:
            submit = not self.config.get("no_submit")
        if submit is None:
            submit = True
        if self._ntask > 0:
            self._queue.commit(*args, submit=submit, **kwargs)

    def monitor_job(self, monitor=True, *args, **kwargs):
        if not monitor:
            return False
        if self._queue:
            status = monitor_job(self._queue, *args, **kwargs)
            return True if status is None else status
        return False

# -------------------------------
# No-Cluster Implementation
# -------------------------------
@Queue.register("NoCluster")
class NoCluster(Queue):
    def __init__(self, config, *args, **kwargs):
        super().__init__(config, *args, **kwargs)
        self.type = "none"
        self._queue_file = None
        self._queue = None

    def create(self, *args, label=None, **kwargs):
        if self.config.get("no_write"):
            return

        self._queue = []
        if self.config.get("queue_sub_dir"):
            os.environ["SLURM_SCRATCH_DIR"] = self.config["queue_sub_dir"]
        if os.getenv("SLURM_SCRATCH_DIR") is None:
            os.environ["SLURM_SCRATCH_DIR"] = os.getcwd()

        now = datetime.now()
        self._queue_file = os.path.join(
            os.getenv("SLURM_SCRATCH_DIR"),
            now.strftime("%d-%m-%Y_%H%M"),
            now.strftime("queue_%d-%m-%Y_%H%M.bash")
        )

        now = datetime.now()
        if label is None:
            label = now.strftime("%d-%m-%Y_%H%M")
        self._queue_file = os.path.join(os.getenv("SLURM_SCRATCH_DIR"), label, now.strftime("%d-%m-%Y_%H%M"),
                                        now.strftime(f"queue_%d-%m-%Y_%H%M.{self.type.lower()}"))


        os.makedirs(os.path.dirname(self._queue_file), exist_ok=True)

    def append(self, cmd, outfile=None, errfile=None):
        if self._queue is None:
            return
        self._ntask += 1
        self._queue.append(f"{cmd} > {outfile} 2> {errfile}")

    def commit(self, *args, submit=True, **kwargs):
        if self._queue is None:
            return
        if self._ntask > 0:
            with open(self._queue_file, "w") as f:
                for line in self._queue:
                    f.write(line + "\n")
            splog.info(f"Manual Queue written to {self._queue_file}")

@Queue.register("Scheduler")
class Scheduler(Queue):
    def __init__(self, config, *args, **kwargs):
        super().__init__(config, *args, **kwargs)
        self.type = config.get('scheduler')
        self._queue_file = None
        self._queue = None
        self._queue_config = {}
        if self.type.lower() in ['slurm','pbs']:
            self.template = os.path.join(os.path.dirname(boss_drp.__file__), 'etc','templates',f'{self.type}.j2')
        else:
            self.template = self.type
            self.type = 'Custom'

    def create(self, *args, label = None, **kwargs):
        if self.config.get("no_write"):
            return
        
        self._queue = []

        if self.config.get("queue_sub_dir"):
            os.environ["SLURM_SCRATCH_DIR"] = self.config["queue_sub_dir"]
        if os.getenv("SLURM_SCRATCH_DIR") is None:
            os.environ["SLURM_SCRATCH_DIR"] = os.getcwd()

        now = datetime.now()
        if label is None:
            label = now.strftime("%d-%m-%Y_%H%M")
        self._queue_file = os.path.join(os.getenv("SLURM_SCRATCH_DIR"), label, now.strftime("%d-%m-%Y_%H%M"),
                                        now.strftime(f"queue_%d-%m-%Y_%H%M.{self.type.lower()}"))
        
        os.makedirs(os.path.dirname(self._queue_file), exist_ok=True)

        self._queue_config = {'label':'_'.join(label.rstrip('/').split(os.sep)), **kwargs}

    def append(self, cmd, outfile=None, errfile=None):
        if self._queue is None:
            return
        self._ntask += 1
        self._queue.append(f"{cmd} > {outfile} 2> {errfile}")

    def commit(self, *args, cmd_file = None, submit=True, **kwargs):
        if self._queue is None:
            return

        pipe = None
        if self._ntask > 0:
            if not parallel_path:
                pipe = {'commands':self._queue, **self._queue_config }
            else:  
                cmd_file = self._queue_file.replace('queue_','run_').replace(f'.{self.type.lower()}','.bash')
                with open(cmd_file, "w") as f:
                    for line in self._queue:
                        f.write(line + "\n")

                pipe = {'commands_file':cmd_file, **self._queue_config }
        if cmd_file is not None:
            pipe = {'commands_file':cmd_file, **self._queue_config }
        if pipe:
            for key in self.config.keys():
                if key  == 'max':
                    continue
                if key not in pipe:
                    pipe[key] = self.config.get(key)
            for key in pipe.keys():
                if pipe[key] == 'None':
                    pipe[key] = None
                
            if pipe['ppn'] is None:
                try:
                    maxppn = int(self.config.get('max.ppn'))
                    pipe['ppn'] = min(maxppn,self._ntask)
                except:
                    pipe['ppn'] = self._ntask

            if pipe['cpus'] is None: 
                pipe['cpus'] = pipe['ppn']
            if pipe['cpus'] > pipe['ppn']:
                pipe['cpus'] = pipe['ppn']

            pipe['log'] = os.path.join(os.path.dirname(self._queue_file),pipe['label'])
            
            with open(self.template) as template:
                with open(self._queue_file, "w", encoding="utf-8") as f:
                    j2_template = Template(template.read())
                    f.write(j2_template.render(**pipe))
            splog.info(f"{self.type.upper()} Batch File written to {self._queue_file}")

@Queue.register("Slurm")
class Slurm(Scheduler):
    def __init__(self, config, *args, **kwargs):
        config.set('scheduler','Slurm')
        super().__init__(config, *args, **kwargs)
        self.type = "Slurm"
        self.template = os.path.join(os.path.dirname(boss_drp.__file__), 'etc','templates','slurm.j2')

@Queue.register("PBS")
class PBS(Scheduler):
    def __init__(self, config, *args, **kwargs):
        config.set('scheduler','PBS')
        super().__init__(config, *args, **kwargs)
        self.type = "PBS"
        self.template = os.path.join(os.path.dirname(boss_drp.__file__), 'etc','templates','pbs.j2')

# -------------------------------
# Usage
# -------------------------------
# Just construct the base class:
# The correct subclass will automatically be returned.
# config = {...}
# q = Queue(config)
# q.create()
# q.append("echo hi")
# q.commit()
