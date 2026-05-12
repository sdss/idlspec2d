import boss_drp
from boss_drp.utils.splog import splog
from sklearn.preprocessing import normalize
from sdsstools.configuration import DEFAULT_PATHS, get_config
import copy
import os
from datetime import datetime, timedelta
import warnings
import ast
import operator
from ruamel.yaml import YAML
import json
from io import StringIO

class QueueConfigError(Exception):
    """This is a custom exception."""
    pass

class QueueConfigWarning(Warning):
    pass

DEFAULT_PATHS = [
    "~/.config/sdss/{name}",
    "~/.config/sdss/{name}/{name}",
    "~/.{name}/{name}",
    f"{os.getenv('BOSS_DRP_DAILY_DIR')}/config"
]


NOT_BOOL = {'general':'*', 'fmjdselect':'*', 'customSettings':['custom_name'], 'custom_name':['schema_file'], 
            'plan':['type'], 'plan.daily':['minscience','dailyplan_logfile'],
            'plan.epoch':['min_epoch_length', 'epochplan_logfile', 'max_epoch_length'],
            'plan.custom':['cartons', 'program', 'catalogids', 'customplan_logfile'],
            'plan.trace':['traceplan_logfile'], 'fibermap':['datamodel'], 'reduce':['map3d', 'docams', 'nitersky'],
            'combine':['minsn2', 'bscore'], 'analyze':'*', 
            'post.fieldmerge':['datamodel_file', 'line_datamodel_file', 'batchwise.ndays'],
            'post':['spCalib', 'healpix'], 'monitor':['pause'], 
            'Summary.batchwise':['ndays', 'backup','n_iter','mjdstart','mjdend','MJD_dir']}

BOOL = {'general':['REMOTE'], 'fmjdselect':['epoch', 'custom', 'dither', 'commissioning','batch_mjd','trace_all_mjds']}

queue_bool = ['no_write','no_submit','exclusive']
def fill_none_with_false(cfg, path=""):
    """
    Recursively replace None with False according to NOT_BOOL / BOOL rules.
    """
    for key, value in cfg.items():
        full_path = f"{path}.{key}" if path else key

        # Recurse into nested dicts
        if isinstance(value, dict):
            fill_none_with_false(value, full_path)
            continue

        # Only act on None values
        if value is not None:
            continue

        if path in queue_bool:
            cfg[path] = False

        rule = NOT_BOOL.get(path, None)
        if rule is None:
            continue

        # Sections marked with '*': only keys in BOOL are set to False
        if rule == '*':
            if key in BOOL.get(path, []):
                cfg[key] = False
            # otherwise leave as None

        # Sections with explicit "do not bool" lists:
        # listed keys stay None, everything else becomes False
        elif isinstance(rule, list):
            if key not in rule:
                cfg[key] = False

        # Sections not mentioned in NOT_BOOL:
        # all None values become False
        else:
            cfg[key] = False

    return cfg
                    
            


def to_td(s):
    h, m, sec = map(int, s.split(':'))
    return timedelta(hours=h, minutes=m, seconds=sec)

class QueueConfig:
    def __init__(self, obj, queue_config, queue_config_file):
        self._obj = obj
        self._queue_config = queue_config
        self._queue_config_file = queue_config_file


    def __getattr__(self, item):
        # pass through all other attribute lookups
        return getattr(self._obj, item)

    def __getitem__(self, key):
        return self._obj[key]

    def __setitem__(self, key, value):
        self._obj[key] = value

    def set(self, name, value, force=False):
        if (value is None) & (not force): 
            return
        elif not force:
            maxval = self._obj.get(f'max.{name}')
        else:
            maxval = None
        if maxval is not None:
            if name in ['wall','walltime']:
                _val = to_td(value)
                _maxval = to_td(maxval)
            else:
                _val = value
                _maxval = int(maxval)
            if _val > _maxval:
                raise QueueConfigError(f"Invalid {name}: {value} > max value ({maxval})")
        self._obj[f'{name}'] = value

    def to_dict(self, label = None):
        bundle = True if self._obj.get('nbundle') is not None else False
        _dict = dict(nodes=str(self._obj.get('nodes')), ppn=str(self._obj.get('ppn')),
                    partition = self._obj.get('partition'), alloc=self._obj.get('alloc'),
                    walltime = self._obj.get('wall'), mem_per_cpu = self._obj.get('mem_per_cpu'),
                    mem = self._obj.get('mem'), nbundle = self._obj.get('nbundle'), bundle=bundle,
                    qos=self._obj.get('qos'), exclusive=self._obj.get('exclusive'),
                    constraint=self._obj.get('constraint'), gres=self._obj.get('gres'))
        for key,val in _dict.items():
            if val == 'None':
                _dict[key] = None
        if label is not None:
            _dict = {'label': label, **_dict}
        return _dict

    def to_str(self):
        queue_str = []
        for key, value in self._obj.items():
            if key in ['max','fallback']:
                continue
            queue_str.append(f"{key}: {value}")
        queue_str = ' \n'.join(queue_str)       
        return queue_str 

    def __deepcopy__(self, memo):
        # deepcopy the underlying object
        copied_obj = copy.deepcopy(self._obj, memo)
        # return a new QueueConfig wrapping the copied object
        return QueueConfig(copied_obj, 
                       self._queue_config, 
                       self._queue_config_file)
    def _fallback(self):
        if self._obj.get('fallback') is None:
            return
        _queue = get_config('queue', allow_user=True, config_file=self._queue_config_file,
                                    config_envvar='BOSS_DRP_QUEUE_CONFIG_PATH')
        avail_config = ','.join(_queue.keys())
        _queue = _queue.get(self._queue_config)
        updates = {}
        for key, var in _queue.items():
            if var != self._obj.get(key):
                updates[key] = self._obj.get(key)

        _queue = get_config('queue', allow_user=True, config_file=self._queue_config_file,
                                    config_envvar='BOSS_DRP_QUEUE_CONFIG_PATH')
        avail_config = ','.join(_queue.keys())
        _queue = _queue.get(self._obj.get('fallback.config'))   

        if _queue is None:
            raise QueueConfigError(f"Missing queue Configuration for {self._obj.get('fallback')} (available: {avail_config})")
        
        self._obj = _queue

        splog.info(updates)

        for key, var in updates.items():
            self.set(key,var)   

        splog.info(self.to_str())    

    def check_fallback(self):
        if self._obj.get('fallback') is None:
            return
        if (self._obj.get('fallback.check') is None) or (self._obj.get('fallback.config') is None):     
            return
        expr = self._obj.get('fallback.check').format(**self.to_dict(''))
        if safe_eval(expr):
            warnings.warn(f"Fallback Config conditions met... Falling back to {self._obj.get('fallback.config')}",
                          QueueConfigWarning)
            self._fallback()
           
class Config:
    def __init__(self):
        self.queue = None
        self._queue = None
        self.pipe = None
        self._queue_config_name = None

        
    # def load(self, config, queue_config, config_file=None, queue_config_file=None):
    def load(self, queue_config=None,  queue_config_file=None, config_name = 'boss_drp', config_file=None,):
        # The pipeline steps will not use this config file, 
        # but could have their own in future if idl is ported to python
        ptest = splog._log.propagate
        splog._log.propagate = False
        try:
            _pipe = get_config(config_name,allow_user=True, config_file=config_file,
                                     config_envvar='BOSS_DRP_PIPE_CONFIG_PATH')
        except:
            _pipe = get_config(config_name,allow_user=True, config_file=config_file,
                                     config_envvar=None)
        self.pipe = _pipe
        splog.info(f'Loaded Pipeline Config from {_pipe._CONFIG_FILE}')

        if queue_config_file is not None:
            os.environ.pop('BOSS_DRP_QUEUE_CONFIG_PATH', None)
        if config_file is not None:
            os.environ.pop('BOSS_DRP_PIPE_CONFIG_PATH', None)
        if queue_config is not None:
            try:
                self._queue = get_config('queue', allow_user=True, config_file=queue_config_file,
                                          config_envvar='BOSS_DRP_QUEUE_CONFIG_PATH')
            except:
                self._queue = get_config('queue', allow_user=True, config_file=queue_config_file,
                                          config_envvar=None)
            avail_config = ','.join(self._queue.keys())
            _queue = self._queue.get(queue_config)
            self._queue_config_name = queue_config

            if _queue is None:
                raise QueueConfigError(f"Missing queue Configuration for {queue_config} (available: {avail_config})")
            
            splog.info(f'Loaded Cluster Queue Config {queue_config} from {self._queue._CONFIG_FILE}')

            self.queue = QueueConfig(_queue, queue_config, queue_config_file)
        splog._log.propagate = ptest

    def __repr__(self):
        return self.__str__()
    
    def __str__(self):
        queue_str = self.queue.to_str() if self.queue is not None else ''
        if self.pipe is not None:
            cfg_str  = "\n ".join(f"{k}: {v}" for k, v in self.pipe['general'].items())
        else:
            cfg_str = ''

        return (cfg_str+
                f"\n{queue_str}\n")

    def full_str(self, stages=None, classes= ['Pipe', 'Queue']):
        parts = []
        if 'Queue' in classes:
            cfgs = {'Pipe':self.pipe, 'Queue': self.queue._obj}
        else:
            cfgs = {'Pipe':self.pipe, 'Queue': None}


        for name in classes:
            cfg = cfgs[name]
            header = [
                '#######################################',
                f'                {name}                ',
                '#######################################'
            ]

            # If stages is provided, filter top-level keys
            if stages is not None:
                if isinstance(stages, str):
                    stages_set = {stages}
                else:
                    stages_set = set(stages)

                # Only keep the keys listed in `stages`
                filtered_cfg = {k: v for k, v in cfg.items() if k in stages_set}
            else:
                filtered_cfg = cfg


            def normalize(obj):
                if isinstance(obj, dict):
                    return {str(k): normalize(v) for k, v in obj.items()}
                elif isinstance(obj, list):
                    return [normalize(v) for v in obj]
                elif isinstance(obj, tuple):
                    return [normalize(v) for v in obj]
                elif hasattr(obj, "item"):  # numpy scalars
                    return obj.item()
                else:
                    return obj
            filtered_cfg = normalize(filtered_cfg)
    
            yaml = YAML()#typ='unsafe')
            yaml.default_flow_style = False  # same as your argument

            stream = StringIO()
            yaml.dump(filtered_cfg, stream)

            yaml_text = stream.getvalue()

            parts.extend(header)
            parts.append(yaml_text)

        return '\n'.join(parts)


    def add_config(self, name):
    #     """Add a copy of self.queue as {name}_queue"""
        if self.queue is None:
            raise QueueConfigError("No queue configuration loaded to copy.")

        if hasattr(self, f'{name}_queue'):
            return

        # Deep-copy the entire QueueConfig instance
        new_queue_config = copy.deepcopy(self.queue)

        setattr(self, f"{name}_queue", new_queue_config)


config = Config()
    

# Safe operators you want to allow
ops = {
    ast.Gt: operator.gt,
    ast.Lt: operator.lt,
    ast.Eq: operator.eq,
    ast.NotEq: operator.ne,
    ast.GtE: operator.ge,
    ast.LtE: operator.le,
}

def safe_eval(expr: str):
    """Evaluate a simple comparison expression safely."""
    tree = ast.parse(expr, mode="eval")

    def _eval(node):
        if isinstance(node, ast.Expression):
            return _eval(node.body)
        if isinstance(node, ast.Compare):
            left = _eval(node.left)
            assert len(node.ops) == 1 and len(node.comparators) == 1
            op = ops[type(node.ops[0])]
            right = _eval(node.comparators[0])
            return op(left, right)
        if isinstance(node, ast.Constant):   # Python 3.8+
            return node.value
        if isinstance(node, ast.Num):        # older Python
            return node.n
        raise ValueError("Unsupported expression:", node)

    return _eval(tree)

def update_key(config, target_key, new_value, path="config", debug = False, skip_section=[]):
    """
    Recursively update all instances of target_key in nested dict/list structures.
    Prints debug info for each update.
    """
    found = False
    if isinstance(config, dict):
        for key, value in config.items():
            if key in skip_section:
                continue
            new_path = f"{path}.{key}"
            if key == target_key:
                if debug:
                    splog.debug(f"{new_path} = {new_value!r}  (type: {type(new_value).__name__})")
                config[key] = new_value
                found = True
            else:
                if update_key(value, target_key, new_value, new_path, debug=debug):
                    found = True

    elif isinstance(config, list):
        for idx, item in enumerate(config):
            new_path = f"{path}[{idx}]"
            if update_key(item, target_key, new_value, new_path, debug=debug):
                found = True

    return found 


def get_inline_comment(container, key):
    """
    Best-effort extraction of an inline comment from a ruamel round-trip node.
    Works for mappings and sequences.
    """
    ca = getattr(container, "ca", None)
    if ca is None or not getattr(ca, "items", None):
        return ""

    item = ca.items.get(key)
    if not item:
        return ""

    # ruamel stores comment tokens in the item metadata; find the first token-like value.
    for part in reversed(item):
        if part is None:
            continue
        if isinstance(part, list):
            for token in part:
                value = getattr(token, "value", None)
                if value:
                    value = value.strip()
                    if value.startswith("#"):
                        return value[1:].strip()
                    return value
        else:
            value = getattr(part, "value", None)
            if value:
                value = value.strip()
                if value.startswith("#"):
                    return value[1:].strip()
                return value

    return ""

def coerce_like(existing, text):
    text_clean = text.strip().lower()

    if isinstance(existing, bool):
        if text_clean in {"1", "true", "yes", "y", "on"}:
            return True
        elif text_clean in {"0", "false", "no", "n", "off"}:
            return False
        else:
            raise ValueError(f"Invalid boolean value: {text}")
    if isinstance(existing, int) and not isinstance(existing, bool):
        return int(text)
    if isinstance(existing, float):
        return float(text)
    return text

def prompt_edit(node, path=()):
    if isinstance(node, dict):
        for key in list(node.keys()):
            value = node[key]
            comment = get_inline_comment(node, key)
            new_path = path + (str(key),)

            if isinstance(value, (dict, list)):
                node[key] = prompt_edit(value, new_path)
                continue

            label = ".".join(new_path)
            suffix = f"  # {comment}" if comment else ""
            prompt = f"{label}{suffix} [{value}]: "

            entered = input(prompt)
            if entered.strip() == "":
                continue

            node[key] = coerce_like(value, entered)

        return node

    if isinstance(node, list):
        for idx, value in enumerate(list(node)):
            comment = get_inline_comment(node, idx)
            new_path = path + (f"[{idx}]",)

            if isinstance(value, (dict, list)):
                node[idx] = prompt_edit(value, new_path)
                continue

            label = "".join(new_path)
            suffix = f"  # {comment}" if comment else ""
            prompt = f"{label}{suffix} [{value}]: "

            entered = input(prompt)
            if entered.strip() == "":
                continue

            node[idx] = coerce_like(value, entered)

        return node

    return node


import click
def show_config_opt(f):
    f = click.option('--show_config', is_flag=True, default=False, help='Show the final config with out excecuting commands')(f)
    return f
def show_config(pipe=True, queue=True):
    if pipe:
        print(config.full_str(classes = ['Pipe']))
        # print('---------Pipeline Config------------')
        # print(json.dumps(config.pipe, indent=2, default=str))

    if queue:
        print(config.full_str(classes = ['Queue']))
        # print('---------Queue Config------------')
        # print(json.dumps(config.queue._obj, indent=2, default=str))
        

