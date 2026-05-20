from tqdm.auto import tqdm


import sys
import time
from tqdm import tqdm
from boss_drp.utils.splog import splog

class progress:
    def __init__(self, iterable=None, total=None, desc=None, log_every=10, print_func=splog.info, **kwargs):
        self.iterable = iterable
        self.total = total
        self.desc = desc or "Progress"
        self.log_every = log_every
        self.is_tty = sys.stderr.isatty()
        self.print_func = print_func

        self._i = 0
        self._start = time.time()
        self._last_iter_n = 0  # tracks auto-iteration increments
        self._last_iter_cnt = 0
        if self.is_tty:
            self._pbar = tqdm(iterable=iterable, total=total, desc=desc, **kwargs)
        else:
            self._pbar = None

    def __iter__(self):
        if self.iterable is None:
            raise TypeError("tqdm iterable is None")

        for item in self.iterable:
            self._i += 1
            self._last_iter_n += 1
            yield item

    def update(self, n=1):
        # If iterable mode is being used, subtract the implicit +1 per iteration
        if self.iterable is not None:
            n = max(0, n - self._last_iter_n)
            self._last_iter_n = 0

        if n <= 0:
            return

        self._i += n

        if self.is_tty:
            if self._pbar:
                self._pbar.update(n)
        else:
            if self._i % self.log_every == 0:
                elapsed = time.time() - self._start
                rate = self._i / elapsed if elapsed > 0 else 0

                msg = f"{self.desc}: {self._i}"
                if self.total:
                    pct = 100 * self._i / self.total
                    msg += f"/{self.total} ({pct:.1f}%)"
                msg += f" | {rate:.2f} it/s"
                self._last_iter_cnt = self._i

                self.print_func(msg)#, file=sys.stderr)

    def close(self):
        if not self.is_tty:
            if self._i > self._last_iter_cnt:
                elapsed = time.time() - self._start
                rate = self._i / elapsed if elapsed > 0 else 0
                msg = f"{self.desc}: {self._i}"
                if self.total:
                    pct = 100 * self._i / self.total
                    msg += f"/{self.total} ({pct:.1f}%)"
                msg += f" | {rate:.2f} it/s"
                self.print_func(msg)#, file=sys.stderr)
        if self._pbar:
            self._pbar.close()