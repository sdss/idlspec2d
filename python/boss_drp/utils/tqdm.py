from tqdm.auto import tqdm


import sys
import time
from tqdm import tqdm
from boss_drp.utils.splog import splog

# def progress(iterable, total=None, desc=None, log_every=None, print_func=splog.info):
#     if sys.stderr.isatty():
#         # Interactive → full tqdm
#         yield from tqdm(iterable, total=total, desc=desc)
#     else:
#         if log_every is None:
#             log_every = max(1, total // 100) if total else 10
#         # Non-interactive → periodic logging
#         start = time.time()
#         for i, item in enumerate(iterable, 1):
#             yield item
#             if i % log_every == 0:
#                 elapsed = time.time() - start
#                 rate = i / elapsed if elapsed > 0 else 0
#                 msg = f"{desc or 'Progress'}: {i}"
#                 if total:
#                     pct = 100 * i / total
#                     msg += f"/{total} ({pct:.1f}%)"
#                 msg += f" | {rate:.2f} it/s"
#                 print_func(msg)#, file=sys.stderr)


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