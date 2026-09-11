import os
import socket
import time

#TODO: replace .lock with .lock.fits.gz or .lock.parquet etc

def lock(file, pause=5, niter=None, logger = print):
    """Acquire a lock using atomic O_CREAT|O_EXCL."""
    lockfile = str(file) + ".lock"
    attempts = 0

    while niter is None or attempts < niter:
        attempts += 1

        try:
            fd = os.open(
                lockfile,
                os.O_WRONLY | os.O_CREAT | os.O_EXCL,
                0o644,
            )

            try:
                with os.fdopen(fd, "w") as f:
                    f.write(f"pid={os.getpid()}\n")
                    f.write(f"host={socket.gethostname()}\n")
                    f.write(f"time={time.time()}\n")
            except Exception:
                # We successfully acquired the lock, but failed
                # while writing its owner information.
                try:
                    os.unlink(lockfile)
                except OSError:
                    pass
                raise

            return True

        except FileExistsError:
            logger(
                f"Lock already acquired ({file}). "
                f"Retrying in {pause} seconds..."
            )
            time.sleep(pause)

        except OSError as e:
            logger(
                f"Error acquiring lock {file}: "
                f"{type(e).__name__}: {e}"
            )
            return False

    return False


def unlock(file, logger = print):
    """Release the lock."""
    lockfile = str(file) + ".lock"

    try:
        os.unlink(lockfile)
        return True

    except FileNotFoundError:
        return True

    except OSError as e:
        logger(
            f"Error releasing lock {file}: "
            f"{type(e).__name__}: {e}"
        )
        return False


