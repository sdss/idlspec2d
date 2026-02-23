from boss_drp.utils.splog import splog
from boss_drp.sos.sos_classes import SOS_config, Consts

import os
import time


class LicenseLock:
    MIN_WAIT_TIME = int(Consts().licensePause) # seconds

    def get_lockFile(self):
        return os.path.join(SOS_config.controlDir,"license.lock")

    def create_or_refresh_lock(self):
        """Create or refresh the lock file's timestamp."""
        try:
            with open(get_lockFile(), "w") as f:
                f.write("")  # Create an empty file or update its timestamp
        except:
            pass

    def is_lock_file_recent(self):
        """Check if the lock file is less than 5 seconds old."""
        if not os.path.exists(get_lockFile()):
            return False
        lock_age = time.time() - os.path.getmtime(get_lockFile())
        return lock_age < self.MIN_WAIT_TIME

    def cleanup(self, logger=splog.info):
        """Remove the lock file if it is stale."""
        try:
            if os.path.exists(get_lockFile()):
                lock_age = time.time() - os.path.getmtime(get_lockFile())
                if lock_age >= self.MIN_WAIT_TIME:
                    os.remove(get_lockFile())
                    logger("Stale license lock file removed.")
        except:
            pass

    def check(self, message=None, logger=splog.info):
        if not SOS_config.pause:
            return
        """Check the lock file timing."""
        try:
            # Wait until the lock file is old enough
            i = 0
            while self.is_lock_file_recent():
                logger("License Lock file is recent. Waiting...")
                time.sleep(1)  # Polling interval
                i+=1
                if i >= 2*self.MIN_WAIT_TIME:
                    logger(f'Paused for 2x{self.MIN_WAIT_TIME}s...Forcing start')
                    break

            # Refresh the lock file and proceed
            self.create_or_refresh_lock()
            #logger("Lock refreshed. Starting process...")
            return
    
        except:
            pass
        finally:
            # Remove the lock file if it becomes stale
            self.cleanup(logger=logger)

licenselock = LicenseLock()
