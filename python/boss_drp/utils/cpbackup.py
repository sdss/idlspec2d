import os
import shutil

def cpbackup(filename):
    """
    Make a backup copy of the specified file by appending '.1', '.2', etc.
    The first unused number is used as a suffix.
    """

    if not os.path.isfile(filename):
        return  # File does not exist; do nothing

    num = 1
    while True:
        backname = f"{filename}.{num}"
        if not os.path.exists(backname):
            break
        num += 1

    shutil.copy2(filename, backname)  # copy2 preserves metadata
