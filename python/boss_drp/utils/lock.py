import os
import time

def lock(file, pause=5, niter=None):
    """Attempt to acquire a file lock by creating a symlink. Retry on failure."""
    i = 0
    while True:
        if niter is not None and i == niter:
            break
        i += 1
        try:
            os.symlink(file, file + '.lock')
            return True
        except FileExistsError:
            print(f"Lock already acquired ({file}). Retrying in {pause} seconds...")
            time.sleep(pause)
        except Exception as e:
            print(f"An error occurred: {e}")
            return False
    return False

def unlock(file):
    """Release the file lock by removing the symlink."""
    try:
        os.unlink(file + '.lock')
        return(True)
    except FileNotFoundError:
        return(True)
    except Exception as e:
        print(f"An error occurred while releasing the lock: {e}")
        return(False)
    return(False)

"""
# Example usage
file = 'path/to/your/file.txt'

if lock(file):
    try:
        # Perform your file operations here
        print("Performing file operations.")
        time.sleep(10)  # Simulate long-running task
    finally:
        unlock(file)
else:
    print("Could not acquire lock. Exiting.")
"""
