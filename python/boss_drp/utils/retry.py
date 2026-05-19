import time
import sys
import traceback

def get_callable_name(func):
    if hasattr(func, "__name__"):
        return func.__name__
    if hasattr(func, "__class__"):
        return func.__class__.__name__
    return str(func)

def retry(func, *args, retries=3, delay=5, exceptions=(Exception,), logger=print, noerr=False, **kwargs):
    """
    Retries a function call with specified retries and delay on failure.

    Parameters:
    - func (callable): The function to execute.
    - retries (int): The number of times to retry. Default is 3.
    - delay (int or float): Delay between retries in seconds. Default is 1 second.
    - exceptions (tuple): Tuple of exceptions to catch and retry. Default is all exceptions. [ex (ValueError,)]
    - splog (callable): the logger function in use (defaults to print)
    - *args: Positional arguments for the function.
    - **kwargs: Keyword arguments for the function.

    Returns:
    - The result of the function if successful.
    
    Raises
    - The exception if the function fails after the given number of retries.
    """
    attempt = 0

    func_name = get_callable_name(func)

    while attempt < retries:
        try:
            return func(*args, **kwargs)
        except exceptions as e:
            exctype, value, tb = sys.exc_info()
            if tb is not None:
                formatted_tb = "".join(traceback.format_exception(exctype, value, tb))
                logger("Exception caught:\n" + formatted_tb)
            attempt += 1
            if attempt >= retries:
                if not noerr:
                    raise
                logger(f"{func_name} failed: {e}.")
                return
            logger(f"{func_name} failed: {e}. Retrying in {delay} seconds...")
            time.sleep(delay)

# Example usage:

#def sample_function(x, y):
#    if x == 0:
#        raise ValueError("x cannot be zero!")
#    return y / x

#try:
#    result = retry(sample_function, retries=5, delay=2, exceptions=(ValueError,), x=0, y=10)
#    print(f"Success: {result}")
#except Exception as e:
#    print(f"Final failure: {e}")
