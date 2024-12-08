import time
from datetime import datetime

def monitor_job(logger, queue1, pause = 300, jobname = '', return_status=False):
    percomp1 = 0
    q1done=False
    while percomp1 < 100:
        if queue1 is not None and not q1done:
            if queue1.get_job_status() is None:
                logger.info(f'Failure in slurm queue for {jobname}')
            t_percomp1 = queue1.get_percent_complete() if not q1done else 100
            if t_percomp1 != percomp1:
                percomp1 = t_percomp1
                logger.info(f'{jobname} {percomp1}% complete at {datetime.today().ctime()}')
        elif not q1done:
            percomp1 = 100
            status = f'{jobname} not submitted at {datetime.today().ctime()}'
            logger.info(status)
        if percomp1 == 100 and not q1done:
            q1done=True
            status = f'Complete {jobname} '
            logger.info(status)
        else:
            time.sleep(pause)
    if return_status:
        return logger, status
    return logger
