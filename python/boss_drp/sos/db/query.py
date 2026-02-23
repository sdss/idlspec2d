from sdssdb.peewee.sdss5db import opsdb, targetdb
from datetime import datetime, timedelta
import os
from astropy.time import Time
import numpy as np

b_camera = 'b2' if os.getenv('OBSERVATORY') == 'LCO' else 'b1'
r_camera = 'r2' if os.getenv('OBSERVATORY') == 'LCO' else 'r1'

def robodamus(sjd = None):
    pred = opsdb.PredSnr
    cf = opsdb.CameraFrame
    exp = opsdb.Exposure
    design = targetdb.Design
    cfg = opsdb.Configuration
    queue = opsdb.Queue

    r1_db = opsdb.Camera.get(label=r_camera)
    b1_db = opsdb.Camera.get(label=b_camera)

    if not sjd:
        use_time = datetime.now() - timedelta(hours=3)
        conditions_pred = [pred.gfa_date_obs > use_time]
        conditions_sos  = [exp.start_time >= use_time]
    else:
        sjd_start = Time(sjd - 0.3, format="mjd").datetime
        sjd_end   = Time(sjd + 0.7, format="mjd").datetime
        use_time = sjd_start
        conditions_pred = [pred.gfa_date_obs > use_time]
        conditions_sos  = [exp.start_time >= use_time]
        conditions_pred += [pred.gfa_date_obs >= sjd_start, pred.gfa_date_obs < sjd_end,]
        conditions_sos += [exp.start_time >= sjd_start, exp.start_time < sjd_end, ]



    predQuery = pred.select(pred.camera_pk,
                            pred.gfa_date_obs,
                            pred.pred_value,
                            pred.design_mode_label)\
                     .where(*conditions_pred)\
                     .order_by(pred.gfa_date_obs.asc()).dicts()

    sosQuery = cf.select(cf.sn2, cf.camera, exp.start_time,
                         design.design_mode)\
                 .join(exp)\
                 .join(cfg)\
                 .join(design, on=(design.design_id == cfg.design_id))\
                 .where(*conditions_sos)\
                .order_by(exp.start_time.asc()).dicts()

    b_pred = [[p["gfa_date_obs"],
               p["pred_value"],
               p["design_mode_label"]]
               for p in predQuery if p["camera_pk"] == b1_db.pk]
    r_pred = [[p["gfa_date_obs"],
               p["pred_value"],
               p["design_mode_label"]]
               for p in predQuery if p["camera_pk"] == r1_db.pk]

    b_sos = [[s["start_time"] + timedelta(minutes=7),
              s["sn2"],
              s["design_mode"]]
              for s in sosQuery if s["camera"] == b1_db.pk]
    r_sos = [[s["start_time"] + timedelta(minutes=7),
              s["sn2"],
              s["design_mode"]]
              for s in sosQuery if s["camera"] == r1_db.pk]

    current = queue.select(design.design_mode)\
                   .join(design)\
                   .where(queue.position == -1).dicts()
    
    if len(current) == 1:
        mode = current[0]["design_mode"]
        last_time = max([b[0] for b in b_sos])
        b_sos.append([last_time + timedelta(minutes=15),
                      None,
                      mode])

    b_pred_trim = list()
    r_pred_trim = list()

    b_colors = list()
    r_colors = list()
    sos_shapes = list()
    for b in b_sos:
        mode = b[2].replace("_no_apogee_skies", "")
        if "bright" in mode:
            b_colors.append("grey")
            r_colors.append("grey")
            sos_shapes.append("x")
            # tom doesn't predict bright, but show something
            mode = "dark_faint"
        elif "plane" in mode:
            b_colors.append("blue")
            r_colors.append("red")
            sos_shapes.append("star-triangle-up")
        else:
            b_colors.append("blue")
            r_colors.append("red")
            sos_shapes.append("circle")
        start_time = b[0] - timedelta(minutes=7)
        end_time = b[0] + timedelta(minutes=7)
        b_times = np.array([b[0] for b in b_pred])
        b_modes = np.array([b[2].strip() for b in b_pred])

        w_time = np.logical_and(
                    b_times > start_time,
                    b_times < end_time,
        )

        w_b = np.where(np.logical_and(
                    b_modes == mode,
                    w_time
        ))

        b_pred_trim.extend(
            [b_pred[i] for i in w_b[0]]
        )

        r_times = np.array([r[0] for r in r_pred])
        r_modes = np.array([r[2].strip() for r in r_pred])

        w_time = np.logical_and(
                    r_times > start_time,
                    r_times < end_time,
        )

        w_r = np.where(np.logical_and(
                    r_modes == mode,
                    w_time
        ))

        r_pred_trim.extend(
            [r_pred[i] for i in w_r[0]]
        )

    outformat = "%Y-%m-%d %H:%M:%S"
    output = {
        "b_sn" : [b[1] for b in b_pred_trim],
        "r_sn" : [r[1] for r in r_pred_trim],
        "b_times" : [b[0].strftime(outformat) for b in b_pred_trim],
        "r_times" : [r[0].strftime(outformat) for r in r_pred_trim],
        "b_sn_sos" : [b[1] for b in b_sos],
        "r_sn_sos" : [r[1] for r in r_sos],
        "b_times_sos" : [b[0].strftime(outformat) for b in b_sos],
        "r_times_sos" : [r[0].strftime(outformat) for r in r_sos],
        "b_colors": b_colors,
        "r_colors": r_colors,
        "sos_shapes": sos_shapes,
    }

    return output