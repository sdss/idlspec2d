#!/usr/bin/env python3

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

from datetime import datetime
from robodamus import robodamus
import os

PLOTLY_TO_MPL_MARKER = {
    "circle": "o",
    "x": "x",
    "star-triangle-up": "^",   # closest common equivalent
}

def plot_robodamus(output, filename="robodamus.png"):
    outformat = "%Y-%m-%d %H:%M:%S"

    # Convert time strings back to datetime
    b_times = [datetime.strptime(t, outformat) for t in output["b_times"]]
    r_times = [datetime.strptime(t, outformat) for t in output["r_times"]]

    b_times_sos = [datetime.strptime(t, outformat) for t in output["b_times_sos"]]
    r_times_sos = [datetime.strptime(t, outformat) for t in output["r_times_sos"]]

    fig, (ax_b, ax_r) = plt.subplots(
        2, 1, figsize=(10, 7), sharex=True
    )

    # ---- B camera subplot ----
    ax_b.plot(
        b_times,
        output["b_sn"],
        label="B predicted",
        color="blue",alpha=.5
    )

    for t, sn, c, m in zip(
        b_times_sos,
        output["b_sn_sos"],
        output["b_colors"],
        output["sos_shapes"],
    ):
        if sn is not None:
            marker = PLOTLY_TO_MPL_MARKER.get(m, "o")
            if marker == "x":
                ax_b.scatter(t, sn, color=c, marker=marker, s=80)
            else:
                ax_b.scatter(t, sn, color=c, marker=marker, s=80, edgecolors="black",linewidths=0.8,)

    ax_b.set_ylabel("B S/N")
    ax_b.legend()
    ax_b.grid(True)

    # ---- R camera subplot ----
    ax_r.plot(
        r_times,
        output["r_sn"],
        label="R predicted",
        color="red",alpha=.5
    )

    for t, sn, c, m in zip(
        r_times_sos,
        output["r_sn_sos"],
        output["r_colors"],
        output["sos_shapes"],
    ):
        if sn is not None:
            marker = PLOTLY_TO_MPL_MARKER.get(m, "o")
            if marker == "x":
                ax_r.scatter(t, sn, color=c, marker=marker, s=80)
            else:
                ax_r.scatter(t, sn, color=c, marker=marker, s=80,edgecolors="black",linewidths=0.8,)

    ax_r.set_ylabel("R S/N")
    ax_r.set_xlabel("Time")
    ax_r.legend()
    ax_r.grid(True)


    legend_elements = [
      Line2D([0], [0], color='blue', label='B predicted'),
      Line2D([0], [0], marker='o', color='w', label='Dark',
           markerfacecolor='blue', markersize=8),
      Line2D([0], [0], marker='^', color='w', label='Plane',
           markerfacecolor='blue', markersize=8),
      Line2D([0], [0], marker='x', color='grey', label='Bright',
           markersize=8),
     ]

    ax_b.legend(handles=legend_elements, loc="upper left", fontsize="xx-small")


    legend_elements = [
      Line2D([0], [0], color='red', label='R predicted'),
      Line2D([0], [0], marker='o', color='w', label='Dark',
           markerfacecolor='red', markersize=8),
      Line2D([0], [0], marker='^', color='w', label='Plane',
           markerfacecolor='red', markersize=8),
      Line2D([0], [0], marker='x', color='grey', label='Bright',
           markersize=8),
     ]
    ax_r.legend(handles=legend_elements, loc="upper left",fontsize="xx-small")

    fig.text(0.99,0.01,'Generated: '+datetime.now().strftime("%Y-%m-%d %H:%M:%S"),ha='right',va='bottom',fontsize=6,color='gray')
    fig.text(0.01,0.99,'Bright predictions use dark faint predicitons',ha='left',va='top',fontsize=6,color='gray')
    fig.suptitle("Robodamus S/N Prediction vs Observed")
    fig.autofmt_xdate()

    plt.tight_layout(rect=[0, 0, 1, 0.96])
    plt.savefig(filename, dpi=150)
    plt.close(fig)

def _get_sjd():
    return Time.now().mjd + 0.3

if __name__ == "__main__":
    import argparse
    from astropy.time import Time
    import numpy as np
    parser = argparse.ArgumentParser(description='Plot Robodamus')
    parser.add_argument('--mjd','-m', help='SJD of reduction', type=float, default=Time.now().mjd + 0.3)
    args = parser.parse_args()
    if args.mjd < 0:
       args.mjd = _get_sjd()+args.mjd
    args.mjd = int(np.floor(args.mjd))
    print(args.mjd)
    output = robodamus(sjd=args.mjd)
    plot_robodamus(output, os.path.join('/data/boss/sos/tests/robodamus',f"robodamus_sn_{int(np.floor(args.mjd))}.png"))