import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.transforms as transforms
import matplotlib as mpl
from astropy.coordinates import Angle, get_body, EarthLocation
from astropy.table import Table
from astropy.time import Time
import astropy.units as u


def make_jatakam_plot(
    test_date_utc_time,
    location,
    date_str,
    nakṣatra_extent,
    rāśi_extent,
    ayanāṃśa,
    tithi,
    tithi_name,
    vāra,
    nakṣatra,
    pāda,
    yoga,
    karaṇa,
    moon_lambda,
    sun_lambda,
    lagna_lambda,
    rāśi_names_file="rashi_names.tex",
    plotfile="jatakam_at_test_time.pdf",
):
    ## Sanskrit typesetting using XeLaTeX ##
    mpl.use("pgf")

    ## TeX preamble ##
    preamble = r"""\usepackage{fontspec}
               \usepackage{polyglossia}
               \usepackage[T1]{fontenc}
               \usepackage{tikz}
               \setmainlanguage{english}
               \setotherlanguages{sanskrit}
               \newfontfamily\devanagarifont[Script=Devanagari]{Sanskrit 2003}
               \XeTeXgenerateactualtext 1
               \newcommand{\sam}[1]{\begin{sanskrit}#1\end{sanskrit}}"""

    params = {
        "font.family": "serif",
        "text.usetex": True,
        "pgf.rcfonts": False,
        "pgf.texsystem": "xelatex",
        "pgf.preamble": preamble,
    }
    mpl.rcParams.update(params)

    fig, ax = plt.subplots()

    ax.vlines([0, 2, 6, 8], ymin=0, ymax=8, color="b")
    ax.vlines([4, 4], ymin=[0, 6], ymax=[2, 8], color="b")
    ax.hlines([0, 2, 6, 8], xmin=0, xmax=8, color="b")
    ax.hlines([4, 4], xmin=[0, 6], xmax=[2, 8], color="b")

    coord_for_rāśi = [
        (3, 7),
        (5, 7),
        (7, 7),
        (7, 5),
        (7, 3),
        (7, 1),
        (5, 1),
        (3, 1),
        (1, 1),
        (1, 3),
        (1, 5),
        (1, 7),
    ]
    rāśi_tab = Table.read(rāśi_names_file, format="latex")

    rāśi_names = rāśi_tab["names"].data

    print(rāśi_names)
    for i in range(len(rāśi_names)):
        ax.text(
            coord_for_rāśi[i][0],
            coord_for_rāśi[i][1],
            "test",
            ha="center",
            va="center",
            color="b",
        )

    ax.set_aspect("equal")
    ax.set_title(r"\sam{{ मम}}")
    ax.axis("off")
    plt.savefig(plotfile)

    return 0


make_jatakam_plot(
    0,
    0,
    0,
    0,
    30,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
)
