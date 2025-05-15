import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.transforms as transforms
import matplotlib as mpl
from astropy import coordinates as coor
from astropy.coordinates import get_body, EarthLocation, SkyCoord
from astropy.table import Table
from astropy.time import Time
import astropy.units as u
from mpl_toolkits.axes_grid1 import make_axes_locatable
from adjustText import adjust_text
from indic_transliteration import sanscript
from indic_transliteration.sanscript import transliterate
from skyfield import almanac
from skyfield.api import load

#### Colors for grahas ####
# Sun - Surya - Red - circle
# Moon - Chandra - White - semicircle
# Mars - Mangala - Red - triangle
# Mercury - Budha - Green - bow (dhanusha)
# Jupiter - Guru - Yellow - lotus
# Venus - Shukra - White - square
# Saturn - Shani - Black - snake
# Rahu - Blue - crocodile
# Ketu - Grey - sword


def calc_rahu_ketu_pos(test_date_utc_time, location):
    rahu_speed = 232238 / 1577916450 * 360  # deg per day

    ts = load.timescale()
    eph = load("de421.bsp")

    t0 = ts.tt_jd(test_date_utc_time.jd - 15)
    t1 = ts.tt_jd(test_date_utc_time.jd + 15)

    t, y = almanac.find_discrete(t0, t1, almanac.moon_nodes(eph))

    sel_rahu = np.where(y == 1)

    time_of_transit_of_moon_at_rahu = Time(
        t[sel_rahu].utc_iso(delimiter="T", places=0)[0][:-1], format="isot", scale="tt"
    )

    lon_of_moon_at_time_of_transit = get_body(
        "moon", time=time_of_transit_of_moon_at_rahu, location=location
    ).geocentrictrueecliptic.lon.value

    delta_t = test_date_utc_time - time_of_transit_of_moon_at_rahu

    rahu_lambda = (lon_of_moon_at_time_of_transit - delta_t.value * rahu_speed) % 360
    ketu_lambda = (rahu_lambda + 180) % 360

    return rahu_lambda, ketu_lambda


def plot_moon_phase(day, drawing_origin, radius, fig, ax):
    center = np.array([0, 0])
    trans = fig.dpi_scale_trans + transforms.ScaledTranslation(
        drawing_origin[0], drawing_origin[1], ax.transData
    )
    if 0 < day <= 7.5:
        circle = mpatches.Circle(
            center,
            radius,
            ec="darkslategrey",
            fc="w",
            lw=0.5,
            transform=trans,
            clip_on=False,
            zorder=10,
        )
        ellipse = mpatches.Ellipse(
            center,
            2 * radius * np.sin(np.pi / 2 - day * 12 * np.pi / 180),
            2 * radius,
            fc="darkslategrey",
            ec=None,
            transform=trans,
            clip_on=False,
            zorder=10,
        )
        wedge = mpatches.Wedge(
            center,
            radius,
            90,
            270,
            fc="darkslategrey",
            ec=None,
            transform=trans,
            clip_on=False,
            zorder=10,
        )
        ax.add_patch(circle)
        ax.add_patch(wedge)
        ax.add_patch(ellipse)

    if 7.5 < day <= 15:
        circle = mpatches.Circle(
            center,
            radius,
            ec="darkslategrey",
            fc="w",
            lw=0.5,
            transform=trans,
            clip_on=False,
            zorder=10,
        )
        ellipse = mpatches.Ellipse(
            center,
            2 * radius * np.sin(np.pi / 2 - day * 12 * np.pi / 180),
            2 * radius,
            fc="w",
            ec=None,
            transform=trans,
            clip_on=False,
            zorder=10,
        )
        wedge = mpatches.Wedge(
            center,
            radius,
            90,
            270,
            fc="darkslategrey",
            ec=None,
            transform=trans,
            clip_on=False,
            zorder=10,
        )
        ax.add_patch(circle)
        ax.add_patch(wedge)
        ax.add_patch(ellipse)

    if 15 < day <= 22.5:
        circle = mpatches.Circle(
            center,
            radius,
            ec="darkslategrey",
            fc="darkslategrey",
            lw=0.5,
            transform=trans,
            clip_on=False,
            zorder=10,
        )
        ellipse = mpatches.Ellipse(
            center,
            2 * radius * np.sin(np.pi / 2 - (day - 15) * 12 * np.pi / 180),
            2 * radius,
            fc="w",
            ec=None,
            transform=trans,
            clip_on=False,
            zorder=10,
        )
        wedge = mpatches.Wedge(
            center,
            radius,
            90,
            270,
            fc="w",
            ec=None,
            transform=trans,
            clip_on=False,
            zorder=10,
        )
        ax.add_patch(circle)
        ax.add_patch(wedge)
        ax.add_patch(ellipse)

    if 22.5 < day <= 30:
        circle = mpatches.Circle(
            center,
            radius,
            ec="darkslategrey",
            fc="darkslategrey",
            lw=0.5,
            transform=trans,
            clip_on=False,
            zorder=10,
        )
        ellipse = mpatches.Ellipse(
            center,
            2 * radius * np.sin(np.pi / 2 - (day - 15) * 12 * np.pi / 180),
            2 * radius,
            fc="darkslategrey",
            ec=None,
            transform=trans,
            clip_on=False,
            zorder=10,
        )
        wedge = mpatches.Wedge(
            center,
            radius,
            90,
            270,
            fc="w",
            ec=None,
            transform=trans,
            clip_on=False,
            zorder=10,
        )
        ax.add_patch(circle)
        ax.add_patch(wedge)
        ax.add_patch(ellipse)

    return 0


def plot_sun(drawing_origin, radius, fig, ax):
    center = np.array([0, 0])
    trans = fig.dpi_scale_trans + transforms.ScaledTranslation(
        drawing_origin[0], drawing_origin[1], ax.transData
    )

    circle = mpatches.Circle(
        center,
        radius,
        ec="orange",
        fc="yellow",
        transform=trans,
        clip_on=False,
        zorder=10,
    )
    ax.add_patch(circle)

    list_of_patches = []
    for theta in np.arange(0, 2 * np.pi, np.pi / 6):
        triangle = mpatches.RegularPolygon(
            [1.2 * radius * np.cos(theta), 1.2 * radius * np.sin(theta)],
            3,
            radius=0.2 * radius,
            ec="yellow",
            fc="orange",
            orientation=theta + np.pi / 6,
            transform=trans,
            clip_on=False,
            linewidth=0.5,
            zorder=10,
        )
        ax.add_patch(triangle)

    return 0


def plot_inner_graha_phase(graha, angle_to_sun, drawing_origin, radius, fig, ax):
    center = np.array([0, 0])
    trans = fig.dpi_scale_trans + transforms.ScaledTranslation(
        drawing_origin[0], drawing_origin[1], ax.transData
    )

    if graha == "budha":
        bright_side_color = "green"
        dark_side_color = "darkgreen"
    elif graha == "shukra":
        bright_side_color = "white"
        dark_side_color = "grey"
    else:
        print(
            "That is not an inner graha! Use the outer graha function to plot ", graha
        )
        return

    if 0 < (angle_to_sun / 12) <= 7.5:
        circle = mpatches.Circle(
            center,
            radius,
            ec=dark_side_color,
            fc=bright_side_color,
            lw=0.5,
            transform=trans,
            clip_on=False,
            zorder=10,
        )
        ellipse = mpatches.Ellipse(
            center,
            2 * radius * np.sin(np.pi / 2 - (angle_to_sun / 12) * 12 * np.pi / 180),
            2 * radius,
            fc=dark_side_color,
            ec=None,
            transform=trans,
            clip_on=False,
            zorder=10,
        )
        wedge = mpatches.Wedge(
            center,
            radius,
            90,
            270,
            fc=dark_side_color,
            ec=None,
            transform=trans,
            clip_on=False,
            zorder=10,
        )
        ax.add_patch(circle)
        ax.add_patch(wedge)
        ax.add_patch(ellipse)

    if 7.5 < (angle_to_sun / 12) <= 15:
        circle = mpatches.Circle(
            center,
            radius,
            ec=dark_side_color,
            fc=bright_side_color,
            lw=0.5,
            transform=trans,
            clip_on=False,
            zorder=10,
        )
        ellipse = mpatches.Ellipse(
            center,
            2 * radius * np.sin(np.pi / 2 - (angle_to_sun / 12) * 12 * np.pi / 180),
            2 * radius,
            fc=bright_side_color,
            ec=None,
            transform=trans,
            clip_on=False,
            zorder=10,
        )
        wedge = mpatches.Wedge(
            center,
            radius,
            90,
            270,
            fc=dark_side_color,
            ec=None,
            transform=trans,
            clip_on=False,
            zorder=10,
        )
        ax.add_patch(circle)
        ax.add_patch(wedge)
        ax.add_patch(ellipse)

    if 15 < (angle_to_sun / 12) <= 22.5:
        circle = mpatches.Circle(
            center,
            radius,
            ec=dark_side_color,
            fc=bright_side_color,
            lw=0.5,
            transform=trans,
            clip_on=False,
            zorder=10,
        )
        ellipse = mpatches.Ellipse(
            center,
            2
            * radius
            * np.sin(np.pi / 2 - ((angle_to_sun / 12) - 15) * 12 * np.pi / 180),
            2 * radius,
            fc=bright_side_color,
            ec=None,
            transform=trans,
            clip_on=False,
            zorder=10,
        )
        wedge = mpatches.Wedge(
            center,
            radius,
            270,
            90,
            fc=dark_side_color,
            ec=None,
            transform=trans,
            clip_on=False,
            zorder=10,
        )
        ax.add_patch(circle)
        ax.add_patch(wedge)
        ax.add_patch(ellipse)

    if 22.5 < (angle_to_sun / 12) <= 30:
        circle = mpatches.Circle(
            center,
            radius,
            ec=dark_side_color,
            fc=bright_side_color,
            lw=0.5,
            transform=trans,
            clip_on=False,
            zorder=10,
        )
        ellipse = mpatches.Ellipse(
            center,
            2
            * radius
            * np.sin(np.pi / 2 - ((angle_to_sun / 12) - 15) * 12 * np.pi / 180),
            2 * radius,
            fc=dark_side_color,
            ec=None,
            transform=trans,
            clip_on=False,
            zorder=10,
        )
        wedge = mpatches.Wedge(
            center,
            radius,
            270,
            90,
            fc=dark_side_color,
            ec=None,
            transform=trans,
            clip_on=False,
            zorder=10,
        )
        ax.add_patch(circle)
        ax.add_patch(wedge)
        ax.add_patch(ellipse)

    return 0


def plot_outer_graha(graha, drawing_origin, radius, fig, ax):
    center = np.array([0, 0])
    trans = fig.dpi_scale_trans + transforms.ScaledTranslation(
        drawing_origin[0], drawing_origin[1], ax.transData
    )

    if graha == "mangala":
        mcolor = "red"
    elif graha == "guru":
        mcolor = "yellow"
    elif graha == "shani":
        mcolor = "black"
    else:
        print("That is not an outer graha! Use the inner graha function to plot", graha)
        return

    circle = mpatches.Circle(
        center, radius, ec="w", fc=mcolor, lw=0.5, transform=trans, clip_on=False
    )
    ax.add_patch(circle)

    return 0


def plot_rahu(drawing_origin, scale, fig, ax):
    center = np.array([0, 0])
    trans = transforms.Affine2D().scale(scale) + transforms.ScaledTranslation(
        drawing_origin[0], drawing_origin[1], ax.transData
    )

    from matplotlib.path import Path

    verts = [
        (-1, 1),  # Vert1
        (1.0, 1.0),  # Vert2
        (1.0, -1.0),  # Pivot1
        (0, -1.0),  # Vert3
        (-1.0, -1.0),  # Pivot2
        (-1.0, 1.0),
        (-1.0, 1.0),
    ]

    codes = [
        Path.MOVETO,
        Path.LINETO,
        Path.CURVE3,
        Path.LINETO,
        Path.CURVE3,
        Path.LINETO,
        Path.CLOSEPOLY,
    ]

    path = Path(verts, codes)

    patch = mpatches.PathPatch(
        path,
        facecolor="blue",
        edgecolor="lightblue",
        lw=0.5 * scale / 5,
        transform=trans,
        clip_on=False,
    )
    ax.add_patch(patch)

    return 0


def plot_ketu(drawing_origin, scale, fig, ax):
    center = np.array([0, 0])
    trans = transforms.Affine2D().scale(scale) + transforms.ScaledTranslation(
        drawing_origin[0], drawing_origin[1], ax.transData
    )

    from matplotlib.path import Path

    verts = [
        (-1, 1),  # Vert1
        (1.0, 1.0),  # Vert2
        (1.0, -1.0),  # Vert3
        (-1.0, -1.0),  # Vert4
        (0, 0),  # Vert5
        (-1.0, 1.0),  # Vert6
        (-1.0, 1.0),
    ]

    codes = [
        Path.MOVETO,
        Path.LINETO,
        Path.LINETO,
        Path.LINETO,
        Path.LINETO,
        Path.LINETO,
        Path.CLOSEPOLY,
    ]

    path = Path(verts, codes)

    patch = mpatches.PathPatch(
        path,
        facecolor="dimgrey",
        edgecolor="lightgrey",
        lw=0.5 * scale / 5,
        transform=trans,
        clip_on=False,
    )
    ax.add_patch(patch)

    return 0


def translit_str(some_str, language):
    lang_dict = {
        "Devanagari": "sanskrit",
        "Grantha": "grantha",
        "Kannada": "kannada",
        "Malayalam": "malayalam",
        "Tamil": "tamil",
        "Telugu": "telugu",
    }
    lang_translit_obj_dict = {
        "Devanagari": sanscript.DEVANAGARI,
        "Grantha": sanscript.GRANTHA,
        "Kannada": sanscript.KANNADA,
        "Malayalam": sanscript.MALAYALAM,
        "Tamil": sanscript.TAMIL,
        "Telugu": sanscript.TELUGU,
    }
    if language == "Devanagari":
        translit = some_str
    else:
        translit = transliterate(
            some_str,
            lang_translit_obj_dict["Devanagari"],
            lang_translit_obj_dict[language],
        )

    s = (
        r"\begin{"
        + lang_dict[language]
        + r"}"
        + translit
        + r"\end{"
        + lang_dict[language]
        + r"}"
    )
    return s


## TeX preamble ##
preamble = r"""\usepackage{fontspec}
           \usepackage{polyglossia}
           \usepackage[T1]{fontenc}
           \usepackage{tikz}
           \usepackage{tgpagella}
           \setmainlanguage{english}
           \setotherlanguages{sanskrit}
           \newfontfamily\devanagarifont[Script=Devanagari]{Sanskrit 2003}
           \setotherlanguage{grantha}
           \newfontfamily\granthafont[Script=Grantha]{Noto Serif Grantha}
           \setotherlanguage{kannada}
           \newfontfamily\kannadafont[Script=Kannada]{Noto Serif Kannada}
           \setotherlanguage{malayalam}
           \newfontfamily\malayalamfont[Script=Tamil]{Noto Serif Malayalam}
           \setotherlanguage{tamil}
           \newfontfamily\tamilfont[Script=Tamil]{Noto Serif Tamil}
           \setotherlanguage{telugu}
           \newfontfamily\telugufont[Script=Telugu]{Noto Serif Telugu} 
           \XeTeXgenerateactualtext 1
           """


def make_circle_plot(
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
    language="Devanagari",
    nakṣatra_names_file="nakshatra_names.tex",
    rāśi_names_file="rashi_names.tex",
    graha_names_file="graha_names.tex",
    plotfile="panchanga_at_test_time.pdf",
):
    ## Sanskrit typesetting using XeLaTeX ##
    mpl.use("pgf")

    params = {
        "font.family": "serif",
        "text.usetex": True,
        "pgf.rcfonts": False,
        "pgf.texsystem": "xelatex",
        "pgf.preamble": preamble,
    }

    mpl.rcParams.update(params)

    if language in ["Tamil", "Grantha", "Malayalam"]:
        tamil_font_size_param = {"font.size": 6}
        mpl.rcParams.update(tamil_font_size_param)

    plt.style.use("dark_background")

    observing_location = EarthLocation.of_address(location)

    fig, ax = plt.subplots(subplot_kw={"projection": "polar"}, figsize=(7, 5))

    ## Setup circle plot with nakṣatra names and rāśi names

    nakṣatra_edges = np.linspace(ayanāṃśa, ayanāṃśa + 26 * nakṣatra_extent, 27) % 360
    nakṣatra_centers = nakṣatra_edges + nakṣatra_extent / 2
    nakṣatra_edge_labels = []

    nakṣatra_names_tab = Table.read(nakṣatra_names_file, format="latex")
    nakṣatra_names = nakṣatra_names_tab["names"].data

    for i in range(len(nakṣatra_edges)):
        coord = nakṣatra_edges[i]
        if coord > 360:
            coord = nakṣatra_edges[i] - 360
        deg_coord = int(coord)
        min_coord = 60 * (coord % 1)
        nakṣatra_edge_labels.append(
            r"%d$^{{\circ}}$\,%d$^{{\prime}}$" % (deg_coord, min_coord)
        )
        ax.plot(np.deg2rad([coord, coord]), [1, 2.5], "grey", alpha=0.3)

        rotation_val = nakṣatra_centers[i]
        color_val = "grey"
        alpha_val = 0.3

        if 90 < nakṣatra_centers[i] < 265:
            rotation_val = nakṣatra_centers[i] + 180

        if nakṣatra_names[i] == nakṣatra:
            color_val = "white"
            alpha_val = 1.0

        ax.text(
            np.deg2rad(nakṣatra_centers[i]),
            1.5,
            translit_str(nakṣatra_names[i], language),
            rotation=rotation_val,
            alpha=alpha_val,
            color=color_val,
            ha="center",
            va="center",
        )

    rāśi_edges = np.linspace(ayanāṃśa, ayanāṃśa + (rāśi_extent * 11), 12)
    rāśi_centers = rāśi_edges + rāśi_extent / 2
    rāśi_names_tab = Table.read(rāśi_names_file, format="latex")
    rāśi_names = rāśi_names_tab["names"].data

    for i in range(len(rāśi_edges)):
        ax.plot(np.deg2rad([rāśi_edges[i], rāśi_edges[i]]), [0, 1], "grey", alpha=0.3)

        rotation_val = rāśi_centers[i]
        color_val = "grey"
        alpha_val = 0.3

        if 90 < rāśi_centers[i] < 265:
            rotation_val = rāśi_centers[i] + 180

        if (
            rāśi_edges[np.floor((sun_lambda - ayanāṃśa) / 30).astype(int)]
            == rāśi_edges[i]
        ):
            color_val = "white"
            alpha_val = 1
        ax.text(
            np.deg2rad(rāśi_centers[i]),
            0.65,
            translit_str(rāśi_names[i], language),
            rotation=rotation_val,
            alpha=alpha_val,
            color=color_val,
            ha="center",
            va="center",
        )

    ### Plotting grahas
    ### Plot Sun
    plot_sun(np.array([np.deg2rad(sun_lambda), 0.9]), 0.08, fig, ax)

    ## Plot moon with phase ##
    plot_moon_phase(tithi, np.array([np.deg2rad(moon_lambda), 1.9]), 0.08, fig, ax)

    ### Inner grahas ###

    mercury_lambda = get_body(
        "mercury", test_date_utc_time, location=observing_location
    ).geocentrictrueecliptic.lon.value
    venus_lambda = get_body(
        "venus", test_date_utc_time, location=observing_location
    ).geocentrictrueecliptic.lon.value

    if sun_lambda > mercury_lambda:
        mercury_angle_to_sun = 360 - (sun_lambda - mercury_lambda)
    else:
        mercury_angle_to_sun = mercury_lambda - sun_lambda

    if sun_lambda > venus_lambda:
        venus_angle_to_sun = 360 - (sun_lambda - venus_lambda)
    else:
        venus_angle_to_sun = venus_lambda - sun_lambda

    plot_inner_graha_phase(
        "budha", mercury_angle_to_sun, [np.deg2rad(mercury_lambda), 1.9], 0.05, fig, ax
    )
    plot_inner_graha_phase(
        "shukra", venus_angle_to_sun, [np.deg2rad(venus_lambda), 1.9], 0.05, fig, ax
    )

    ### Outer grahas ###

    mars_lambda = get_body(
        "mars", test_date_utc_time, location=observing_location
    ).geocentrictrueecliptic.lon.value
    jupiter_lambda = get_body(
        "jupiter", test_date_utc_time, location=observing_location
    ).geocentrictrueecliptic.lon.value
    saturn_lambda = get_body(
        "saturn", test_date_utc_time, location=observing_location
    ).geocentrictrueecliptic.lon.value

    plot_outer_graha("mangala", [np.deg2rad(mars_lambda), 1.9], 0.05, fig, ax)
    plot_outer_graha("guru", [np.deg2rad(jupiter_lambda), 1.9], 0.05, fig, ax)
    plot_outer_graha("shani", [np.deg2rad(saturn_lambda), 1.9], 0.05, fig, ax)

    ##  rahu and ketu  ##

    rahu_lambda, ketu_lambda = calc_rahu_ketu_pos(
        test_date_utc_time, observing_location
    )

    rahu_lambda, ketu_lambda = rahu_lambda, ketu_lambda

    plot_rahu([np.deg2rad(rahu_lambda), 1.9], 5, fig, ax)
    plot_ketu([np.deg2rad(ketu_lambda), 1.9], 5, fig, ax)

    ## lagna ##

    ax.text(
        np.deg2rad(lagna_lambda),
        0.9,
        translit_str("लग्न", language),
        color="yellow",
        rotation=lagna_lambda + 90,
        ha="center",
        va="center",
    )

    ################ legend for grahas #####################
    graha_legend_ax = fig.add_axes([1 - 0.15, 0.1, 0.1, 0.8], polar=False)
    graha_legend_ax.axis("off")
    graha_tab = Table.read(graha_names_file, format="latex")

    graha_names = graha_tab["names"].data

    for graha_name_id in range(len(graha_names[:-1])):
        graha_legend_ax.text(
            0.5,
            1 - (graha_name_id / 10) - 0.1,
            translit_str(graha_names[graha_name_id], language),
            va="center",
        )

    plot_sun((0.3, 0.9), 0.08, fig, graha_legend_ax)
    plot_moon_phase(tithi, (0.3, 0.8), 0.08, fig, graha_legend_ax)

    plot_inner_graha_phase(
        "budha", mercury_angle_to_sun, (0.3, 0.7), 0.05, fig, graha_legend_ax
    )
    plot_inner_graha_phase(
        "shukra", venus_angle_to_sun, (0.3, 0.6), 0.05, fig, graha_legend_ax
    )
    plot_outer_graha("mangala", (0.3, 0.5), 0.05, fig, graha_legend_ax)
    plot_outer_graha("guru", (0.3, 0.4), 0.05, fig, graha_legend_ax)
    plot_outer_graha("shani", (0.3, 0.3), 0.05, fig, graha_legend_ax)

    plot_rahu((0.3, 0.2), 5, fig, graha_legend_ax)
    plot_ketu((0.3, 0.1), 5, fig, graha_legend_ax)

    ################ legend for panchanga #####################
    pañcāṅga_ax = fig.add_axes([0.05, 0.1, 0.18, 0.8], polar=False)
    pañcāṅga_ax.axis("off")

    pañcāṅga_ax.text(
        0.0,
        0.8,
        location,
        va="center",
    )

    pañcāṅga_ax.text(0.0, 0.75, date_str.split(" ")[0], va="center")

    pañcāṅga_ax.text(0.0, 0.7, date_str.split(" ")[1], va="center")
    pañcāṅga_ax.text(
        0.2,
        0.5,
        translit_str("पञ्चाङ्गः ", language),
        bbox=dict(facecolor="none", edgecolor="white"),
        ha="center",
        va="center",
    )

    pañcāṅga_ax.text(0.0, 0.4, translit_str(tithi_name, language), va="center")

    pañcāṅga_ax.text(0.0, 0.35, translit_str(vāra, language), va="center")

    pañcāṅga_ax.text(
        0.0,
        0.3,
        translit_str(nakṣatra, language)
        + r"\hspace{5pt}"
        + translit_str(pāda, language),
        va="center",
    )

    pañcāṅga_ax.text(0.0, 0.25, translit_str(yoga, language), va="center")
    pañcāṅga_ax.text(0.0, 0.2, translit_str(karaṇa, language), va="center")

    ax.set_title(translit_str("ॐ ", language), fontsize=20)

    # ---------------------------------------

    ax.set_xticks(np.deg2rad(nakṣatra_edges))
    ax.set_xticklabels(nakṣatra_edge_labels, fontsize=5)
    ax.set_rmax(2)
    ax.set_yticks([1])
    ax.set_yticklabels([""])
    ax.grid(axis="x")

    plt.savefig(plotfile)

    return 0


def make_jatakam_plot(
    test_date_utc_time,
    location,
    date_str,
    rāśi_extent,
    ayanāṃśa,
    tithi_name,
    vāra,
    nakṣatra,
    pāda,
    yoga,
    karaṇa,
    moon_lambda,
    sun_lambda,
    lagna_lambda,
    language="Devanagari",
    rāśi_names_file="rashi_names.tex",
    plotfile="jatakam_at_test_time.pdf",
):
    plt.style.use("default")
    ## Sanskrit typesetting using XeLaTeX ##
    mpl.use("pgf")

    params = {
        "font.family": "serif",
        "text.usetex": True,
        "pgf.rcfonts": False,
        "pgf.texsystem": "xelatex",
        "pgf.preamble": preamble,
    }

    mpl.rcParams.update(params)

    if language in ["Tamil", "Grantha", "Malayalam"]:
        tamil_font_size_param = {"font.size": 6}
        mpl.rcParams.update(tamil_font_size_param)

    fig, ax = plt.subplots()
    divider = make_axes_locatable(ax)
    pañcāṅga_ax = divider.append_axes("left", 0.5, sharey=ax)

    fig.suptitle(
        translit_str(r" ॐ ", language),
        color="b",
    )

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
    ### Inner grahas ###
    observing_location = EarthLocation.of_address(location)
    mercury_lambda = get_body(
        "mercury", test_date_utc_time, location=observing_location
    ).geocentrictrueecliptic.lon.value
    venus_lambda = get_body(
        "venus", test_date_utc_time, location=observing_location
    ).geocentrictrueecliptic.lon.value

    ### Outer grahas ###

    mars_lambda = get_body(
        "mars", test_date_utc_time, location=observing_location
    ).geocentrictrueecliptic.lon.value
    jupiter_lambda = get_body(
        "jupiter", test_date_utc_time, location=observing_location
    ).geocentrictrueecliptic.lon.value
    saturn_lambda = get_body(
        "saturn", test_date_utc_time, location=observing_location
    ).geocentrictrueecliptic.lon.value

    ##  rahu and ketu  ##

    rahu_lambda, ketu_lambda = calc_rahu_ketu_pos(
        test_date_utc_time, observing_location
    )

    rahu_lambda, ketu_lambda = rahu_lambda, ketu_lambda

    lst_of_grahas_pos = [
        sun_lambda,
        moon_lambda,
        mercury_lambda,
        venus_lambda,
        mars_lambda,
        jupiter_lambda,
        saturn_lambda,
        rahu_lambda,
        ketu_lambda,
        lagna_lambda,
    ]
    graha_tab = Table.read("graha_names.tex", format="latex")

    graha_names = graha_tab["names"].data

    texts = []

    for graha_id in range(len(graha_names)):
        sel_rāśi_id = np.floor(
            (lst_of_grahas_pos[graha_id] - ayanāṃśa) / rāśi_extent
        ).astype(int)
        t = ax.text(
            coord_for_rāśi[sel_rāśi_id][0],
            coord_for_rāśi[sel_rāśi_id][1],
            translit_str(graha_names[graha_id], language),
            color="b",
            ha="center",
            # va="center",
        )
        texts.append(t)

    adjust_text(texts, ax=ax)

    ################ legend for panchanga #####################
    pañcāṅga_ax.text(
        0.2,
        7.5,
        location,
        color="b",
        va="center",
    )

    pañcāṅga_ax.text(
        0.2,
        7.0,
        date_str.split(" ")[0],
        color="b",
        va="center",
    )

    pañcāṅga_ax.text(
        0.2,
        6.5,
        date_str.split(" ")[1],
        color="b",
        va="center",
    )

    pañcāṅga_ax.text(
        1,
        4.0,
        translit_str("पञ्चाङ्ः", language),
        bbox=dict(facecolor="none", edgecolor="b"),
        color="b",
        va="center",
        ha="center",
    )

    pañcāṅga_ax.text(
        0.2,
        2.5,
        translit_str(tithi_name, language),
        color="b",
        va="center",
    )

    pañcāṅga_ax.text(
        0.2,
        2.0,
        translit_str(vāra, language),
        color="b",
        va="center",
    )

    pañcāṅga_ax.text(
        0.2,
        1.5,
        translit_str(nakṣatra, language)
        + r"\hspace{5pt}"
        + translit_str(pāda, language),
        color="b",
        va="center",
    )

    pañcāṅga_ax.text(
        0.2,
        1.0,
        translit_str(yoga, language),
        color="b",
        va="center",
    )
    pañcāṅga_ax.text(
        0.2,
        0.5,
        translit_str(karaṇa, language),
        color="b",
        va="center",
    )

    ax.set_aspect("equal")
    # pañcāṅga_ax.set_aspect("equal")
    pañcāṅga_ax.axis("off")
    ax.axis("off")
    plt.savefig(plotfile, bbox_inches="tight")

    return 0


def make_sky_plot(
    test_date_utc_time,
    location,
    date_str,
    nakṣatra_extent,
    rāśi_extent,
    ayanāṃśa,
    tithi,
    tithi_name,
    language,
    nakṣatra_names_file="nakshatra_names.tex",
    rāśi_names_file="rashi_names.tex",
):
    ## Sanskrit typesetting using XeLaTeX ##
    mpl.use("pgf")

    params = {
        "font.family": "serif",
        "text.usetex": True,
        "pgf.rcfonts": False,
        "pgf.texsystem": "xelatex",
        "pgf.preamble": preamble,
    }

    mpl.rcParams.update(params)

    if language in ["Tamil", "Grantha", "Malayalam"]:
        tamil_font_size_param = {"font.size": 6}
        mpl.rcParams.update(tamil_font_size_param)

    plt.style.use("dark_background")

    observing_location = EarthLocation.of_address(location)

    fig = plt.figure(facecolor="midnightblue", figsize=(12, 5))

    fig.suptitle(
        translit_str(r" ॐ ", language),
    )
    ax_above = fig.add_subplot(121, projection="polar")
    ax_above.set_facecolor("midnightblue")
    ax_above.set_theta_zero_location("N")
    ax_above.set_theta_direction(-1)
    ax_above.grid(True, which="major", linestyle="dotted")

    degree_sign = r"$^{\circ}$"

    # For positively-increasing range (e.g., range(1, 90, 15)),
    # labels go from middle to outside.
    r_labels = [
        "90" + degree_sign,
        "",
        "60" + degree_sign,
        "",
        "30" + degree_sign,
        "",
        "0" + degree_sign + " Alt.",
    ]

    theta_labels = []
    for chunk in range(0, 7):
        label_angle = (0.0 * u.deg * (1 / u.deg)) + (chunk * 45.0)
        while label_angle >= 360.0:
            label_angle -= 360.0
        if chunk == 0:
            theta_labels.append("N " + "\n" + str(label_angle) + degree_sign + " Az")
        elif chunk == 2:
            theta_labels.append("E" + "\n" + str(label_angle) + degree_sign)
        elif chunk == 4:
            theta_labels.append("S" + "\n" + str(label_angle) + degree_sign)
        elif chunk == 6:
            theta_labels.append("W" + "\n" + str(label_angle) + degree_sign)
        else:
            theta_labels.append(str(label_angle) + degree_sign)
    theta_labels.append("")

    ax_below = fig.add_subplot(122, projection="polar")
    ax_below.set_facecolor("black")
    ax_below.set_theta_zero_location("N")
    ax_below.grid(True, which="major", linestyle="dotted")

    # For positively-increasing range (e.g., range(1, 90, 15)),
    # labels go from middle to outside.
    r_labels = [
        "90" + degree_sign,
        "",
        "60" + degree_sign,
        "",
        "30" + degree_sign,
        "",
        "0" + degree_sign + " Alt.",
    ]

    theta_labels = []
    for chunk in range(0, 7):
        label_angle = (0.0 * u.deg * (1 / u.deg)) + (chunk * 45.0)
        while label_angle >= 360.0:
            label_angle -= 360.0
        if chunk == 0:
            theta_labels.append("N " + "\n" + str(label_angle) + degree_sign + " Az")
        elif chunk == 2:
            theta_labels.append("E" + "\n" + str(label_angle) + degree_sign)
        elif chunk == 4:
            theta_labels.append("S" + "\n" + str(label_angle) + degree_sign)
        elif chunk == 6:
            theta_labels.append("W" + "\n" + str(label_angle) + degree_sign)
        else:
            theta_labels.append(str(label_angle) + degree_sign)
    theta_labels.append("")

    aa_frame = coor.AltAz(obstime=test_date_utc_time, location=observing_location)
    # Plot ecliptic
    """
    coord = SkyCoord(
        lat=np.repeat(0, 100) * u.deg,
        lon=np.linspace(0, 360, 100) * u.deg,
        frame="geocentrictrueecliptic",
    )
    
    az = coord.transform_to(aa_frame).az.rad
    alt = coord.transform_to(aa_frame).alt.deg

    sel_above = np.where(alt > 0)

    ax_above.plot(
        az[sel_above],
        90 - alt[sel_above],
        color="w",
        linestyle="none",
        marker="o",
        markersize=2,
        zorder=0,
    )

    sel_below = np.where(alt < 0)
    ax_below.plot(
        az[sel_below],
        90 + alt[sel_below],
        color="w",
        linestyle="none",
        marker="o",
        markersize=2,
        zorder=0,
    )
    """
    # Plot nakṣatra grid
    nakṣatra_edges = np.linspace(ayanāṃśa, ayanāṃśa + 26 * nakṣatra_extent, 27) % 360
    nakṣatra_centers = nakṣatra_edges + nakṣatra_extent / 2
    nakṣatra_names_tab = Table.read(nakṣatra_names_file, format="latex")
    nakṣatra_names = nakṣatra_names_tab["names"].data

    for edge in nakṣatra_edges:
        edge_coord = SkyCoord(
            lat=np.linspace(0, 30, 50) * u.deg,
            lon=np.repeat(edge, 50) * u.deg,
            frame="geocentrictrueecliptic",
        )
        edge_az = edge_coord.transform_to(aa_frame).az.rad
        edge_alt = edge_coord.transform_to(aa_frame).alt.deg

        edge_sel_above = np.where(edge_alt > 0)

        ax_above.plot(
            edge_az[edge_sel_above],
            90 - edge_alt[edge_sel_above],
            color="w",
            linewidth=1,
            zorder=0,
            alpha=0.5,
        )

        edge_sel_below = np.where(edge_alt < 0)
        ax_below.plot(
            edge_az[edge_sel_below],
            90 + edge_alt[edge_sel_below],
            color="w",
            linewidth=1,
            zorder=0,
            alpha=0.5,
        )

    for id, center in enumerate(nakṣatra_centers):
        center_coord = SkyCoord(
            lat=10 * u.deg, lon=center * u.deg, frame="geocentrictrueecliptic"
        )

        center_az = center_coord.transform_to(aa_frame).az.rad
        center_alt = center_coord.transform_to(aa_frame).alt.deg

        rot_val = 90

        if center_alt > 0:
            ax_above.text(
                center_az,
                90 - center_alt,
                translit_str(nakṣatra_names[id], language),
                rotation=rot_val,
                ha="center",
                va="center",
                fontsize=7,
                alpha=0.5,
            )

        else:
            ax_below.text(
                center_az,
                90 + center_alt,
                translit_str(nakṣatra_names[id], language),
                rotation=rot_val,
                ha="center",
                va="center",
                fontsize=7,
                alpha=0.5,
            )

    # Plot rāśī grid
    rāśi_edges = np.linspace(ayanāṃśa, ayanāṃśa + (rāśi_extent * 11), 12)
    rāśi_centers = rāśi_edges + rāśi_extent / 2
    rāśi_names_tab = Table.read(rāśi_names_file, format="latex")
    rāśi_names = rāśi_names_tab["names"].data

    for edge in rāśi_edges:
        edge_coord = SkyCoord(
            lat=np.linspace(-30, 0, 50) * u.deg,
            lon=np.repeat(edge, 50) * u.deg,
            frame="geocentrictrueecliptic",
        )
        edge_az = edge_coord.transform_to(aa_frame).az.rad
        edge_alt = edge_coord.transform_to(aa_frame).alt.deg

        edge_sel_above = np.where(edge_alt > 0)

        ax_above.plot(
            edge_az[edge_sel_above],
            90 - edge_alt[edge_sel_above],
            color="w",
            linewidth=1,
            zorder=0,
            alpha=0.5,
        )

        edge_sel_below = np.where(edge_alt < 0)
        ax_below.plot(
            edge_az[edge_sel_below],
            90 + edge_alt[edge_sel_below],
            color="w",
            linewidth=1,
            zorder=0,
            alpha=0.5,
        )

    for id, center in enumerate(rāśi_centers):
        center_coord = SkyCoord(
            lat=-10 * u.deg, lon=center * u.deg, frame="geocentrictrueecliptic"
        )

        center_az = center_coord.transform_to(aa_frame).az.rad
        center_alt = center_coord.transform_to(aa_frame).alt.deg

        rot_val = 90

        if center_alt > 0:
            ax_above.text(
                center_az,
                90 - center_alt,
                translit_str(rāśi_names[id], language),
                rotation=rot_val,
                ha="center",
                va="center",
                fontsize=7,
                alpha=0.5,
            )

        else:
            ax_below.text(
                center_az,
                90 + center_alt,
                translit_str(rāśi_names[id], language),
                rotation=rot_val,
                ha="center",
                va="center",
                fontsize=7,
                alpha=0.5,
            )

    # Plot grahas
    sun = get_body("sun", test_date_utc_time, location=observing_location).transform_to(
        aa_frame
    )
    if sun.alt.deg > 0:
        plot_sun((sun.az.rad, 90 - sun.alt.deg), 0.08, fig, ax_above)
    else:
        plot_sun((sun.az.rad, 90 + sun.alt.deg), 0.08, fig, ax_below)

    moon = get_body(
        "moon", test_date_utc_time, location=observing_location
    ).transform_to(aa_frame)

    if moon.alt.deg > 0:
        plot_moon_phase(
            30 - tithi, (moon.az.rad, 90 - moon.alt.deg), 0.08, fig, ax_above
        )
    else:
        plot_moon_phase(
            30 - tithi, (moon.az.rad, 90 + moon.alt.deg), 0.08, fig, ax_below
        )
    ### Inner grahas ###
    sun = get_body("sun", test_date_utc_time, location=observing_location)

    mercury = get_body("mercury", test_date_utc_time, location=observing_location)
    venus = get_body("venus", test_date_utc_time, location=observing_location)

    if sun.geocentrictrueecliptic.lon.deg > mercury.geocentrictrueecliptic.lon.deg:
        mercury_angle_to_sun = 360 - (
            sun.geocentrictrueecliptic.lon.deg - mercury.geocentrictrueecliptic.lon.deg
        )
    else:
        mercury_angle_to_sun = (
            mercury.geocentrictrueecliptic.lon.deg - sun.geocentrictrueecliptic.lon.deg
        )

    if sun.geocentrictrueecliptic.lon.deg > venus.geocentrictrueecliptic.lon.deg:
        venus_angle_to_sun = 360 - (
            sun.geocentrictrueecliptic.lon.deg - venus.geocentrictrueecliptic.lon.deg
        )
    else:
        venus_angle_to_sun = (
            venus.geocentrictrueecliptic.lon.deg - sun.geocentrictrueecliptic.lon.deg
        )

    if mercury.transform_to(aa_frame).alt.deg > 0:
        plot_inner_graha_phase(
            "budha",
            360 - mercury_angle_to_sun,
            [
                mercury.transform_to(aa_frame).az.rad,
                90 - mercury.transform_to(aa_frame).alt.deg,
            ],
            0.05,
            fig,
            ax_above,
        )
    else:
        plot_inner_graha_phase(
            "budha",
            360 - mercury_angle_to_sun,
            [
                mercury.transform_to(aa_frame).az.rad,
                90 + mercury.transform_to(aa_frame).alt.deg,
            ],
            0.05,
            fig,
            ax_below,
        )

    if venus.transform_to(aa_frame).alt.deg > 0:
        plot_inner_graha_phase(
            "shukra",
            360 - venus_angle_to_sun,
            [
                venus.transform_to(aa_frame).az.rad,
                90 - venus.transform_to(aa_frame).alt.deg,
            ],
            0.05,
            fig,
            ax_above,
        )
    else:
        plot_inner_graha_phase(
            "shukra",
            360 - venus_angle_to_sun,
            [
                venus.transform_to(aa_frame).az.rad,
                90 + venus.transform_to(aa_frame).alt.deg,
            ],
            0.05,
            fig,
            ax_below,
        )

    # Outer grahas

    mars = get_body(
        "mars", test_date_utc_time, location=observing_location
    ).transform_to(aa_frame)
    jupiter = get_body(
        "jupiter", test_date_utc_time, location=observing_location
    ).transform_to(aa_frame)
    saturn = get_body(
        "saturn", test_date_utc_time, location=observing_location
    ).transform_to(aa_frame)

    if mars.alt.deg > 0:
        plot_outer_graha(
            "mangala", [mars.az.rad, 90 - mars.alt.deg], 0.05, fig, ax_above
        )
    else:
        plot_outer_graha(
            "mangala", [mars.az.rad, 90 + mars.alt.deg], 0.05, fig, ax_below
        )

    if jupiter.alt.deg > 0:
        plot_outer_graha(
            "guru", [jupiter.az.rad, 90 - jupiter.alt.deg], 0.05, fig, ax_above
        )
    else:
        plot_outer_graha(
            "guru", [jupiter.az.rad, 90 + jupiter.alt.deg], 0.05, fig, ax_below
        )

    if saturn.alt.deg > 0:
        plot_outer_graha(
            "shani", [saturn.az.rad, 90 - saturn.alt.deg], 0.05, fig, ax_above
        )
    else:
        plot_outer_graha(
            "shani", [saturn.az.rad, 90 + saturn.alt.deg], 0.05, fig, ax_below
        )
    ax_above.plot(0, 0, marker="$\star$", color="white")
    # Set ticks and labels.

    alpha_val = 0.2

    ax_above.set_rgrids(range(1, 106, 15), r_labels, angle=-45, alpha=alpha_val)
    ax_above.set_thetagrids(range(0, 360, 45), theta_labels, alpha=alpha_val)
    ax_above.grid(alpha=alpha_val)
    ax_above.set_title("Above horizon")

    ax_below.set_rgrids(range(1, 106, 15), r_labels, angle=-45, alpha=alpha_val)
    ax_below.set_thetagrids(range(0, 360, 45), theta_labels, alpha=alpha_val)
    ax_below.grid(alpha=alpha_val)
    ax_below.set_title("Below horizon")

    plt.savefig("sky_plot_at_test_time.pdf", bbox_inches="tight")

    return fig, ax_above, ax_below
