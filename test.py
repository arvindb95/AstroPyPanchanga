import matplotlib.pyplot as plt
import matplotlib as mpl

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
           \setotherlanguage{kannada}
           \newfontfamily\kannadafont[Script=Kannada]{Gubbi}
           \newcommand{\sam}[1]{\begin{sanskrit}#1\end{sanskrit}}
           \newcommand{\kan}[1]{\begin{kannada}#1\end{kannada}}"""

params = {
    "font.family": "serif",
    "text.usetex": True,
    "pgf.rcfonts": False,
    "pgf.texsystem": "xelatex",
    "pgf.preamble": preamble,
}
mpl.rcParams.update(params)

fig, ax = plt.subplots()

ax.text(0.5, 0.5, r"\kan{ॐ }  ", fontsize=20)

plt.savefig("test.pdf")
