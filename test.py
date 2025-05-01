import matplotlib.pyplot as plt
import matplotlib as mpl
from indic_transliteration import sanscript
from indic_transliteration.sanscript import SchemeMap, SCHEMES, transliterate

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

set_language = "kannada"

text = "ॐ"

print(text)

final_text = transliterate(text, sanscript.DEVANAGARI, sanscript.KANNADA)

ax.text(
    0.5,
    0.5,
    r"\begin{%s} %s \end{%s}" % (set_language, final_text, set_language),
    fontsize=20,
)

plt.savefig("test.pdf")
