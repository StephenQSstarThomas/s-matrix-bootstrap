"""LaTeX version of the comparison document for the authors of arXiv:2309.12402; compiled with tectonic.

Usage: python scripts/sdp/authors_comparison_tex.py --root RESULTS_ROOT --docdir docs/reproduction_2309
Writes REPRODUCTION_COMPARISON_EN.tex next to the figures/ folder and compiles REPRODUCTION_COMPARISON_EN.pdf.
"""
from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from authors_comparison_en import gather  # noqa: E402


def f(x, d=4):
    """Number for LaTeX: scientific notation below 1e-2, math-mode minus, '--' for missing."""
    if x is None:
        return "--"
    ax = abs(x)
    if ax != 0 and ax < 1e-2:
        m, e = f"{x:.{max(d-1,2)}e}".split("e")
        return f"\\ensuremath{{{float(m):.{max(d-2,2)}g}\\times10^{{{int(e)}}}}}"
    t = f"{x:.{d}g}"
    return f"\\ensuremath{{{t}}}" if t.startswith("-") else t


def pct(a, b):
    return f"{(a / b - 1) * 100:+.1f}\\,\\%"


def esc(s):
    return s.replace("%", "\\%").replace("_", "\\_").replace("&", "\\&")


def main(argv=None):
    p = argparse.ArgumentParser(); p.add_argument("--root", required=True); p.add_argument("--docdir", required=True)
    a = p.parse_args(argv); g = gather(a.root); D = Path(a.docdir)
    c1, c2, c4, c5, c6, c7, c8, face = g["c1rows"], g["c2"], g["c4rows"], g["c5"], g["c6"], g["c7"], g["c8rows"], g["face"]
    tip, mreg, srmom = g["tip"], g["mreg"], g["srmom"]
    ra, rf, rd = c5.get("SRa", {}), c5.get("SRa_ffs0", {}), c5.get("SRd_ffs0", {})
    def rng(name, fn):
        row = face.get(name, {}).get(fn)
        if not row: return "--"
        lo, hi = row["range_lo"], row["range_hi"]
        return f"$[{f(lo,4)},\\,{f(hi,4)}]$" if lo is not None else f"$\\le {f(hi,4)}$"
    s2min, s2max = min(r["rms20_deg"] for r in c7.values()), max(r["rms20_deg"] for r in c7.values())
    s0min, s0max = min(r["rms00_deg"] for r in c7.values()), max(r["rms00_deg"] for r in c7.values())
    ladder = ", ".join(f"{v:.4f}" for v in c2.get("monotone_in_eps", {}).get("ends_desc_eps", []))
    mreg_x = pct(mreg["verification"]["f00_3"], tip["verification"]["f00_3"]) if mreg and tip else "--"
    mreg_dm = f(1000 * (mreg["observables"]["P1"]["crossing_90_GeV"] - tip["observables"]["P1"]["crossing_90_GeV"]), 2) if mreg and tip else "--"
    mreg_de = f(mreg["observables"]["P1"]["min_eta_below_1p2GeV"] - tip["observables"]["P1"]["min_eta_below_1p2GeV"], 2) if mreg and tip else "--"
    srm = f(srmom["verification"]["functional"]["value"], 5) if srmom else "--"

    tex = r"""\documentclass[11pt,a4paper]{article}
\usepackage[margin=21mm,top=20mm,bottom=22mm]{geometry}
\usepackage{lmodern}
\usepackage[T1]{fontenc}
\usepackage{microtype}
\usepackage{amsmath,amssymb}
\usepackage{booktabs,tabularx,array,longtable}
\usepackage{graphicx,xcolor,float}
\usepackage[font=small,labelformat=empty,skip=4pt]{caption}
\usepackage{enumitem}
\usepackage[hidelinks]{hyperref}
\usepackage{titlesec}
\definecolor{ink}{HTML}{1B2430}\definecolor{quote}{HTML}{5B6673}\definecolor{rule}{HTML}{9AA4AF}
\titleformat{\section}{\large\bfseries\color{ink}}{\thesection}{0.7em}{}
\titleformat{\subsection}{\normalsize\bfseries\color{ink}}{\thesubsection}{0.6em}{}
\titlespacing*{\section}{0pt}{16pt}{6pt}\titlespacing*{\subsection}{0pt}{12pt}{4pt}
\newcommand{\paper}[1]{{\color{quote}\itshape #1}}
\newcolumntype{Y}{>{\raggedright\arraybackslash}X}
\setlength{\parskip}{4pt}\setlength{\parindent}{0pt}\setlength{\emergencystretch}{2.5em}
\renewcommand{\arraystretch}{1.15}
\begin{document}
{\LARGE\bfseries Reproduction of arXiv:2309.12402 with SDPB}\\[4pt]
{\large Comparison against the paper's figures}\\[6pt]
{\color{quote}Prepared for Yifei He and Martin Kruczenski. 14 September 2026.\\
Code and receipts: branch \texttt{sdpb-2309-regularised}, folder \texttt{docs/reproduction\_2309}, of\\ \url{https://github.com/StephenQSstarThomas/s-matrix-bootstrap}}
\vspace{6pt}\hrule\vspace{8pt}

\section{What was done}
We reproduced the finite S-matrix / form-factor bootstrap of \emph{Bootstrapping gauge theories} (arXiv:2309.12402 v3) from the equations of Sections~2 and~3 of the paper and compared the result with each figure. The optimiser is SDPB 3.1.0 (192-bit arithmetic, duality gap $10^{-6}$, primal and dual error $10^{-10}$); no other solver was used. The discretised operators---the conformal map and grid (3.58)--(3.62), the cot kernel (3.67), the angular projections of (2.7), the current kernels (3.68)--(3.75)---were derived from the paper and checked independently: a Mathematica script recomputes the formulas at $M=50$; an independent reconstruction of the (2.7) projections agrees with the production operators to $5\times10^{-15}$ over all 3876 coefficients; every accepted solution is re-verified against all original constraints in Arb interval arithmetic. The comparison thresholds for each figure were fixed before the corresponding runs (files \texttt{receipts/GATE\_PREREGISTRATION.md}, \texttt{receipts/GATE\_LOG.md}).

The paper leaves five numerical items unstated. They decide the outcome for Figures~8--10, so we list them with the readings we ran.

\begin{tabularx}{\textwidth}{@{}>{\raggedright\arraybackslash}p{0.27\textwidth}>{\raggedright\arraybackslash}p{0.27\textwidth}Y@{}}
\toprule
Item & Paper text & Readings run \\
\midrule
Norm of the sum-rule tolerance (3.73), $\epsilon^{SR}=2\times10^{-3}$ & \paper{with some norm} & per-moment box on the raw moments (main line); per-wave $L^2$ ball (the structure of your 2403 code); relative 10\,\% per moment \\
Evaluation point of the factor between $F$ and $\mathcal F$ in (3.75) & \paper{which we evaluate at $s=s_0$} (in the estimate of $\epsilon^{FF}$) & factor at each node $s_i>s_0$ (main line); factor frozen at $s_0$ \\
Regularisation of the double spectral density & not mentioned in 2309; 2103.11484 Section~3 and your 2403 code use an $M$-bound & $|\rho_{ij}|\le M_{\rm reg}$ with $M_{\rm reg}$ fixed by an in-model rule (omitted-wave unitarity and $L$-stability): $10^2$ for the chiral and UV stages, $10^3$ for the pure stage; $\ell^2$, $\ell^4$ and $M_{\rm reg}\times10$ controls \\
Amplitude at a boundary point & \paper{only points at the boundary have partial waves associated with them} & the solver's optimal point, plus a diagnostic of the whole near-optimal face (Section~4) \\
Chiral norm in (3.64) & \paper{with some norm} & one combined 8-dimensional $L^2$ norm (as in your 2403 code); two separate 4-dimensional norms as a control \\
\bottomrule
\end{tabularx}

Everything else is taken as printed: $M=50$, $\phi_i=(i-\tfrac12)\pi/M$, $\nu_0=0$, $L=10$ waves per isospin, $s_0=(1.2\,\text{GeV})^2$, $\alpha_s=0.4$, $m_u=4$\,MeV, $m_d=7.3$\,MeV, the condensates (2.54), $m_\pi=140$\,MeV, $f_\pi=92$\,MeV, the printed sum-rule numbers (2.56) times $s_0^{\,n+2}$, $\epsilon^\chi=2\times10^{-3}$, $\epsilon^{FF}=6\times10^{-5}$, chiral points $s=\tfrac12,1,\tfrac32,2$, moments $n=0,1$ (S0) and $n=-1,0$ (P1).

\section{Summary}
\begin{tabularx}{\textwidth}{@{}l>{\raggedright\arraybackslash}p{0.24\textwidth}Y>{\raggedright\arraybackslash}p{0.22\textwidth}@{}}
\toprule
Figure & Paper statement & Our result & Agreement \\
\midrule
"""
    tex += (f"Fig.~3 & pure-unitarity region, $M=50$, $L=10$ & $+x$ end {f(c1['f00_max']['ours'],6)} vs {f(c1['f00_max']['paper'],6)} ({pct(c1['f00_max']['ours'],c1['f00_max']['paper'])}); $-x$ end {f(c1['f00_min']['ours'],5)} vs {f(c1['f00_min']['paper'],5)}; $f_1^1$ range $[{f(c1['f11_min']['ours'],4)},{f(c1['f11_max']['ours'],4)}]$ vs $[{f(c1['f11_min']['paper'],4)},{f(c1['f11_max']['paper'],4)}]$ & $+x$ end agrees; our region lies inside yours on the $-x$ and $\\pm y$ sides \\\\\n"
            f"Fig.~4 & chiral constraints collapse the region onto $f_1^1=-f_0^0/15$ & $\\epsilon^\\chi=2\\times10^{{-3}}$: $+x$ end {f(c2['x_end_eps002']['ours'],6)} vs {f(c2['x_end_eps002']['paper'],6)} ({pct(c2['x_end_eps002']['ours'],c2['x_end_eps002']['paper'])}); $x_{{\\rm ref}}$ section width {f(c2['xref_width_eps002']['width'],4)} vs {f(c2['xref_width_eps002']['paper'],4)}; six-tolerance ladder monotone & agrees \\\\\n"
            f"Fig.~5 & subthreshold waves nearly linear; the S0 chiral zero moves with $\\epsilon^\\chi$ & RMS $\\le 6.4\\,\\%$ of $f_0^0(3)$ against the digitised curves at $\\epsilon^\\chi=2,4,6\\times10^{{-3}}$; S0 zero at $s=0.426$, $0.293$, absent & agrees \\\\\n"
            f"Fig.~7 & chiral-only phases: S0, S2 agree with experiment, P1 has no $\\rho$ & RMS against the digitised curves S0 {f(c4['S0']['rms_deg'],3)}$^\\circ$, S2 {f(c4['S2']['rms_deg'],3)}$^\\circ$, P1 {f(c4['P1']['rms_deg'],3)}$^\\circ$; no P1 crossing below 1.2\\,GeV & agrees \\\\\n"
            f"Fig.~8 & with the sum rules the upper boundary shrinks much more than the lower & upper shrink / lower rise at $x_{{\\rm ref}}$: {f(ra['upper shrink at x_ref']['ours'],3)} / {f(ra['lower rise at x_ref']['ours'],3)} (node factor), {f(rf['upper shrink at x_ref']['ours'],3)} / {f(rf['lower rise at x_ref']['ours'],3)} (frozen factor); paper $2.32\\times10^{{-4}}$ / $3.7\\times10^{{-5}}$. $+x$ end {f(ra['UV +x end']['ours'],5)} / {f(rf['UV +x end']['ours'],5)} vs 0.08112 & the asymmetry does not appear in any feasible reading \\\\\n"
            f"Fig.~9 & P1 phase shift with the $\\rho$; 90$^\\circ$ crossing about 6\\,\\% above 770\\,MeV & crossings: tip {f(c6['tip']['crossing_MeV'],4)}\\,MeV ($|S|_{{\\min}}=0.585$), mid {f(c6['mid']['crossing_MeV'],4)}\\,MeV ($|S|_{{\\min}}=0.077$), ref: $|S|$ dips to 0.24 at 0.73\\,GeV (708\\,MeV on the $+180^\\circ$ branch); paper 813--827\\,MeV & a $\\rho$-like structure appears at every point once the UV constraints are added; it sits 40--110\\,MeV lower than yours and is strongly inelastic \\\\\n"
            f"Fig.~10 & S0 and S2 coincide at low energy and spread above & S2: RMS {f(s2min,3)}--{f(s2max,3)}$^\\circ$, endpoints $-28$ to $-32^\\circ$ (paper $-18$ to $-36^\\circ$); S0: RMS {f(s0min,3)}--{f(s0max,3)}$^\\circ$, endpoint at 1.196\\,GeV 181--192$^\\circ$ (paper 87--106$^\\circ$) & S2 agrees; S0 rises faster than yours above 0.6\\,GeV \\\\\n"
            f"Fig.~11 & good convergence in $L$; the $\\rho$ peak shifts with $M$ & P1 crossing at $M=50$, $L=8/10/12$: {f(c8[(50,8)]['P1_crossing_MeV'],4)} / {f(c8[(50,10)]['P1_crossing_MeV'],4)} / {f(c8[(50,12)]['P1_crossing_MeV'],4)}\\,MeV; $L=10$, $M=45/50/60$: {f(c8[(45,10)]['P1_crossing_MeV'],4)} / {f(c8[(50,10)]['P1_crossing_MeV'],4)} / {f(c8[(60,10)]['P1_crossing_MeV'],4)}\\,MeV; S0 at 1\\,GeV 150--165$^\\circ$ (paper about 95$^\\circ$) & $L$-stability agrees; the $M$ shift is 37\\,MeV; the S0 level differs \\\\\n")
    tex += r"""\bottomrule
\end{tabularx}

\section{Figure-by-figure comparison}
In each figure the left (or upper) panel is the paper's figure rendered from the arXiv source; the right (or lower) panel shows our accepted solutions with the digitised paper curve or boundary overlaid. Paper statements are quoted in \paper{grey italics}.

\subsection{Fig.~3, pure unitarity region}
\paper{``we have used $M=50$ in (3.61) for discretization and imposed unitarity for 10 partial waves per isospin.''}

"""
    tex += (f"Ours: 24 support directions at 15$^\\circ$ steps, all accepted, $M=50$, $L=10$, $M_{{\\rm reg}}=10^3$ (the six first omitted waves stay within $|S|\\le1.02$ at every node and the $L=8/10/12$ endpoints agree to 2\\,\\%). Extrema against the digitised figure: $f_0^0$ max {f(c1['f00_max']['ours'],6)} vs {f(c1['f00_max']['paper'],6)} ({pct(c1['f00_max']['ours'],c1['f00_max']['paper'])}), $f_0^0$ min {f(c1['f00_min']['ours'],5)} vs {f(c1['f00_min']['paper'],5)}, $f_1^1$ min {f(c1['f11_min']['ours'],4)} vs {f(c1['f11_min']['paper'],4)}, $f_1^1$ max {f(c1['f11_max']['ours'],4)} vs {f(c1['f11_max']['paper'],4)}. The $+x$ end is insensitive to $M_{{\\rm reg}}$ (2.2116 / 2.2245 / 2.2305 at $10^2$ / $10^3$ / $10^4$). The three other extrema lie inside your region by 15--23\\,\\%; these are the directions where the amplitudes are large and the $|\\rho_{{ij}}|$ bound is active. We read this as a measurement of the unstated regulariser rather than of the model; an $M_{{\\rm reg}}=10^4$ control for these three directions is running.\n\n"
            "\\begin{figure}[H]\\centering\\includegraphics[width=\\textwidth]{figures/fig3.png}\\caption{Fig.~3: paper (left) and our 24 support points with the digitised boundary (right).}\\end{figure}\n\n"
            "\\subsection{Fig.~4, chiral constraints}\n\\paper{``restricted by the chiral constraints (3.64) with tolerances $6\\times10^{-3}$, $4\\times10^{-3}$, $2\\times10^{-3}$, $1\\times10^{-3}$, $6\\times10^{-4}$, $2\\times10^{-4}$ (from the outer shape inward) \\dots\\ with some norm''}\n\n"
            f"Ours: one combined 8-dimensional $L^2$ norm on the four-point residuals (the packaging of your 2403 code), $M_{{\\rm reg}}=10^2$. At $\\epsilon^\\chi=2\\times10^{{-3}}$ the $+x$ end is {f(c2['x_end_eps002']['ours'],6)} (paper {f(c2['x_end_eps002']['paper'],6)}, {pct(c2['x_end_eps002']['ours'],c2['x_end_eps002']['paper'])}) and the vertical section at $x_{{\\rm ref}}$ has width {f(c2['xref_width_eps002']['width'],4)} (paper {f(c2['xref_width_eps002']['paper'],4)}). The $+x$ ends of the six tolerances decrease monotonically: {ladder} for $\\epsilon^\\chi=6\\times10^{{-3}}\\dots2\\times10^{{-4}}$. The two-norm control moves the $+x$ end by less than the ladder spacing.\n\n"
            "\\begin{figure}[H]\\centering\\includegraphics[width=\\textwidth]{figures/fig4.png}\\caption{Fig.~4: paper (left); our six $+x$ ends and the $\\epsilon^\\chi=2\\times10^{-3}$ sections with the digitised $\\epsilon^\\chi=2\\times10^{-3}$ boundary (right).}\\end{figure}\n\n"
            "\\subsection{Fig.~5, subthreshold partial waves}\n\\paper{``For the larger tolerances the partial waves are not approximately linear \\dots\\ the chiral zero of the amplitude disappears for the blue points. \\dots\\ The value $\\epsilon^\\chi=0.002$ \\dots\\ allows the physical value of $f_\\pi$.''}\n\n"
            "Ours: the curves are evaluated from the Arb coefficients of the accepted leaves on $0.05\\le s\\le3.95$ at the upper-branch $x_{\\rm ref}$ representative of each tolerance. The RMS against your digitised curves is at most 6.4\\,\\% of $f_0^0(3)$ (threshold 8\\,\\%); the S0 zero sits at $s=0.426$ ($2\\times10^{-3}$), $0.293$ ($4\\times10^{-3}$) and is absent at $6\\times10^{-3}$, as in the figure.\n\n"
            "\\begin{figure}[H]\\centering\\includegraphics[width=\\textwidth]{figures/fig5.png}\\caption{Fig.~5: the paper's three panels (top) and ours (bottom).}\\end{figure}\n\n"
            "\\subsection{Fig.~7, chiral-only phase shifts}\n\\paper{``We choose a point closest to the black dot \\dots\\ The partial waves at the magenta point and other nearby agree very well with experimental values for the S0 and S2 waves but not for the P1.''}\n\n"
            f"Ours: the representative is fixed by a rule set before any phase was read (nearest upper-branch endpoint to the black dot among the vertical sections at $x_{{\\rm ref}}+\\{{-0.002,-0.001,0,+0.001,+0.002\\}}$). RMS against your digitised curves: S0 {f(c4['S0']['rms_deg'],3)}$^\\circ$, S2 {f(c4['S2']['rms_deg'],3)}$^\\circ$, P1 {f(c4['P1']['rms_deg'],3)}$^\\circ$; P1 has no 90$^\\circ$ crossing below 1.2\\,GeV.\n\n"
            "\\begin{figure}[H]\\centering\\includegraphics[width=\\textwidth]{figures/fig7.png}\\caption{Fig.~7: paper (top) and ours (bottom).}\\end{figure}\n\n"
            "\\subsection{Fig.~8, region after the sum rules and the form-factor bounds}\n\\paper{``the plots are produced by taking $\\epsilon^{SR}=2\\times10^{-3}$ in (3.73) and $\\epsilon^{FF}=6\\times10^{-5}$''}; caption: \\paper{``The upper boundary shrinks notable but the lower not so much.''}\n\n"
            "Ours, three readings of $\\epsilon^{SR}$ and two of the (3.75) factor; all values at the $x_{\\rm ref}$ section, paper values from the digitised figure:\n\n"
            "\\begin{tabularx}{\\textwidth}{@{}Yrrrr@{}}\\toprule Reading & UV $+x$ end & upper shrink & lower rise & ratio \\\\\\midrule\n"
            f"raw box $2\\times10^{{-3}}$ per moment, factor at each node & {f(ra['UV +x end']['ours'],5)} & {f(ra['upper shrink at x_ref']['ours'],3)} & {f(ra['lower rise at x_ref']['ours'],3)} & {f(ra['upper/lower ratio']['ours'],3)} \\\\\n"
            f"raw box, factor frozen at $s_0$ & {f(rf['UV +x end']['ours'],5)} & {f(rf['upper shrink at x_ref']['ours'],3)} & {f(rf['lower rise at x_ref']['ours'],3)} & {f(rf['upper/lower ratio']['ours'],3)} \\\\\n"
            f"per-wave $L^2$ ball $2\\times10^{{-3}}$, factor frozen at $s_0$ & {f(rd['UV +x end']['ours'],5)} & {f(rd['upper shrink at x_ref']['ours'],3)} & {f(rd['lower rise at x_ref']['ours'],3)} & {f(rd['upper/lower ratio']['ours'],3)} \\\\\n"
            "paper (digitised) & 0.08112 & $2.32\\times10^{-4}$ & $3.7\\times10^{-5}$ & 6.3 \\\\\\bottomrule\\end{tabularx}\n\n"
            f"The relative reading (10\\,\\% of each moment) is infeasible in this model. We proved this by removing the S0 $n=0$ box, keeping the other three 10\\,\\% boxes and every other constraint, and minimising that moment: the minimum is {srm}, which is $1.57$ times the QCD value, so no point of the feasible set lies within 10\\,\\% of it. In the raw-box reading the accepted solutions sit at the upper edge of three of the four boxes (S0 $n=0$ at $+84\\,\\%$ of its target, S0 $n=1$ at $+1.8\\,\\%$, P1 $n=0$ at $+1.4\\,\\%$): the sum rules bind, in the direction of larger spectral integrals.\n\n"
            "\\begin{figure}[H]\\centering\\includegraphics[width=\\textwidth]{figures/fig8.png}\\caption{Fig.~8: paper (left); our chiral and chiral+UV boundary points with the digitised boundaries and a zoom on the $x_{\\rm ref}$ region (right).}\\end{figure}\n\n"
            "\\subsection{Fig.~9, P1 phase shift}\n\\paper{``The resonance energy where the phase shift crosses $\\pi/2$ is slightly shifted from the real world data on the mass of the rho at 770 MeV by roughly 6\\,\\%.''} Digitised crossings: 0.827, 0.824, 0.813\\,GeV.\n\n"
            f"Ours, at the three representatives of the same frozen rule (tip $=$ $+x$ end; ref $=$ nearest upper-branch point to the black dot; mid $=$ its neighbour towards the tip): tip crossing {f(c6['tip']['crossing_MeV'],4)}\\,MeV with $|S_{{P1}}|=0.585$ at 0.79\\,GeV; mid {f(c6['mid']['crossing_MeV'],4)}\\,MeV with $|S_{{P1}}|$ down to 0.077; ref: $|S_{{P1}}|$ falls to 0.24 at 0.73\\,GeV, and the node values admit two phase branches (the nearest-node lift swings to $-47^\\circ$; the $+180^\\circ$ branch crosses 90$^\\circ$ at 708\\,MeV). With the chiral constraints alone (Fig.~7) there is no crossing below 1.2\\,GeV, so the resonance is produced by the sum rules and the form-factor bounds, as in the paper; its position is 40--110\\,MeV below yours and the wave is far from elastic there. The paper does not show $|S|$; we add it in the second panel.\n\n"
            "\\begin{figure}[H]\\centering\\includegraphics[width=\\textwidth]{figures/fig9.png}\\caption{Fig.~9: paper (left) and our P1 phases at the representatives with the three digitised curves (right).}\\end{figure}\n"
            "\\begin{figure}[H]\\centering\\includegraphics[width=0.62\\textwidth]{figures/fig9_eta.png}\\caption{$|S_{P1}|$ at the same points (not shown in the paper).}\\end{figure}\n\n"
            "\\subsection{Fig.~10, S0 and S2 phase shifts}\n\\paper{``the three bootstrap (red, pink, light pink) curves \\dots\\ coincide at low energy but they spread at high energy.''}\n\n"
            f"Ours: S2 agrees with your curves (RMS {f(s2min,3)}--{f(s2max,3)}$^\\circ$; endpoints $-28$ to $-32^\\circ$ at 1.196\\,GeV against your $-18$ to $-36^\\circ$). S0 agrees below about 0.6\\,GeV and then rises faster: it crosses 90$^\\circ$ at 0.69\\,GeV and reaches 181--192$^\\circ$ at 1.196\\,GeV, where your curves are at 87--106$^\\circ$. With the chiral constraints alone our S0 matches your Fig.~7 to 0.4$^\\circ$, so the difference is introduced by the UV constraints.\n\n"
            "\\begin{figure}[H]\\centering\\includegraphics[width=\\textwidth]{figures/fig10.png}\\caption{Fig.~10: paper (top) and ours (bottom).}\\end{figure}\n\n"
            "\\subsection{Fig.~11, dependence on $M$ and $L$}\n\\paper{``For $M=50$, we have plotted the results with $L=8$, $10$, $12$. In all three partial waves there is a reasonably good convergence. \\dots\\ the P1 phase shifts show a $\\rho$ resonance with the peak slightly shifted depending on $M$.''}\n\n"
            "\\begin{tabularx}{\\textwidth}{@{}lrrrr@{}}\\toprule $(M,L)$ & UV $+x$ end & P1 crossing (MeV) & $\\min|S_{P1}|$ below 1.2\\,GeV & S0 at 1\\,GeV ($^\\circ$) \\\\\\midrule\n")
    for ml in ((50, 8), (50, 10), (50, 12), (45, 10), (60, 10)):
        r = c8[ml]; tex += f"({ml[0]}, {ml[1]}) & {f(r['x_tip'],5)} & {f(r['P1_crossing_MeV'],4)} & {f(r['min_eta_P1'],3)} & {f(r['S0_at_1GeV_deg'],3)} \\\\\n"
    tex += ("\\bottomrule\\end{tabularx}\n\nThe $L$ dependence is 1.5\\,MeV in the crossing and 0.3\\,\\% in the endpoint. The $M$ dependence is 37\\,MeV, with $M=45$ and $M=60$ both below $M=50$; in your figure the ordering is $M=60<M=45<M=50$. The S0 level at 1\\,GeV is 150--165$^\\circ$ at every configuration against about 95$^\\circ$ in your figure.\n\n"
            "\\begin{figure}[H]\\centering\\includegraphics[width=\\textwidth]{figures/fig11.png}\\caption{Fig.~11: paper (top) and ours (bottom).}\\end{figure}\n\n"
            "\\section{Why the fine features of Figs.~8--10 do not follow from the stated problem}\n"
            "Two diagnostics were run on the accepted UV solutions. Both are exact up to the solver tolerance and were registered before the runs.\n\n"
            "\\textbf{Near-optimal face.} At a boundary point with support direction $d$ and optimal value $v^*$ we add the constraint $d\\cdot(f_0^0,f_1^1)\\ge v^*-2\\times10^{-6}$ (twice the duality-gap threshold), keep the section, and maximise or minimise one node functional that is linear in the primal variables: $1-\\operatorname{Re}S(s_i)$, $\\operatorname{Im}S(s_i)$, $\\operatorname{Im}F(s_i)$. The range of the functional over this slab shows how far the amplitude at that boundary point is determined by the extremal problem. Your Fig.~9 phases at the nodes 0.792\\,GeV (68.7--76.2$^\\circ$) and 0.864\\,GeV (112.8--124.1$^\\circ$) have $\\cos2\\delta<0$, so any amplitude with those phases has $1-\\operatorname{Re}S>1$ whatever $|S|$ is.\n\n"
            "\\begin{tabularx}{\\textwidth}{@{}>{\\raggedright\\arraybackslash}p{0.17\\textwidth}>{\\raggedright\\arraybackslash}p{0.2\\textwidth}>{\\raggedright\\arraybackslash}p{0.17\\textwidth}Y@{}}\\toprule Representative & Functional & Range over the face & Consequence \\\\\\midrule\n"
            f"tip ($+x$ end) & $1-\\operatorname{{Re}}S_{{P1}}$ (0.792\\,GeV) & {rng('uv_SRa_tip','ImKH_P1_38')} & pinned to $\\pm0.03$ around the tip's own value 1.470; below the 1.72--1.85 of your curves at $|S|=1$ \\\\\n"
            f"tip & $1-\\operatorname{{Re}}S_{{P1}}$ (0.864\\,GeV) & {rng('uv_SRa_tip','ImKH_P1_39')} & below 1: your 112.8--124.1$^\\circ$ excluded for every $|S|$ \\\\\n"
            f"tip & $\\operatorname{{Im}}S_{{P1}}$ (0.792\\,GeV) & {rng('uv_SRa_tip','ImS_P1_38')} & negative on the whole face: every amplitude there is past the resonance at 0.79\\,GeV \\\\\n"
            f"$x_{{\\rm ref}}$ upper (ref) & $1-\\operatorname{{Re}}S_{{P1}}$ (0.792\\,GeV) & {rng('uv_SRa_section_hi','ImKH_P1_38')} & wide (undetermined), yet below 1 everywhere: your 68.7--76.2$^\\circ$ excluded for every $|S|$ \\\\\n"
            f"$x_{{\\rm ref}}$ upper & $1-\\operatorname{{Re}}S_{{P1}}$ (0.864\\,GeV) & {rng('uv_SRa_section_hi','ImKH_P1_39')} & below 1: excluded for every $|S|$ \\\\\n"
            f"$x_{{\\rm ref}}$ upper & $\\operatorname{{Im}}F_1$ (0.792\\,GeV) & {rng('uv_SRa_section_hi','ImF_P1_38')} & the vector form factor is essentially free at this point \\\\\n"
            f"tip & $1-\\operatorname{{Re}}S_{{S0}}$ (0.792\\,GeV) & {rng('uv_SRa_tip','ImKH_S0_38')} & the S0 phase is fixed at 111$^\\circ$ on the face; only $|S|$ varies (0.79--0.99) \\\\\\bottomrule\\end{{tabularx}}\n\n"
            "The tip face is rigid in the P1 observables and does not contain your P1 shape. The $x_{\\rm ref}$ face is degenerate---the premise that a boundary point carries one set of partial waves does not hold there---and it still excludes your P1 phase at both nodes. Choosing a different point on these faces, or a different solver, cannot reproduce Fig.~9 from this constraint set.\n\n"
            f"\\textbf{{Regulariser scale.}} Raising $M_{{\\rm reg}}$ from $10^2$ to $10^3$ moves the UV $+x$ end by {mreg_x}, the P1 crossing by {mreg_dm}\\,MeV and $\\min|S_{{P1}}|$ by {mreg_de}; the Fig.~8 asymmetry does not appear at either scale.\n\n"
            "\\begin{figure}[H]\\centering\\includegraphics[width=0.8\\textwidth]{figures/face_ranges.png}\\caption{Ranges of the node functionals over the near-optimal faces (bars) against the floor $1-\\operatorname{Re}S=1$ required by the paper's phases (dashed).}\\end{figure}\n\n"
            "\\section{Questions}\nThe following items would let us close the gap or confirm that it is real. We have your 2403 and 2505 codes and can see how they treat these points; we are asking specifically about the 2309 runs.\n"
            "\\begin{enumerate}[leftmargin=*,itemsep=2pt]\n"
            "\\item Which norm was used in (3.73), and is the QCD value compared with the raw moment or with the normalised one (divided by $s_0^{\\,n+2}$)? The raw per-moment box is the only reading with $\\epsilon^{SR}=2\\times10^{-3}$ that is feasible in our implementation.\n"
            "\\item In (3.75), is the factor between $F$ and $\\mathcal F$ evaluated at each node above $s_0$ or frozen at $s_0$?\n"
            "\\item Was a bound on the double spectral density $\\rho_{ij}$ (as in 2103.11484 Section~3, or the $\\ell^4$ bound of your 2403 code) applied in the 2309 runs, and at what scale? Without any bound the finite problem is ill-posed at 192-bit precision; with our in-model rule the $+x$ ends agree with your Figs.~3 and~4 to 0.4\\,\\%.\n"
            "\\item Were the Fig.~9--10 amplitudes the solver's optimal point of the support problem, or was a saturation / Watson step (as in your later code) already used to select the amplitude at each of the three points?\n"
            "\\item Could you share the $(f_0^0(3),f_1^1(3))$ coordinates and $|S_{P1}|$ of the red, pink and light pink points, and the value of $s_0$ used for Figs.~8--11?\n"
            "\\end{enumerate}\n\n"
            "\\section{Receipts}\n\\texttt{receipts/C1\\_RESULT.json} \\dots\\ \\texttt{C8\\_RESULT.json} (per-figure verdicts on the pre-registered rules), \\texttt{receipts/FACE\\_RESULT.json} (face diagnostic), \\texttt{receipts/UV\\_REPRESENTATIVE\\_SELECTION.json} and \\texttt{IR\\_REPRESENTATIVE\\_SELECTION.json} (frozen representative choices), \\texttt{receipts/GATE\\_LOG.md} (time-stamped log of every run and decision). Each accepted leaf keeps its PMP, the SDPB output and the Arb verification (\\texttt{report.json}); these are available on request (several GB).\n"
            "\\end{document}\n")
    D = D.resolve(); tex_path = D / "REPRODUCTION_COMPARISON_EN.tex"; tex_path.write_text(tex)
    r = subprocess.run(["tectonic", "-o", str(D), str(tex_path)], capture_output=True, text=True, cwd=str(D))
    print(r.stdout[-600:], r.stderr[-1500:])
    print("pdf", (D / "REPRODUCTION_COMPARISON_EN.pdf").exists())
    return r.returncode


if __name__ == "__main__":
    raise SystemExit(main())
