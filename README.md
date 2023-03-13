# ANN203 (2023)

This repo contains the homeworks for the
[`ANN203`](https://synapses.ensta-paris.fr/catalogue/2022-2023/ue/7679/ann203-methodes-numeriques-matricielles-avancees-analyse-et-experimentation)
course. 
## Getting started
 
 * Install [Julia]((https://julialang.org/downloads/). I suggest using
   [Juliaup](https://github.com/JuliaLang/juliaup) for an easier version
   control.
 * Install [Pluto](https://github.com/fonsp/Pluto.jl)
 * Open a terminal and enter the *Julia REPL* by typing `julia`
 * Enter the following lines to open a browser with *Pluto.jl*
 ```julia
using Pluto
Pluto.run()
 ```

You should see something like this on your web browser:

![GitHub Logo](/pluto_home_page.png)

The last step is to load the homework assignment. Copy and paste the following
URL on the **Open from file** options and press `Open`:

<https://github.com/maltezfaria/ANN203/blob/main/ANN203_tp0.jl> 

That should load the notebook for the `TP0` and start running some code.

## To `Julia` from `Matlab`

If you are familiar with `Matlab`, it should not be too difficult to get started
with `Julia`. A lot of the syntax looks familiar, and `Julia`'s standard library
comes with a `LinearAlgebra` package implementing many of the things that you
love about `Matlab` (e.g. the *magical* backslash `\` operator). Here is a [list
of noteworthy
differences](https://docs.julialang.org/en/v1/manual/noteworthy-differences/)
that may help you get started translating code.
 

## Workflows and `Pluto.jl`

As you will notice in your first assignment, the homeworks for this class are
available both in the form of a *static* HTML file (which you can open in any
browser), and in the form of a `Pluto` notebook.

[`Pluto`](https://github.com/fonsp/Pluto.jl) provides a notebook-like
environment for combining `Julia` with pretty text that can easily be rendered
on a website, and you can think of it as a web-based *workflow* for developing
code (similar to e.g. *Jupyter* notebooks if you have heard of those). While in
general I don't recommend using `Pluto` for coding anything long or complex,
most of the code you will write in this class is simple enough that using
`Pluto` as a *web-based*
[IDE](https://en.wikipedia.org/wiki/Integrated_development_environment) may be a
viable option. My suggestion is to give it a try to see if it fits you.

Should you find yourself hating `Pluto.jl` (and missing `matlab`), I suggest you
try [Visual Studio Code](https://code.visualstudio.com) with its [julia
extension](https://code.visualstudio.com/docs/languages/julia) to have an
*out-of-the-box* experience as similar to `matlab` as possible. Keep in mind
that I still expect you to submit your code as a `Pluto` notebook, but you can
always just copy-and-past the code to `Pluto` after you get it to work.