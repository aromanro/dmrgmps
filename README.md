# dmrgmps
A Density Matrix Renormalization Group program with Matrix Product State and Matrix Product Operator formalism

> [!NOTE]
> This project is vibed into existence, starting with the old dmrg project.
> It's not in the state as it was generated, I looked briefly over the code and changed/removed/added quite a bit.
> Still have to check all the details... so I'm not 100% it's ok.

As already mentioned on the blog, there is a more modern method to deal with those blocks... a method I already used in some other projects, like [iTEBD](https://github.com/aromanro/TEBD) or [QCSim](https://github.com/aromanro/QCSim) (more specifically, the MPS and MPO simulators).

I used Claude Opus 4.8, the old [dmrg project](https://github.com/aromanro/dmrg) as a starting point (I simply copied and renamed some things before starting to port it), then I used a prompt explaining what I needed, asking the AI to look at the dmrg project and also pointing out the MPS implementation from the [QCSim project](https://github.com/aromanro/QCSim) (I explicitely named all the relevant files).
It wasn't perfect at first, it had bugs which were caught and fixed by the tests I asked to be implemented, then still nothing was displayed in the chart because the operators were not transformed properly.
I also used GPT 5.5 to check the code and I also checked the code myself, although not all details. Some of the tests are added/modified by me, I also had to remove or change portions of the old code and I had to change quite a bit of the new code because there was a lot of switching between tensors and matrices... some of the computations that could be done with tensor contractions were done by slicing and operating with matrices.
Anyway, in the end everything looks ok, it seems to get (almost) the same results as the old code, in some cases even better results.
All the tests pass.
With this I basically used all the tokens for one month from my GitHub Copilot Pro+ subscription, so apparently I stressed AI quite a bit to obtain this.

This program displays the results for the whole chain, unlike the original one. It has three methods: one-site, one site with subspace expansion and two-sites.
It's also faster, as expected.


Description on Computational Physics blog (this is for the old project, implemented the 'old' way, without MPS and MPO): https://compphys.go.ro/density-matrix-renormalization-group/

The program can be compiled with Visual C++ on Windows (native, with mfc). 
It uses Eigen library for matrices.

Currenly only Heisenberg chains are implemented (both spin 1/2 and 1).

### TOOLS

The project compiles on Windows with Visual Studio 2026 (the code can be compiled with older versions starting with VS 2015, but it's currently maintained with VS 2026 and C++ 17 or higher).

### LIBRARIES

Besides mfc and other typical VC++ runtime libraries, the program uses GDI+ for drawing.
Eigen 5 is also used. It worked with an older version as well, but I recommend using ver 5 because it has some fixed issues. 

### PROGRAM IN ACTION

This is for the old project, the new one displays the whole chain, otherwise the look is identical (except some names changed from dmrg to dmrgmps).
There is also a combo box in the Options dialog box that allows selecting one of the three methods.

[![Program video](https://img.youtube.com/vi/LHnebbE_XP4/0.jpg)](https://youtu.be/LHnebbE_XP4)

