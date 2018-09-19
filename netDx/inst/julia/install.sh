#!/bin/bash
julia -e "Pkg.rm("ModMashup")"
julia -e 'Pkg.clone("https://github.com/BaderLab/ModMashup_private.jl.git")'
