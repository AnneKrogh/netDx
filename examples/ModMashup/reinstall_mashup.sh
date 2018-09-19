#!/bin/bash
cd ~/software/ModMashup.jl
git commit -m "debugging" .
git push
julia julia_upMash.jl
