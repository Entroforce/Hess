#!/bin/sh
git submodule init &&
git submodule update &&
cp indigo/indigo_version.h indigo/Indigo/api/c/indigo/src/
cd hess-preparation && make CONF=Release && cd .. &&
cd hess-empirical-docking && make CONF=Release && cd ..
