#!/bin/bash
workdir=~/work/
PATH=/bin:/usr/bin:/sbin:/usr/sbin; export PATH
export LANG=C
cd $workdir
git add .
git commit -m 'daily update'
git push origin master