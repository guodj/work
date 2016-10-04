#!/bin/bash
basedir=~/backup
PATH=/bin:/usr/bin:/sbin:/usr/sbin; export PATH
export LANG=C
basefile=$basedir/week.$(date +%Y-%m-%d).tar.gz
rmfile=$basedir/week.$(date --date='7 days ago' +%Y-%m-%d).tar.gz
[ ! -d "$basedir" ] && mkdir $basedir
cd ~
tar -zpc -f $basefile work document desktop /data/output/
rm -f $rmfile
