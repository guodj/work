#!/bin/bash
basedir=/home/gdj/backup
PATH=/bin:/usr/bin:/sbin:/usr/sbin; export PATH
export LANG=C
basefile=$basedir/day.$(date +%Y-%m-%d).tar.gz
rmfile=$basedir/day.$(date --date='2 days ago' +%Y-%m-%d).tar.gz
[ ! -d "$basedir" ] && mkdir $basedir
cd /home/gdj
tar -zpc -f $basefile\
    .vimrc .vim .viminfo\
    .bashrc\
    .config/matplotlib/matplotlibrc\
    .muttrc .mutt .msmtprc .getmail .mail .mailcap\
    .backupday.sh .backupweek.sh .backupwork.sh .backupblog.sh\
rm -f $rmfile
