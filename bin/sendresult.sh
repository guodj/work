#!/bin/bash
touch ~/tmp/logfile
stime=$(date "+%F %T")
if [ "$($1)" ]; then
    etime=$(date "+%F %T")
    echo "Great, '$1' ran successfully!" >> ~/tmp/logfile
    echo "Start time: $stime" >> ~/tmp/logfile
    echo "Stop time: $etime" >> ~/tmp/logfile
    mutt -s "Good news" guodj90@qq.com < ~/tmp/logfile
else
    echo "'$1' has something wrong." | mutt -s "Bad news" guodj90@qq.com
fi
rm ~/tmp/logfile
