#! /bin/sh

/home/scott/rsync-2.4.3/rsync -ar --rsh=ssh --rsync-path=/p/rsync/v2_4_3/rsync \
      --log-format="`startapo.sh /data/spectro/%f`" \
      "scott@sdsshost:/data/spectro/5*" /data/spectro/
/home/scott/rsync-2.4.3/rsync -ar --rsh=ssh --rsync-path=/p/rsync/v2_4_3/rsync \
      --log-format="/astrolog/%f" "scott@sdsshost:/astrolog/5*" /astrolog/   



