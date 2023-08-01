#!/bin/bash

if [ ${distance} = "0" ]
then
scil_filter_tractogram.py ${tractogram} ${basename}__${out_extension}.trk \
    --filtering_list ${filtering_list} ${extract_masks} -f \
    --display_count  > ${basename}__${out_extension}.txt;
else
scil_filter_tractogram.py ${tractogram} ${basename}__${out_extension}.trk \
    --filtering_list ${filtering_list} ${extract_masks} -f \
     --overwrite_distance both_ends include ${distance} --overwrite_distance either_end include ${distance} --display_count  > ${basename}__${out_extension}.txt;
fi

if [ ${keep} = "true" ]
then
  scil_streamlines_math.py difference ${tractogram} \
                                      ${basename}__${out_extension}.trk \
                                      ${sid}__${remaining_extension}.trk ;
fi
 