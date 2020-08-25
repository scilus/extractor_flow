#!/bin/bash

scil_filter_tractogram.py ${tractogram} ${basename}__${out_extension}.trk \
    --filtering_list ${filtering_list} -f --display_count > ${basename}__${out_extension}.txt;

if ${params.debug}
then
  scil_streamlines_math.py difference ${tractogram} \
                                      ${basename}__${out_extension}.trk \
                                      ${sid}__${remaining_extension}.trk ;
fi
