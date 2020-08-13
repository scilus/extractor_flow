#!/bin/bash

scil_filter_tractogram.py ${tractogram} ${basename}_${out_extension}.trk \
    --filtering_list ${filtering_list} -f;

if ${params.keep}
then
  scil_streamlines_math.py difference ${tractogram} \
                                      ${basename}_${out_extension}.trk \
                                      ${sid}_${remaining_extension}.trk ;
  scil_count_streamlines.py ${basename}_${out_extension}.trk > ${basename}_${out_extension}.txt;
  scil_count_streamlines.py ${sid}_${remaining_extension}.trk > ${sid}_${remaining_extension}.txt;
fi
