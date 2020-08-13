#!/bin/bash

scil_filter_tractogram.py ${tractogram} ${basename}_${out_extension}.trk \
  --drawn_roi ${atlas} ${mode} ${criteria} -f --display_count > ${basename}_${out_extension}.txt;

if ${params.keep}
then
  scil_streamlines_math.py difference ${tractogram} \
                                      ${basename}_${out_extension}.trk \
                                      ${sid}_${remaining_extension}.trk;
fi
