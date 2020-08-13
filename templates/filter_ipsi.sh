#!/bin/bash

scil_filter_tractogram.py ${tractogram} ${basename}_all_brainstem_either_end_CC_ipsi_${side}.trk \
    --drawn_roi ${atlasROI1} either_end include \
    --drawn_roi ${atlasROI2} any include -f;
