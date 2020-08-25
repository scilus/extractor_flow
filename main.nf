#!/usr/bin/env nextflow

params.root = false
params.help = false
params.debug = true

if(params.help) {
    usage = file("$baseDir/USAGE")

    bindings = ["rois_folder":"$params.rois_folder",
                "filtering_lists_folder": "$params.filtering_lists_folder"]

    engine = new groovy.text.SimpleTemplateEngine()
    template = engine.createTemplate(usage.text).make(bindings)
    print template.toString()
    return
    }

if (params.root){
    log.info "Input: $params.root"
    root = file(params.root)
    in_tractogram = Channel
        .fromFilePairs("$root/**/*.trk",
                       size:1,
                       maxDepth:1,
                       flat: true) {it.parent.name}
}
else {
    error "Error ~ Please use --root for the input data."
}

sides = params.sides?.tokenize(',')

process Remove_Out_JHU {
    cpus 1
    tag = "Remove out of JHU"

    input:
      set sid, file(tractogram) from in_tractogram

    output:
      set sid, "${sid}__wb_in_JHU.trk" into wb_for_rm_crossing_gyri
      file "${sid}__wb_in_JHU.txt" optional true
      set sid, "${sid}__wb_out_JHU.trk" optional true
      file "${sid}__wb_out_JHU.txt" optional true

    script:
    atlas=params.rois_folder+params.atlas.JHU
    mode="all"
    criteria="include"
    out_extension='wb_in_JHU'
    remaining_extension='wb_out_JHU'
    basename="${sid}"

    template "filter_with_atlas.sh"
}

process Remove_Crossing_Gyri {
  cpus 1
  tag = "Remove crossing the gyri limits of the JHU template"

  input:
    set sid, file(tractogram) from wb_for_rm_crossing_gyri

  output:
   set sid, "${sid}__wb_rm_crossing_gyri.trk" into wb_for_pruning
   file "${sid}__wb_rm_crossing_gyri.txt" optional true
   set sid, "${sid}__wb_crossing_gyri.trk" optional true
   file "${sid}__wb_crossing_gyri.txt" optional true

   script:
   atlas=params.rois_folder+params.atlas.shell_limits
   mode="any"
   criteria="exclude"
   out_extension='wb_rm_crossing_gyri'
   remaining_extension='wb_crossing_gyri'
   basename="${sid}"

   template "filter_with_atlas.sh"
}

process Pruning {
  cpus 1
  tag = "Pruning min ${params.min_streaminline_lenght}mm"

  input:
    set sid, file(tractogram) from wb_for_pruning

  output:
    set sid, "${sid}__wb_min${params.min_streaminline_lenght}.trk" into wb_for_rmloop
    file "${sid}__wb_min${params.min_streaminline_lenght}.txt" optional true
    set sid, "${sid}__wb_max${params.min_streaminline_lenght}.trk" optional true
    file "${sid}__wb_max${params.min_streaminline_lenght}.txt" optional true

  script:

    """
    scil_filter_streamlines_by_length.py ${tractogram} --minL ${params.min_streaminline_lenght} \
                                         --maxL ${params.max_streaminline_lenght} ${sid}__wb_min${params.min_streaminline_lenght}.trk \
                                         -f --display_counts > ${sid}__wb_min${params.min_streaminline_lenght}.txt
    if ${params.debug}
    then
      scil_streamlines_math.py difference ${tractogram} ${sid}__wb_min${params.min_streaminline_lenght}.trk ${sid}__wb_max${params.min_streaminline_lenght}.trk -f
    fi
    """
}

process Remove_Loops {
  cpus 1
  tag = "Remove loop"

  input:
    set sid, file(wb_min20) from wb_for_rmloop

  output:
    set sid, "${sid}__wb_min20_noloop.trk" into wb_for_rm_end_in_cc_dwm
    set sid, "${sid}__wb_loops.trk" optional true
    file "${sid}__wb_min20_noloop.txt" optional true
    file "${sid}__wb_loops.txt" optional true

  script:

    """
    scil_detect_streamlines_loops.py ${wb_min20} ${sid}__wb_min20_noloop.trk \
                                     -a ${params.loop_angle_threshold}  \
                                     --looping_tractogram ${sid}__wb_loops.trk \
                                     -f
    if ${params.debug};
    then
      scil_count_streamlines.py ${sid}__wb_loops.trk > ${sid}__wb_loops.txt
      scil_count_streamlines.py ${sid}__wb_min20_noloop.trk > ${sid}__wb_min20_noloop.txt
    fi
    """
}

process Removing_End_In_CC_DWM {
  cpus 1
  tag = "Remove end in CC DWM"

  input:
    set sid, file(wb_min20_noloop) from wb_for_rm_end_in_cc_dwm

  output:
    set sid, "${sid}__wb_clean01.trk" into wb_for_cerebellum_split
    file "${sid}__wb_clean01.txt" optional true

  script:
  """
  scil_filter_tractogram.py ${wb_min20_noloop} ${sid}__wb_clean01.trk \
                    --drawn_roi ${params.rois_folder}${params.atlas.cc} either_end exclude \
                    --drawn_roi ${params.rois_folder}${params.atlas.dwm} either_end exclude \
                    -f --display_count > ${sid}__wb_clean01.txt
  """
}

process Split_Either_End_In_Cerebellum {
  cpus 1
  tag = "Split either end in cerebellum"

  input:
    set sid, file(tractogram) from wb_for_cerebellum_split

  output:
    set sid, "${sid}__wb_clean01_nocereb.trk" into not_end_in_cerebellum
    set sid, "${sid}__all_cereb.trk" into endInCerebellum
    file "${sid}__all_cereb.txt" optional true
    file "${sid}__wb_clean01_nocereb.txt" optional true

  script:
  atlas=params.rois_folder+params.atlas.cerebellum
  mode=params.mode.either_end
  criteria=params.criteria.exclude
  out_extension='wb_clean01_nocereb'
  remaining_extension='all_cereb'
  basename="${sid}"

  template "filter_with_atlas.sh"
}

process Split_Either_End_In_Brainstem {
  cpus 1
  tag = "Split either end in brainstem"

  input:
    set sid, file(tractogram) from not_end_in_cerebellum

  output:
    set sid, "${sid}__wb_clean02.trk" into not_end_in_brainstem
    set sid, "${sid}__all_brainstem.trk" into endInBrainstem
    file "${sid}__wb_clean02.txt" optional true
    file "${sid}__all_brainstem.txt" optional true

  script:
    atlas=params.rois_folder+params.atlas.brainstem
    mode=params.mode.either_end
    criteria=params.criteria.exclude
    out_extension='wb_clean02'
    remaining_extension='all_brainstem'
    basename="${sid}"

    template "filter_with_atlas.sh"
}

process Remove_GM {
  cpus 1
  tag = "Brain Remove GM"

  input:
    set sid, file(tractogram) from not_end_in_brainstem

  output:
    set sid, "${sid}__wb_either_CGM_SWM.trk" into endInCGMSWI
    set sid, "${sid}__no_CGM_SWM.trk" into endNotInCGMSWI
    file "${sid}__wb_either_CGM_SWM.txt" optional true
    file "${sid}__no_CGM_SWM.txt" optional true

  script:
    atlas=params.rois_folder+params.atlas.CGM_SWM
    mode=params.mode.either_end
    criteria=params.criteria.include
    out_extension='wb_either_CGM_SWM'
    remaining_extension='no_CGM_SWM'
    basename="${sid}"

    template "filter_with_atlas.sh"
}

process Split_CC {
  cpus 1
  tag = "split CC"

  input:
    set sid, file(tractogram) from endInCGMSWI

  output:
    set sid, "${sid}__tmp_CC.trk" into inCC4BG, inCC4Other
    set sid, "${sid}__wb_either_CGM_SWM__noCC.trk" into notInCC
    file "${sid}__wb_either_CGM_SWM_noCC.txt" optional true
    file "${sid}__tmp_CC.txt" optional true

  script:
  atlas=params.rois_folder+params.atlas.midline
  mode=params.mode.any
  criteria=params.criteria.exclude
  out_extension="noCC"
  remaining_extension='tmp_CC'
  basename=tractogram.getSimpleName()

  template "filter_with_atlas.sh"
}

process Split_CC_BG {
  cpus 1
  tag = "Split Basal Ganglia"

  input:
    set sid, file(tractogram) from inCC4BG
    each side from sides

  output:
    set sid, "${sid}__contra_BG_${side}.trk" into inCCBG
    file "${sid}__contra_BG_${side}.txt" optional true

  script:
  atlas=params.rois_folder+params.atlas.subcortical+"_${side}.nii.gz"
  mode=params.mode.either_end
  criteria=params.criteria.include
  out_extension="contra_BG_" + "${side}"
  remaining_extension=params.ext.notUsed
  basename="${sid}"

  template "filter_with_atlas.sh"
}

process Remove_Unplausible_CC {
  cpus 1
  tag = "Remove unplausible CC"

  input:
    set sid, file(tractogram) from inCC4Other

  output:
    set sid, "${sid}__CC_Cx.trk" into ccCleaned
    set sid, "${sid}__CC_lost.trk" optional true into CC_lost
    file "${sid}__CC_Cx.txt" optional true
    file "${sid}__CC_lost.txt" optional true


  script:
  """
  scil_filter_tractogram.py ${tractogram} ${sid}__CC_Cx.trk \
    --drawn_roi ${params.rois_folder}${params.atlas.allsubcortical} either_end exclude \
    --drawn_roi ${params.rois_folder}${params.atlas.brainstemINF} any exclude \
    --drawn_roi ${params.rois_folder}${params.atlas.ic} any exclude \
    --drawn_roi ${params.rois_folder}${params.atlas.allThal} any exclude \
    --drawn_roi ${params.rois_folder}${params.atlas.midline} either_end exclude \
    --drawn_roi ${params.rois_folder}${params.atlas.allL} both_ends exclude \
    --drawn_roi ${params.rois_folder}${params.atlas.allR} both_ends exclude -f;
  if ${params.debug}
  then
    scil_streamlines_math.py difference ${tractogram} ${sid}__CC_Cx.trk ${sid}__CC_lost.trk
    scil_count_streamlines.py ${sid}__CC_lost.trk > ${sid}__CC_lost.txt
    scil_count_streamlines.py ${sid}__CC_Cx.trk > ${sid}__CC_Cx.txt
  fi
  """
}

process Split_NoCC_Asso_BG {
  cpus 1
  tag = "Split not CC in asso BG and not BG"

  input:
    set sid, file(tractogram) from notInCC

  output:
    set sid, "${sid}__all_subcortical_from_CGM_SWM_noCC_f.trk" into assoBG
    file "${sid}__all_subcortical_from_CGM_SWM_noCC_f.txt" optional true
    set sid, "${sid}__asso_noBG.trk" into assoNoBG
    file "${sid}__asso_noBG.txt" optional true

  script:
  atlas=params.rois_folder+params.atlas.allsubcortical
  mode=params.mode.either_end
  criteria=params.criteria.include
  out_extension="all_subcortical_from_CGM_SWM_noCC_f"
  remaining_extension='asso_noBG'
  basename="${sid}"

  template "filter_with_atlas.sh"
}

process Split_Asso_In_Hemi {
  cpus 1
  tag = "Split asso in hemispheres"

  input:
    set sid, file(tractogram) from assoNoBG
    each side from sides

  output:
    set sid, val(side), "${sid}__asso_${side}.trk" into asso
    set sid, "${sid}__asso_${side}_lost.trk" optional true into assoLost
    file "${sid}__asso_${side}.txt" optional true
    file "${sid}__asso_${side}_lost.txt" optional true

  script:
  if (side=='L')
    atlas=params.rois_folder+params.atlas.all+"_R.nii.gz"
  else
    atlas=params.rois_folder+params.atlas.all+"_L.nii.gz"

  mode=params.mode.any
  criteria=params.criteria.exclude
  out_extension="asso_${side}"
  remaining_extension="asso_${side}_lost"
  basename="${sid}"

  template "filter_with_atlas.sh"
}

assoLost.groupTuple().map{it.flatten().toList()}.set{assoLostGroup}

process Extract_Lost {
  cpus 1
  tag = "Extract asso lost"

  input:
    set sid, file(trk1), file(trk2) from assoLostGroup

  output:
    file "${sid}__asso_lost.txt"
    set sid, "${sid}__asso_lost.trk" into allassoLost


  script:
  """
  scil_streamlines_math.py intersection ${trk1} \
                                        ${trk2} \
                                        ${sid}__asso_lost.trk
  scil_count_streamlines.py ${sid}__asso_lost.trk > ${sid}__asso_lost.txt
  """
}

process Split_UShape_CGM_Asso {
  cpus 1
  tag = "Extracting U-shaped and streamlines restricted to Cortical GM and removing them from asso"

  input:
    set sid, val(side), file(tractogram) from asso

  output:
    set sid, "${sid}__asso_only_in_CGM_${side}.trk" into assoCGM
    set sid, "${sid}__asso_Ushape_${side}.trk" into assoUShape
    set sid, val(side), "${sid}__asso_f_${side}.trk" into assof
    file "${sid}__asso_only_in_CGM_${side}.txt" optional true
    file "${sid}__asso_Ushape_${side}.txt" optional true
    file "${sid}__asso_f_${side}.txt" optional true

  script:
  """
    scil_filter_tractogram.py ${tractogram} ${sid}__tmp1_${side}.trk \
      --drawn_roi ${params.rois_folder}${params.atlas.CGM}_${side}.nii.gz ${params.mode.both_ends} include -f

    scil_streamlines_math.py difference ${tractogram} ${sid}__tmp1_${side}.trk \
                             ${sid}__asso_SWM_${side}.trk -f

    scil_filter_tractogram.py ${sid}__tmp1_${side}.trk ${sid}__asso_only_in_CGM_${side}.trk \
      --drawn_roi ${params.rois_folder}${params.atlas.SWM}_${side}.nii.gz ${params.mode.any} exclude

    scil_streamlines_math.py difference ${sid}__tmp1_${side}.trk ${sid}__asso_only_in_CGM_${side}.trk \
                                 ${sid}__tmp2_${side}.trk -f

    scil_filter_tractogram.py ${sid}__tmp2_${side}.trk ${sid}__asso_Ushape_${side}.trk \
      --drawn_roi ${params.rois_folder}${params.atlas.DWM}_${side}.nii.gz ${params.mode.any} exclude

    scil_streamlines_math.py difference ${sid}__tmp2_${side}.trk ${sid}__asso_Ushape_${side}.trk \
                               ${sid}__asso_DWM_${side}.trk -f

    scil_streamlines_math.py concatenate ${sid}__asso_DWM_${side}.trk ${sid}__asso_SWM_${side}.trk ${sid}__asso_f_${side}.trk -f

    if ${params.debug}
    then
      scil_count_streamlines.py ${sid}__asso_only_in_CGM_${side}.trk > ${sid}__asso_only_in_CGM_${side}.txt
      scil_count_streamlines.py ${sid}__asso_Ushape_${side}.trk > ${sid}__asso_Ushape_${side}.txt
      scil_count_streamlines.py ${sid}__asso_f_${side}.trk > ${sid}__asso_f_${side}.txt
    fi
  """
}

process Remove_Unplausible_Long_Range_Asso {
  cpus 1
  tag = "Extracting unplausible long-range association streamlines passing through subcortical structures (Cd, Put, GP, Thal, Amyg)"

  input:
    set sid, val(side), file(tractogram) from assof

  output:
    set sid, val(side), "${sid}__asso_all_intra_inter_${side}.trk" into assoAllIntraInter
    set sid, "${sid}__asso_lost2_${side}.trk" into assoLost2
    file "${sid}__asso_all_intra_inter_${side}.txt" optional true
    file "${sid}__asso_lost2_${side}.txt" optional true

  script:
  atlas=params.rois_folder+params.atlas.allsubcortical
  mode="any"
  criteria="exclude"
  out_extension="asso_all_intra_inter_${side}"
  remaining_extension="asso_lost2_${side}"
  basename="${sid}"

  template "filter_with_atlas.sh"
}

inCCBG.groupTuple().map{it.flatten().toList()}.set{inCCBG_List}
assoUShape.groupTuple().map{it.flatten().toList()}.set{assoUShape_list}

assoAllIntraInter.into{assoAllIntraInter_for_filtering;
                       assoAllIntraInter_plausible}

assoAllIntraInter_plausible.groupTuple().map{it.flatten().toList()}.set{assoAllIntraInter_list}


ccCleaned.into{ccCleanedPlausible;
               CC_for_homotopic}

endInBrainstem.into{brainstemForPlausible; inBrainstem}

endInCerebellum.into{endInCerebellumForPlausible; inCerebellum}



assoCGM.groupTuple().map{it.flatten().toList()}.set{assoCGM_list}
assoLost2.groupTuple().map{it.flatten().toList()}.set{assoLost2_list}


/*
  Cerebellum
*/

process Cerebellum_Remove_GM {

  cpus 1
  tag = "Cereb - Remove GM"

  input:
    set sid, file(tractogram) from inCerebellum

  output:
    set sid, "${sid}__all_cereb_nocx.trk" into ccEndNotInCGMSWI
    set sid, "${sid}__cereb_bin_01.trk" optional true into ccEndInCGMSWI
    file "${sid}__all_cereb_nocx.txt" optional true
    file "${sid}__cereb_bin_01.txt" optional true

  script:
    atlas=params.rois_folder+params.atlas.CGM_SWM
    mode=params.mode.either_end
    criteria=params.criteria.exclude
    out_extension='all_cereb_nocx'
    remaining_extension='cereb_bin_01'
    basename="${sid}"

    template "filter_with_atlas.sh"
}

process Cerebellum_Remove_Cerebellar_Hemisphere {
  cpus 1
  tag = "Cereb - Remove GM"

  input:
    set sid, file(tractogram) from ccEndNotInCGMSWI

  output:
    set sid, "${sid}__all_cereb_nocx_f.trk" into noCXf
    set sid, "${sid}__cereb_bin_02.trk" optional true into bin02
    file "${sid}__all_cereb_nocx_f.txt" optional true
    file "${sid}__cereb_bin_02.txt" optional true

  script:
    atlas=params.rois_folder+params.atlas.cerebellum_GWM
    mode=params.mode.either_end
    criteria=params.criteria.include
    out_extension='all_cereb_nocx_f'
    remaining_extension='cereb_bin_02'
    basename="${sid}"

    template "filter_with_atlas.sh"
}

process Cerebellum_Split {
  cpus 1
  tag = "Cereb - Split"

  input:
    set sid, file(tractogram) from noCXf

  output:
    set sid, "${sid}__all_cereb_nocx_in_cereb.trk" into noCX_inCereb
    set sid, "${sid}__all_cereb_nocx_out_cereb.trk" into noCX_outCereb
    file "${sid}__all_cereb_nocx_in_cereb.txt" optional true
    file "${sid}__all_cereb_nocx_out_cereb.txt" optional true

  script:
  atlas=params.rois_folder+params.atlas.cerebellum
  mode=params.mode.both_ends
  criteria=params.criteria.include
  out_extension='all_cereb_nocx_in_cereb'
  remaining_extension='all_cereb_nocx_out_cereb'
  basename="${sid}"

  template "filter_with_atlas.sh"
}

process Cerebellum_No_Loop{
  cpus 1
  tag = "Cereb - no loop"

  input:
    set sid, file(tractogram) from noCX_inCereb

  output:
    set sid, "${sid}__all_cereb_nocx_in_cereb_f.trk" into noC
    set sid, "${sid}__cereb_bin_03.trk" into bin03
    file "${sid}__all_cereb_nocx_in_cereb_f.txt" optional true

  script:
  """
    scil_filter_tractogram.py ${tractogram} ${sid}__all_cereb_nocx_in_cereb_f.trk \
      --drawn_roi ${params.rois_folder}${params.atlas.medulla} any exclude \
      --drawn_roi ${params.rois_folder}${params.atlas.midbrainNoSCP} any exclude \
      -f \
      --display_count > ${sid}__all_cereb_nocx_in_cereb_f.txt

    scil_streamlines_math.py difference ${tractogram} ${sid}__all_cereb_nocx_in_cereb_f.trk ${sid}__cereb_bin_03.trk
  """
}

process Cerebellum_In_Medulla {
  cpus 1
  tag = "Cereb - In Medulla"

  input:
    set sid, file(tractogram) from noCX_outCereb

  output:
    set sid, "${sid}__all_cereb_nocx_out_cereb_in_medulla.trk" into inMedulla
    set sid, "${sid}__all_cereb_nocx_out_cereb_no_medulla.trk" into noMedulla
    file "${sid}__all_cereb_nocx_out_cereb_in_medulla.txt" optional true
    file "${sid}__all_cereb_nocx_out_cereb_no_medulla.txt" optional true

  script:
  """
  scil_filter_tractogram.py ${tractogram} ${sid}__all_cereb_nocx_out_cereb_in_medulla.trk \
    --drawn_roi ${params.rois_folder}${params.atlas.medulla} either_end include \
    --drawn_roi ${params.rois_folder}${params.atlas.thal_CP_MidBrain} any exclude \
    --drawn_roi ${params.rois_folder}${params.atlas.MCPant} any exclude \
    --drawn_roi ${params.rois_folder}${params.atlas.CGM_SWM_DWM} any exclude -f

  scil_streamlines_math.py difference ${tractogram} ${sid}__all_cereb_nocx_out_cereb_in_medulla.trk ${sid}__all_cereb_nocx_out_cereb_no_medulla.trk

  if ${params.debug}
  then
    scil_count_streamlines.py ${sid}__all_cereb_nocx_out_cereb_in_medulla.trk > ${sid}__all_cereb_nocx_out_cereb_in_medulla.txt
    scil_count_streamlines.py ${sid}__all_cereb_nocx_out_cereb_no_medulla.trk > ${sid}__all_cereb_nocx_out_cereb_no_medulla.txt
  fi
  """
}

process Cerebellum_In_Pons{
  cpus 1
  tag = "Cereb - In Pons"

  input:
    set sid, file(tractogram) from noMedulla

  output:
    set sid, "${sid}__all_cereb_nocx_out_cereb_in_pons.trk" into inPons
    set sid, "${sid}__all_cereb_nocx_out_cereb_no_pons.trk" into noPons
    file "${sid}__all_cereb_nocx_out_cereb_in_pons.txt" optional true
    file "${sid}__all_cereb_nocx_out_cereb_no_pons.txt" optional true


  script:
  """
  scil_filter_tractogram.py ${tractogram} ${sid}__all_cereb_nocx_out_cereb_in_pons.trk \
    --drawn_roi ${params.rois_folder}${params.atlas.pons} either_end include \
    --drawn_roi ${params.rois_folder}${params.atlas.medulla} any exclude \
    --drawn_roi ${params.rois_folder}${params.atlas.thal_CP_MidBrain} any exclude \
    --drawn_roi ${params.rois_folder}${params.atlas.CGM_SWM_DWM} any exclude \
    --drawn_roi ${params.rois_folder}${params.atlas.MCPant} any include -f

  scil_streamlines_math.py difference ${tractogram} ${sid}__all_cereb_nocx_out_cereb_in_pons.trk ${sid}__all_cereb_nocx_out_cereb_no_pons.trk -f

  if ${params.debug}
  then
    scil_count_streamlines.py ${sid}__all_cereb_nocx_out_cereb_in_pons.trk > ${sid}__all_cereb_nocx_out_cereb_in_pons.txt
    scil_count_streamlines.py ${sid}__all_cereb_nocx_out_cereb_no_pons.trk > ${sid}__all_cereb_nocx_out_cereb_no_pons.txt
  fi
  """
}

process Cerebellum_In_Mid_Brain{
  cpus 1
  tag = "Cereb - In the midbrain"

  input:
    set sid, file(tractogram) from noPons

  output:
    set sid, "${sid}__all_cereb_nocx_out_cereb_in_midbrain.trk" into inMidBrain
    set sid, "${sid}__all_cereb_nocx_out_cereb_no_midbrain.trk" into noMidBrain
    file "${sid}__all_cereb_nocx_out_cereb_in_midbrain.txt" optional true
    file "${sid}__all_cereb_nocx_out_cereb_no_midbrain.txt" optional true


  script:
  """
    scil_filter_tractogram.py ${tractogram} ${sid}__all_cereb_nocx_out_cereb_in_midbrain.trk \
      --drawn_roi ${params.rois_folder}${params.atlas.midbrain} either_end include \
      --drawn_roi ${params.rois_folder}${params.atlas.thal_CP_Medulla} any exclude \
      --drawn_roi ${params.rois_folder}${params.atlas.MCPant} any exclude \
      --drawn_roi ${params.rois_folder}${params.atlas.CGM_SWM_DWM} any exclude -f

    scil_streamlines_math.py difference ${tractogram} ${sid}__all_cereb_nocx_out_cereb_in_midbrain.trk ${sid}__all_cereb_nocx_out_cereb_no_midbrain.trk -f

    if ${params.debug}
    then
      scil_count_streamlines.py ${sid}__all_cereb_nocx_out_cereb_in_midbrain.trk > ${sid}__all_cereb_nocx_out_cereb_in_midbrain.txt
      scil_count_streamlines.py ${sid}__all_cereb_nocx_out_cereb_no_midbrain.trk > ${sid}__all_cereb_nocx_out_cereb_no_midbrain.txt
    fi
  """
}

process Cerebellum_In_Red_Nucleus_Thalamus{
  cpus 1
  tag = "Cereb - One termination in the red nucleus and/or the thalamus"

  input:
    set sid, file(tractogram) from noMidBrain

  output:
    set sid, "${sid}__all_cereb_nocx_out_cereb_in_redN_and_thalamus.trk" into inRedNThal
    set sid, "${sid}__cereb_bin_04.trk" into cereb_bin_04
    file "${sid}__all_cereb_nocx_out_cereb_in_redN_and_thalamus.txt" optional true
    file "${sid}__cereb_bin_04.txt" optional true

  script:
  """
  scil_filter_tractogram.py ${tractogram} ${sid}__all_cereb_nocx_out_cereb_in_redN_and_thalamus.trk \
    --drawn_roi ${params.rois_folder}${params.atlas.BG_IC_CP_MCPant_Medulla} any exclude \
    --drawn_roi ${params.rois_folder}${params.atlas.CGM_SWM_DWM} any exclude -f

  scil_streamlines_math.py difference ${tractogram} ${sid}__all_cereb_nocx_out_cereb_in_redN_and_thalamus.trk ${sid}__cereb_bin_04.trk -f

  if ${params.debug}
  then
    scil_count_streamlines.py ${sid}__all_cereb_nocx_out_cereb_in_redN_and_thalamus.trk > ${sid}__all_cereb_nocx_out_cereb_in_redN_and_thalamus.txt
    scil_count_streamlines.py ${sid}__cereb_bin_04.trk > ${sid}__cereb_bin_04.txt
  fi
  """
}

/*
  Brainstem
*/

process Brainstem_End_In_Brainstem {

  cpus 1
  tag = "Brainstem - Split both_ends and either_end"

  input:
    set sid, file(tractogram) from inBrainstem

  output:
    set sid, "${sid}__all_brainstem_both_ends.trk" into bothEndInBrainstem
    set sid, "${sid}__all_brainstem_either_end.trk" into eitherEndInBrainstem
    file "${sid}__all_brainstem_both_ends.txt" optional true
    file "${sid}__all_brainstem_either_end.txt" optional true

  script:
    atlas=params.rois_folder+params.atlas.brainstem
    mode=params.mode.both_ends
    criteria=params.criteria.include
    out_extension='all_brainstem_both_ends'
    remaining_extension='all_brainstem_either_end'
    basename="${sid}"

    template "filter_with_atlas.sh"
}

process Brainstem_Medulla{

  cpus 1
  tag = "Brainstem - Segmentation Medulla"

  input:
    set sid, file(tractogram) from bothEndInBrainstem

  output:
    set sid, "${sid}__all_brainstem_both_ends_Medulla.trk" into EndInMedulla
    set sid, "${sid}__all_brainstem_both_ends_noMedulla.trk" into NotEndInMedulla
    file "${sid}__all_brainstem_both_ends_Medulla.txt" optional true
    file "${sid}__all_brainstem_both_ends_noMedulla.txt" optional true

  script:
  """
  scil_filter_tractogram.py ${tractogram} ${sid}__all_brainstem_both_ends_Medulla.trk \
    --drawn_roi ${params.rois_folder}${params.atlas.CGM_SWM_subcortical} both_ends exclude \
    --drawn_roi ${params.rois_folder}${params.atlas.midbrain} both_ends exclude \
    --drawn_roi ${params.rois_folder}${params.atlas.CP} both_ends exclude \
    --drawn_roi ${params.rois_folder}${params.atlas.pons} both_ends exclude \
    --drawn_roi ${params.rois_folder}${params.atlas.medulla} both_ends include \
    --drawn_roi ${params.rois_folder}${params.atlas.pons} any exclude -f

  scil_filter_tractogram.py ${tractogram} ${sid}__all_brainstem_both_ends_noMedulla.trk \
    --drawn_roi ${params.rois_folder}${params.atlas.CGM_SWM_subcortical} both_ends exclude \
    --drawn_roi ${params.rois_folder}${params.atlas.midbrain} both_ends exclude \
    --drawn_roi ${params.rois_folder}${params.atlas.CP} both_ends exclude \
    --drawn_roi ${params.rois_folder}${params.atlas.pons} both_ends exclude \
    --drawn_roi ${params.rois_folder}${params.atlas.medulla} both_ends exclude -f
  """
}

process Brainstem_noCP{

  cpus 1
  tag = "Brainstem - No  ends in cereberal peduncles"

  input:
    set sid, file(tractogram) from NotEndInMedulla

  output:
    set sid, "${sid}__all_brainstem_both_ends_noMedulla_noCP.trk" into brainstemNoCP
    set sid, "${sid}__brainstem_bin_02.trk" into brainstemBin02
    file "${sid}__all_brainstem_both_ends_noMedulla_noCP.txt" optional true
    file "${sid}__brainstem_bin_02.txt" optional true

  when:
    params.debug

  script:
    atlas=params.rois_folder+params.atlas.CP
    mode=params.mode.either_end
    criteria=params.criteria.exclude
    out_extension='all_brainstem_both_ends_noMedulla_noCP'
    remaining_extension='brainstem_bin_02'
    basename="${sid}"

    template "filter_with_atlas.sh"
}

/*
  CC Homotopic
*/

cc_homotopic_pairs=params.cc_homotopic_pairs?.tokenize(',')

process CC_Homotopic {
  cpus 1
  tag = "CC Homotopic"

  input:
    set sid, file(tractogram) from CC_for_homotopic
    each pair from cc_homotopic_pairs

  output:
    set sid, "${sid}__cc_homotopic_${pair}.trk" into CC_Homotopic

  script:
    File ori_file = new File(params.filtering_lists_folder+"CC_homo_${pair}_filtering_list_f.txt")
    File new_file = File.createTempFile("tmp_filtering_list",".tmp")
    for (line in ori_file.readLines()) {
      new_file << line.replace("drawn_roi /JHU_template_GIN_dil/", "drawn_roi "+params.rois_folder+"/")+"\n"}

    filtering_list=new_file.absolutePath
    out_extension="cc_homotopic_${pair}"
    remaining_extension="garbage_${pair}"
    basename="${sid}"

    template "filter_with_list.sh"
}

/*
  ASSO
*/

asso_lists=params.asso_lists?.tokenize(',')

process asso_part1 {
  cpus 1
  tag = "Asso filtering to get asso_SLS->SLF+AF and asso_ILS->IFOF+UF"

  input:
    set sid, val(side), file(tractogram) from assoAllIntraInter_for_filtering
    each asso_list from asso_lists

  output:
    set sid, val(side), "${sid}__asso_${asso_list}_${side}.trk" into assoAllIntraInter_filtered
    set sid, "${sid}__asso_lost_${asso_list}_${side}.trk" into assoLost3
    file "${sid}__asso_${asso_list}_${side}.txt" optional true
    file "${sid}__asso_lost_${asso_list}_${side}.txt" optional true

  script:

  File ori_file = new File(params.filtering_lists_folder+"ASSO_${asso_list}_${side}_filtering_list.txt")
  File new_file = File.createTempFile("tmp_filtering_list",".tmp")
  for (line in ori_file.readLines()) {
    new_file << line.replace("drawn_roi /JHU_template_GIN_dil/", "drawn_roi "+params.rois_folder+"/")+"\n"}


  filtering_list=new_file.absolutePath
  out_extension="asso_${asso_list}_${side}"
  remaining_extension="asso_lost_${asso_list}_${side}"
  basename="${sid}"

  template "filter_with_list.sh"
}
