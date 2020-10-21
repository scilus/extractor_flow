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

in_tractogram.into{for_remove_out_JHU; in_tractogram_for_get_unplausible}

sides = params.sides?.tokenize(',')
Channel.from(sides).into{sides_ipsi;
                         sides_split_CC_BG;
                         sides_split_asso_in_hemi;
                         sides_split_BG_Thal;
                         sides_split_BG_Put}

process Remove_Out_JHU {
    cpus 1
    tag = "Remove out of JHU"

    input:
      set sid, file(tractogram) from for_remove_out_JHU

    output:
      set sid, "${sid}__wb_in_JHU.trk" into wb_for_rm_crossing_gyri
      file "${sid}__wb_in_JHU.txt" optional true
      set sid, "${sid}__wb_out_JHU.trk" optional true
      file "${sid}__wb_out_JHU.txt" optional true

    script:
    atlas=params.rois_folder+params.atlas.JHU
    mode=params.mode.both_ends
    criteria=params.criteria.include
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

  script:

    """
    scil_detect_streamlines_loops.py ${wb_min20} ${sid}__wb_min20_noloop.trk \
                                     -a ${params.loop_angle_threshold}  \
                                     --looping_tractogram ${sid}__wb_loops.trk \
                                     --display_counts > ${sid}__wb_min20_noloop.txt \
                                     -f
    """
}

process Removing_End_In_CC_DWM {
  cpus 1
  tag = "Remove end in CC DWM"

  input:
    set sid, file(wb_min20_noloop) from wb_for_rm_end_in_cc_dwm

  output:
    set sid, "${sid}__wb_clean01.trk" into wb_for_extract_end_in_cerebellum, for_extract_unplausible
    file "${sid}__wb_clean01.txt" optional true

  script:
  """
  scil_filter_tractogram.py ${wb_min20_noloop} ${sid}__wb_clean01.trk \
                    --drawn_roi ${params.rois_folder}${params.atlas.cc} either_end exclude \
                    --drawn_roi ${params.rois_folder}${params.atlas.dwm} either_end exclude \
                    -f --display_count > ${sid}__wb_clean01.txt
  """
}

in_tractogram_for_get_unplausible.join(for_extract_unplausible).set{unplausible_streamlines}

process extract_unplausible{
  cpus 1
  tag = 'Extract unplausible'

  input:
    set sid, file(tractogram1), file(tractogram2) from unplausible_streamlines

  output:
    set sid, "${sid}_unplausible_streamlines.trk"

  script:
  """
  scil_streamlines_math.py difference ${tractogram1} \
                                      ${tractogram2} \
                                      ${sid}_unplausible_streamlines.trk;
  """
}

process extract_end_in_cerebellum {
  cpus 1
  tag = "Extract either end in cerebellum"

  input:
    set sid, file(tractogram) from wb_for_extract_end_in_cerebellum

  output:
    set sid, "${sid}__wb_clean01_nocereb.trk" into wb_for_extract_end_in_brainstem
    set sid, "${sid}__all_cerebellum.trk" into cerebellum_for_rm_cortex_gm
    file "${sid}__all_cerebellum.txt" optional true
    file "${sid}__wb_clean01_nocereb.txt" optional true

  script:
  atlas=params.rois_folder+params.atlas.cerebellum
  mode=params.mode.either_end
  criteria=params.criteria.exclude
  out_extension='wb_clean01_nocereb'
  remaining_extension='all_cerebellum'
  basename="${sid}"

  template "filter_with_atlas.sh"
}

/*
  Cerebellum
*/

process Cerebellum_Remove_Cortex_GM {

  cpus 1
  tag = "Cereb - Remove End In Cortex GM"

  input:
    set sid, file(tractogram) from cerebellum_for_rm_cortex_gm

  output:
    set sid, "${sid}__all_cereb_nocx.trk" into cerebellum_for_end_in_GWM
    set sid, "${sid}__cereb_bin_01.trk" optional true into cerebellum_bin_01
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

process Cerebellum_Remove_Not_End_In_GWM {
  cpus 1
  tag = "Cereb - Remove not end in cerebellum GM"

  input:
    set sid, file(tractogram) from cerebellum_for_end_in_GWM

  output:
    set sid, "${sid}__all_cereb_nocx_f.trk" into cerebellum_for_split_both_end_in_cerebellum
    set sid, "${sid}__cereb_bin_02.trk" optional true into cerebellum_bin_02
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
  tag = "Cereb - Split - Both end in cerebellum"

  input:
    set sid, file(tractogram) from cerebellum_for_split_both_end_in_cerebellum

  output:
    set sid, "${sid}__all_cereb_nocx_in_cereb.trk" into cerebellum_for_remove_loops
    set sid, "${sid}__all_cereb_nocx_out_cereb.trk" into cerebellum_for_in_medulla
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
    set sid, file(tractogram) from cerebellum_for_remove_loops

  output:
    set sid, "${sid}__all_cereb_nocx_in_cereb_f.trk" into cerebellum_for_merge_plausible
    set sid, "${sid}__cereb_bin_03.trk" into cerebellum_bin_03
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
    set sid, file(tractogram) from cerebellum_for_in_medulla

  output:
    set sid, "${sid}__all_cereb_nocx_out_cereb_in_medulla.trk" into cerebellum_in_medulla_for_merge_plausible
    set sid, "${sid}__all_cereb_nocx_out_cereb_no_medulla.trk" into cerebellum_for_in_pons
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
    set sid, file(tractogram) from cerebellum_for_in_pons

  output:
    set sid, "${sid}__all_cereb_nocx_out_cereb_in_pons.trk" into cerebellum_in_pons_for_merge_plausible
    set sid, "${sid}__all_cereb_nocx_out_cereb_no_pons.trk" into cerebellum_for_in_mid_brain
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
    set sid, file(tractogram) from cerebellum_for_in_mid_brain

  output:
      set sid, "${sid}__all_cereb_nocx_out_cereb_in_midbrain.trk" into cerebellum_in_midbrain_for_merge_plausible
    set sid, "${sid}__all_cereb_nocx_out_cereb_no_midbrain.trk" into cerebellum_for_in_red_nuclei
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
    set sid, file(tractogram) from cerebellum_for_in_red_nuclei

  output:
    set sid, "${sid}__all_cereb_nocx_out_cereb_in_redN_and_thalamus.trk" into cerebellum_in_redN_for_merge_plausible
    set sid, "${sid}__cereb_bin_04.trk" into cerebellum_bin_04
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

cerebellum_for_merge_plausible.join(cerebellum_in_medulla_for_merge_plausible).join(cerebellum_in_midbrain_for_merge_plausible).join(cerebellum_in_pons_for_merge_plausible).join(cerebellum_in_redN_for_merge_plausible).set{cerebellum_plausible}
process Cerebellum_merge_plausible{
 cpus 1
 tag = "Merge plausible Cerebellum"

 input:
   set sid, file(trk01), file(trk02), file(trk03), file(trk04), file(trk05) from cerebellum_plausible

 output:
  set sid, "${sid}__cerebellum_plausible.trk"
  file "${sid}__all_cereb_nocx_out_cereb_in_redN_and_thalamus.txt" optional true

 script:
 """
 scil_streamlines_math.py lazy_concatenate ${trk01} \
  ${trk02} \
  ${trk03} \
  ${trk04} \
  ${trk05} \
  ${sid}__cerebellum_plausible.trk -f

  if ${params.debug}
  then
    scil_count_streamlines.py ${sid}__cerebellum_plausible.trk > ${sid}__cerebellum_plausible.txt
  fi
 """
}

cerebellum_bin_01.join(cerebellum_bin_02).join(cerebellum_bin_03).join(cerebellum_bin_04).set{cerebellum_unplausible}
process Cerebellum_merge_unplausible {
 cpus 1
 tag = "Merge unplausible Cerebellum"

 input:
  set sid, file(trk01), file(trk02), file(trk03), file(trk04) from cerebellum_unplausible

 output:
  set sid, "${sid}__cerebellum_unplausible.trk"
  file "${sid}__cerebellum_unplausible.txt" optional true

 script:
 """
 scil_streamlines_math.py lazy_concatenate ${trk01} \
   ${trk02} \
   ${trk03} \
   ${trk04} ${sid}__cerebellum_unplausible.trk -f

   if ${params.debug}
   then
    scil_count_streamlines.py ${sid}__cerebellum_unplausible.trk > ${sid}__cerebellum_unplausible.txt
  fi
 """
}

/*
  END Cerebellum
*/


process Split_Either_End_In_Brainstem {
  cpus 1
  tag = "Split either end in brainstem"

  input:
    set sid, file(tractogram) from wb_for_extract_end_in_brainstem

  output:
    set sid, "${sid}__wb_clean02.trk" into wb_for_split_end_in_CGMSWI
    set sid, "${sid}__all_brainstem.trk" into brainstem_for_both_end_in_brainstem
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

/*
  Brainstem
*/

process Brainstem_End_In_Brainstem {

  cpus 1
  tag = "Brainstem - Split both_ends and either_end"

  input:
    set sid, file(tractogram) from brainstem_for_both_end_in_brainstem

  output:
    set sid, "${sid}__all_brainstem_both_ends.trk" into brainstem_for_medulla_segmentation
    set sid, "${sid}__all_brainstem_either_end.trk" into brainstem_for_through_cc
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

/*
  BOTH END Brainstem
*/

process brainstem_be_medulla_segmentation{

  cpus 1
  tag = "Brainstem - Both End - Segmentation Medulla"

  input:
    set sid, file(tractogram) from brainstem_for_medulla_segmentation

  output:
    set sid, "${sid}__all_brainstem_both_ends_Medulla.trk" into brainstem_be_medulla_for_plausible
    set sid, "${sid}__all_brainstem_both_ends_no_medulla.trk" into brainstem_for_noCP
    file "${sid}__all_brainstem_both_ends_Medulla.txt" optional true
    file "${sid}__all_brainstem_both_ends_no_medulla.txt" optional true

  script:
  """
  scil_filter_tractogram.py ${tractogram} ${sid}__all_brainstem_both_ends_Medulla.trk \
    --drawn_roi ${params.rois_folder}${params.atlas.CGM_SWM_subcortical} both_ends exclude \
    --drawn_roi ${params.rois_folder}${params.atlas.midbrain} both_ends exclude \
    --drawn_roi ${params.rois_folder}${params.atlas.CP} both_ends exclude \
    --drawn_roi ${params.rois_folder}${params.atlas.pons} both_ends exclude \
    --drawn_roi ${params.rois_folder}${params.atlas.medulla} both_ends include \
    --drawn_roi ${params.rois_folder}${params.atlas.pons} any exclude -f

  scil_filter_tractogram.py ${tractogram} ${sid}__all_brainstem_both_ends_no_medulla.trk \
    --drawn_roi ${params.rois_folder}${params.atlas.CGM_SWM_subcortical} both_ends exclude \
    --drawn_roi ${params.rois_folder}${params.atlas.midbrain} both_ends exclude \
    --drawn_roi ${params.rois_folder}${params.atlas.CP} both_ends exclude \
    --drawn_roi ${params.rois_folder}${params.atlas.pons} both_ends exclude \
    --drawn_roi ${params.rois_folder}${params.atlas.medulla} both_ends exclude -f
  """
}

process brainstem_be_noCP{

  cpus 1
  tag = "Brainstem - Both End - No ends in cereberal peduncles"

  input:
    set sid, file(tractogram) from brainstem_for_noCP

  output:
    set sid, "${sid}__all_brainstem_both_ends_no_medulla_noCP.trk" into brainstem_for_split_midbrain_pons
    set sid, "${sid}__brainstem_bin_02.trk" into brainstem_bin_03
    file "${sid}__all_brainstem_both_ends_no_medulla_noCP.txt" optional true
    file "${sid}__brainstem_bin_02.txt" optional true

  script:
    atlas=params.rois_folder+params.atlas.CP
    mode=params.mode.either_end
    criteria=params.criteria.exclude
    out_extension='all_brainstem_both_ends_no_medulla_noCP'
    remaining_extension='brainstem_bin_02'
    basename="${sid}"

    template "filter_with_atlas.sh"
}

process brainstem_be_split_midbrain_pons{

  cpus 1
  tag = "Brainstem - Both End - Split MidBrain and Pons"

  input:
    set sid, file(tractogram) from brainstem_for_split_midbrain_pons

  output:
    set sid, "${sid}__all_brainstem_both_ends_noMedulla_noCP_Midbrain.trk" into brainstem_be_midbrain_for_plausible
    set sid, "${sid}__all_brainstem_both_ends_noMedulla_noCP_Pons.trk" into brainstem_be_pons_for_plausible
    file "${sid}__all_brainstem_both_ends_noMedulla_noCP_Midbrain.txt" optional true
    file "${sid}__all_brainstem_both_ends_noMedulla_noCP_Pons.txt" optional true

  script:
    atlas=params.rois_folder+params.atlas.midbrain
    mode=params.mode.either_end
    criteria=params.criteria.include
    out_extension='all_brainstem_both_ends_noMedulla_noCP_Midbrain'
    remaining_extension='all_brainstem_both_ends_noMedulla_noCP_Pons'
    basename="${sid}"

    template "filter_with_atlas.sh"
}

brainstem_be_medulla_for_plausible.join(brainstem_be_midbrain_for_plausible).join(brainstem_be_pons_for_plausible).set{plausible_brainstem}

process brainstem_be_merge_plausible{
  cpus 1
  tag = "Brainstem - Both End - Concatenating plausible streamlines Brainstem"

  input:
  set sid, file(trk01), file(trk02), file(trk03) from plausible_brainstem

  output:
  set sid, "${sid}__brainstem_be_plausible_tracks.trk" into brainstem_merge_plausible_01
  file "${sid}__brainstem_be_plausible_tracks.txt" optional true

  script:
  """
  scil_streamlines_math.py concatenate ${trk01} \
                                       ${trk02} \
                                       ${trk03} \
                                       ${sid}__brainstem_be_plausible_tracks.trk -f
  if ${params.debug}
  then
    scil_count_streamlines.py ${sid}__brainstem_be_plausible_tracks.trk > ${sid}__brainstem_be_plausible_tracks.txt
  fi
  """
}

/*
  END BOTH END Brainstem
*/

/*
  EITHER END Brainstem
*/

process brainstem_ee_through_CC{
  cpus 1
  tag = "Brainstem - either end - through CC"

  input:
    set sid, file(tractogram) from brainstem_for_through_cc

  output:
    set sid, "${sid}__all_brainstem_either_end_CC.trk" into brainstem_for_ipsi
    set sid, "${sid}__all_brainstem_either_end_noCC.trk" into brainstem_for_merge_ipsi_noCC_01
    file "${sid}__all_brainstem_either_end_CC.txt" optional true
    file "${sid}__all_brainstem_either_end_noCC.txt" optional true

  script:
    atlas=params.rois_folder+params.atlas.cc
    mode=params.mode.any
    criteria=params.criteria.include
    out_extension='all_brainstem_either_end_CC'
    remaining_extension='all_brainstem_either_end_noCC'
    basename="${sid}"

  template "filter_with_atlas.sh"
}

process brainstem_ee_ipsi{
  cpus 1
  tag = "Brainstem - either end - IPSI"

  input:
    set sid, file(tractogram) from brainstem_for_ipsi
    each side from sides_ipsi

  output:
    set sid, "${sid}__all_brainstem_either_end_CC_ipsi_${side}.trk" into ipsi

  script:
  """
  if [["$side"=='L']];
  then
    roi1=${params.rois_folder}${params.atlas.CGM_SWM_side}_L.nii.gz
    roi2=${params.rois_folder}${params.atlas.ccR}
  else
    roi1=${params.rois_folder}${params.atlas.CGM_SWM_side}_R.nii.gz
    roi2=${params.rois_folder}${params.atlas.ccL}
  fi

  scil_filter_tractogram.py ${tractogram} ${sid}__all_brainstem_either_end_CC_ipsi_${side}.trk \
    --drawn_roi \${roi1} either_end include \
    --drawn_roi \${roi2} any exclude -f;
  """
}

ipsi.groupTuple().map{it.flatten().toList()}.set{brainstem_for_merge_ipsi_noCC_02}
brainstem_for_merge_ipsi_noCC_01.join(brainstem_for_merge_ipsi_noCC_02).set{brainstem_for_merge_ipsi_noCC}

process brainstem_ee_merge_ipsi_and_noCC{
  cpus 1
  tags = 'Brainstem - Merge IPSI and either end noCC'

  input:
  set sid, file(trk01), file(trk02), file(trk03) from brainstem_for_merge_ipsi_noCC

  output:
  set sid, "${sid}__all_brainstem_either_end_tmp.trk" into brainstem_for_noBG_noCP_noThal

  script:
  """
  scil_streamlines_math.py concatenate ${trk01} \
                                       ${trk02} \
                                       ${trk03} \
                                       ${sid}__all_brainstem_either_end_tmp.trk
  """
}

process brainstem_ee_noBG_noCP_noThal{
  cpus 1
  tags = 'Brainstem - either end - remove BG, CP, Thal'

  input:
  set sid, file(tractogram) from brainstem_for_noBG_noCP_noThal

  output:
  set sid, "${sid}__all_brainstem_either_end_tmp_noBG_noCP_noThal.trk" into brainstem_for_ee_cgmswm, brainstem_for_red_nucleus
  set sid, "${sid}__all_brainstem_either_tmp_Thal.trk" into brainstem_for_thalamus

  script:
  """
  scil_filter_tractogram.py ${tractogram} ${sid}__all_brainstem_either_end_tmp_noBG_noCP_noThal.trk \
    --drawn_roi ${params.rois_folder}${params.atlas.Put_Cd_GP_Amyg_SNr} ${params.mode.either_end} exclude  \
    --drawn_roi ${params.rois_folder}${params.atlas.CP} ${params.mode.either_end} exclude  \
    --drawn_roi ${params.rois_folder}${params.atlas.allThal} ${params.mode.either_end} exclude -f

  scil_filter_tractogram.py ${tractogram} ${sid}__all_brainstem_either_tmp_Thal.trk \
    --drawn_roi ${params.rois_folder}${params.atlas.Put_Cd_GP_Amyg_SNr} ${params.mode.either_end} exclude  \
    --drawn_roi ${params.rois_folder}${params.atlas.CP} ${params.mode.either_end} exclude  \
    --drawn_roi ${params.rois_folder}${params.atlas.allThal} ${params.mode.either_end} include -f
  """
}

process brainstem_ee_thalamus{
  cpus 1
  tags = 'Brainstem - either end - Thalamus'

  input:
  set sid, file(tractogram) from brainstem_for_thalamus

  output:
  set sid, "${sid}__all_brainstem_thalamus.trk" into brainstem_ee_all_thalamus_for_plausible

  script:
  """
  scil_filter_tractogram.py ${tractogram} ${sid}__tmp.trk \
    --drawn_roi ${params.rois_folder}${params.atlas.CP} ${params.mode.any} exclude -f

  scil_detect_streamlines_loops.py ${sid}__tmp.trk ${sid}__all_brainstem_thalamus.trk -a 270 -f
  """
}

process brainstem_ee_end_CGMSWM{
  cpus 1
  tag = 'Brainstem - either end - Either end In CGM-SWM'

  input:
    set sid, file(tractogram) from brainstem_for_ee_cgmswm

  output:
    set sid, "${sid}__all_brainstem_either_end_tmp_CGM_SWM.trk" into brainstem_for_split_thal

  script:
  """
    scil_filter_tractogram.py ${tractogram} ${sid}__all_brainstem_either_end_tmp_CGM_SWM.trk \
      --drawn_roi ${params.rois_folder}${params.atlas.CGM_SWM} ${params.mode.either_end} include -f
  """
}

process brainstem_ee_red_nucleus{
  cpus 1
  tag = 'Brainstem - Either end - End in Red Nucleus'

  input:
    set sid, file(tractogram) from brainstem_for_red_nucleus

  output:
    set sid, "${sid}__all_brainstem_red_nucleus.trk" into brainstem_ee_red_nucleus_for_plausible

  script:
  """
    scil_filter_tractogram.py ${tractogram} ${sid}__all_brainstem_either_end_tmp_RedN.trk \
      --drawn_roi ${params.rois_folder}${params.atlas.CGM_SWM} ${params.mode.either_end} exclude -f
    scil_outlier_rejection.py ${sid}__all_brainstem_either_end_tmp_RedN.trk \
                              ${sid}__all_brainstem_red_nucleus.trk \
                              --remaining_bundle ${sid}__all_brainstem_either_end_bin_05.trk \
                              --alpha 0.4 -f
  """
}

process brainstem_ee_end_CGMSWM_split_thal{
  cpus 1
  tag = 'Brainstem - End in CGMSWI no Thal'

  input:
    set sid, file(tractogram) from brainstem_for_split_thal

  output:
    set sid, "${sid}__all_brainstem_either_end_tmp_CGM_SWM_tmp_noThal.trk" into brainstem_for_merge_and_split_pons_01
    set sid, "${sid}__all_brainstem_either_end_tmp_CGM_SWM_tmp_Thal.trk" into brainstem_for_split_IC


  script:
  """
    scil_filter_tractogram.py ${tractogram} ${sid}__all_brainstem_either_end_tmp_CGM_SWM_tmp_noThal.trk \
      --drawn_roi ${params.rois_folder}${params.atlas.allThal} ${params.mode.any} exclude -f
    scil_streamlines_math.py difference ${tractogram} ${sid}__all_brainstem_either_end_tmp_CGM_SWM_tmp_noThal.trk \
      ${sid}__all_brainstem_either_end_tmp_CGM_SWM_tmp_Thal.trk -f
  """
}

process brainstem_ee_end_CGMSWM_Thal_split_IC{
  cpus 1
  tag = 'Brainstem - Either End - End in CGMSWI Thal split IC'

  input:
    set sid, file(tractogram) from brainstem_for_split_IC

  output:
    set sid, "${sid}__all_brainstem_either_end_tmp_CGM_SWM_tmp_Thal_IC.trk" into brainstem_for_split_CP
    set sid, "${sid}__all_brainstem_either_end_tmp_CGM_SWM_tmp_Thal_noIC.trk" into brainstem_for_3rd_order_PoC_01

  script:
  """
    scil_filter_tractogram.py ${tractogram} ${sid}__all_brainstem_either_end_tmp_CGM_SWM_tmp_Thal_IC.trk \
      --drawn_roi ${params.rois_folder}${params.atlas.ic} ${params.mode.any} include -f
    scil_streamlines_math.py difference ${tractogram} ${sid}__all_brainstem_either_end_tmp_CGM_SWM_tmp_Thal_IC.trk \
      ${sid}__all_brainstem_either_end_tmp_CGM_SWM_tmp_Thal_noIC.trk -f
  """
}

process brainstem_ee_end_CGMSWM_Thal_IC_split_CP{
  cpus 1
  tag = 'Brainstem - Either End - End in CGMSWI Thal IC split CP'

  input:
    set sid, file(tractogram) from brainstem_for_split_CP

  output:
    set sid, "${sid}__all_brainstem_either_end_tmp_CGM_SWM_tmp_Thal_IC_CP.trk" into brainstem_for_merge_and_split_pons_02
    set sid, "${sid}__all_brainstem_either_end_tmp_CGM_SWM_tmp_Thal_IC_noCP.trk" into brainstem_for_3rd_order_PoC_02

  script:
  """
    scil_filter_tractogram.py ${tractogram} ${sid}__all_brainstem_either_end_tmp_CGM_SWM_tmp_Thal_IC_CP.trk \
      --drawn_roi ${params.rois_folder}${params.atlas.CP} ${params.mode.any} include -f
    scil_streamlines_math.py difference ${tractogram} ${sid}__all_brainstem_either_end_tmp_CGM_SWM_tmp_Thal_IC_CP.trk \
      ${sid}__all_brainstem_either_end_tmp_CGM_SWM_tmp_Thal_IC_noCP.trk -f
  """
}

brainstem_for_3rd_order_PoC_01.join(brainstem_for_3rd_order_PoC_02).set{brainstem_for_3rd_order_PoC}

process brainstem_ee_3rd_order_PoC{
  cpus 1
  tag = 'Brainstem - Either End - 3rd_order_PoC'

  input:
    set sid, file(trk01), file(trk02) from brainstem_for_3rd_order_PoC

  output:
    set sid, "${sid}__all_brainstem_3rd_order_PoC.trk" into brainstem_ee_3rd_order_PoC_for_plausible

  script:
  """
    scil_streamlines_math.py lazy_concatenate ${trk01} ${trk02} ${sid}__tmp1.trk
    scil_filter_tractogram.py ${sid}__tmp1.trk ${sid}__tmp2.trk \
      --drawn_roi ${params.rois_folder}${params.atlas.PoCGWM_all} ${params.mode.either_end} include \
      --drawn_roi ${params.rois_folder}${params.atlas.medulla} ${params.mode.either_end} include -f
    scil_outlier_rejection.py ${sid}__tmp2.trk ${sid}__all_brainstem_3rd_order_PoC.trk --alpha 0.4
  """
}

brainstem_for_merge_and_split_pons_01.join(brainstem_for_merge_and_split_pons_02).set{brainstem_for_merge_and_split_pons}
process brainstem_ee_merge_split_pons{
  cpus 1
  tag = 'Brainstem - Either End - Merge noThal and Thal CP IC - Split Pons'

  input:
    set sid, file(trk01), file(trk02) from brainstem_for_merge_and_split_pons

  output:
    set sid, "${sid}__all_brainstem_either_end_tmp_CGM_SWM_tmp_noThal_eePons.trk" into brainstem_for_corticopontine_frontal
    set sid, "${sid}__all_brainstem_either_end_tmp_CGM_SWM_tmp_noThal_no_eePons.trk" into brainstem_for_merge_not_ee_frontal_not_ee_pons_02

  script:
  """
    scil_streamlines_math.py concatenate ${trk01} ${trk02} ${sid}__tmp01.trk
    scil_filter_tractogram.py ${sid}__tmp01.trk ${sid}__all_brainstem_either_end_tmp_CGM_SWM_tmp_noThal_eePons.trk \
      --drawn_roi ${params.rois_folder}${params.atlas.pons} ${params.mode.any} include -f
    scil_streamlines_math.py difference ${sid}__tmp01.trk ${sid}__all_brainstem_either_end_tmp_CGM_SWM_tmp_noThal_eePons.trk \
        ${sid}__all_brainstem_either_end_tmp_CGM_SWM_tmp_noThal_no_eePons.trk -f
  """
}

process brainstem_ee_corticopontic_frontal{
  cpus 1
  tag = 'Brainstem - End in CGMSWI Thal IC noCP'

  input:
    set sid, file(tractogram) from brainstem_for_corticopontine_frontal

  output:
    set sid, "${sid}__all_brainstem_corticopontine_frontal.trk" into brainstem_corticopontine_frontal_for_plausible
    set sid, "${sid}__all_brainstem_ee_tmp_CGM_SWM_tmp_noThal_ee_Pons_no_ee_frontal.trk" into brainstem_for_merge_not_ee_frontal_not_ee_pons_01

  script:
  """
    scil_filter_tractogram.py ${tractogram} ${sid}__tmp01.trk \
      --drawn_roi ${params.rois_folder}${params.atlas.frontal} ${params.mode.either_end} include -f
    scil_streamlines_math.py difference ${tractogram} ${sid}__tmp01.trk \
     ${sid}__all_brainstem_ee_tmp_CGM_SWM_tmp_noThal_ee_Pons_no_ee_frontal.trk -f

    scil_detect_streamlines_loops.py ${sid}__tmp01.trk ${sid}__tmp02.trk -a 240 -f
    scil_outlier_rejection.py ${sid}__tmp01.trk ${sid}__all_brainstem_corticopontine_frontal.trk --alpha 0.45 -f
  """
}

brainstem_for_merge_not_ee_frontal_not_ee_pons_01.join(brainstem_for_merge_not_ee_frontal_not_ee_pons_02).set{brainstem_for_merge_not_ee_frontal_not_ee_pons_split_pons}
process brainstem_ee_merge_not_ee_frontal_not_ee_pons_split_pons {
  cpus 1
  tag = 'Brainstem - Either End - Merge noThal and Thal CP IC - Split Pons'

  input:
   set sid, file(trk01), file(trk02) from brainstem_for_merge_not_ee_frontal_not_ee_pons_split_pons

  output:
   set sid, "${sid}__all_brainstem_either_end_tmp_CGM_SWM_tmp_noThal_ee_Pons.trk" into brainstem_for_corticopontine_parietotemporoccipital
   set sid, "${sid}__all_brainstem_either_end_tmp_CGM_SWM_tmp_noThal_no_ee_Pons.trk" into brainstem_for_merge_ee_pons_no_ee_parietooccipital_01

 script:
 """
  scil_streamlines_math.py concatenate ${trk01} ${trk02} ${sid}__tmp01.trk

  scil_filter_tractogram.py ${sid}__tmp01.trk ${sid}__all_brainstem_either_end_tmp_CGM_SWM_tmp_noThal_ee_Pons.trk \
    --drawn_roi ${params.rois_folder}${params.atlas.pons} ${params.mode.either_end} include -f

  scil_streamlines_math.py difference ${sid}__tmp01.trk ${sid}__all_brainstem_either_end_tmp_CGM_SWM_tmp_noThal_ee_Pons.trk \
      ${sid}__all_brainstem_either_end_tmp_CGM_SWM_tmp_noThal_no_ee_Pons.trk -f
 """
}

process brainstem_ee_corticopontic_parietotemporooccipital{
  cpus 1
  tag = 'Brainstem - End in CGMSWI Thal IC noCP'

  input:
    set sid, file(tractogram) from brainstem_for_corticopontine_parietotemporoccipital

  output:
    set sid, "${sid}__all_brainstem_corticopontine_parietotemporooccipital.trk" into brainstem_ee_corticopontine_parietotemporooccipital_for_plausible
    set sid, "${sid}__all_brainstem_not_ee_parieto_occipital.trk" into brainstem_for_merge_ee_pons_no_ee_parietooccipital_02

  script:
  """
    scil_filter_tractogram.py ${tractogram} ${sid}__tmp_01.trk \
      --drawn_roi ${params.rois_folder}${params.atlas.parietotemporooccipital} ${params.mode.either_end} include -f
    scil_streamlines_math.py difference ${tractogram} ${sid}__tmp_01.trk \
      ${sid}__all_brainstem_not_ee_parieto_occipital.trk -f
    scil_detect_streamlines_loops.py ${sid}__tmp_01.trk ${sid}__tmp_02.trk -a 240 -f
    scil_outlier_rejection.py ${sid}__tmp_01.trk "${sid}__all_brainstem_corticopontine_parietotemporooccipital.trk" --alpha 0.45 -f
  """
}

brainstem_for_merge_ee_pons_no_ee_parietooccipital_01.join(brainstem_for_merge_ee_pons_no_ee_parietooccipital_02).set{brainstem_for_merge_ee_pons_no_ee_parietooccipital}
process brainstem_ee_merge_no_ee_pons_no_ee_parieto_temporo_occipital_split_fronto_parietal {
 cpus 1
 tag = 'Brainstem - Either End - '

 input:
  set sid, file(trk01), file(trk02) from brainstem_for_merge_ee_pons_no_ee_parietooccipital

 output:
  set sid, "${sid}__all_brainstem_no_ee_fronto_parietal.trk" into brainstem_for_corticotectal_01
  set sid, "${sid}__all_brainstem_ee_fronto_parietal.trk" into brainstem_for_pyramidal

 script:
 """
  scil_streamlines_math.py concatenate ${trk01} ${trk02} ${sid}__tmp_01.trk -f
  scil_filter_tractogram.py ${sid}__tmp_01.trk  ${sid}__all_brainstem_ee_fronto_parietal.trk \
    --drawn_roi ${params.rois_folder}${params.atlas.frontoparietal} ${params.mode.either_end} include -f
  scil_streamlines_math.py difference ${sid}__tmp_01.trk ${sid}__all_brainstem_ee_fronto_parietal.trk \
    ${sid}__all_brainstem_no_ee_fronto_parietal.trk -f
 """
}

process brainstem_ee_pyramidal{
  cpus 1

  tag = 'Brainstem - End in FrontoParietal'

  input:
    set sid, file(tractogram) from brainstem_for_pyramidal

  output:
    set sid, "${sid}__all_brainstem_pyramidal.trk" into brainstem_pyramidal_for_merge
    set sid, "${sid}__all_brainstem_not_ee_medulla.trk" into brainstem_for_corticotectal_02

  script:
  """
    scil_filter_tractogram.py ${tractogram} ${sid}__tmp_01.trk \
      --drawn_roi ${params.rois_folder}${params.atlas.medulla} ${params.mode.either_end} include -f

    scil_streamlines_math.py difference ${tractogram} ${sid}__tmp_01.trk \
      ${sid}__all_brainstem_not_ee_medulla.trk -f

    scil_detect_streamlines_loops.py ${sid}__tmp_01.trk ${sid}__tmp_02.trk -a 240
    scil_outlier_rejection.py ${sid}__tmp_02.trk  ${sid}__all_brainstem_pyramidal.trk --alpha 0.6 -f
  """
}

brainstem_for_corticotectal_01.join(brainstem_for_corticotectal_02).set{brainstem_for_corticotectal}
process brainstem_ee_corticotectal{
  cpus 1
  tag = 'Brainstem - Corticotectal'

  input:
    set sid, file(trk01), file(trk02) from brainstem_for_corticotectal

  output:
    set sid, "${sid}__all_brainstem_corticotectal.trk" into brainstem_corticotectal_for_plausible

  script:
  """
    scil_streamlines_math.py concatenate ${trk01} ${trk02} ${sid}__tmp_01.trk
    scil_filter_tractogram.py ${sid}__tmp_01.trk ${sid}__tmp_02.trk \
      --drawn_roi ${params.rois_folder}${params.atlas.midbrain} ${params.mode.either_end} include -f
    scil_detect_streamlines_loops.py ${sid}__tmp_02.trk ${sid}__tmp_noloop.trk -a 240
    scil_outlier_rejection.py ${sid}__tmp_noloop.trk ${sid}__all_brainstem_corticotectal.trk --alpha 0.5
  """
}

brainstem_ee_all_thalamus_for_plausible.join(brainstem_ee_red_nucleus_for_plausible).join(brainstem_ee_3rd_order_PoC_for_plausible).join(brainstem_corticopontine_frontal_for_plausible).join(brainstem_pyramidal_for_merge).join(brainstem_ee_corticopontine_parietotemporooccipital_for_plausible).join(brainstem_corticotectal_for_plausible).set{brainstem_ee_merge_plausible}
process brainstem_ee_merge_plausible{
  cpus 1
  tag = "Brainstem - Either End - Concatenating plausible streamlines Brainstem"

  input:
  set sid, file(trk01), file(trk02), file(trk03), file(trk04), file(trk05), file(trk06), file(trk07) from brainstem_ee_merge_plausible

  output:
  set sid, "${sid}__brainstem_either_ends_plausible_tracks.trk" into brainstem_merge_plausible_02
  file "${sid}__brainstem_either_ends_plausible_tracks.txt" optional true

  script:
  """
  scil_streamlines_math.py concatenate ${trk01} \
                                       ${trk02} \
                                       ${trk03} \
                                       ${trk04} \
                                       ${trk05} \
                                       ${trk06} \
                                       ${trk07} \
                                       ${sid}__brainstem_either_ends_plausible_tracks.trk -f
  if ${params.debug}
  then
    scil_count_streamlines.py ${sid}__brainstem_either_ends_plausible_tracks.trk > ${sid}__brainstem_either_ends_plausible_tracks.txt
  fi
  """
}

/*
  END EITHER END Brainstem
*/

brainstem_merge_plausible_01.join(brainstem_merge_plausible_02).set{brainstem_merge_plausible}
process brainstem_merge_plausible{
  cpus 1
  tag = "Brainstem - Concatenating plausible streamlines"

  input:
  set sid, file(trk01), file(trk02) from brainstem_merge_plausible

  output:
  set sid, "${sid}__brainstem_plausible_tracks.trk"
  file "${sid}__brainstem_plausible_tracks.txt" optional true

  script:
  """
  scil_streamlines_math.py concatenate ${trk01} \
                                       ${trk02} \
                                       ${sid}__brainstem_plausible_tracks.trk -f
  if ${params.debug}
  then
    scil_count_streamlines.py ${sid}__brainstem_plausible_tracks.trk > ${sid}__brainstem_plausible_tracks.txt
  fi
  """
}

/*
  END Brainstem
*/

process Remove_out_of_CGM_DWM {
  cpus 1
  tag = "Brain - Either end in CGM SWM"

  input:
    set sid, file(tractogram) from wb_for_split_end_in_CGMSWI

  output:
    set sid, "${sid}__wb_either_CGM_SWM.trk" into wb_for_extract_all_commissural
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

process Extracting_all_commissural {
  cpus 1
  tag = "Extract all commissural"

  input:
    set sid, file(tractogram) from wb_for_extract_all_commissural

  output:
    set sid, "${sid}__tmp_CC.trk" into cc_for_ee_BG, cc_for_remove_unplausible
    set sid, "${sid}__wb_either_CGM_SWM__noCC.trk" into no_cc_for_split_asso_BG
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
    set sid, file(tractogram) from cc_for_ee_BG
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

process First_cc_cleaning {
  cpus 1
  tag = "First_cc_cleaning"

  input:
    set sid, file(tractogram) from cc_for_remove_unplausible

  output:
    set sid, "${sid}__CC_Cx.trk" into cc_for_merge_plausible_01
    set sid, "${sid}__CC_lost.trk" optional true
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

process Split_no_CC_Asso_and_BG {
  cpus 1
  tag = "Split not CC in asso BG and not BG"

  input:
    set sid, file(tractogram) from no_cc_for_split_asso_BG

  output:
    set sid, "${sid}__all_subcortical_from_CGM_SWM_noCC_f.trk" into asso_BG_for_split_Thal, asso_BG_for_split_Put
    file "${sid}__all_subcortical_from_CGM_SWM_noCC_f.txt" optional true
    set sid, "${sid}__asso_noBG.trk" into asso_noBG_for_split_hemi
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

bg_list=params.bg_lists?.tokenize(',')
Channel.from(bg_list).into{bg_thal_list;
                           bg_put_list}

process Split_BG_Thal {
  cpus 1
  tag = "Split BG Thal"

  input:
    set sid, file(tractogram) from asso_BG_for_split_Thal
    each list from bg_thal_list
    each side from sides_split_BG_Thal

  output:
    set sid, "${sid}__BG_ipsi_Thal_${list}_${side}.trk" into BG_ipsi_Thal_for_merge

  script:
    filtering_list=params.filtering_lists_folder+"BG_ipsi_Thal_${list}_${side}_f.txt"
    out_extension="BG_ipsi_Thal_${list}_${side}"
    remaining_extension="garbage_BG_ipsi_Thal_${list}_${side}"
    basename="${sid}"

    template "filter_with_list.sh"
}

BG_ipsi_Thal_for_merge.groupTuple().map{it}.set{BG_ipsi_Thal_list_for_merge}

process Merge_BG_Thal{
  cpus 1
  tag = "Merge BG Thal"
  input:
    set sid, file(tractogram) from BG_ipsi_Thal_list_for_merge
  output:
    set sid, "${sid}__BG_ipsi_Thal_all.trk"

  script:
  """
  scil_streamlines_math.py concatenate ${tractogram} ${sid}__BG_ipsi_Thal_all.trk -f
  """
}

process Split_BG_Put {
  cpus 1
  tag = "Split BG Put"

  input:
    set sid, file(tractogram) from asso_BG_for_split_Put
    each list from bg_put_list
    each side from sides_split_BG_Put

  output:
    set sid, "${sid}__BG_ipsi_Put_${list}_${side}.trk" into BG_ipsi_Put_for_merge

  script:
    filtering_list=params.filtering_lists_folder+"BG_ipsi_Put_${list}_${side}_f.txt"
    out_extension="BG_ipsi_Put_${list}_${side}"
    remaining_extension="garbage_BG_ipsi_Put_${list}_${side}"
    basename="${sid}"

    template "filter_with_list.sh"
}

BG_ipsi_Put_for_merge.groupTuple().map{it}.set{BG_ipsi_Put_list_for_merge}

process Merge_BG_Put{
  cpus 1
  tag = "Merge BG Put"
  input:
    set sid, file(tractogram) from BG_ipsi_Put_list_for_merge
  output:
    set sid, "${sid}__BG_ipsi_Put_all.trk"

  script:
  """
  scil_streamlines_math.py concatenate ${tractogram} ${sid}__BG_ipsi_Put_all.trk -f
  """
}

process Split_Asso_In_Hemi {
  cpus 1
  tag = "Split asso in hemispheres"

  input:
    set sid, file(tractogram) from asso_noBG_for_split_hemi
    each side from sides

  output:
    set sid, val(side), "${sid}__asso_${side}.trk" into asso_for_extract_u_shape
    set sid, "${sid}__asso_${side}_lost.trk" optional true
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

process Split_UShape_CGM_Asso {
  cpus 1
  tag = "Extracting U-shaped and streamlines restricted to Cortical GM and removing them from asso"

  input:
    set sid, val(side), file(tractogram) from asso_for_extract_u_shape

  output:
    set sid, "${sid}__asso_only_in_CGM_${side}.trk" into assoCGM
    set sid, "${sid}__asso_Ushape_${side}.trk" into assoUShape
    set sid, val(side), "${sid}__asso_f_${side}.trk" into asso_for_remove_long_range
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
    set sid, val(side), file(tractogram) from asso_for_remove_long_range

  output:
    set sid, val(side), "${sid}__asso_all_intra_inter_${side}.trk" into asso_all_intra_inter
    set sid, "${sid}__asso_lost2_${side}.trk"
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

asso_all_intra_inter.into{asso_all_intra_inter_for_filtering;
                       asso_all_intra_inter_plausible}

asso_all_intra_inter_plausible.groupTuple().map{it.flatten().toList()}.set{asso_all_intra_inter_list}


cc_for_merge_plausible_01.into{ccCleanedPlausible; CC_for_homotopic}

assoCGM.groupTuple().map{it.flatten().toList()}.set{assoCGM_list}


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
      filtering_list=params.filtering_lists_folder+"CC_homo_${pair}_filtering_list_f.txt"
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
    set sid, val(side), file(tractogram) from asso_all_intra_inter_for_filtering
    each asso_list from asso_lists

  output:
    set sid, val(side), "${sid}__asso_${asso_list}_${side}.trk" into asso_all_intra_inter_filtered
    set sid, "${sid}__asso_lost_${asso_list}_${side}.trk"
    file "${sid}__asso_${asso_list}_${side}.txt" optional true
    file "${sid}__asso_lost_${asso_list}_${side}.txt" optional true

  script:
    filtering_list=params.filtering_lists_folder+"ASSO_${asso_list}_${side}_filtering_list.txt"
    out_extension="asso_${asso_list}_${side}"
    remaining_extension="asso_lost_${asso_list}_${side}"
    basename="${sid}"

    template "filter_with_list.sh"
}
