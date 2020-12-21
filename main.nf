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
                         sides_split_BG_Thal;
                         sides_split_BG_Put;
                         sides_split_BG_Caud;
                         sides_split_asso_ventral_f_o_f_p;
                         sides_split_asso_ventral_f_t;
                         sides_split_asso_ventral;
                         sides_split_asso_dorsal_f_o_f_t;
                         sides_split_asso_dorsal_f_p;
                         sides_split_asso_dorsal;
                         sides_split_asso_p_o;
                         sides_split_asso_o_t;
                         sides_split_asso_p_t;
                         sides_split_asso_ins;
                         sides_split_asso_cing;
                         sides_split_asso_frontal_be;
                         sides_split_asso_frontal_ee;
                         sides_split_asso_occipital_be;
                         sides_split_asso_occipital_ee;
                         sides_split_asso_parietal_be;
                         sides_split_asso_parietal_ee;
                         sides_split_asso_temporal_be;
                         sides_split_asso_temporal_ee;}

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
    set sid, "${sid}__wb_no_In_CC_DWM.trk" into unplausible_for_fornix
    file "${sid}__wb_clean01.txt" optional true
    file "${sid}__wb_no_In_CC_DWM.txt" optional true

  script:
  """
  scil_filter_tractogram.py ${wb_min20_noloop} ${sid}__wb_clean01.trk \
                    --drawn_roi ${params.rois_folder}${params.atlas.cc} either_end exclude \
                    --drawn_roi ${params.rois_folder}${params.atlas.dwm} either_end exclude \
                    -f --display_count > ${sid}__wb_clean01.txt
  scil_streamlines_math.py difference ${wb_min20_noloop} ${sid}__wb_clean01.trk ${sid}__wb_no_In_CC_DWM.trk -f
  scil_count_streamlines.py {sid}__wb_no_In_CC_DWM.trk > ${sid}__wb_no_In_CC_DWM.txt
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
    set sid, "${sid}__wb_either_CGM_SWM_noCC.trk" into no_cc_for_split_asso_BG
    file "${sid}__wb_either_CGM_SWM_noCC.txt" optional true
    file "${sid}__tmp_CC.txt" optional true

  script:
  atlas=params.rois_folder+params.atlas.midline
  mode=params.mode.any
  criteria=params.criteria.exclude
  out_extension="wb_either_CGM_SWM_noCC"
  remaining_extension='tmp_CC'
  basename="${sid}"

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
    scil_streamlines_math.py difference ${tractogram} ${sid}__CC_Cx.trk ${sid}__CC_lost.trk # -CC_BG
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
    set sid, "${sid}__all_subcortical_from_CGM_SWM_noCC_f.trk" into asso_BG_for_split_Thal, asso_BG_for_split_Put, asso_BG_for_split_Caud
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

bg_caud_list=params.bg_caud_lists?.tokenize(',')

process Split_BG_Caud {
  cpus 1
  tag = "Split BG Caud"

  input:
    set sid, file(tractogram) from asso_BG_for_split_Caud
    each list from bg_caud_list
    each side from sides_split_BG_Caud

  output:
    set sid, "${sid}__BG_ipsi_Caud_${list}_${side}.trk" into BG_ipsi_Caud_for_merge

  script:
    filtering_list=params.filtering_lists_folder+"BG_ipsi_Caud_${list}_${side}_f.txt"
    out_extension="BG_ipsi_Caud_${list}_${side}"
    remaining_extension="garbage_BG_ipsi_Caud_${list}_${side}"
    basename="${sid}"

    template "filter_with_list.sh"
}

BG_ipsi_Caud_for_merge.groupTuple().map{it}.set{BG_ipsi_Caud_list_for_merge}

process Merge_BG_Caud{
  cpus 1
  tag = "Merge BG Caud"
  input:
    set sid, file(tractogram) from BG_ipsi_Caud_list_for_merge

  output:
    set sid, "${sid}__BG_ipsi_Caud_all.trk"

  script:
  """
  scil_streamlines_math.py concatenate ${tractogram} ${sid}__BG_ipsi_Caud_all.trk -f
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
    set sid, val(side), "${sid}__asso_only_in_CGM_${side}.trk" into assoCGM
    set sid, val(side), "${sid}__asso_Ushape_${side}.trk" into assoUShape
    set sid, val(side), "${sid}__asso_Ushape_${side}_u.trk" into assoUShapef
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

    track_vis ${sid}__asso_Ushape_${side}.trk -u 0.5 1 -nr -o ${sid}__asso_Ushape_${side}_u.trk

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

asso_all_intra_inter.into{asso_all_intra_inter_for_ventral_f_o_f_p_filtering;
                          asso_all_intra_inter_for_ventral_f_t_filtering;
                          asso_all_intra_inter_for_dorsal_f_o_f_t_filtering;
                          asso_all_intra_inter_for_dorsal_f_p_filtering;
                          asso_all_intra_inter_for_p_o_filtering;
                          asso_all_intra_inter_for_p_t_filtering;
                          asso_all_intra_inter_for_o_t_filtering;
                          asso_all_intra_inter_for_ins_filtering;
                          asso_all_intra_inter_for_cing_filtering;
                          asso_all_intra_inter_for_be_frontal_filtering;
                          asso_all_intra_inter_for_ee_frontal_filtering;
                          asso_all_intra_inter_for_be_occipital_filtering;
                          asso_all_intra_inter_for_ee_occipital_filtering;
                          asso_all_intra_inter_for_be_parietal_filtering;
                          asso_all_intra_inter_for_ee_parietal_filtering;
                          asso_all_intra_inter_for_be_temporal_filtering;
                          asso_all_intra_inter_for_ee_temporal_filtering
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
    set sid, "${sid}__cc_homotopic_${pair}.trk" into CC_Homotopic_for_merge

    script:
      filtering_list=params.filtering_lists_folder+"CC_homo_${pair}_filtering_list_f.txt"
      out_extension="cc_homotopic_${pair}"
      remaining_extension="garbage_${pair}"
      basename="${sid}"

      template "filter_with_list.sh"
}

CC_Homotopic_for_merge.groupTuple().map{it}.set{CC_Homotopic_list_for_merge}

process CC_Homotopic_merge {
  cpus 1
  tag = "Merge CC Homotopic"

input:
  set sid, file(tractogram) from CC_Homotopic_list_for_merge

output:
  set sid, "${sid}__CC_homo.trk"

script:
  """
  scil_streamlines_math.py concatenate ${tractogram} ${sid}__CC_homo.trk
  """
}

/*
  ASSO VENTRAL
*/

asso_ventral_f_t_list=params.asso_ventral_f_t_lists?.tokenize(',')

process asso_ventral_f_t {
  cpus 1
  tag = "Asso Ventral filtering"

  input:
    set sid, val(side), file(tractogram) from asso_all_intra_inter_for_ventral_f_t_filtering
    each asso_list from asso_ventral_f_t_list

  output:
    set sid, val(side), "${sid}__asso_F_${asso_list}_${side}.trk" into asso_all_intra_inter_ventral_f_t_for_merge
    set sid, "${sid}__asso_lost_${asso_list}_${side}.trk"
    file "${sid}__asso_F_${asso_list}_${side}.txt" optional true
    file "${sid}__asso_lost_${asso_list}_${side}.txt" optional true

  script:
    filtering_list=params.filtering_lists_folder+"ASSO_F_${asso_list}_${side}_filtering_list.txt"
    out_extension="asso_F_${asso_list}_${side}"
    remaining_extension="asso_lost_${asso_list}_${side}"
    basename="${sid}"

    template "filter_with_list.sh"
}

asso_all_intra_inter_ventral_f_t_for_merge.groupTuple(by:[0,1]).map{it}.set{asso_all_intra_inter_ventral_f_t_list_for_merge}

process merge_asso_ventral_f_t {
  cpus 1
  tag = "Merge Asso Ventral F_T"

  input:
    set sid, val(side), file(tractogram) from asso_all_intra_inter_ventral_f_t_list_for_merge

  output:
    set sid, val(side), "${sid}__asso_F_T_ventral_f_${side}.trk" into asso_all_intra_inter_ventral_all_f_t_for_merge

  script:
  """
  scil_streamlines_math.py concatenate ${tractogram} ${sid}__asso_F_T_ventral_f_${side}.trk -f
  """
}

asso_ventral_f_o_f_p_lists=params.asso_ventral_f_o_f_p_lists?.tokenize(',')

process asso_ventral_f_o_f_p {
  cpus 1
  tag = "Asso Ventral filtering"

  input:
    set sid, val(side), file(tractogram) from asso_all_intra_inter_for_ventral_f_o_f_p_filtering
    each asso_list from asso_ventral_f_o_f_p_lists

  output:
    set sid, val(side), "${sid}__asso_F_${asso_list}_${side}.trk" into asso_all_intra_inter_ventral_all_f_o_f_p_for_merge
    set sid, "${sid}__asso_lost_${asso_list}_${side}.trk"
    file "${sid}__asso_${asso_list}_${side}.txt" optional true
    file "${sid}__asso_lost_${asso_list}_${side}.txt" optional true

  script:
    filtering_list=params.filtering_lists_folder+"ASSO_F_${asso_list}_${side}_filtering_list.txt"
    out_extension="asso_F_${asso_list}_${side}"
    remaining_extension="asso_lost_${asso_list}_${side}"
    basename="${sid}"

    template "filter_with_list.sh"
}

asso_all_intra_inter_ventral_all_f_t_for_merge.groupTuple(by:[0,1]).join(asso_all_intra_inter_ventral_all_f_o_f_p_for_merge.groupTuple(by:[0,1]), by:[0,1]).map{it.flatten().toList()}.set{asso_all_intra_inter_ventral_all_for_merge}

process merge_asso_ventral {
  cpus 1
  tag = "Merge Asso Ventral"

  input:
    set sid, val(side), file(trk01), file(trk02), file(trk03) from asso_all_intra_inter_ventral_all_for_merge

  output:
    set sid, val(side), "${sid}__asso_all_ventral_f_${side}.trk"

  script:
  """
  scil_streamlines_math.py concatenate ${trk01} ${trk02} ${trk03} ${sid}__asso_all_ventral_f_${side}.trk -f
  """
}

/*
  ASSO DORSAL
*/

asso_dorsal_f_p_lists=params.asso_dorsal_f_p_lists?.tokenize(',')

process asso_dorsal_f_p {
  cpus 1
  tag = "Asso Dorsal F_P filtering"

  input:
    set sid, val(side), file(tractogram) from asso_all_intra_inter_for_dorsal_f_p_filtering
    each asso_list from asso_dorsal_f_p_lists

  output:
    set sid, val(side), "${sid}__asso_${asso_list}_${side}.trk" into asso_all_intra_inter_dorsal_f_p_for_merge
    set sid, "${sid}__asso_lost_${asso_list}_${side}.trk"
    file "${sid}__asso_${asso_list}_${side}.txt" optional true
    file "${sid}__asso_lost_${asso_list}_${side}.txt" optional true

  script:
    filtering_list=params.filtering_lists_folder+"ASSO_F_${asso_list}_${side}_filtering_list.txt"
    out_extension="asso_${asso_list}_${side}"
    remaining_extension="asso_lost_${asso_list}_${side}"
    basename="${sid}"

    template "filter_with_list.sh"
}

asso_all_intra_inter_dorsal_f_p_for_merge.groupTuple(by:[0,1]).map{it}.set{asso_all_intra_inter_dorsal_f_p_list_for_merge}

process merge_asso_dorsal_f_p {
  cpus 1
  tag = "Merge Asso Dorsal F_P"

  input:
    set sid, val(side), file(tractogram) from asso_all_intra_inter_dorsal_f_p_list_for_merge

  output:
    set sid, val(side), "${sid}__asso_F_P_dorsal_f_${side}.trk" into asso_all_intra_inter_dorsal_all_f_p_for_merge

  script:
  """
  scil_streamlines_math.py union ${tractogram} ${sid}__asso_F_P_dorsal_f_${side}.trk -f
  """
}

asso_dorsal_f_o_f_t_list=params.asso_dorsal_f_o_f_t_lists?.tokenize(',')

process asso_dorsal_f_o_f_t {
  cpus 1
  tag = "Asso Dorsal F_O F_T filtering"

  input:
    set sid, val(side), file(tractogram) from asso_all_intra_inter_for_dorsal_f_o_f_t_filtering
    each asso_list from asso_dorsal_f_o_f_t_list

  output:
    set sid, val(side), "${sid}__asso_${asso_list}_${side}.trk" into asso_all_intra_inter_dorsal_all_f_o_f_t_for_merge
    set sid, "${sid}__asso_lost_${asso_list}_${side}.trk"
    file "${sid}__asso_${asso_list}_${side}.txt" optional true
    file "${sid}__asso_lost_${asso_list}_${side}.txt" optional true

  script:
    filtering_list=params.filtering_lists_folder+"ASSO_F_${asso_list}_${side}_filtering_list.txt"
    out_extension="asso_${asso_list}_${side}"
    remaining_extension="asso_lost_${asso_list}_${side}"
    basename="${sid}"

    template "filter_with_list.sh"
}

asso_all_intra_inter_dorsal_all_f_p_for_merge.groupTuple(by:[0,1]).join(asso_all_intra_inter_dorsal_all_f_o_f_t_for_merge.groupTuple(by:[0,1]), by:[0,1]).map{it.flatten().toList()}.set{asso_all_intra_inter_dorsal_all_for_merge}

process merge_asso_dorsal {
  cpus 1
  tag = "Merge Asso Dorsal"

  input:
    set sid, val(side), file(trk01), file(trk02), file(trk03) from asso_all_intra_inter_dorsal_all_for_merge

  output:
    set sid, val(side), "${sid}__asso_all_dorsal_f_${side}.trk"

  script:
  """
  scil_streamlines_math.py concatenate ${trk01} ${trk02} ${trk03} ${sid}__asso_all_dorsal_f_${side}.trk -f
  """
}

/*
  ASSO P_O
*/

asso_p_o_list=params.asso_p_o_lists?.tokenize(',')

process asso_p_o {
  cpus 1
  tag = "Asso P_O filtering"

  input:
    set sid, val(side), file(tractogram) from asso_all_intra_inter_for_p_o_filtering
    each asso_list from asso_p_o_list

  output:
    set sid, val(side), "${sid}__asso_${asso_list}_${side}.trk" into asso_intra_inter_p_o_for_merge
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

asso_intra_inter_p_o_for_merge.groupTuple(by:[0,1]).map{it}.set{asso_intra_inter_p_o_list_for_merge}

process merge_p_o {
  cpus 1
  tag = "Merge Asso P_O"

  input:
    set sid, val(side), file(tractogram) from asso_intra_inter_p_o_list_for_merge

  output:
    set sid, val(side), "${sid}__asso_all_P_O_f_${side}.trk"

  script:
  """
  scil_streamlines_math.py concatenate ${tractogram} ${sid}__asso_all_P_O_f_${side}.trk -f
  """
}

/*
  ASSO P_T
*/

asso_p_t_list=params.asso_p_t_lists?.tokenize(',')

process asso_p_t {
  cpus 1
  tag = "Asso P_T filtering"

  input:
    set sid, val(side), file(tractogram) from asso_all_intra_inter_for_p_t_filtering
    each asso_list from asso_p_t_list

  output:
    set sid, val(side), "${sid}__asso_${asso_list}_${side}.trk" into asso_intra_inter_p_t_for_merge
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

asso_intra_inter_p_t_for_merge.groupTuple(by:[0,1]).map{it}.set{asso_intra_inter_p_t_list_for_merge}

process merge_p_t {
  cpus 1
  tag = "Merge Asso P_T"

  input:
    set sid, val(side), file(tractogram) from asso_intra_inter_p_t_list_for_merge

  output:
    set sid, val(side), "${sid}__asso_all_P_T_f_${side}.trk"

  script:
  """
  scil_streamlines_math.py concatenate ${tractogram} ${sid}__asso_all_P_T_f_${side}.trk -f
  """
}

/*
  ASSO O_T
*/

asso_o_t_list=params.asso_o_t_lists?.tokenize(',')

process asso_o_t {
  cpus 1
  tag = "Asso O_T filtering"

  input:
    set sid, val(side), file(tractogram) from asso_all_intra_inter_for_o_t_filtering
    each asso_list from asso_o_t_list

  output:
    set sid, val(side), "${sid}__asso_${asso_list}_${side}.trk" into asso_intra_inter_o_t_for_merge
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

asso_intra_inter_o_t_for_merge.groupTuple(by:[0,1]).map{it}.set{asso_intra_inter_o_t_list_for_merge}

process merge_o_t {
  cpus 1
  tag = "Merge Asso O_T"

  input:
    set sid, val(side), file(tractogram) from asso_intra_inter_o_t_list_for_merge

  output:
    set sid, val(side), "${sid}__asso_all_O_T_f_${side}.trk"

  script:
  """
  scil_streamlines_math.py concatenate ${tractogram} ${sid}__asso_all_O_T_f_${side}.trk -f
  """
}

/*
  ASSO Ins
*/

asso_ins_list=params.asso_ins_lists?.tokenize(',')

process asso_ins {
  cpus 1
  tag = "Asso Ins filtering"

  input:
    set sid, val(side), file(tractogram) from asso_all_intra_inter_for_ins_filtering
    each asso_list from asso_ins_list

  output:
    set sid, val(side), "${sid}__asso_${asso_list}_${side}.trk" into asso_intra_inter_ins_for_merge
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

asso_intra_inter_ins_for_merge.groupTuple(by:[0,1]).map{it}.set{asso_intra_inter_ins_list_for_merge}

process merge_ins {
  cpus 1
  tag = "Merge Asso Ins"

  input:
    set sid, val(side), file(tractogram) from asso_intra_inter_ins_list_for_merge

  output:
    set sid, val(side), "${sid}__asso_all_Ins_f_${side}.trk"

  script:
  """
  scil_streamlines_math.py concatenate ${tractogram} ${sid}__asso_all_Ins_f_${side}.trk -f
  """
}

/*
  ASSO CING
*/

process asso_Cing {
  cpus 1
  tag = "Asso Cing filtering"

  input:
    set sid, val(side), file(tractogram) from asso_all_intra_inter_for_cing_filtering

  output:
    set sid, val(side), "${sid}__asso_all_Cing_${side}.trk"
    set sid, "${sid}__asso_lost_Cing_${side}.trk"
    file "${sid}__asso_all_Cing_${side}.txt" optional true
    file "${sid}__asso_lost_Cing_${side}.txt" optional true

  script:
    filtering_list=params.filtering_lists_folder+"ASSO_Cing_${side}_filtering_list.txt"
    out_extension="asso_all_Cing_${side}"
    remaining_extension="asso_lost_Cing_${side}"
    basename="${sid}"

    template "filter_with_list.sh"
}

/*
 BE ASSO FRONTAL: extracting all streamlines with both ends in a frontal gyrus (U-shape > 20 mm)
*/

asso_frontal_be_list=params.asso_frontal_be_lists?.tokenize(',')
process asso_be_frontal_gyrus {
  cpus 1
  tag = "Extract be frontal gyrus"

  input:
    set sid, val(side), file(tractogram) from asso_all_intra_inter_for_be_frontal_filtering
    each gyrus from asso_frontal_be_list

  output:
    set sid, val(side), val(gyrus), "${sid}_asso_intra_be_frontal_${gyrus}_${side}_u.trk" into asso_frontal_be_for_merge

  script:
  """
  scil_filter_tractogram.py ${tractogram} tmp.trk --filtering_list /filtering_lists/ASSO_be_${gyrus}_${side}_filtering_list.txt -f
  track_vis tmp.trk -u 0.5 1 -nr -o ${sid}_asso_intra_be_frontal_${gyrus}_${side}_u.trk
  """
}

asso_frontal_be_for_merge.groupTuple(by:[0,1]).map{it}.set{asso_frontal_be_list_for_merge}
process merge_asso_be_frontal_gyrus{
  cpus 1
  tag = "Merge asso be frontal gyrus"

  input:
    set sid, val(side),  val(gyrus), file(tractogram) from asso_frontal_be_list_for_merge

  output:
    set sid, val(side), "${sid}_asso_all_intraF_f_${side}_u.trk"

  script:
  """
    scil_streamlines_math.py lazy_concatenate ${tractogram} ${sid}_asso_all_intraF_f_${side}_u.trk -f
  """
}

/*
 EE ASSO FRONTAL: extracting all streamlines with either ends in a frontal gyrus (U-shape > 20 mm)
*/

asso_frontal_ee_list = Channel.from(['SFG_MFG', 70],['SFG_IFG', 70], ['SFG_PrCG', 90], ['SFG_FrOrbG', 70], ['MFG_IFG', 70], ['MFG_PrCG', 110], ['MFG_FrOrbG', 60], ['IFG_PrCG', 110], ['IFG_FrOrbG', 60])
asso_all_intra_inter_for_ee_frontal_filtering.combine(asso_frontal_ee_list).set{asso_frontal_ee_for_extract}
process asso_ee_frontal_gyrus {
  cpus 1
  tag = "Extract be frontal gyrus"

  input:
    set sid, val(side), file(tractogram), val(gyrus), val(max_length) from asso_frontal_ee_for_extract

  output:
    set sid, val(side), val(gyrus), "${sid}_asso_intra_ee_frontal_${gyrus}_${side}.trk" into asso_frontal_ee_for_merge

  script:
  """
  scil_filter_tractogram.py ${tractogram} tmp_01.trk --filtering_list /filtering_lists/ASSO_ee_${gyrus}_${side}_filtering_list.txt -f
  scil_filter_streamlines_by_length.py tmp_01.trk tmp_02.trk --maxL ${max_length} -f
  track_vis tmp_02.trk -u 0.5 1 -nr -o ${sid}_asso_intra_ee_frontal_${gyrus}_${side}.trk
  """
}

asso_frontal_ee_for_merge.groupTuple(by:[0,1]).map{it}.set{asso_frontal_ee_list_for_merge}
process merge_asso_ee_frontal_gyrus{
  cpus 1
  tag = "Merge asso ee frontal gyrus"

  input:
    set sid, val(side),  val(gyrus), file(tractogram) from asso_frontal_ee_list_for_merge

  output:
    set sid, val(side), "${sid}_asso_all_intraF_f_${side}_u.trk"

  script:
  """
    scil_streamlines_math.py lazy_concatenate ${tractogram} ${sid}_asso_all_intraF_f_${side}_u.trk -f
  """
}

/*
 BE ASSO OCCIPITAL: extracting all streamlines with both ends in a occipital gyrus (U-shape > 20 mm)
*/

asso_occipital_be_list=params.asso_occipital_be_lists?.tokenize(',')
process asso_be_occipital_gyrus {
  cpus 1
  tag = "Extract be occipital gyrus"

  input:
    set sid, val(side), file(tractogram) from asso_all_intra_inter_for_be_occipital_filtering
    each gyrus from asso_occipital_be_list

  output:
    set sid, val(side), val(gyrus), "${sid}_asso_intra_be_occipital_${gyrus}_${side}_u.trk" into asso_occipital_be_for_merge

  script:
  """
  scil_filter_tractogram.py ${tractogram} tmp.trk --filtering_list /filtering_lists/ASSO_be_${gyrus}_${side}_filtering_list.txt -f
  track_vis tmp.trk -u 0.5 1 -nr -o ${sid}_asso_intra_be_occipital_${gyrus}_${side}_u.trk
  """
}

asso_occipital_be_for_merge.groupTuple(by:[0,1]).map{it}.set{asso_occipital_be_list_for_merge}
process merge_asso_be_occipital_gyrus{
  cpus 1
  tag = "Merge asso be occipital gyrus"

  input:
    set sid, val(side),  val(gyrus), file(tractogram) from asso_occipital_be_list_for_merge

  output:
    set sid, val(side), "${sid}_asso_all_intraO_f_${side}_u.trk"

  script:
  """
    scil_streamlines_math.py lazy_concatenate ${tractogram} ${sid}_asso_all_intraO_f_${side}_u.trk -f
  """
}

/*
 EE ASSO OCCIPITAL: extracting all streamlines with either ends in a occipital gyrus (U-shape > 20 mm)
*/

asso_occipital_ee_list = Channel.from(['MOG_SOG', 60],['MOG_IOG', 50], ['MOG_CuG', 60], ['SOG_CuG', 30], ['CuG_LG', 60])
asso_all_intra_inter_for_ee_occipital_filtering.combine(asso_occipital_ee_list).set{asso_occipital_ee_for_extract}
process asso_ee_occipital_gyrus {
  cpus 1
  tag = "Extract be occipital gyrus"

  input:
    set sid, val(side), file(tractogram), val(gyrus), val(max_length) from asso_occipital_ee_for_extract

  output:
    set sid, val(side), val(gyrus), "${sid}_asso_intra_ee_occipital_${gyrus}_${side}.trk" into asso_occipital_ee_for_merge

  script:
  """
  scil_filter_tractogram.py ${tractogram} tmp_01.trk --filtering_list /filtering_lists/ASSO_ee_${gyrus}_${side}_filtering_list.txt -f
  scil_filter_streamlines_by_length.py tmp_01.trk tmp_02.trk --maxL ${max_length} -f
  track_vis tmp_02.trk -u 0.5 1 -nr -o ${sid}_asso_intra_ee_occipital_${gyrus}_${side}.trk
  """
}

asso_occipital_ee_for_merge.groupTuple(by:[0,1]).map{it}.set{asso_occipital_ee_list_for_merge}
process merge_asso_ee_occipital_gyrus{
  cpus 1
  tag = "Merge asso ee occipital gyrus"

  input:
    set sid, val(side),  val(gyrus), file(tractogram) from asso_occipital_ee_list_for_merge

  output:
    set sid, val(side), "${sid}_asso_all_intraO_f_${side}_u.trk"

  script:
  """
    scil_streamlines_math.py lazy_concatenate ${tractogram} ${sid}_asso_all_intraO_f_${side}_u.trk -f
  """
}

/*
 BE ASSO PARIETAL: extracting all streamlines with both ends in a parietal gyrus (U-shape > 20 mm)
*/

asso_parietal_be_list=params.asso_parietal_be_lists?.tokenize(',')
process asso_be_parietal_gyrus {
  cpus 1
  tag = "Extract be occipital gyrus"

  input:
    set sid, val(side), file(tractogram) from asso_all_intra_inter_for_be_parietal_filtering
    each gyrus from asso_parietal_be_list

  output:
    set sid, val(side), val(gyrus), "${sid}_asso_intra_be_parietal_${gyrus}_${side}_u.trk" into asso_parietal_be_for_merge

  script:
  """
  scil_filter_tractogram.py ${tractogram} tmp.trk --filtering_list /filtering_lists/ASSO_be_${gyrus}_${side}_filtering_list.txt -f
  track_vis tmp.trk -u 0.5 1 -nr -o ${sid}_asso_intra_be_parietal_${gyrus}_${side}_u.trk
  """
}

asso_parietal_be_for_merge.groupTuple(by:[0,1]).map{it}.set{asso_parietal_be_list_for_merge}
process merge_asso_be_parietal_gyrus{
  cpus 1
  tag = "Merge asso be parietal gyrus"

  input:
    set sid, val(side),  val(gyrus), file(tractogram) from asso_parietal_be_list_for_merge

  output:
    set sid, val(side), "${sid}_asso_all_intraP_f_${side}_u.trk"

  script:
  """
    scil_streamlines_math.py lazy_concatenate ${tractogram} ${sid}_asso_all_intraP_f_${side}_u.trk -f
  """
}

/*
 EE ASSO PARIETAL: extracting all streamlines with either ends in a parietal gyrus (U-shape > 20 mm)
*/

asso_parietal_ee_list = Channel.from(['SPG_PoCG', 50], ['SPG_AG', 80], ['SPG_SMG', 70], ['SPG_PrCuG', 50], ['AG_PoCG', 10000], ['AG_SMG', 90], ['AG_PrCuG', 90] , ['SMG_PoCG', 60], ['SMG_PrCuG',100], ['PoCG_PrCuG', 80])
asso_all_intra_inter_for_ee_parietal_filtering.combine(asso_parietal_ee_list).set{asso_parietal_ee_for_extract}
process asso_ee_parietal_gyrus {
  cpus 1
  tag = "Extract be parietal gyrus"

  input:
    set sid, val(side), file(tractogram), val(gyrus), val(max_length) from asso_parietal_ee_for_extract

  output:
    set sid, val(side), val(gyrus), "${sid}_asso_intra_ee_parietal_${gyrus}_${side}.trk" into asso_parietal_ee_for_merge

  script:
  """
  scil_filter_tractogram.py ${tractogram} tmp_01.trk --filtering_list /filtering_lists/ASSO_ee_${gyrus}_${side}_filtering_list.txt -f
  scil_filter_streamlines_by_length.py tmp_01.trk tmp_02.trk --maxL ${max_length} -f
  track_vis tmp_02.trk -u 0.5 1 -nr -o ${sid}_asso_intra_ee_parietal_${gyrus}_${side}.trk
  """
}

asso_parietal_ee_for_merge.groupTuple(by:[0,1]).map{it}.set{asso_parietal_ee_list_for_merge}
process merge_asso_ee_parietal_gyrus{
  cpus 1
  tag = "Merge asso ee parietal gyrus"

  input:
    set sid, val(side),  val(gyrus), file(tractogram) from asso_parietal_ee_list_for_merge

  output:
    set sid, val(side), "${sid}_asso_all_intraP_f_${side}.trk"

  script:
  """
    scil_streamlines_math.py lazy_concatenate ${tractogram} ${sid}_asso_all_intraP_f_${side}.trk -f
  """
}

/*
 BE ASSO TEMPORAL: extracting all streamlines with both ends in a temporal gyrus and merge (U-shape > 20 mm)
*/

asso_temporal_be_list=params.asso_temporal_be_lists?.tokenize(',')
process asso_be_temporal_gyrus {
  cpus 1
  tag = "Extract be temporal gyrus"

  input:
    set sid, val(side), file(tractogram) from asso_all_intra_inter_for_be_temporal_filtering
    each gyrus from asso_temporal_be_list

  output:
    set sid, val(side), val(gyrus), "${sid}_asso_intra_be_temporal_${gyrus}_${side}_u.trk" into asso_temporal_be_for_merge

  script:
  """
  scil_filter_tractogram.py ${tractogram} tmp.trk --filtering_list /filtering_lists/ASSO_be_${gyrus}_${side}_filtering_list.txt -f
  track_vis tmp.trk -u 0.5 1 -nr -o ${sid}_asso_intra_be_temporal_${gyrus}_${side}_u.trk
  """
}

asso_temporal_be_for_merge.groupTuple(by:[0,1]).map{it}.set{asso_temporal_be_list_for_merge}
process merge_asso_be_temporal_gyrus{
  cpus 1
  tag = "Merge asso be temporal gyrus"

  input:
    set sid, val(side),  val(gyrus), file(tractogram) from asso_temporal_be_list_for_merge

  output:
    set sid, val(side), "${sid}_asso_all_intraT_f_${side}_u.trk"

  script:
  """
    scil_streamlines_math.py lazy_concatenate ${tractogram} ${sid}_asso_all_intraT_f_${side}_u.trk -f
  """
}

/*
 EE ASSO TEMPORAL: extracting all streamlines with either ends in a temporal gyrus and merge (U-shape > 20 mm)
*/

asso_temporal_ee_list = Channel.from(['STG_MTG', 60], ['STG_ITG',80], ['STG_Tpole',110], ['MTG_ITG',60], ['MTG_Tpole', 100000], ['ITG_Tpole', 60])
asso_all_intra_inter_for_ee_temporal_filtering.combine(asso_temporal_ee_list).set{asso_temporal_ee_for_extract}
process asso_ee_temporal_gyrus {
  cpus 1
  tag = "Extract be temporal gyrus"

  input:
    set sid, val(side), file(tractogram), val(gyrus), val(max_length) from asso_temporal_ee_for_extract

  output:
    set sid, val(side), val(gyrus), "${sid}_asso_intra_ee_temporal_${gyrus}_${side}.trk" into asso_temporal_ee_for_merge

  script:
  """
  scil_filter_tractogram.py ${tractogram} tmp_01.trk --filtering_list /filtering_lists/ASSO_ee_${gyrus}_${side}_filtering_list.txt -f
  scil_filter_streamlines_by_length.py tmp_01.trk tmp_02.trk --maxL ${max_length} -f
  track_vis tmp_02.trk -u 0.5 1 -nr -o ${sid}_asso_intra_ee_temporal_${gyrus}_${side}.trk
  """
}

asso_temporal_ee_for_merge.groupTuple(by:[0,1]).map{it}.set{asso_temporal_ee_list_for_merge}
process merge_asso_ee_temporal_gyrus{
  cpus 1
  tag = "Merge asso ee temporal gyrus"

  input:
    set sid, val(side),  val(gyrus), file(tractogram) from asso_temporal_ee_list_for_merge

  output:
    set sid, val(side), "${sid}_asso_all_intraT_f_${side}.trk"

  script:
  """
    scil_streamlines_math.py lazy_concatenate ${tractogram} ${sid}_asso_all_intraT_f_${side}.trk -f
  """
}
