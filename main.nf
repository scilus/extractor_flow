#!/usr/bin/env nextflow

params.input = false
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

if (params.input){
    log.info "Input: $params.input"
    root = file(params.input)
    in_tractogram = Channel
        .fromFilePairs("$root/**/*.trk",
                       size:1,
                       maxDepth:1,
                       flat: true) {it.parent.name}
}
else {
    error "Error ~ Please use --input for the input data."
}

in_tractogram.into{for_remove_invalid_streamlines; in_tractogram_for_extract_first_unplausible; in_tractogram_for_extract_unplausible}

sides = params.sides?.tokenize(',')
Channel.from(sides).into{sides_ipsi;
                         sides_split_CC_BG;
                         sides_split_BG_Thal;
                         sides_split_BG_Put;
                         sides_split_BG_Caud;
                         side_corticopontineF;
                         side_corticopontinePOT;
                         side_cst}

process Remove_invalid_streamlines {
    cpus 1

    input:
      set sid, file(tractogram) from for_remove_invalid_streamlines

    output:
      set sid, "${sid}__rm_invalid_streamlines.trk" into for_remove_out_not_JHU

    script:
    """
      scil_remove_invalid_streamlines.py ${tractogram} ${sid}__rm_invalid_streamlines.trk --cut_invalid --remove_single_point -f
    """
}

process Remove_out_not_JHU {
    cpus 1

    input:
      set sid, file(tractogram) from for_remove_out_not_JHU

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

process Remove_crossing_Gyri {
  cpus 1

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

/*
Pruning min ${params.min_streaminline_lenght}mm
*/
process Pruning {
  cpus 1

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

process remove_loops {
  cpus 1

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

process remove_ee_CC_DWM {
  cpus 1

  input:
    set sid, file(wb_min20_noloop) from wb_for_rm_end_in_cc_dwm

  output:
    set sid, "${sid}__wb_clean01.trk" into wb_for_extract_end_in_cerebellum, wb_for_extract_first_unplausible
    set sid, "${sid}__wb_no_In_CC_DWM.trk"
    file "${sid}__wb_clean01.txt" optional true
    file "${sid}__wb_no_In_CC_DWM.txt" optional true

  script:
  """
  scil_filter_tractogram.py ${wb_min20_noloop} ${sid}__wb_clean01.trk \
                    --drawn_roi ${params.rois_folder}${params.atlas.cc} either_end exclude \
                    --drawn_roi ${params.rois_folder}${params.atlas.dwm} either_end exclude \
                    -f --display_count > ${sid}__wb_clean01.txt
  scil_streamlines_math.py difference ${wb_min20_noloop} ${sid}__wb_clean01.trk ${sid}__wb_no_In_CC_DWM.trk -f
  scil_count_streamlines.py ${sid}__wb_no_In_CC_DWM.trk > ${sid}__wb_no_In_CC_DWM.txt
  """
}

in_tractogram_for_extract_first_unplausible.join(wb_for_extract_first_unplausible).set{unplausible_streamlines}
process extract_first_unplausible{
  cpus 1

  input:
    set sid, file(tractogram1), file(tractogram2) from unplausible_streamlines

  output:
    set sid, "${sid}__unplausible_streamlines.trk" into unplausible_for_fornix

  script:
  """
  scil_streamlines_math.py difference ${tractogram1} \
                                      ${tractogram2} \
                                      ${sid}__unplausible_streamlines.trk;
  """
}

process extract_fornix{
  cpus 1

  input:
    set sid, file(tractogram) from unplausible_for_fornix

  output:
    set sid, "${sid}__fornix_f.trk" into fornix_for_trk_plausible, fornix_for_rename
    file "${sid}__fornix_f.txt"

  script:
    filtering_list=params.filtering_lists_folder+"fornix_filtering_list.txt"
    out_extension="fornix_f"
    remaining_extension="unplausible_wo_fornix"
    basename="${sid}"

    template "filter_with_list.sh"
}


process extract_ee_cerebellum {
  cpus 1

  input:
    set sid, file(tractogram) from wb_for_extract_end_in_cerebellum

  output:
    set sid, "${sid}__wb_clean01_nocereb.trk" into wb_for_extract_end_in_brainstem
    set sid, "${sid}__all_cerebellum.trk" into ee_cerebellum_for_extract_plausible
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


process extract_plausible_cerebellum {
  cpus 1

  input:
    set sid, file(tractogram) from ee_cerebellum_for_extract_plausible

  output:
    set sid, "${sid}__all_cerebellum_plausibles.trk" into cerebellum_for_trk_plausible, cerebellum_for_rename
    file "${sid}__all_in_cerebellum_nocx_nocerebwm.trk"
    file "${sid}__all_in_cerebellum_in_Medulla.trk"
    file "${sid}__all_in_cerebellum_in_Pons.trk"
    file "${sid}__all_in_cerebellum_in_Midbrain.trk"
    file "${sid}__all_in_cerebellum_in_redN_and_Thal.trk"

  script:
  """
  scil_filter_tractogram.py ${tractogram} ${sid}__all_in_cerebellum_nocx_nocerebwm.trk --filtering_list /filtering_lists/list_cerebellum/cerebellum_nocx_in_cereb_filtering_list.txt -f
  scil_filter_tractogram.py ${tractogram} ${sid}__all_in_cerebellum_in_Medulla.trk --filtering_list /filtering_lists/list_cerebellum/cerebellum_in_medulla_filtering_list.txt -f
  scil_filter_tractogram.py ${tractogram} ${sid}__all_in_cerebellum_in_Pons.trk --filtering_list /filtering_lists/list_cerebellum/cerebellum_in_pons_filtering_list.txt -f
  scil_filter_tractogram.py ${tractogram} ${sid}__all_in_cerebellum_in_Midbrain.trk --filtering_list /filtering_lists/list_cerebellum/cerebellum_in_midbrain_filtering_list.txt -f
  scil_filter_tractogram.py ${tractogram} ${sid}__all_in_cerebellum_in_redN_and_Thal.trk --filtering_list /filtering_lists/list_cerebellum/cerebellum_in_rednucleus_and_thalamus_filtering_list.txt -f
  scil_streamlines_math.py concatenate ${sid}__all_in_*.trk ${sid}__all_cerebellum_plausibles.trk -f
  scil_streamlines_math.py difference ${sid}__all_cerebellum.trk ${sid}__all_cerebellum_plausibles.trk ${sid}__all_cerebellum_unplausibles.trk -f
  """
}


/*
  END Cerebellum
*/

process extract_ee_brainstem {
  cpus 1

  input:
    set sid, file(tractogram) from wb_for_extract_end_in_brainstem

  output:
    set sid, "${sid}__wb_clean02.trk" into wb_for_split_end_in_CGMSWI
    set sid, "${sid}__all_brainstem.trk" into all_brainstem_for_extract_plausible
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

process extract_plausible_brainstem {
  cpus 1

  input:
    set sid, file(tractogram) from all_brainstem_for_extract_plausible
  output:
    set sid, "${sid}__all_brainstem_plausibles.trk" into brainstem_for_trk_plausible, brainstem_for_rename
    file "${sid}__all_brainstem_unplausibles.trk"
    file "${sid}__be_midbrain.trk"
    file "${sid}__be_medulla.trk"
    file "${sid}__be_pons.trk"
    file "${sid}__ee_thalamus.trk"
    file "${sid}__ee_red_nucleus.trk"
    set sid, "${sid}__ee_fronto_pontine.trk" into brainstem_corticopontine_frontal_for_rename
    set sid, "${sid}__ee_parietotemporooccipital_pontine.trk" into brainstem_ee_corticopontine_parietotemporooccipital_for_rename
    set sid, "${sid}__ee_pyramidal.trk" into brainstem_pyramidal_for_rename
    file "${sid}__ee_cortico_tectal.trk"

  script:
  """
  # Extract be midbrain
  scil_filter_tractogram.py ${sid}__all_brainstem.trk ${sid}__be_midbrain.trk --filtering_list /filtering_lists/list_brainstem/brainstem_be_midbrain_filtering_list.txt -f
  # Extract be medulla
  scil_filter_tractogram.py ${sid}__all_brainstem.trk ${sid}__be_medulla.trk --filtering_list /filtering_lists/list_brainstem/brainstem_be_medulla_filtering_list.txt -f
  # Extract be pons
  scil_filter_tractogram.py ${sid}__all_brainstem.trk ${sid}__be_pons.trk --filtering_list /filtering_lists/list_brainstem/brainstem_be_pons_filtering_list.txt -f

  # Extract ee thalamus
  scil_filter_tractogram.py ${sid}__all_brainstem.trk ${sid}__ee_thalamus.trk --filtering_list /filtering_lists/list_brainstem/brainstem_ee_thalamus_filtering_list.txt -f
  # Extract ee red_nucleus
  scil_filter_tractogram.py ${sid}__all_brainstem.trk ${sid}__ee_red_nucleus.trk --filtering_list /filtering_lists/list_brainstem/brainstem_ee_red_nucleus_filtering_list.txt -f

  # Prepartion for fronto-pontine, parietotemporooccipito-pontine, pyramidal, cortico-tectal
  scil_filter_tractogram.py ${sid}__all_brainstem.trk ${sid}__ee_tmp_01.trk --filtering_list /filtering_lists/list_brainstem/brainstem_ee_tmp_01_filtering_list.txt -f
  scil_filter_tractogram.py ${sid}__all_brainstem.trk ${sid}__ee_tmp_02.trk --filtering_list /filtering_lists/list_brainstem/brainstem_ee_tmp_02_filtering_list.txt -f
  scil_streamlines_math.py concatenate ${sid}__ee_tmp_01.trk ${sid}__ee_tmp_02.trk ${sid}__ee_tmp_03.trk -f

  # Extract ee Fronto-pontine
  scil_filter_tractogram.py ${sid}__ee_tmp_03.trk ${sid}__ee_fronto_pontine.trk --filtering_list /filtering_lists/list_brainstem/brainstem_ee_F_pontine_filtering_list.txt -f
  # Extract ee ParietoTemporooccipital pontine
  scil_filter_tractogram.py ${sid}__ee_tmp_03.trk ${sid}__ee_parietotemporooccipital_pontine.trk --filtering_list /filtering_lists/list_brainstem/brainstem_ee_PTO_pontine_filtering_list.txt -f
  # Extract ee Pyramidal
  scil_filter_tractogram.py ${sid}__ee_tmp_03.trk ${sid}__ee_pyramidal.trk --filtering_list /filtering_lists/list_brainstem/brainstem_ee_pyramidal_filtering_list.txt -f
  # Extract ee Tectal
  scil_filter_tractogram.py ${sid}__ee_tmp_03.trk ${sid}__ee_cortico_tectal.trk --filtering_list /filtering_lists/list_brainstem/brainstem_ee_cortico_tectal_filtering_list.txt -f
  scil_filter_streamlines_by_length.py ${sid}__ee_cortico_tectal.trk ${sid}__ee_cortico_tectal.trk --maxL 100 -f

  rm -f ${sid}__*tmp_*.trk

  scil_streamlines_math.py concatenate ${sid}__be_*.trk ${sid}__ee_*.trk ${sid}__all_brainstem_plausibles.trk -f
  scil_streamlines_math.py difference ${sid}__all_brainstem.trk ${sid}__all_brainstem_plausibles.trk ${sid}__all_brainstem_unplausibles.trk -f
  """
}

/*
Brain - Either end in CGM SWM
*/

process remove_out_of_CGM_DWM {
  cpus 1

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

process extract_all_commissural {
  cpus 1

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

process split_CC_BG {
  cpus 1

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
  remaining_extension="notUsed"
  basename="${sid}"

  template "filter_with_atlas.sh"
}

process first_cc_cleaning {
  cpus 1

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

/*
Split not CC in asso BG and not BG
*/

process Split_no_CC_Asso_and_BG {
  cpus 1

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
/*
BG THAL
*/
process split_BG_Thal {
  cpus 1

  input:
    set sid, file(tractogram) from asso_BG_for_split_Thal
    each list from bg_thal_list
    each side from sides_split_BG_Thal

  output:
    set sid, "${sid}__BG_ipsi_Thal_${list}_${side}.trk" into BG_ipsi_Thal_for_merge
    set sid, val(side), "${sid}__BG_ipsi_Thal_${list}_${side}.trk" into BG_ipsi_Thal_split_for_merge
    set sid, val(side), val(list), "${sid}__BG_ipsi_Thal_${list}_${side}.trk" into BG_ipsi_Thal_for_filter_CuGWM, BG_ipsi_Thal_for_filter_LGWM

  script:
    filtering_list=params.filtering_lists_folder+"BG_ipsi_Thal_${list}_${side}_f.txt"
    out_extension="BG_ipsi_Thal_${list}_${side}"
    remaining_extension="garbage_BG_ipsi_Thal_${list}_${side}"
    basename="${sid}"

    template "filter_with_list.sh"
}

BG_ipsi_Thal_split_for_merge.groupTuple(by:[0,1]).set{BG_ipsi_Thal_for_rename}

BG_ipsi_Thal_for_filter_CuGWM.filter{it[2]=='CuGWM'}.set{CuGWM_for_combine}
BG_ipsi_Thal_for_filter_LGWM.filter{it[2]=='LGWM'}.set{LGWM_for_combine}
CuGWM_for_combine.concat(LGWM_for_combine).groupTuple(by:[0,1]).set{optic_radiation_for_rename}

BG_ipsi_Thal_for_merge.groupTuple().map{it}.set{BG_ipsi_Thal_list_for_merge}

process merge_BG_Thal{
  cpus 1

  input:
    set sid, file(tractogram) from BG_ipsi_Thal_list_for_merge

  output:
    set sid, "${sid}__BG_ipsi_Thal_all.trk" into BG_ipsi_Thal_for_trk_plausible

  script:
  """
  scil_streamlines_math.py concatenate ${tractogram} ${sid}__BG_ipsi_Thal_all.trk -f
  """
}

/*
BG PUT
*/
process split_BG_Put {
  cpus 1

  input:
    set sid, file(tractogram) from asso_BG_for_split_Put
    each list from bg_put_list
    each side from sides_split_BG_Put

  output:
    set sid, "${sid}__BG_ipsi_Put_${list}_${side}.trk" into BG_ipsi_Put_for_merge
    set sid, val(side), "${sid}__BG_ipsi_Put_${list}_${side}.trk" into BG_ipsi_Put_for_rename

  script:
    filtering_list=params.filtering_lists_folder+"BG_ipsi_Put_${list}_${side}_f.txt"
    out_extension="BG_ipsi_Put_${list}_${side}"
    remaining_extension="garbage_BG_ipsi_Put_${list}_${side}"
    basename="${sid}"

    template "filter_with_list.sh"
}

BG_ipsi_Put_for_merge.groupTuple().map{it}.set{BG_ipsi_Put_list_for_merge}

process merge_BG_Put{
  cpus 1

  input:
    set sid, file(tractogram) from BG_ipsi_Put_list_for_merge

  output:
    set sid, "${sid}__BG_ipsi_Put_all.trk" into BG_ipsi_Put_for_trk_plausible

  script:
  """
  scil_streamlines_math.py concatenate ${tractogram} ${sid}__BG_ipsi_Put_all.trk -f
  """
}

/*
BG CAUD
*/
bg_caud_list=params.bg_caud_lists?.tokenize(',')
process split_BG_Caud {
  cpus 1

  input:
    set sid, file(tractogram) from asso_BG_for_split_Caud
    each list from bg_caud_list
    each side from sides_split_BG_Caud

  output:
    set sid, "${sid}__BG_ipsi_Caud_${list}_${side}.trk" into BG_ipsi_Caud_for_merge
    set sid, val(side), "${sid}__BG_ipsi_Caud_${list}_${side}.trk" into BG_ipsi_Caud_for_rename

  script:
    filtering_list=params.filtering_lists_folder+"BG_ipsi_Caud_${list}_${side}_f.txt"
    out_extension="BG_ipsi_Caud_${list}_${side}"
    remaining_extension="garbage_BG_ipsi_Caud_${list}_${side}"
    basename="${sid}"

    template "filter_with_list.sh"
}

BG_ipsi_Caud_for_merge.groupTuple().map{it}.set{BG_ipsi_Caud_list_for_merge}

process merge_BG_Caud{
  cpus 1

  input:
    set sid, file(tractogram) from BG_ipsi_Caud_list_for_merge

  output:
    set sid, "${sid}__BG_ipsi_Caud_all.trk" into BG_ipsi_Caud_for_trk_plausible

  script:
  """
  scil_streamlines_math.py concatenate ${tractogram} ${sid}__BG_ipsi_Caud_all.trk -f
  """
}

process split_asso_in_hemi {
  cpus 1

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

/*
Extracting U-shaped and streamlines restricted to Cortical GM and removing them from asso
*/

process split_ushape_CGM_asso {
  cpus 1

  input:
    set sid, val(side), file(tractogram) from asso_for_extract_u_shape

  output:
    set sid, val(side), "${sid}__asso_only_in_CGM_${side}.trk" into assoCGM
    set sid, val(side), "${sid}__asso_Ushape_${side}.trk" into assoUShape
    set sid, "${sid}__asso_Ushape_${side}_u.trk" into asso_u_shape_for_trk_plausible
    set sid, val(side), "${sid}__asso_Ushape_${side}_u.trk" into asso_u_shape_for_rename

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

/*
Extracting unplausible long-range association streamlines passing through subcortical structures (Cd, Put, GP, Thal, Amyg)
*/

process Remove_Unplausible_Long_Range_Asso {
  cpus 1

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

/*
inCCBG.groupTuple().map{it.flatten().toList()}.set{inCCBG_List}
assoUShape.groupTuple().map{it.flatten().toList()}.set{assoUShape_list}
*/

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

  input:
    set sid, file(tractogram) from CC_for_homotopic
    each pair from cc_homotopic_pairs

  output:
    set sid, "${sid}__cc_homotopic_${pair}.trk" into CC_Homotopic_for_merge
    set sid, val(pair), "${sid}__cc_homotopic_${pair}.trk" into CC_Homotopic_for_filter_AGWM, CC_Homotopic_for_filter_CingGWM, CC_Homotopic_for_filter_CuGWM, CC_Homotopic_for_filter_FuGWM, CC_Homotopic_for_filter_Hippo, CC_Homotopic_for_filter_IFGWM, CC_Homotopic_for_filter_Ins, CC_Homotopic_for_filter_IOGWM, CC_Homotopic_for_filter_ITGWM, CC_Homotopic_for_filter_LFOGWM, CC_Homotopic_for_filter_LGWM, CC_Homotopic_for_filter_MFGWM, CC_Homotopic_for_filter_MFOGWM, CC_Homotopic_for_filter_MOGWM, CC_Homotopic_for_filter_MTGWM, CC_Homotopic_for_filter_PHG, CC_Homotopic_for_filter_PoCGWM, CC_Homotopic_for_filter_PrCGWM, CC_Homotopic_for_filter_PrCuGWM, CC_Homotopic_for_filter_RGGWM, CC_Homotopic_for_filter_SFGWM, CC_Homotopic_for_filter_SMGWM, CC_Homotopic_for_filter_SOGWM, CC_Homotopic_for_filter_SPGWM, CC_Homotopic_for_filter_STGWM, CC_Homotopic_for_filter_T_pole_gwm

  script:
    filtering_list=params.filtering_lists_folder+"CC_homo_${pair}_filtering_list_f.txt"
    out_extension="cc_homotopic_${pair}"
    remaining_extension="garbage_${pair}"
    basename="${sid}"

    template "filter_with_list.sh"
}

/*
Filter + Concat frontal
*/
CC_Homotopic_for_filter_IFGWM.filter{it[1]=='IFGWM'}.set{CC_IFGWM_for_combine_frontal}
CC_Homotopic_for_filter_SFGWM.filter{it[1]=='SFGWM'}.set{CC_SFGWM_for_combine_frontal}
CC_Homotopic_for_filter_MFGWM.filter{it[1]=='MFGWM'}.set{CC_MFGWM_for_combine_frontal}
CC_Homotopic_for_filter_MFOGWM.filter{it[1]=='MFOGWM'}.set{CC_MFOGWM_for_combine_frontal}
CC_Homotopic_for_filter_LFOGWM.filter{it[1]=='LFOGWM'}.set{CC_LFOGWM_for_combine_frontal}
CC_Homotopic_for_filter_PrCGWM.filter{it[1]=='PrCGWM'}.set{CC_PrCGWM_for_combine_frontal}
CC_Homotopic_for_filter_RGGWM.filter{it[1]=='RGGWM'}.set{CC_RGGWM_for_combine_frontal}

CC_IFGWM_for_combine_frontal.concat(CC_SFGWM_for_combine_frontal).concat(CC_MFGWM_for_combine_frontal).concat(CC_MFOGWM_for_combine_frontal).concat(CC_LFOGWM_for_combine_frontal).concat(CC_PrCGWM_for_combine_frontal).concat(CC_RGGWM_for_combine_frontal).groupTuple(by:0).set{CC_Homotopic_frontal_for_rename}

/*
Filter + Concat occipital
*/
CC_Homotopic_for_filter_SOGWM.filter{it[1]=='SOGWM'}.set{CC_SOGWM_for_combine_occipital}
CC_Homotopic_for_filter_MOGWM.filter{it[1]=='MOGWM'}.set{CC_MOGWM_for_combine_occipital}
CC_Homotopic_for_filter_IOGWM.filter{it[1]=='IOGWM'}.set{CC_IOGWM_for_combine_occipital}
CC_Homotopic_for_filter_CuGWM.filter{it[1]=='CuGWM'}.set{CC_CuGWM_for_combine_occipital}
CC_Homotopic_for_filter_LGWM.filter{it[1]=='LGWM'}.set{CC_LGWM_for_combine_occipital}

CC_SOGWM_for_combine_occipital.concat(CC_MOGWM_for_combine_occipital).concat(CC_IOGWM_for_combine_occipital).concat(CC_CuGWM_for_combine_occipital).concat(CC_LGWM_for_combine_occipital).groupTuple(by:0).set{CC_Homotopic_occipital_for_rename}

/*
Filter + Concat temporal
*/
CC_Homotopic_for_filter_STGWM.filter{it[1]=='STGWM'}.set{CC_STGWM_for_combine_temporal}
CC_Homotopic_for_filter_T_pole_gwm.filter{it[1]=='T_pole_gwm'}.set{CC_T_pole_gwm_for_combine_temporal}
CC_Homotopic_for_filter_MTGWM.filter{it[1]=='MTGWM'}.set{CC_MTGWM_for_combine_temporal}
CC_Homotopic_for_filter_ITGWM.filter{it[1]=='ITGWM'}.set{CC_ITGWM_for_combine_temporal}
CC_Homotopic_for_filter_PHG.filter{it[1]=='PHG'}.set{CC_PHG_for_combine_temporal}
CC_Homotopic_for_filter_Hippo.filter{it[1]=='Hippo'}.set{CC_Hippo_for_combine_temporal}
CC_Homotopic_for_filter_FuGWM.filter{it[1]=='FuGWM'}.set{CC_FuGWM_for_combine_temporal}

CC_STGWM_for_combine_temporal.concat(CC_T_pole_gwm_for_combine_temporal).concat(CC_MTGWM_for_combine_temporal).concat(CC_ITGWM_for_combine_temporal).concat(CC_PHG_for_combine_temporal).concat(CC_Hippo_for_combine_temporal).concat(CC_FuGWM_for_combine_temporal).groupTuple(by:0).set{CC_Homotopic_temporal_for_rename}

/*
Filter + Concat parietal
*/
CC_Homotopic_for_filter_SPGWM.filter{it[1]=='SPGWM'}.set{CC_SPGWM_for_combine_parietal}
CC_Homotopic_for_filter_SMGWM.filter{it[1]=='SMGWM'}.set{CC_SMGWM_for_combine_parietal}
CC_Homotopic_for_filter_PrCuGWM.filter{it[1]=='PrCuGWM'}.set{CC_PrCuGWM_for_combine_parietal}
CC_Homotopic_for_filter_PoCGWM.filter{it[1]=='PoCGWM'}.set{CC_PoCGWM_for_combine_parietal}
CC_Homotopic_for_filter_AGWM.filter{it[1]=='AGWM'}.set{CC_AGWM_for_combine_parietal}

CC_SPGWM_for_combine_parietal.concat(CC_SMGWM_for_combine_parietal).concat(CC_PrCuGWM_for_combine_parietal).concat(CC_PoCGWM_for_combine_parietal).concat(CC_AGWM_for_combine_parietal).groupTuple(by:0).set{CC_Homotopic_parietal_for_rename}


/*
Filter CC Cingulum
*/
CC_Homotopic_for_filter_CingGWM.filter{it[1]=='CingGWM'}.set{CC_Homotopic_cingulum_for_rename}

/*
Filter CC Ins
*/
CC_Homotopic_for_filter_Ins.filter{it[1]=='Ins'}.set{CC_Homotopic_insular_for_rename}


/*
MERGE CC_Homotopic
*/
CC_Homotopic_for_merge.groupTuple().map{it}.set{CC_Homotopic_list_for_merge}

process CC_Homotopic_merge {
  cpus 1

input:
  set sid, file(tractogram) from CC_Homotopic_list_for_merge

output:
  set sid, "${sid}__CC_homo.trk" into CC_homo_for_trk_plausible, CC_homo_for_renaming

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

  input:
    set sid, val(side), file(trk01), file(trk02), file(trk03) from asso_all_intra_inter_ventral_all_for_merge

  output:
    set sid, "${sid}__asso_all_ventral_f_${side}.trk" into asso_all_ventral_for_trk_plausible
    set sid, val(side), "${sid}__asso_all_ventral_f_${side}.trk" into asso_all_ventral_for_split_ifof_uf
  script:
  """
  scil_streamlines_math.py concatenate ${trk01} ${trk02} ${trk03} ${sid}__asso_all_ventral_f_${side}.trk -f
  """
}

process split_asso_ventral_ifof_uf {
  cpus 1

  input:
    set sid, val(side), file(tractogram) from asso_all_ventral_for_split_ifof_uf

  output:
    set sid, val(side), "${sid}__asso_IFOF_f_${side}.trk" into asso_IFOF_for_rename
    set sid, val(side), "${sid}__asso_UF_f_${side}.trk" into asso_UF_for_rename

  script:
  """
  scil_filter_tractogram.py ${tractogram} ${sid}__asso_IFOF_f_${side}.trk  \
                            --drawn_roi ${params.rois_folder}${params.atlas.unc_no_ifof}${side}.nii.gz any include -f

  scil_filter_tractogram.py ${tractogram} ${sid}__asso_UF_f_${side}.trk \
                            --drawn_roi ${params.rois_folder}${params.atlas.unc_no_ifof}${side}.nii.gz any exclude -f
  """
}

/*
  ASSO DORSAL
*/

asso_dorsal_f_p_lists=params.asso_dorsal_f_p_lists?.tokenize(',')

process asso_dorsal_f_p {
  cpus 1

  input:
    set sid, val(side), file(tractogram) from asso_all_intra_inter_for_dorsal_f_p_filtering
    each asso_list from asso_dorsal_f_p_lists

  output:
    set sid, val(side), "${sid}__asso_${asso_list}_${side}.trk" into asso_all_intra_inter_dorsal_f_p_for_merge
    set sid, val(side), val(asso_list), "${sid}__asso_${asso_list}_${side}.trk" into asso_all_intra_inter_dorsal_f_p_for_rename
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

  input:
    set sid, val(side), file(tractogram) from asso_all_intra_inter_for_dorsal_f_o_f_t_filtering
    each asso_list from asso_dorsal_f_o_f_t_list

  output:
    set sid, val(side), "${sid}__asso_${asso_list}_${side}.trk" into asso_all_intra_inter_dorsal_all_f_o_f_t_for_merge
    set sid, val(side), val(asso_list), "${sid}__asso_${asso_list}_${side}.trk" into asso_all_intra_inter_dorsal_all_f_T_for_filter, asso_all_intra_inter_dorsal_all_f_O_for_filter
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

asso_all_intra_inter_dorsal_all_f_T_for_filter.filter{it[2]=='T_dorsal_f'}.set{asso_all_intra_inter_dorsal_all_f_T_for_rename}
asso_all_intra_inter_dorsal_all_f_O_for_filter.filter{it[2]=='O_dorsal_f'}.set{asso_all_intra_inter_dorsal_all_f_O_for_rename}

asso_all_intra_inter_dorsal_all_f_p_for_merge.groupTuple(by:[0,1]).join(asso_all_intra_inter_dorsal_all_f_o_f_t_for_merge.groupTuple(by:[0,1]), by:[0,1]).map{it.flatten().toList()}.set{asso_all_intra_inter_dorsal_all_for_merge}

process merge_asso_dorsal {
  cpus 1

  input:
    set sid, val(side), file(trk01), file(trk02), file(trk03) from asso_all_intra_inter_dorsal_all_for_merge

  output:
    set sid, "${sid}__asso_all_dorsal_f_${side}.trk" into asso_all_dorsal_for_trk_plausible

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

  input:
    set sid, val(side), file(tractogram) from asso_intra_inter_p_o_list_for_merge

  output:
    set sid, "${sid}__asso_all_P_O_f_${side}.trk" into all_P_O_for_trk_plausible

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

  input:
    set sid, val(side), file(tractogram) from asso_intra_inter_p_t_list_for_merge

  output:
    set sid, "${sid}__asso_all_P_T_f_${side}.trk" into all_P_T_for_trk_plausible

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

  input:
    set sid, val(side), file(tractogram) from asso_intra_inter_o_t_list_for_merge

  output:
    set sid, "${sid}__asso_all_O_T_f_${side}.trk" into all_O_T_for_trk_plausible
    set sid, val(side), "${sid}__asso_all_O_T_f_${side}.trk" into all_O_T_for_rename

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

  input:
    set sid, val(side), file(tractogram) from asso_intra_inter_ins_list_for_merge

  output:
    set sid, "${sid}__asso_all_Ins_f_${side}.trk" into Ins_for_trk_plausible

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

  input:
    set sid, val(side), file(tractogram) from asso_all_intra_inter_for_cing_filtering

  output:
    set sid, "${sid}__asso_all_Cing_${side}.trk" into Cing_for_trk_plausible
    set sid, val(side), "${sid}__asso_all_Cing_${side}.trk" into Cing_for_rename
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

  input:
    set sid, val(side),  val(gyrus), file(tractogram) from asso_frontal_be_list_for_merge

  output:
    set sid, "${sid}_asso_all_intraF_be_f_${side}_u.trk" into asso_all_intraF_be_for_trk_plausible

  script:
  """
    scil_streamlines_math.py lazy_concatenate ${tractogram} ${sid}_asso_all_intraF_be_f_${side}_u.trk -f
  """
}

/*
 EE ASSO FRONTAL: extracting all streamlines with either ends in a frontal gyrus (U-shape > 20 mm)
*/

asso_frontal_ee_list = Channel.from(['SFG_MFG', 70],['SFG_IFG', 70], ['SFG_PrCG', 90], ['SFG_FrOrbG', 70], ['MFG_IFG', 70], ['MFG_PrCG', 110], ['MFG_FrOrbG', 60], ['IFG_PrCG', 110], ['IFG_FrOrbG', 60])
asso_all_intra_inter_for_ee_frontal_filtering.combine(asso_frontal_ee_list).set{asso_frontal_ee_for_extract}
process asso_ee_frontal_gyrus {
  cpus 1

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

  input:
    set sid, val(side),  val(gyrus), file(tractogram) from asso_frontal_ee_list_for_merge

  output:
    set sid, "${sid}_asso_all_intraF_ee_f_${side}_u.trk" into asso_all_intraF_ee_for_trk_plausible

  script:
  """
    scil_streamlines_math.py lazy_concatenate ${tractogram} ${sid}_asso_all_intraF_ee_f_${side}_u.trk -f
  """
}

/*
 BE ASSO OCCIPITAL: extracting all streamlines with both ends in a occipital gyrus (U-shape > 20 mm)
*/

asso_occipital_be_list=params.asso_occipital_be_lists?.tokenize(',')
process asso_be_occipital_gyrus {
  cpus 1

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

  input:
    set sid, val(side),  val(gyrus), file(tractogram) from asso_occipital_be_list_for_merge

  output:
    set sid, "${sid}_asso_all_intraO_be_f_${side}_u.trk" into asso_all_intraO_be_for_trk_plausible

  script:
  """
    scil_streamlines_math.py lazy_concatenate ${tractogram} ${sid}_asso_all_intraO_be_f_${side}_u.trk -f
  """
}

/*
 EE ASSO OCCIPITAL: extracting all streamlines with either ends in a occipital gyrus (U-shape > 20 mm)
*/

asso_occipital_ee_list = Channel.from(['MOG_SOG', 60],['MOG_IOG', 50], ['MOG_CuG', 60], ['SOG_CuG', 30], ['CuG_LG', 60])
asso_all_intra_inter_for_ee_occipital_filtering.combine(asso_occipital_ee_list).set{asso_occipital_ee_for_extract}
process asso_ee_occipital_gyrus {
  cpus 1

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

  input:
    set sid, val(side),  val(gyrus), file(tractogram) from asso_occipital_ee_list_for_merge

  output:
    set sid, "${sid}_asso_all_intraO_ee_f_${side}_u.trk" into asso_all_intraO_ee_for_trk_plausible

  script:
  """
    scil_streamlines_math.py lazy_concatenate ${tractogram} ${sid}_asso_all_intraO_ee_f_${side}_u.trk -f
  """
}

/*
 BE ASSO PARIETAL: extracting all streamlines with both ends in a parietal gyrus (U-shape > 20 mm)
*/

asso_parietal_be_list=params.asso_parietal_be_lists?.tokenize(',')
process asso_be_parietal_gyrus {
  cpus 1

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

  input:
    set sid, val(side),  val(gyrus), file(tractogram) from asso_parietal_be_list_for_merge

  output:
    set sid, "${sid}_asso_all_intraP_be_f_${side}_u.trk" into asso_all_intraP_be_for_trk_plausible

  script:
  """
    scil_streamlines_math.py lazy_concatenate ${tractogram} ${sid}_asso_all_intraP_be_f_${side}_u.trk -f
  """
}

/*
 EE ASSO PARIETAL: extracting all streamlines with either ends in a parietal gyrus (U-shape > 20 mm)
*/

asso_parietal_ee_list = Channel.from(['SPG_PoCG', 50], ['SPG_AG', 80], ['SPG_SMG', 70], ['SPG_PrCuG', 50], ['AG_PoCG', 10000], ['AG_SMG', 90], ['AG_PrCuG', 90] , ['SMG_PoCG', 60], ['SMG_PrCuG',100], ['PoCG_PrCuG', 80])
asso_all_intra_inter_for_ee_parietal_filtering.combine(asso_parietal_ee_list).set{asso_parietal_ee_for_extract}
process asso_ee_parietal_gyrus {
  cpus 1

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

  input:
    set sid, val(side),  val(gyrus), file(tractogram) from asso_parietal_ee_list_for_merge

  output:
    set sid, "${sid}_asso_all_intraP_ee_f_${side}.trk" into asso_all_intraP_ee_for_trk_plausible

  script:
  """
    scil_streamlines_math.py lazy_concatenate ${tractogram} ${sid}_asso_all_intraP_ee_f_${side}.trk -f
  """
}

/*
 BE ASSO TEMPORAL: extracting all streamlines with both ends in a temporal gyrus and merge (U-shape > 20 mm)
*/

asso_temporal_be_list=params.asso_temporal_be_lists?.tokenize(',')
process asso_be_temporal_gyrus {
  cpus 1

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

  input:
    set sid, val(side),  val(gyrus), file(tractogram) from asso_temporal_be_list_for_merge

  output:
    set sid, "${sid}_asso_all_intraT_be_f_${side}_u.trk" into asso_all_intraT_be_for_trk_plausible

  script:
  """
    scil_streamlines_math.py lazy_concatenate ${tractogram} ${sid}_asso_all_intraT_be_f_${side}_u.trk -f
  """
}

/*
 EE ASSO TEMPORAL: extracting all streamlines with either ends in a temporal gyrus and merge (U-shape > 20 mm)
*/

asso_temporal_ee_list = Channel.from(['STG_MTG', 60], ['STG_ITG',80], ['STG_Tpole',110], ['MTG_ITG',60], ['MTG_Tpole', 100000], ['ITG_Tpole', 60])
asso_all_intra_inter_for_ee_temporal_filtering.combine(asso_temporal_ee_list).set{asso_temporal_ee_for_extract}
process asso_ee_temporal_gyrus {
  cpus 1

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

  input:
    set sid, val(side),  val(gyrus), file(tractogram) from asso_temporal_ee_list_for_merge

  output:
    set sid, "${sid}_asso_all_intraT_ee_f_${side}.trk" into asso_all_intraT_ee_for_trk_plausible

  script:
  """
    scil_streamlines_math.py lazy_concatenate ${tractogram} ${sid}_asso_all_intraT_ee_f_${side}.trk -f
  """
}

fornix_for_trk_plausible.concat(cerebellum_for_trk_plausible,brainstem_for_trk_plausible,BG_ipsi_Thal_for_trk_plausible,BG_ipsi_Put_for_trk_plausible,BG_ipsi_Caud_for_trk_plausible,asso_u_shape_for_trk_plausible,CC_homo_for_trk_plausible,asso_all_dorsal_for_trk_plausible,asso_all_ventral_for_trk_plausible,all_P_O_for_trk_plausible,all_P_T_for_trk_plausible,all_O_T_for_trk_plausible,Ins_for_trk_plausible,Cing_for_trk_plausible,asso_all_intraF_be_for_trk_plausible,asso_all_intraF_ee_for_trk_plausible,asso_all_intraP_be_for_trk_plausible,asso_all_intraP_ee_for_trk_plausible,asso_all_intraO_be_for_trk_plausible,asso_all_intraO_ee_for_trk_plausible,asso_all_intraT_be_for_trk_plausible,asso_all_intraT_ee_for_trk_plausible).groupTuple(by: 0).set{merge_trk_plausible}

process merge_trk_plausible{
  cpus 1

  input:
    set sid, file(tractogram) from merge_trk_plausible

  output:
    set sid, "${sid}__plausible.trk" into plausible_for_extract_unplausible

  script:
  """
    scil_streamlines_math.py lazy_concatenate ${tractogram} ${sid}__plausible.trk -f
  """
}

in_tractogram_for_extract_unplausible.join(plausible_for_extract_unplausible).set{for_trk_unplausible}

process extract_trk_unplausible{
  cpus 1

  input:
    set sid, file(trk01), file(trk02) from for_trk_unplausible
  output:
    set sid, "${sid}__unplausible.trk"

  script:
  """
    scil_streamlines_math.py difference ${trk01} ${trk02} ${sid}__unplausible.trk -f
  """
}

/*
RENAME CC CC_Homotopic
*/
process rename_cc_homotopic {
  cpus 1

  input:
    set sid, val(list), file(trk01) from CC_Homotopic_frontal_for_rename
    set sid, val(list), file(trk02) from CC_Homotopic_occipital_for_rename
    set sid, val(list), file(trk03) from CC_Homotopic_temporal_for_rename
    set sid, val(list), file(trk04) from CC_Homotopic_parietal_for_rename
    set sid, val(list), file(trk05) from CC_Homotopic_insular_for_rename
    set sid, val(list), file(trk06) from CC_Homotopic_cingulum_for_rename
  output:
    set sid, "${sid}__cc_homotopic_frontal.trk"
    set sid, "${sid}__cc_homotopic_occipital.trk"
    set sid, "${sid}__cc_homotopic_temporal.trk"
    set sid, "${sid}__cc_homotopic_parietal.trk"
    set sid, "${sid}__cc_homotopic_insular.trk"
    set sid, "${sid}__cc_homotopic_cingulum.trk"

  script:
  """
  scil_streamlines_math.py lazy_concatenate ${trk01} "${sid}__cc_homotopic_frontal.trk" -f
  scil_streamlines_math.py lazy_concatenate ${trk02} "${sid}__cc_homotopic_occipital.trk" -f
  scil_streamlines_math.py lazy_concatenate ${trk03} "${sid}__cc_homotopic_temporal.trk" -f
  scil_streamlines_math.py lazy_concatenate ${trk04} "${sid}__cc_homotopic_parietal.trk" -f
  cp ${trk05} ${sid}__cc_homotopic_insular.trk -f
  cp ${trk06} ${sid}__cc_homotopic_cingulum.trk -f
  """
}

/*
RENAME CORTICO_STRIATE
*/
BG_ipsi_Caud_for_rename.concat(BG_ipsi_Put_for_rename).groupTuple(by:[0,1]).set{cortico_striate_for_rename}
process rename_cortico_striate {
  cpus 1

  input:
    set sid, val(side), file(tractogram) from cortico_striate_for_rename

  output:
    set sid, "${sid}__corticostriatal_${side}.trk"

  script:
  """
    scil_streamlines_math.py lazy_concatenate ${tractogram} ${sid}__corticostriatal_${side}.trk -f
  """
}

/*
RENAME Corona radiata
*/
process rename_coronaradiata {
  cpus 1

  input:
    set sid, val(side), file(tractogram) from BG_ipsi_Thal_for_rename

  output:
    set sid, val(side), "${sid}__coronaradiata_${side}.trk"

  script:
  """
    scil_streamlines_math.py lazy_concatenate ${tractogram} ${sid}__coronaradiata_${side}.trk -f
  """
}

/*
RENAME OPTICAL RADIATION
*/
process rename_optical_radiation {
  cpus 1

  input:
    set sid, val(side), val(list), file(tractogram) from optic_radiation_for_rename

  output:
    set sid, "${sid}__optical_radiation_${side}.trk"

  script:
  """
    scil_streamlines_math.py lazy_concatenate ${tractogram} ${sid}__optical_radiation_${side}.trk -f
  """
}

/*
RENAME U SHAPE
*/
process rename_ushape {
  cpus 1

  input:
    set sid, val(side), file(tractogram) from asso_u_shape_for_rename

  output:
    set sid, "${sid}__ushape_${side}.trk"

  script:
  """
    cp ${tractogram} ${sid}__ushape_${side}.trk
  """
}

/*
RENAME CING
*/
process rename_cing {
  cpus 1

  input:
    set sid, val(side), file(tractogram) from Cing_for_rename

  output:
    set sid, "${sid}__cing_${side}.trk"

  script:
  """
    cp ${tractogram} ${sid}__cing_${side}.trk
  """
}

/*
RENAME SLF
*/
asso_all_intra_inter_dorsal_all_f_O_for_rename.concat(asso_all_intra_inter_dorsal_f_p_for_rename).groupTuple(by:[0,1]).set{slf_for_rename}
process rename_slf {
  cpus 1

  input:
    set sid, val(side), val(list), file(tractogram) from slf_for_rename

  output:
    set sid, "${sid}__slf_${side}.trk"

  script:
  """
    scil_streamlines_math.py lazy_concatenate ${tractogram} ${sid}__slf_${side}.trk -f
  """
}

/*
RENAME AF
*/
process rename_af {
  cpus 1

  input:
    set sid, val(side), asso_list, file(tractogram) from asso_all_intra_inter_dorsal_all_f_T_for_rename

  output:
    set sid, "${sid}__af_${side}.trk"

  script:
  """
    cp ${tractogram} ${sid}__af_${side}.trk -f
  """
}

/*
RENAME Cortico-pontine_F
*/
process rename_corticopontine_F {
  cpus 1

  input:
    set sid, file(tractogram) from brainstem_corticopontine_frontal_for_rename
    each side from side_corticopontineF
  output:
    set sid, "${sid}__corticopontine_frontal_${side}.trk"

  script:
  """
  scil_filter_tractogram.py ${tractogram} ${sid}__corticopontine_frontal_${side}.trk --drawn_roi ${params.rois_folder}${params.atlas.frontal_side}${side}.nii.gz either_end include -f
  """
}
/*
RENAME cortico-pontine_POT
*/
process rename_corticopontine_POT {
  cpus 1

  input:
    set sid, file(tractogram) from brainstem_ee_corticopontine_parietotemporooccipital_for_rename
    each side from side_corticopontinePOT

  output:
    set sid, "${sid}__corticopontine_POT_${side}.trk"

  script:
  """
    scil_filter_tractogram.py ${tractogram} ${sid}__corticopontine_POT_${side}.trk --drawn_roi ${params.rois_folder}${params.atlas.parieto_temporo_occipital_side}${side}.nii.gz either_end include -f
  """
}

/*
RENAME Pyramidal tract (CST)
*/
process rename_cst {
  cpus 1

  input:
    set sid, file(tractogram) from brainstem_pyramidal_for_rename
    each side from side_cst

  output:
    set sid, "${sid}__cst_${side}.trk"

  script:
  """
    scil_filter_tractogram.py ${tractogram} ${sid}__cst_${side}.trk --drawn_roi ${params.rois_folder}${params.atlas.fronto_parietal_side}${side}.nii.gz either_end include -f
  """
}

/*
RENAME fornix
*/
process rename_fornix {
  cpus 1

  input:
    set sid, file(tractogram) from fornix_for_rename

  output:
    set sid, "${sid}__fornix.trk"

  script:
  """
    cp ${tractogram} ${sid}__fornix.trk -f
  """
}

/*
RENAME IFOF
*/
process rename_ifof {
  cpus 1

  input:
    set sid, val(side), file(tractogram) from asso_IFOF_for_rename

  output:
    set sid, "${sid}__ifof_${side}.trk"

  script:
  """
    cp ${tractogram} ${sid}__ifof_${side}.trk -f
  """
}

/*
RENAME UF
*/
process rename_uf {
  cpus 1

  input:
    set sid, val(side), file(tractogram) from asso_UF_for_rename

  output:
    set sid, "${sid}__uf_${side}.trk"

  script:
  """
    cp ${tractogram} ${sid}__uf_${side}.trk -f
  """
}

/*
RENAME ILF
*/
process rename_ilf {
  cpus 1

  input:
    set sid, val(side), file(tractogram) from all_O_T_for_rename

  output:
    set sid, "${sid}__ilf_${side}.trk"

  script:
  """
    cp ${tractogram} ${sid}__ilf_${side}.trk -f
  """
}

/*
RENAME BRAINSTEM
*/
process rename_brainstem {
  cpus 1

  input:
    set sid, file(tractogram) from brainstem_for_rename

  output:
    set sid, "${sid}__brainstem.trk"

  script:
  """
    cp ${tractogram} ${sid}__brainstem.trk -f
  """
}

/*
RENAME CEREBELLUM
*/
process rename_cerebellum {
  cpus 1

  input:
    set sid, file(tractogram) from cerebellum_for_rename

  output:
    set sid, "${sid}__cerebellum.trk"

  script:
  """
    cp ${tractogram} ${sid}__cerebellum.trk -f
  """
}
