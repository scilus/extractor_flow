#!/usr/bin/env nextflow

params.input = false
params.help = false
params.debug = true


if(params.help) {
    usage = file("$baseDir/USAGE")
    cpu_count = Runtime.runtime.availableProcessors()
    bindings = ["rois_folder":"$params.rois_folder",
                "FLF": "$params.FLF",
                "run_bet":"$params.run_bet",
                "orig":"$params.orig",
                "extended":"$params.extended",
                "keep_intermediate_steps":"$params.keep_intermediate_steps",
                "quick_registration": "$params.quick_registration",
                "cpu_count":"$cpu_count",
                "processes_bet_register_t1":"$params.processes_bet_register_t1",
                "processes_major_filtering":"$params.processes_major_filtering"]  

    engine = new groovy.text.SimpleTemplateEngine()
    template = engine.createTemplate(usage.text).make(bindings)
    print template.toString()
    return
    }

log.info "Extractor_flow pipeline"
log.info "==================="
log.info "Start time: $workflow.start"
log.info ""

workflow.onComplete {
    log.info "Pipeline completed at: $workflow.complete"
    log.info "Execution status: ${ workflow.success ? 'OK' : 'failed' }"
    log.info "Execution duration: $workflow.duration"
}

if (!params.keep_intermediate_steps) {
  log.info "Warning: You won't be able to resume your processing if you don't use the option --keep_intermediate_steps"
  log.info ""
}

if (params.input){
    log.info "Input: $params.input"
    root = file(params.input)
    in_tractogram = Channel
        .fromFilePairs("$root/**/*.trk",
                       size:1,
                       maxDepth:1,
                       flat: true) {it.parent.name}

    in_tractogram.into{check_trks;
                       in_tractogram_for_unplausible;
                       in_tractogram_for_transformation;
                       in_tractogram_for_mix}


    Channel
    .fromPath("$root/**/*_t1.nii.gz",
              maxDepth:1)
             .map{[it.parent.name, it]}
             .into{t1s_for_register;
                   t1s_for_register_back;
                   t1s_for_copy_to_orig;
                   check_t1s;
                   t1s_empty}
}
else {
    error "Error ~ Please use --input for the input data."
}

check_trks.count().into{check_subjects_number; number_subj_for_null_check}
check_t1s.count().into{number_t1s_for_compare; number_t1s_check_with_orig}

number_subj_for_null_check
.subscribe{a -> if (a == 0)
    error "Error ~ No subjects found. Please check the naming convention, your --input path."}

check_subjects_number
  .concat(number_t1s_for_compare)
  .toList()
  .subscribe{a, b -> if (a != b && b > 0)
  error "Error ~ Some subjects have a T1w and others don't.\n" +
        "Please be sure to have the same acquisitions for all subjects."}

if (params.orig){
    number_t1s_check_with_orig
      .subscribe{a -> if (a == 0)
      error "Error ~ You cannot use --orig without having any T1w in the orig space."}
}

sides = params.sides?.tokenize(',')
Channel.from(sides).into{sides_ipsi;
                         sides_split_CC_BG;
                         sides_split_BG_Thal;
                         sides_split_BG_Put;
                         sides_split_BG_Caud;
                         side_corticopontineF;
                         side_corticopontinePOT;
                         side_cst}

/* BEGINNING TRANSFO */

process Register_T1 {
    publishDir = params.final_output_mni_space
    cpus params.processes_bet_register_t1

    input:
    set sid, file(t1) from t1s_for_register

    output:
    set sid, "${sid}__output0GenericAffine.mat", "${sid}__output1InverseWarp.nii.gz", "${sid}__output1Warp.nii.gz" into transformation_for_trk
    file "${sid}__t1_${params.template_space}.nii.gz"
    file "${sid}__t1_bet_mask.nii.gz" optional true
    file "${sid}__t1_bet.nii.gz" optional true

    script:
    if (params.run_bet){
    """
        export ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS=1
        export OMP_NUM_THREADS=1
        export OPENBLAS_NUM_THREADS=1
        export ANTS_RANDOM_SEED=1234

        antsBrainExtraction.sh -d 3 -a $t1 -e $params.template_t1/t1_template.nii.gz\
            -o bet/ -m $params.template_t1/t1_brain_probability_map.nii.gz -u 0
        scil_image_math.py convert bet/BrainExtractionMask.nii.gz ${sid}__t1_bet_mask.nii.gz --data_type uint8
        scil_image_math.py multiplication $t1 ${sid}__t1_bet_mask.nii.gz ${sid}__t1_bet.nii.gz

        ${params.registration_script} -d 3 -m ${sid}__t1_bet.nii.gz -f ${params.rois_folder}${params.atlas.template} -n ${task.cpus} -o "${sid}__output" -t s
        mv ${sid}__outputWarped.nii.gz ${sid}__t1_${params.template_space}.nii.gz
    """
    }
    else{
    """
        export ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS=1
        export OMP_NUM_THREADS=1
        export OPENBLAS_NUM_THREADS=1
        export ANTS_RANDOM_SEED=1234

        ${params.registration_script} -d 3 -m ${t1} -f ${params.rois_folder}${params.atlas.template} -n ${task.cpus} -o "${sid}__output" -t s
        mv ${sid}__outputWarped.nii.gz ${sid}__t1_${params.template_space}.nii.gz
    """
    }
}

transformation_for_trk.into{transformation_for_trk_registration;
                            transformation_for_join_with_t1}
if (params.orig) {
    t1s_for_register_back
        .cross(transformation_for_join_with_t1)
        .map { [ it[0][0], it[0][1], it[1][1], it[1][2], it[1][3]] }
        .into{transformation_and_t1_for_transformation_to_orig;
             transformation_and_t1_for_transformation_to_orig_bundles}
}

transformation_for_trk_registration
    .cross(in_tractogram_for_transformation)
    .map { [ it[0][0], it[0][1], it[0][2], it[0][3], it[1][1] ] }
    .set{trk_and_template_for_transformation_to_template}


process Transform_TRK {
    publishDir = params.final_output_mni_space
    cpus 1

    input:
    set sid, file(transfo), file(inv_deformation), file(deformation), file(trk) from trk_and_template_for_transformation_to_template

    output:
    set sid, "${trk.getSimpleName()}_${params.template_space}.trk" into transformed_for_remove_out_not_JHU, transformed_for_unplausible

    script:
    """
    scil_apply_transform_to_tractogram.py $trk ${params.rois_folder}${params.atlas.template} ${transfo} ${trk.getSimpleName()}_${params.template_space}.trk --remove_invalid --inverse --in_deformation ${inv_deformation}
    """
}

/* END TRANSFO */

for_remove_invalid_streamlines = Channel.empty()
if (t1s_empty.count().get()==0){
  for_remove_invalid_streamlines = for_remove_invalid_streamlines.mix(in_tractogram_for_mix)
  in_tractogram_for_unplausible.into{trk_for_extract_first_unplausible; trk_for_extract_unplausible}
}
else{
  transformed_for_unplausible.into{trk_for_extract_first_unplausible; trk_for_extract_unplausible}
}

process Remove_invalid_streamlines {
    cpus 1

    input:
      set sid, file(tractogram) from for_remove_invalid_streamlines

    output:
      set sid, "${sid}__rm_invalid_streamlines.trk" into rm_invalid_for_remove_out_not_JHU
      file "${sid}__t1_mni_space.nii.gz"

    script:
    """
      scil_remove_invalid_streamlines.py ${tractogram} ${sid}__rm_invalid_streamlines.trk --cut_invalid --remove_single_point -f
      cp ${params.rois_folder}${params.atlas.template} ${sid}__t1_mni_space.nii.gz
    """
}

rm_invalid_for_remove_out_not_JHU.mix(transformed_for_remove_out_not_JHU).set{for_major_filtering}


process Major_filtering {
    cpus params.processes_major_filtering

    input:
      set sid, file(tractogram) from for_major_filtering

    output:
      set sid, "${sid}__wb_clean01.trk" into wb_for_extract_end_in_cerebellum, wb_for_extract_first_unplausible

    script:
    """
      scil_filter_tractogram_anatomically.py ${tractogram} \
        ${params.rois_folder}${params.atlas.JHU_8} \
        ${sid} \
        --minL ${params.min_streaminline_lenght} \
        --maxL ${params.max_streaminline_lenght} \
        -a ${params.loop_angle_threshold} \
        --csf_bin ${params.rois_folder}${params.atlas.shell_limits} \
        --processes ${params.processes_major_filtering}\
        -f
      mv ${sid}/${tractogram.getSimpleName()}_filtered.trk ${sid}__wb_clean01.trk
    """
}

trk_for_extract_first_unplausible.join(wb_for_extract_first_unplausible).set{unplausible_streamlines}


process Extract_first_unplausible{
  cpus 1

  input:
    set sid, file(tractogram1), file(tractogram2) from unplausible_streamlines

  output:
    set sid, "${sid}__unplausible_streamlines.trk" into unplausible_for_fornix

  script:
  """
  scil_tractogram_math.py difference ${tractogram1} \
                                      ${tractogram2} \
                                      ${sid}__unplausible_streamlines.trk;
  """
}


process Extract_fornix{
  cpus 1

  input:
    set sid, file(tractogram) from unplausible_for_fornix

  output:
    set sid, "${sid}__fornix_f.trk" into fornix_for_trk_plausible, fornix_for_rename
    file "${sid}__fornix_f.txt"
    file "${sid}__unplausible_streamlines_wo_fornix.trk" optional true

  script:
    filtering_list=params.FLF+"fx.txt"
    out_extension="fornix_f"
    remaining_extension="unplausible_streamlines_wo_fornix"
    basename="${sid}"
    keep="$params.keep_intermediate_steps"
    extract_masks=""

    template "filter_with_list.sh"
}


process Extract_ee_cerebellum {
  cpus 1

  input:
    set sid, file(tractogram) from wb_for_extract_end_in_cerebellum

  output:
    set sid, "${sid}__wb_clean01_nocereb.trk" into wb_for_extract_end_in_brainstem
    set sid, "${sid}__all_cerebellum.trk" into ee_cerebellum_for_extract_plausible
    file "${sid}__all_cerebellum.txt" optional true
    file "${sid}__wb_clean01_nocereb.txt" optional true

  script:
  filtering_list=params.FLF+"out_cerebellum.txt"
  out_extension="wb_clean01_nocereb"
  remaining_extension="all_cerebellum"
  basename="${sid}"
  keep=true
  extract_masks=""

  template "filter_with_list.sh"
}


process Extract_plausible_cerebellum {
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
  scil_filter_tractogram.py ${tractogram} ${sid}__tmp_in_cerebellum.trk\
        --filtering_list ${params.FLF}in_cerebellum.txt -f
  scil_filter_tractogram.py ${sid}__tmp_in_cerebellum.trk ${sid}__all_in_cerebellum_nocx_nocerebwm.trk\
        --filtering_list ${params.FLF}cerebellum_nocx_in_cereb.txt -f
  scil_filter_tractogram.py ${sid}__tmp_in_cerebellum.trk ${sid}__all_in_cerebellum_in_Medulla.trk\
        --filtering_list ${params.FLF}cerebellum_in_medulla.txt -f
  scil_filter_tractogram.py ${sid}__tmp_in_cerebellum.trk ${sid}__all_in_cerebellum_in_Pons.trk\
        --filtering_list ${params.FLF}cerebellum_in_pons.txt -f
  scil_filter_tractogram.py ${sid}__tmp_in_cerebellum.trk ${sid}__all_in_cerebellum_in_Midbrain.trk\
        --filtering_list ${params.FLF}cerebellum_in_midbrain.txt -f
  scil_filter_tractogram.py ${sid}__tmp_in_cerebellum.trk ${sid}__all_in_cerebellum_in_redN_and_Thal.trk\
        --filtering_list ${params.FLF}cerebellum_in_rednucleus_and_thalamus.txt -f
  scil_tractogram_math.py union ${sid}__all_in_*.trk ${sid}__all_cerebellum_plausibles.trk -f
  """
}

/*
  END Cerebellum
*/

process Extract_ee_brainstem {
  cpus 1

  input:
    set sid, file(tractogram) from wb_for_extract_end_in_brainstem

  output:
    set sid, "${sid}__wb_clean02.trk" into wb_for_split_end_in_CGMSWI
    set sid, "${sid}__all_brainstem.trk" into all_brainstem_for_extract_plausible
    file "${sid}__wb_clean02.txt" optional true
    file "${sid}__all_brainstem.txt" optional true

  script:
    filtering_list=params.FLF+"out_brainstem.txt"
    out_extension="wb_clean02"
    remaining_extension="all_brainstem"
    basename="${sid}"
    keep=true
    extract_masks=""

    template "filter_with_list.sh"
}

/*
  Brainstem
*/

process Extract_plausible_brainstem {
  cpus 1

  input:
    set sid, file(tractogram) from all_brainstem_for_extract_plausible
  output:
    set sid, "${sid}__all_brainstem_plausibles.trk" into brainstem_for_trk_plausible, brainstem_for_rename
    file "${sid}__all_brainstem_unplausibles.trk" optional true
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
  scil_filter_tractogram.py ${sid}__all_brainstem.trk ${sid}__be_midbrain.trk\
      --filtering_list ${params.FLF}brainstem_be_midbrain.txt -f
  # Extract be medulla
  scil_filter_tractogram.py ${sid}__all_brainstem.trk ${sid}__be_medulla.trk\
      --filtering_list ${params.FLF}brainstem_be_medulla.txt -f
  # Extract be pons
  scil_filter_tractogram.py ${sid}__all_brainstem.trk ${sid}__be_pons.trk\
      --filtering_list ${params.FLF}brainstem_be_pons.txt -f

  # Extract ee thalamus
  scil_filter_tractogram.py ${sid}__all_brainstem.trk ${sid}__ee_thalamus.trk\
      --filtering_list ${params.FLF}brainstem_ee_thalamus.txt -f
  # Extract ee red_nucleus
  scil_filter_tractogram.py ${sid}__all_brainstem.trk ${sid}__ee_red_nucleus.trk\
      --filtering_list ${params.FLF}brainstem_ee_red_nucleus.txt -f

  # Prepartion for fronto-pontine, parietotemporooccipito-pontine, pyramidal, cortico-tectal
  scil_filter_tractogram.py ${sid}__all_brainstem.trk ${sid}__ee_tmp_01.trk\
      --filtering_list ${params.FLF}brainstem_ee_tmp_01.txt -f
  scil_filter_tractogram.py ${sid}__all_brainstem.trk ${sid}__ee_tmp_02.trk\
      --filtering_list ${params.FLF}brainstem_ee_tmp_02.txt -f

  scil_tractogram_math.py union ${sid}__ee_tmp_01.trk ${sid}__ee_tmp_02.trk\
      ${sid}__ee_tmp_03.trk -f

  # Extract ee Fronto-pontine
  scil_filter_tractogram.py ${sid}__ee_tmp_03.trk ${sid}__ee_fronto_pontine.trk\
      --filtering_list ${params.FLF}brainstem_ee_F_pontine.txt -f

  # Extract ee ParietoTemporooccipital pontine
  scil_filter_tractogram.py ${sid}__ee_tmp_03.trk ${sid}__ee_parietotemporooccipital_pontine.trk\
      --filtering_list ${params.FLF}brainstem_ee_PTO_pontine.txt -f

  # Extract ee Pyramidal
  scil_filter_tractogram.py ${sid}__ee_tmp_03.trk ${sid}__ee_pyramidal.trk\
      --filtering_list ${params.FLF}brainstem_ee_pyramidal.txt -f

  # Extract ee Tectal
  scil_filter_tractogram.py ${sid}__ee_tmp_03.trk ${sid}__ee_cortico_tectal.trk\
      --filtering_list ${params.FLF}brainstem_ee_cortico_tectal.txt -f
  scil_filter_streamlines_by_length.py ${sid}__ee_cortico_tectal.trk ${sid}__ee_cortico_tectal.trk --maxL 100 -f

  rm -f ${sid}__*tmp_*.trk

  scil_tractogram_math.py union ${sid}__be_*.trk ${sid}__ee_*.trk ${sid}__all_brainstem_plausibles.trk -f

  if ${params.keep_intermediate_steps}
  then
    scil_tractogram_math.py difference ${sid}__all_brainstem.trk ${sid}__all_brainstem_plausibles.trk ${sid}__all_brainstem_unplausibles.trk -f
  fi
  """
}

/*
Brain - Either end in CGM SWM
*/

process Remove_out_of_CGM_DWM {
  cpus 1

  input:
    set sid, file(tractogram) from wb_for_split_end_in_CGMSWI

  output:
    set sid, "${sid}__wb_either_CGM_SWM.trk" into wb_for_extract_all_commissural
    set sid, "${sid}__no_CGM_SWM.trk" optional true
    file "${sid}__wb_either_CGM_SWM.txt" optional true
    file "${sid}__no_CGM_SWM.txt" optional true

  script:
    filtering_list=params.FLF+"ee_CGM_SWM.txt"
    out_extension="wb_either_CGM_SWM"
    remaining_extension="no_CGM_SWM"
    basename="${sid}"
    keep="$params.keep_intermediate_steps"
    extract_masks=""

    template "filter_with_list.sh"
}

process Extract_all_commissural {
  cpus 1

  input:
    set sid, file(tractogram) from wb_for_extract_all_commissural

  output:
    set sid, "${sid}__tmp_CC.trk" into cc_for_extract_CC_Cx, cc_for_extract_AC_Cx, cc_for_extract_CC_BG, cc_tmp_for_commissural
    set sid, "${sid}__wb_either_CGM_SWM_noCC.trk" into no_cc_for_split_asso_BG
    file "${sid}__wb_either_CGM_SWM_noCC.txt" optional true
    file "${sid}__tmp_CC.txt" optional true

  script:
  filtering_list=params.FLF+"commissural.txt"
  out_extension="wb_either_CGM_SWM_noCC"
  remaining_extension="tmp_CC"
  basename="${sid}"
  keep=true
  extract_masks=""

  template "filter_with_list.sh"
}


process Extract_plausible_CC_Cx {
  cpus 1

  input:
    set sid, file(tractogram) from cc_for_extract_CC_Cx

  output:
    set sid, "${sid}__in_CC_Cx_f.trk" into cc_for_merge_plausible_01
    file "mask_atlas_roi_*.nii.gz" optional true

  script:
    filtering_list=params.FLF+"CC_Cx.txt"
    out_extension="in_CC_Cx_f"
    remaining_extension="garbage"
    basename="${sid}"
    keep="$params.keep_intermediate_steps"
    extract_masks=""

    template "filter_with_list.sh"
}

process Extract_plausible_AC_Cx {
  cpus 1

  input:
    set sid, file(tractogram) from cc_for_extract_AC_Cx

  output:
    set sid, "${sid}__in_AC_Cx_f.trk" into accx_for_trk_plausible, accx_for_rename, accx_for_commissural

  script:
    filtering_list=params.FLF+"AC_Cx.txt"
    out_extension="in_AC_Cx_f"
    remaining_extension="garbage"
    basename="${sid}"
    keep="$params.keep_intermediate_steps"
    extract_masks=""

      template "filter_with_list.sh"
}

process Extract_plausible_CC_BG {
  cpus 1

  input:
    set sid, file(tractogram) from cc_for_extract_CC_BG

  output:
    set sid, "${sid}__in_CC_BG_f.trk" into ccbg_for_trk_plausible, ccbg_for_commissural

  script:
    filtering_list=params.FLF+"CC_BG.txt"
    out_extension="in_CC_BG_f"
    remaining_extension="garbage"
    basename="${sid}"
    keep="$params.keep_intermediate_steps"
    extract_masks=""

    template "filter_with_list.sh"
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
  filtering_list=params.FLF+"all_BG.txt"
  out_extension="all_subcortical_from_CGM_SWM_noCC_f"
  remaining_extension="asso_noBG"
  basename="${sid}"
  keep=true
  extract_masks=""

  template "filter_with_list.sh"
}

bg_list=params.bg_lists?.tokenize(',')
Channel.from(bg_list).into{bg_thal_list;
                           bg_put_list}
/*
BG THAL
*/
process Split_BG_Thal {
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
    filtering_list=params.FLF+"BG_ipsi_Thal_${list}_${side}.txt"
    out_extension="BG_ipsi_Thal_${list}_${side}"
    remaining_extension="garbage_BG_ipsi_Thal_${list}_${side}"
    basename="${sid}"
    keep="$params.keep_intermediate_steps"
    extract_masks=""

    template "filter_with_list.sh"
}

BG_ipsi_Thal_split_for_merge.groupTuple(by:[0,1]).set{BG_ipsi_Thal_for_rename}

BG_ipsi_Thal_for_filter_CuGWM.filter{it[2]=='CuGWM'}.set{CuGWM_for_combine}
BG_ipsi_Thal_for_filter_LGWM.filter{it[2]=='LGWM'}.set{LGWM_for_combine}
CuGWM_for_combine.concat(LGWM_for_combine).groupTuple(by:[0,1]).set{optic_radiation_for_rename}

BG_ipsi_Thal_for_merge.groupTuple().map{it}.set{BG_ipsi_Thal_list_for_merge}

process Merge_BG_Thal{
  cpus 1

  input:
    set sid, file(tractogram) from BG_ipsi_Thal_list_for_merge

  output:
    set sid, "${sid}__BG_ipsi_Thal_all.trk" into BG_ipsi_Thal_for_trk_plausible

  script:
  """
  scil_tractogram_math.py union ${tractogram} ${sid}__BG_ipsi_Thal_all.trk -f
  """
}

/*
BG PUT
*/
process Split_BG_Put {
  cpus 1

  input:
    set sid, file(tractogram) from asso_BG_for_split_Put
    each list from bg_put_list
    each side from sides_split_BG_Put

  output:
    set sid, "${sid}__BG_ipsi_Put_${list}_${side}.trk" into BG_ipsi_Put_for_merge
    set sid, val(side), "${sid}__BG_ipsi_Put_${list}_${side}.trk" into BG_ipsi_Put_for_rename

  script:
    filtering_list=params.FLF+"BG_ipsi_Put_${list}_${side}.txt"
    out_extension="BG_ipsi_Put_${list}_${side}"
    remaining_extension="garbage_BG_ipsi_Put_${list}_${side}"
    basename="${sid}"
    keep="true"
    extract_masks=""

    template "filter_with_list.sh"
}

BG_ipsi_Put_for_merge.groupTuple().map{it}.set{BG_ipsi_Put_list_for_merge}

process Merge_BG_Put{
  cpus 1

  input:
    set sid, file(tractogram) from BG_ipsi_Put_list_for_merge

  output:
    set sid, "${sid}__BG_ipsi_Put_all.trk" into BG_ipsi_Put_for_trk_plausible

  script:
  """
  scil_tractogram_math.py union ${tractogram} ${sid}__BG_ipsi_Put_all.trk -f
  """
}

/*
BG CAUD
*/
bg_caud_list=params.bg_caud_lists?.tokenize(',')
process Split_BG_Caud {
  cpus 1

  input:
    set sid, file(tractogram) from asso_BG_for_split_Caud
    each list from bg_caud_list
    each side from sides_split_BG_Caud

  output:
    set sid, "${sid}__BG_ipsi_Caud_${list}_${side}.trk" into BG_ipsi_Caud_for_merge
    set sid, val(side), "${sid}__BG_ipsi_Caud_${list}_${side}.trk" into BG_ipsi_Caud_for_rename

  script:
    filtering_list=params.FLF+"BG_ipsi_Caud_${list}_${side}.txt"
    out_extension="BG_ipsi_Caud_${list}_${side}"
    remaining_extension="garbage_BG_ipsi_Caud_${list}_${side}"
    basename="${sid}"
    keep="true"
    extract_masks=""

    template "filter_with_list.sh"
}

BG_ipsi_Caud_for_merge.groupTuple().map{it}.set{BG_ipsi_Caud_list_for_merge}

process Merge_BG_Caud{
  cpus 1

  input:
    set sid, file(tractogram) from BG_ipsi_Caud_list_for_merge

  output:
    set sid, "${sid}__BG_ipsi_Caud_all.trk" into BG_ipsi_Caud_for_trk_plausible

  script:
  """
  scil_tractogram_math.py union ${tractogram} ${sid}__BG_ipsi_Caud_all.trk -f
  """
}

process Split_asso_in_hemi {
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
  """
  scil_filter_tractogram.py ${tractogram} ${sid}__asso_L.trk\
   --filtering_list ${params.FLF}asso_L.txt -f
   scil_filter_tractogram.py ${tractogram} ${sid}__asso_R.trk\
    --filtering_list ${params.FLF}asso_R.txt -f

  scil_tractogram_math.py difference ${tractogram} ${sid}__asso_L.trk ${sid}__asso_R.trk ${sid}__asso_lost.trk -f
  """
}

/*
Extracting U-shaped and streamlines restricted to Cortical GM and removing them from asso
*/

process Split_ushape_CGM_asso {
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
      --filtering_list ${params.FLF}all_in_CGM_${side}.txt -f

    scil_tractogram_math.py difference ${tractogram} ${sid}__tmp1_${side}.trk \
                             ${sid}__asso_SWM_${side}.trk -f

    scil_filter_tractogram.py ${sid}__tmp1_${side}.trk ${sid}__asso_only_in_CGM_${side}.trk \
      --filtering_list ${params.FLF}not_in_SWM_${side}.txt

    scil_tractogram_math.py difference ${sid}__tmp1_${side}.trk ${sid}__asso_only_in_CGM_${side}.trk \
                                 ${sid}__tmp2_${side}.trk -f

    scil_filter_tractogram.py ${sid}__tmp2_${side}.trk ${sid}__asso_Ushape_${side}.trk \
      --filtering_list ${params.FLF}not_in_DWM_${side}.txt

    scil_extract_ushape.py ${sid}__asso_Ushape_${side}.trk --minU 0.5 --maxU 1 ${sid}__asso_Ushape_${side}_u.trk -f

    scil_tractogram_math.py difference ${sid}__tmp2_${side}.trk ${sid}__asso_Ushape_${side}.trk \
                               ${sid}__asso_DWM_${side}.trk -f

    scil_tractogram_math.py union ${sid}__asso_DWM_${side}.trk ${sid}__asso_SWM_${side}.trk ${sid}__asso_f_${side}.trk -f

    if ${params.keep_intermediate_steps}
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
    set sid, "${sid}__asso_lost2_${side}.trk" optional true
    file "${sid}__asso_all_intra_inter_${side}.txt" optional true
    file "${sid}__asso_lost2_${side}.txt" optional true

  script:
    filtering_list=params.FLF+"not_in_BG.txt"
    out_extension="asso_all_intra_inter_${side}"
    remaining_extension="asso_lost2_${side}"
    basename="${sid}"
    keep="$params.keep_intermediate_steps"
    extract_masks=""

    template "filter_with_list.sh"
}

/*
inCCBG.groupTuple().map{it.flatten().toList()}.set{inCCBG_List}
assoUShape.groupTuple().map{it.flatten().toList()}.set{assoUShape_list}
*/

asso_all_intra_inter.into{asso_all_intra_inter_for_ventral_filtering;
                          asso_all_intra_inter_for_dorsal_f_p_filtering;
                          asso_all_intra_inter_for_dorsal_f_o_f_t_filtering;
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
    filtering_list=params.FLF+"CC_homo_${pair}.txt"
    out_extension="cc_homotopic_${pair}"
    remaining_extension="garbage_${pair}"
    basename="${sid}"
    keep="$params.keep_intermediate_steps"
    extract_masks=""

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
  set sid, "${sid}__CC_homo.trk" into CC_homo_for_trk_plausible, CC_homo_for_renaming, cc_homo_for_commissural

script:
  """
  scil_tractogram_math.py union ${tractogram} ${sid}__CC_homo.trk
  """
}

/*
  COMMISSURAL
*/

cc_tmp_for_commissural.join(accx_for_commissural).join(ccbg_for_commissural).join(cc_homo_for_commissural).into{all_cc_for_commissural;toto}
toto.println()

process CC_all_commissural {
  cpus 1

  input:
    set sid, file(tmp_cc), file(accx), file(ccbg), file(cc_homo) from all_cc_for_commissural

  output:
    set sid, "${sid}__plausible_commissural.trk" into plausible_commissural_for_register_to_orig
    file "${sid}__unplausible_commissural.trk" optional true

  script:
  """
    scil_tractogram_math.py union ${accx} ${ccbg} ${cc_homo} ${sid}__plausible_commissural.trk -f

    if ${params.keep_intermediate_steps}
    then
      scil_tractogram_math.py difference ${tmp_cc} ${sid}__plausible_commissural.trk ${sid}__unplausible_commissural.trk -f
    fi
  """
}

/*
  ASSO VENTRAL
*/

asso_ventral_lists=params.asso_ventral_lists?.tokenize(',')

process Asso_ventral {
  cpus 1

  input:
    set sid, val(side), file(tractogram) from asso_all_intra_inter_for_ventral_filtering
    each asso_list from asso_ventral_lists

  output:
    set sid, val(side), "${sid}__asso_F_${asso_list}_ventral_f_${side}.trk" into asso_all_intra_inter_ventral_for_merge

  script:
    filtering_list=params.FLF+"ASSO_F_${asso_list}_ventral_${side}.txt"
    out_extension="asso_F_${asso_list}_ventral_f_${side}"
    remaining_extension="asso_lost_${side}"
    basename="${sid}"
    keep=false
    extract_masks=""

    template "filter_with_list.sh"
}

asso_all_intra_inter_ventral_for_merge.groupTuple(by:[0,1]).map{it.flatten().toList()}.set{asso_all_intra_inter_ventral_all_for_merge}

process Merge_asso_ventral {
  cpus 1

  input:
    set sid, val(side), file(trk01), file(trk02), file(trk03) from asso_all_intra_inter_ventral_all_for_merge

  output:
    set sid, "${sid}__asso_all_ventral_f_${side}.trk" into asso_all_ventral_for_trk_plausible
    set sid, val(side), "${sid}__asso_all_ventral_f_${side}.trk" into asso_all_ventral_for_split_ifof_uf

  script:
  """
  scil_tractogram_math.py union ${trk01} ${trk02} ${trk03} ${sid}__asso_all_ventral_f_${side}.trk -f
  """
}

process Split_asso_ventral_ifof_uf {
  cpus 1

  input:
    set sid, val(side), file(tractogram) from asso_all_ventral_for_split_ifof_uf

  output:
    set sid, val(side), "${sid}__asso_IFOF_f_${side}.trk" into asso_IFOF_for_rename
    set sid, val(side), "${sid}__asso_UF_f_${side}.trk" into asso_UF_for_rename

  script:
    filtering_list=params.FLF+"split_IFOF_UF_${side}.txt"
    out_extension="asso_IFOF_f_${side}"
    remaining_extension="asso_UF_f_${side}"
    basename="${sid}"
    keep=true
    extract_masks=""

    template "filter_with_list.sh"
}

/*
  ASSO DORSAL
*/

asso_dorsal_f_p_lists=params.asso_dorsal_f_p_lists?.tokenize(',')

process Asso_dorsal_f_p {
  cpus 1

  input:
    set sid, val(side), file(tractogram) from asso_all_intra_inter_for_dorsal_f_p_filtering
    each asso_list from asso_dorsal_f_p_lists

  output:
    set sid, val(side), "${sid}__asso_${asso_list}_${side}.trk" into asso_all_intra_inter_dorsal_f_p_for_merge
    set sid, val(side), val(asso_list), "${sid}__asso_${asso_list}_${side}.trk" into asso_all_intra_inter_dorsal_f_p_for_rename
    set sid, "${sid}__asso_lost_${asso_list}_${side}.trk" optional true
    file "${sid}__asso_${asso_list}_${side}.txt" optional true
    file "${sid}__asso_lost_${asso_list}_${side}.txt" optional true

  script:
    filtering_list=params.FLF+"ASSO_${asso_list}_${side}.txt"
    out_extension="asso_${asso_list}_${side}"
    remaining_extension="asso_lost_${asso_list}_${side}"
    basename="${sid}"
    keep="$params.keep_intermediate_steps"
    extract_masks=""

    template "filter_with_list.sh"
}

asso_all_intra_inter_dorsal_f_p_for_merge.groupTuple(by:[0,1]).map{it}.set{asso_all_intra_inter_dorsal_f_p_list_for_merge}

process Merge_asso_dorsal_f_p {
  cpus 1

  input:
    set sid, val(side), file(tractogram) from asso_all_intra_inter_dorsal_f_p_list_for_merge

  output:
    set sid, val(side), "${sid}__asso_F_P_dorsal_f_${side}.trk" into asso_all_intra_inter_dorsal_all_f_p_for_merge

  script:
  """
  scil_tractogram_math.py union ${tractogram} ${sid}__asso_F_P_dorsal_f_${side}.trk -f
  """
}

asso_dorsal_f_o_f_t_list=params.asso_dorsal_f_o_f_t_lists?.tokenize(',')

process Asso_dorsal_f_o_f_t {
  cpus 1

  input:
    set sid, val(side), file(tractogram) from asso_all_intra_inter_for_dorsal_f_o_f_t_filtering
    each asso_list from asso_dorsal_f_o_f_t_list

  output:
    set sid, val(side), "${sid}__asso_${asso_list}_${side}.trk" into asso_all_intra_inter_dorsal_all_f_o_f_t_for_merge
    set sid, val(side), val(asso_list), "${sid}__asso_${asso_list}_${side}.trk" into asso_all_intra_inter_dorsal_all_f_T_for_filter, asso_all_intra_inter_dorsal_all_f_O_for_filter
    set sid, "${sid}__asso_lost_${asso_list}_${side}.trk" optional true
    file "${sid}__asso_${asso_list}_${side}.txt" optional true
    file "${sid}__asso_lost_${asso_list}_${side}.txt" optional true

  script:
    filtering_list=params.FLF+"ASSO_${asso_list}_${side}.txt"
    out_extension="asso_${asso_list}_${side}"
    remaining_extension="asso_lost_${asso_list}_${side}"
    basename="${sid}"
    keep="$params.keep_intermediate_steps"
    extract_masks=""

    template "filter_with_list.sh"
}

asso_all_intra_inter_dorsal_all_f_T_for_filter.filter{it[2]=='F_T_dorsal'}.set{asso_all_intra_inter_dorsal_all_f_T_for_rename}
asso_all_intra_inter_dorsal_all_f_O_for_filter.filter{it[2]=='F_O_dorsal'}.set{asso_all_intra_inter_dorsal_all_f_O_for_rename}

asso_all_intra_inter_dorsal_all_f_p_for_merge.groupTuple(by:[0,1]).join(asso_all_intra_inter_dorsal_all_f_o_f_t_for_merge.groupTuple(by:[0,1]), by:[0,1]).map{it.flatten().toList()}.set{asso_all_intra_inter_dorsal_all_for_merge}

process Merge_asso_dorsal {
  cpus 1

  input:
    set sid, val(side), file(trk01), file(trk02), file(trk03) from asso_all_intra_inter_dorsal_all_for_merge

  output:
    set sid, "${sid}__asso_all_dorsal_f_${side}.trk" into asso_all_dorsal_for_trk_plausible

  script:
  """
  scil_tractogram_math.py union ${trk01} ${trk02} ${trk03} ${sid}__asso_all_dorsal_f_${side}.trk -f
  """
}

/*
  ASSO P_O
*/

asso_p_o_list=params.asso_p_o_lists?.tokenize(',')

process Asso_p_o {
  cpus 1

  input:
    set sid, val(side), file(tractogram) from asso_all_intra_inter_for_p_o_filtering
    each asso_list from asso_p_o_list

  output:
    set sid, val(side), "${sid}__asso_${asso_list}_${side}.trk" into asso_intra_inter_p_o_for_merge
    set sid, "${sid}__asso_lost_${asso_list}_${side}.trk" optional true
    file "${sid}__asso_${asso_list}_${side}.txt" optional true
    file "${sid}__asso_lost_${asso_list}_${side}.txt" optional true

  script:
    filtering_list=params.FLF+"ASSO_${asso_list}_${side}.txt"
    out_extension="asso_${asso_list}_${side}"
    remaining_extension="asso_lost_${asso_list}_${side}"
    basename="${sid}"
    keep="$params.keep_intermediate_steps"
    extract_masks=""

    template "filter_with_list.sh"
}

asso_intra_inter_p_o_for_merge.groupTuple(by:[0,1]).map{it}.set{asso_intra_inter_p_o_list_for_merge}

process Merge_p_o {
  cpus 1

  input:
    set sid, val(side), file(tractogram) from asso_intra_inter_p_o_list_for_merge

  output:
    set sid, "${sid}__asso_all_P_O_f_${side}.trk" into all_P_O_for_trk_plausible

  script:
  """
  scil_tractogram_math.py union ${tractogram} ${sid}__asso_all_P_O_f_${side}.trk -f
  """
}

/*
  ASSO P_T
*/

asso_p_t_list=params.asso_p_t_lists?.tokenize(',')

process Asso_p_t {
  cpus 1

  input:
    set sid, val(side), file(tractogram) from asso_all_intra_inter_for_p_t_filtering
    each asso_list from asso_p_t_list

  output:
    set sid, val(side), "${sid}__asso_${asso_list}_${side}.trk" into asso_intra_inter_p_t_for_merge
    set sid, "${sid}__asso_lost_${asso_list}_${side}.trk" optional true
    file "${sid}__asso_${asso_list}_${side}.txt" optional true
    file "${sid}__asso_lost_${asso_list}_${side}.txt" optional true

  script:
    filtering_list=params.FLF+"ASSO_${asso_list}_${side}.txt"
    out_extension="asso_${asso_list}_${side}"
    remaining_extension="asso_lost_${asso_list}_${side}"
    basename="${sid}"
    keep="$params.keep_intermediate_steps"
    extract_masks=""

    template "filter_with_list.sh"
}

asso_intra_inter_p_t_for_merge.groupTuple(by:[0,1]).map{it}.set{asso_intra_inter_p_t_list_for_merge}

process Merge_p_t {
  cpus 1

  input:
    set sid, val(side), file(tractogram) from asso_intra_inter_p_t_list_for_merge

  output:
    set sid, "${sid}__asso_all_P_T_f_${side}.trk" into all_P_T_for_trk_plausible

  script:
  """
  scil_tractogram_math.py union ${tractogram} ${sid}__asso_all_P_T_f_${side}.trk -f
  """
}

/*
  ASSO O_T
*/

asso_o_t_list=params.asso_o_t_lists?.tokenize(',')

process Asso_o_t {
  cpus 1

  input:
    set sid, val(side), file(tractogram) from asso_all_intra_inter_for_o_t_filtering
    each asso_list from asso_o_t_list

  output:
    set sid, val(side), "${sid}__asso_${asso_list}_${side}.trk" into asso_intra_inter_o_t_for_merge
    set sid, "${sid}__asso_lost_${asso_list}_${side}.trk" optional true
    file "${sid}__asso_${asso_list}_${side}.txt" optional true
    file "${sid}__asso_lost_${asso_list}_${side}.txt" optional true

  script:
    filtering_list=params.FLF+"ASSO_${asso_list}_${side}.txt"
    out_extension="asso_${asso_list}_${side}"
    remaining_extension="asso_lost_${asso_list}_${side}"
    basename="${sid}"
    keep="$params.keep_intermediate_steps"
    extract_masks=""

    template "filter_with_list.sh"
}

asso_intra_inter_o_t_for_merge.groupTuple(by:[0,1]).map{it}.set{asso_intra_inter_o_t_list_for_merge}

process Merge_o_t {
  cpus 1

  input:
    set sid, val(side), file(tractogram) from asso_intra_inter_o_t_list_for_merge

  output:
    set sid, "${sid}__asso_all_O_T_f_${side}.trk" into all_O_T_for_trk_plausible
    set sid, val(side), "${sid}__asso_all_O_T_f_${side}.trk" into all_O_T_for_rename

  script:
  """
  scil_tractogram_math.py union ${tractogram} ${sid}__asso_all_O_T_f_${side}.trk -f
  """
}

/*
  ASSO Ins
*/

asso_ins_list=params.asso_ins_lists?.tokenize(',')

process Asso_ins {
  cpus 1

  input:
    set sid, val(side), file(tractogram) from asso_all_intra_inter_for_ins_filtering
    each asso_list from asso_ins_list

  output:
    set sid, val(side), "${sid}__asso_${asso_list}_${side}.trk" into asso_intra_inter_ins_for_merge
    set sid, "${sid}__asso_lost_${asso_list}_${side}.trk" optional true
    file "${sid}__asso_${asso_list}_${side}.txt" optional true
    file "${sid}__asso_lost_${asso_list}_${side}.txt" optional true

  script:
    filtering_list=params.FLF+"ASSO_${asso_list}_${side}.txt"
    out_extension="asso_${asso_list}_${side}"
    remaining_extension="asso_lost_${asso_list}_${side}"
    basename="${sid}"
    keep="$params.keep_intermediate_steps"
    extract_masks=""

    template "filter_with_list.sh"
}

asso_intra_inter_ins_for_merge.groupTuple(by:[0,1]).map{it}.set{asso_intra_inter_ins_list_for_merge}

process Merge_ins {
  cpus 1

  input:
    set sid, val(side), file(tractogram) from asso_intra_inter_ins_list_for_merge

  output:
    set sid, "${sid}__asso_all_Ins_f_${side}.trk" into Ins_for_trk_plausible

  script:
  """
  scil_tractogram_math.py union ${tractogram} ${sid}__asso_all_Ins_f_${side}.trk -f
  """
}

/*
  ASSO CING
*/

process Asso_Cing {
  cpus 1

  input:
    set sid, val(side), file(tractogram) from asso_all_intra_inter_for_cing_filtering

  output:
    set sid, "${sid}__asso_all_Cing_${side}.trk" into Cing_for_trk_plausible
    set sid, val(side), "${sid}__asso_all_Cing_${side}.trk" into Cing_for_rename
    set sid, "${sid}__asso_lost_Cing_${side}.trk" optional true
    file "${sid}__asso_all_Cing_${side}.txt" optional true
    file "${sid}__asso_lost_Cing_${side}.txt" optional true

  script:
    filtering_list=params.FLF+"ASSO_Cing_${side}.txt"
    out_extension="asso_all_Cing_${side}"
    remaining_extension="asso_lost_Cing_${side}"
    basename="${sid}"
    keep="$params.keep_intermediate_steps"
    extract_masks=""

    template "filter_with_list.sh"
}

/*
 BE ASSO FRONTAL: extracting all streamlines with both ends in a frontal gyrus (U-shape > 20 mm)
*/

asso_frontal_be_list=params.asso_frontal_be_lists?.tokenize(',')
process Asso_be_frontal_gyrus {
  cpus 1

  input:
    set sid, val(side), file(tractogram) from asso_all_intra_inter_for_be_frontal_filtering
    each gyrus from asso_frontal_be_list

  output:
    set sid, val(side), val(gyrus), "${sid}_asso_intra_be_frontal_${gyrus}_${side}_u.trk" into asso_frontal_be_for_merge

  script:
  """
  scil_filter_tractogram.py ${tractogram} tmp.trk\
    --filtering_list ${params.FLF}ASSO_be_${gyrus}_${side}.txt -f
  scil_extract_ushape.py tmp.trk --minU 0.5 --maxU 1\
    ${sid}_asso_intra_be_frontal_${gyrus}_${side}_u.trk -f
  """
}

asso_frontal_be_for_merge.groupTuple(by:[0,1]).map{it}.set{asso_frontal_be_list_for_merge}
process Merge_asso_be_frontal_gyrus{
  cpus 1

  input:
    set sid, val(side),  val(gyrus), file(tractogram) from asso_frontal_be_list_for_merge

  output:
    set sid, "${sid}_asso_all_intraF_be_f_${side}_u.trk" into asso_all_intraF_be_for_trk_plausible

  script:
  """
    scil_tractogram_math.py union ${tractogram}\
      ${sid}_asso_all_intraF_be_f_${side}_u.trk -f
  """
}

/*
 EE ASSO FRONTAL: extracting all streamlines with either ends in a frontal gyrus (U-shape > 20 mm)
*/

asso_frontal_ee_list = Channel.from(['SFG_MFG', 70],
                                    ['SFG_IFG', 70],
                                    ['SFG_PrCG', 90],
                                    ['SFG_FrOrbG', 70],
                                    ['MFG_IFG', 70],
                                    ['MFG_PrCG', 110],
                                    ['MFG_FrOrbG', 60],
                                    ['IFG_PrCG', 110],
                                    ['IFG_FrOrbG', 60])
asso_all_intra_inter_for_ee_frontal_filtering.combine(asso_frontal_ee_list).set{asso_frontal_ee_for_extract}
process Asso_ee_frontal_gyrus {
  cpus 1

  input:
    set sid, val(side), file(tractogram), val(gyrus), val(max_length) from asso_frontal_ee_for_extract

  output:
    set sid, val(side), val(gyrus), "${sid}_asso_intra_ee_frontal_${gyrus}_${side}.trk" into asso_frontal_ee_for_merge

  script:
  """
  scil_filter_tractogram.py ${tractogram} tmp_01.trk\
    --filtering_list ${params.FLF}ASSO_ee_${gyrus}_${side}.txt -f
  scil_filter_streamlines_by_length.py tmp_01.trk tmp_02.trk\
    --maxL ${max_length} -f
  scil_extract_ushape.py tmp_02.trk\
    --minU 0.5\
    --maxU 1\
    ${sid}_asso_intra_ee_frontal_${gyrus}_${side}.trk -f
  """
}

asso_frontal_ee_for_merge.groupTuple(by:[0,1]).map{it}.set{asso_frontal_ee_list_for_merge}
process Merge_asso_ee_frontal_gyrus{
  cpus 1

  input:
    set sid, val(side),  val(gyrus), file(tractogram) from asso_frontal_ee_list_for_merge

  output:
    set sid, "${sid}_asso_all_intraF_ee_f_${side}_u.trk" into asso_all_intraF_ee_for_trk_plausible

  script:
  """
    scil_tractogram_math.py union ${tractogram} ${sid}_asso_all_intraF_ee_f_${side}_u.trk -f
  """
}

/*
 BE ASSO OCCIPITAL: extracting all streamlines with both ends in a occipital gyrus (U-shape > 20 mm)
*/

asso_occipital_be_list=params.asso_occipital_be_lists?.tokenize(',')
process Asso_be_occipital_gyrus {
  cpus 1

  input:
    set sid, val(side), file(tractogram) from asso_all_intra_inter_for_be_occipital_filtering
    each gyrus from asso_occipital_be_list

  output:
    set sid, val(side), val(gyrus), "${sid}_asso_intra_be_occipital_${gyrus}_${side}_u.trk" into asso_occipital_be_for_merge

  script:
  """
  scil_filter_tractogram.py ${tractogram} tmp.trk \
    --filtering_list ${params.FLF}ASSO_be_${gyrus}_${side}.txt -f
  scil_extract_ushape.py tmp.trk\
    --minU 0.5\
    --maxU 1\
    ${sid}_asso_intra_be_occipital_${gyrus}_${side}_u.trk -f
  """
}

asso_occipital_be_for_merge.groupTuple(by:[0,1]).map{it}.set{asso_occipital_be_list_for_merge}
process Merge_asso_be_occipital_gyrus{
  cpus 1

  input:
    set sid, val(side),  val(gyrus), file(tractogram) from asso_occipital_be_list_for_merge

  output:
    set sid, "${sid}_asso_all_intraO_be_f_${side}_u.trk" into asso_all_intraO_be_for_trk_plausible

  script:
  """
    scil_tractogram_math.py union ${tractogram} ${sid}_asso_all_intraO_be_f_${side}_u.trk -f
  """
}

/*
 EE ASSO OCCIPITAL: extracting all streamlines with either ends in a occipital gyrus (U-shape > 20 mm)
*/

asso_occipital_ee_list = Channel.from(['MOG_SOG', 60],['MOG_IOG', 50], ['MOG_CuG', 60], ['SOG_CuG', 30], ['CuG_LG', 60])
asso_all_intra_inter_for_ee_occipital_filtering.combine(asso_occipital_ee_list).set{asso_occipital_ee_for_extract}
process Asso_ee_occipital_gyrus {
  cpus 1

  input:
    set sid, val(side), file(tractogram), val(gyrus), val(max_length) from asso_occipital_ee_for_extract

  output:
    set sid, val(side), val(gyrus), "${sid}_asso_intra_ee_occipital_${gyrus}_${side}.trk" into asso_occipital_ee_for_merge

  script:
  """
  scil_filter_tractogram.py ${tractogram} tmp_01.trk\
    --filtering_list ${params.FLF}ASSO_ee_${gyrus}_${side}.txt -f
  scil_filter_streamlines_by_length.py tmp_01.trk tmp_02.trk\
    --maxL ${max_length} -f
  scil_extract_ushape.py tmp_02.trk\
    --minU 0.5\
    --maxU 1\
    ${sid}_asso_intra_ee_occipital_${gyrus}_${side}.trk -f
  """
}

asso_occipital_ee_for_merge.groupTuple(by:[0,1]).map{it}.set{asso_occipital_ee_list_for_merge}
process Merge_asso_ee_occipital_gyrus{
  cpus 1

  input:
    set sid, val(side),  val(gyrus), file(tractogram) from asso_occipital_ee_list_for_merge

  output:
    set sid, "${sid}_asso_all_intraO_ee_f_${side}_u.trk" into asso_all_intraO_ee_for_trk_plausible

  script:
  """
    scil_tractogram_math.py union ${tractogram} ${sid}_asso_all_intraO_ee_f_${side}_u.trk -f
  """
}

/*
 BE ASSO PARIETAL: extracting all streamlines with both ends in a parietal gyrus (U-shape > 20 mm)
*/

asso_parietal_be_list=params.asso_parietal_be_lists?.tokenize(',')
process Asso_be_parietal_gyrus {
  cpus 1

  input:
    set sid, val(side), file(tractogram) from asso_all_intra_inter_for_be_parietal_filtering
    each gyrus from asso_parietal_be_list

  output:
    set sid, val(side), val(gyrus), "${sid}_asso_intra_be_parietal_${gyrus}_${side}_u.trk" into asso_parietal_be_for_merge

  script:
  """
  scil_filter_tractogram.py ${tractogram} tmp.trk\
    --filtering_list ${params.FLF}ASSO_be_${gyrus}_${side}.txt -f
  scil_extract_ushape.py tmp.trk\
    --minU 0.5\
    --maxU 1\
    ${sid}_asso_intra_be_parietal_${gyrus}_${side}_u.trk -f
  """
}

asso_parietal_be_for_merge.groupTuple(by:[0,1]).map{it}.set{asso_parietal_be_list_for_merge}
process Merge_asso_be_parietal_gyrus{
  cpus 1

  input:
    set sid, val(side),  val(gyrus), file(tractogram) from asso_parietal_be_list_for_merge

  output:
    set sid, "${sid}_asso_all_intraP_be_f_${side}_u.trk" into asso_all_intraP_be_for_trk_plausible

  script:
  """
    scil_tractogram_math.py union ${tractogram} ${sid}_asso_all_intraP_be_f_${side}_u.trk -f
  """
}

/*
 EE ASSO PARIETAL: extracting all streamlines with either ends in a parietal gyrus (U-shape > 20 mm)
*/

asso_parietal_ee_list = Channel.from(['SPG_PoCG', 50], ['SPG_AG', 80], ['SPG_SMG', 70], ['SPG_PrCuG', 50], ['AG_PoCG', 10000], ['AG_SMG', 90], ['AG_PrCuG', 90] , ['SMG_PoCG', 60], ['SMG_PrCuG',100], ['PoCG_PrCuG', 80])
asso_all_intra_inter_for_ee_parietal_filtering.combine(asso_parietal_ee_list).set{asso_parietal_ee_for_extract}
process Asso_ee_parietal_gyrus {
  cpus 1

  input:
    set sid, val(side), file(tractogram), val(gyrus), val(max_length) from asso_parietal_ee_for_extract

  output:
    set sid, val(side), val(gyrus), "${sid}_asso_intra_ee_parietal_${gyrus}_${side}.trk" into asso_parietal_ee_for_merge

  script:
  """
  scil_filter_tractogram.py ${tractogram} tmp_01.trk\
    --filtering_list ${params.FLF}ASSO_ee_${gyrus}_${side}.txt -f
  scil_filter_streamlines_by_length.py tmp_01.trk tmp_02.trk\
    --maxL ${max_length} -f
  scil_extract_ushape.py tmp_02.trk\
    --minU 0.5\
    --maxU 1\
    ${sid}_asso_intra_ee_parietal_${gyrus}_${side}.trk -f
  """
}

asso_parietal_ee_for_merge.groupTuple(by:[0,1]).map{it}.set{asso_parietal_ee_list_for_merge}
process Merge_asso_ee_parietal_gyrus{
  cpus 1

  input:
    set sid, val(side),  val(gyrus), file(tractogram) from asso_parietal_ee_list_for_merge

  output:
    set sid, "${sid}_asso_all_intraP_ee_f_${side}.trk" into asso_all_intraP_ee_for_trk_plausible

  script:
  """
    scil_tractogram_math.py union ${tractogram} ${sid}_asso_all_intraP_ee_f_${side}.trk -f
  """
}

/*
 BE ASSO TEMPORAL: extracting all streamlines with both ends in a temporal gyrus and merge (U-shape > 20 mm)
*/

asso_temporal_be_list=params.asso_temporal_be_lists?.tokenize(',')
process Asso_be_temporal_gyrus {
  cpus 1

  input:
    set sid, val(side), file(tractogram) from asso_all_intra_inter_for_be_temporal_filtering
    each gyrus from asso_temporal_be_list

  output:
    set sid, val(side), val(gyrus), "${sid}_asso_intra_be_temporal_${gyrus}_${side}_u.trk" into asso_temporal_be_for_merge

  script:
  """
  scil_filter_tractogram.py ${tractogram} tmp.trk\
    --filtering_list ${params.FLF}ASSO_be_${gyrus}_${side}.txt -f
  scil_extract_ushape.py tmp.trk\
    --minU 0.5\
    --maxU 1\
    ${sid}_asso_intra_be_temporal_${gyrus}_${side}_u.trk -f
  """
}

asso_temporal_be_for_merge.groupTuple(by:[0,1]).map{it}.set{asso_temporal_be_list_for_merge}
process Merge_asso_be_temporal_gyrus{
  cpus 1

  input:
    set sid, val(side),  val(gyrus), file(tractogram) from asso_temporal_be_list_for_merge

  output:
    set sid, "${sid}_asso_all_intraT_be_f_${side}_u.trk" into asso_all_intraT_be_for_trk_plausible

  script:
  """
    scil_tractogram_math.py union ${tractogram} ${sid}_asso_all_intraT_be_f_${side}_u.trk -f
  """
}

/*
 EE ASSO TEMPORAL: extracting all streamlines with either ends in a temporal gyrus and merge (U-shape > 20 mm)
*/

asso_temporal_ee_list = Channel.from(['STG_MTG', 60], ['STG_ITG',80], ['STG_Tpole',110], ['MTG_ITG',60], ['MTG_Tpole', 100000], ['ITG_Tpole', 60])
asso_all_intra_inter_for_ee_temporal_filtering.combine(asso_temporal_ee_list).set{asso_temporal_ee_for_extract}
process Asso_ee_temporal_gyrus {
  cpus 1

  input:
    set sid, val(side), file(tractogram), val(gyrus), val(max_length) from asso_temporal_ee_for_extract

  output:
    set sid, val(side), val(gyrus), "${sid}_asso_intra_ee_temporal_${gyrus}_${side}.trk" into asso_temporal_ee_for_merge

  script:
  """
  scil_filter_tractogram.py ${tractogram} tmp_01.trk\
    --filtering_list ${params.FLF}ASSO_ee_${gyrus}_${side}.txt -f
  scil_filter_streamlines_by_length.py tmp_01.trk tmp_02.trk\
    --maxL ${max_length} -f
  scil_extract_ushape.py tmp_02.trk\
    --minU 0.5\
    --maxU 1\
    ${sid}_asso_intra_ee_temporal_${gyrus}_${side}.trk -f
  """
}

asso_temporal_ee_for_merge.groupTuple(by:[0,1]).map{it}.set{asso_temporal_ee_list_for_merge}
process Merge_asso_ee_temporal_gyrus{
  cpus 1

  input:
    set sid, val(side),  val(gyrus), file(tractogram) from asso_temporal_ee_list_for_merge

  output:
    set sid, "${sid}_asso_all_intraT_ee_f_${side}.trk" into asso_all_intraT_ee_for_trk_plausible

  script:
  """
    scil_tractogram_math.py union ${tractogram} ${sid}_asso_all_intraT_ee_f_${side}.trk -f
  """
}

fornix_for_trk_plausible.concat(cerebellum_for_trk_plausible,brainstem_for_trk_plausible,BG_ipsi_Thal_for_trk_plausible,BG_ipsi_Put_for_trk_plausible,BG_ipsi_Caud_for_trk_plausible,asso_u_shape_for_trk_plausible,CC_homo_for_trk_plausible,asso_all_dorsal_for_trk_plausible,asso_all_ventral_for_trk_plausible,all_P_O_for_trk_plausible,all_P_T_for_trk_plausible,all_O_T_for_trk_plausible,Ins_for_trk_plausible,Cing_for_trk_plausible,asso_all_intraF_be_for_trk_plausible,asso_all_intraF_ee_for_trk_plausible,asso_all_intraP_be_for_trk_plausible,asso_all_intraP_ee_for_trk_plausible,asso_all_intraO_be_for_trk_plausible,asso_all_intraO_ee_for_trk_plausible,asso_all_intraT_be_for_trk_plausible,asso_all_intraT_ee_for_trk_plausible, accx_for_trk_plausible, ccbg_for_trk_plausible).groupTuple(by: 0).set{merge_trk_plausible}

process Merge_trk_plausible{
  publishDir = params.final_output_mni_space
  cpus 1

  input:
    set sid, file(tractogram) from merge_trk_plausible

  output:
    set sid, "${sid}__plausible_${params.template_space}.trk" into plausible_for_extract_unplausible, trk_plausible_for_register_plausible_to_orig

  script:
  """
    scil_tractogram_math.py union ${tractogram} ${sid}__plausible_${params.template_space}.trk -f
  """
}

trk_for_extract_unplausible.join(plausible_for_extract_unplausible).set{for_trk_unplausible}

process Extract_trk_unplausible{
  publishDir = params.final_output_mni_space
  cpus 1

  input:
    set sid, file(trk01), file(trk02) from for_trk_unplausible
  output:
    set sid, "${sid}__unplausible_${params.template_space}.trk" into trk_unplausible_for_register_to_orig

  script:
  """
    scil_tractogram_math.py difference ${trk01} ${trk02} ${sid}__unplausible_${params.template_space}.trk -f
  """
}

/*
RENAME CC CC_Homotopic
*/
process Rename_cc_homotopic {
  publishDir = params.final_output_bundles_mni_space
  cpus 1

  input:
    set sid, val(list), file(trk01) from CC_Homotopic_frontal_for_rename
    set sid, val(list), file(trk02) from CC_Homotopic_occipital_for_rename
    set sid, val(list), file(trk03) from CC_Homotopic_temporal_for_rename
    set sid, val(list), file(trk04) from CC_Homotopic_parietal_for_rename
    set sid, val(list), file(trk05) from CC_Homotopic_insular_for_rename
    set sid, val(list), file(trk06) from CC_Homotopic_cingulum_for_rename
  output:
    set sid, "${sid}__cc_homotopic_frontal_${params.template_space}.trk" into cc_homotopic_frontal_for_register_to_orig
    set sid, "${sid}__cc_homotopic_occipital_${params.template_space}.trk" into cc_homotopic_occipital_for_register_to_orig
    set sid, "${sid}__cc_homotopic_temporal_${params.template_space}.trk" into cc_homotopic_temporal_for_register_to_orig
    set sid, "${sid}__cc_homotopic_parietal_${params.template_space}.trk" into cc_homotopic_parietal_for_register_to_orig
    set sid, "${sid}__cc_homotopic_insular_${params.template_space}.trk" into cc_homotopic_insular_for_register_to_orig
    set sid, "${sid}__cc_homotopic_cingulum_${params.template_space}.trk" into cc_homotopic_cingulum_for_register_to_orig

  when:
    params.extended

  script:
  """
  scil_tractogram_math.py union ${trk01} "${sid}__cc_homotopic_frontal_${params.template_space}.trk" -f
  scil_tractogram_math.py union ${trk02} "${sid}__cc_homotopic_occipital_${params.template_space}.trk" -f
  scil_tractogram_math.py union ${trk03} "${sid}__cc_homotopic_temporal_${params.template_space}.trk" -f
  scil_tractogram_math.py union ${trk04} "${sid}__cc_homotopic_parietal_${params.template_space}.trk" -f
  cp ${trk05} ${sid}__cc_homotopic_insular_${params.template_space}.trk -f
  cp ${trk06} ${sid}__cc_homotopic_cingulum_${params.template_space}.trk -f
  """
}

/*
RENAME CORTICO_STRIATE
*/
BG_ipsi_Caud_for_rename.concat(BG_ipsi_Put_for_rename).groupTuple(by:[0,1]).set{corticostriate_for_rename}
process Rename_cortico_striate {
  publishDir = params.final_output_bundles_mni_space
  cpus 1

  input:
    set sid, val(side), file(tractogram) from corticostriate_for_rename

  output:
    set sid, "${sid}__corticostriatal_${side}_${params.template_space}.trk" into corticostriatal_for_register_to_orig

  when:
    params.extended

  script:
  """
    scil_tractogram_math.py union ${tractogram} ${sid}__corticostriatal_${side}_${params.template_space}.trk -f
  """
}

/*
RENAME Corona radiata
*/
process Rename_coronaradiata {
  publishDir = params.final_output_bundles_mni_space
  cpus 1

  input:
    set sid, val(side), file(tractogram) from BG_ipsi_Thal_for_rename

  output:
    set sid, "${sid}__coronaradiata_${side}_${params.template_space}.trk" into coronaradiata_for_register_to_orig

  when:
    params.extended

  script:
  """
    scil_tractogram_math.py union ${tractogram} ${sid}__coronaradiata_${side}_${params.template_space}.trk -f
  """
}

/*
RENAME OPTICAL RADIATION
*/
process Rename_optical_radiation {
  publishDir = params.final_output_bundles_mni_space
  cpus 1

  input:
    set sid, val(side), val(list), file(tractogram) from optic_radiation_for_rename

  output:
    set sid, "${sid}__optical_radiation_${side}_${params.template_space}.trk" into optical_radiation_for_register_to_orig

  when:
    params.extended

  script:
  """
    scil_tractogram_math.py union ${tractogram} ${sid}__optical_radiation_${side}_${params.template_space}.trk -f
  """
}

/*
RENAME U SHAPE
*/
process Rename_ushape {
  publishDir = params.final_output_bundles_mni_space
  cpus 1

  input:
    set sid, val(side), file(tractogram) from asso_u_shape_for_rename

  output:
    set sid, "${sid}__ushape_${side}_${params.template_space}.trk" into ushape_for_register_to_orig

  when:
    params.extended

  script:
  """
    cp ${tractogram} ${sid}__ushape_${side}_${params.template_space}.trk
  """
}

/*
RENAME CING
*/
process Rename_cing {
  publishDir = params.final_output_bundles_mni_space
  cpus 1

  input:
    set sid, val(side), file(tractogram) from Cing_for_rename

  output:
    set sid, "${sid}__cing_${side}_${params.template_space}.trk" into cing_for_register_to_orig

  when:
    params.extended

  script:
  """
    cp ${tractogram} ${sid}__cing_${side}_${params.template_space}.trk
  """
}

/*
RENAME SLF
*/
asso_all_intra_inter_dorsal_all_f_O_for_rename.concat(asso_all_intra_inter_dorsal_f_p_for_rename).groupTuple(by:[0,1]).set{slf_for_rename}
process Rename_slf {
  publishDir = params.final_output_bundles_mni_space
  cpus 1

  input:
    set sid, val(side), val(list), file(tractogram) from slf_for_rename

  output:
    set sid, "${sid}__slf_${side}_${params.template_space}.trk" into slf_for_register_to_orig

  when:
    params.extended

  script:
  """
    scil_tractogram_math.py union ${tractogram} ${sid}__slf_${side}_${params.template_space}.trk -f
  """
}

/*
RENAME AF
*/
process Rename_af {
  publishDir = params.final_output_bundles_mni_space
  cpus 1

  input:
    set sid, val(side), asso_list, file(tractogram) from asso_all_intra_inter_dorsal_all_f_T_for_rename

  output:
    set sid, "${sid}__af_${side}_${params.template_space}.trk" into af_for_register_to_orig

  when:
    params.extended

  script:
  """
    cp ${tractogram} ${sid}__af_${side}_${params.template_space}.trk -f
  """
}

/*
RENAME Cortico-pontine_F
*/
process Rename_corticopontine_F {
  publishDir = params.final_output_bundles_mni_space
  cpus 1

  input:
    set sid, file(tractogram) from brainstem_corticopontine_frontal_for_rename
    each side from side_corticopontineF
  output:
    set sid, "${sid}__corticopontine_frontal_${side}_${params.template_space}.trk" into corticopontine_frontal_for_register_to_orig

  when:
    params.extended

  script:
    filtering_list=params.FLF+"frontal_${side}.txt"
    out_extension="corticopontine_frontal_${side}_${params.template_space}"
    remaining_extension="lost"
    basename="${sid}"
    keep=false
    extract_masks=""

    template "filter_with_list.sh"
}

/*
RENAME cortico-pontine_POT
*/
process Rename_corticopontine_POT {
  publishDir = params.final_output_bundles_mni_space
  cpus 1

  input:
    set sid, file(tractogram) from brainstem_ee_corticopontine_parietotemporooccipital_for_rename
    each side from side_corticopontinePOT

  output:
    set sid, "${sid}__corticopontine_POT_${side}_${params.template_space}.trk" into corticopontine_POT_for_register_to_orig

  when:
    params.extended

  script:
    filtering_list=params.FLF+"parieto_temporo_occipital_${side}.txt"
    out_extension="corticopontine_POT_${side}_${params.template_space}"
    remaining_extension="lost"
    basename="${sid}"
    keep=false
    extract_masks=""

    template "filter_with_list.sh"
}

/*
RENAME Pyramidal tract (CST)
*/
process Rename_cst {
  publishDir = params.final_output_bundles_mni_space
  cpus 1

  input:
    set sid, file(tractogram) from brainstem_pyramidal_for_rename
    each side from side_cst

  output:
    set sid, "${sid}__cst_${side}_${params.template_space}.trk" into cst_for_register_to_orig

  when:
    params.extended

  script:
    filtering_list=params.FLF+"fronto_parietal_${side}.txt"
    out_extension="cst_${side}_${params.template_space}"
    remaining_extension="lost"
    basename="${sid}"
    keep=false
    extract_masks=""

    template "filter_with_list.sh"
}

/*
RENAME fornix
*/
process Rename_fornix {
  publishDir = params.final_output_bundles_mni_space
  cpus 1

  input:
    set sid, file(tractogram) from fornix_for_rename

  output:
    set sid, "${sid}__fornix_${params.template_space}.trk" into fornix_for_register_to_orig

  when:
    params.extended

  script:
  """
    cp ${tractogram} ${sid}__fornix_${params.template_space}.trk -f
  """
}

/*
RENAME IFOF
*/
process Rename_ifof {
  publishDir = params.final_output_bundles_mni_space
  cpus 1

  input:
    set sid, val(side), file(tractogram) from asso_IFOF_for_rename

  output:
    set sid, "${sid}__ifof_${side}_${params.template_space}.trk" into ifof_for_register_to_orig

  when:
    params.extended

  script:
  """
    cp ${tractogram} ${sid}__ifof_${side}_${params.template_space}.trk -f
  """
}

/*
RENAME UF
*/
process Rename_uf {
  publishDir = params.final_output_bundles_mni_space
  cpus 1

  input:
    set sid, val(side), file(tractogram) from asso_UF_for_rename

  output:
    set sid, "${sid}__uf_${side}_${params.template_space}.trk" into uf_for_register_to_orig

  when:
    params.extended

  script:
  """
    cp ${tractogram} ${sid}__uf_${side}_${params.template_space}.trk -f
  """
}

/*
RENAME ILF
*/
process Rename_ilf {
  publishDir = params.final_output_bundles_mni_space
  cpus 1

  input:
    set sid, val(side), file(tractogram) from all_O_T_for_rename

  output:
    set sid, "${sid}__ilf_${side}_${params.template_space}.trk" into ilf_for_register_to_orig

  when:
    params.extended

  script:
  """
    cp ${tractogram} ${sid}__ilf_${side}_${params.template_space}.trk -f
  """
}

/*
RENAME BRAINSTEM
*/
process Rename_brainstem {
  publishDir = params.final_output_bundles_mni_space
  cpus 1

  input:
    set sid, file(tractogram) from brainstem_for_rename

  output:
    set sid, "${sid}__brainstem_${params.template_space}.trk" into brainstem_for_register_to_orig

  when:
    params.extended

  script:
  """
    cp ${tractogram} ${sid}__brainstem_${params.template_space}.trk -f
  """
}

/*
RENAME CEREBELLUM
*/
process Rename_cerebellum {
  publishDir = params.final_output_bundles_mni_space
  cpus 1

  input:
    set sid, file(tractogram) from cerebellum_for_rename

  output:
    set sid, "${sid}__cerebellum_${params.template_space}.trk" into cerebellum_for_register_to_orig

  when:
    params.extended

  script:
  """
    cp ${tractogram} ${sid}__cerebellum_${params.template_space}.trk -f
  """
}

/*
RENAME AC_CX
*/
process Rename_accx {
  publishDir = params.final_output_bundles_mni_space
  cpus 1

  input:
    set sid, file(tractogram) from accx_for_rename

  output:
    set sid, "${sid}__accx_${params.template_space}.trk" into accx_for_register_to_orig

  when:
    params.extended

  script:
  """
    cp ${tractogram} ${sid}__accx_${params.template_space}.trk -f
  """
}


trks_for_register = Channel.empty()
bundles_for_register = Channel.empty()

if (params.orig){
    if (params.extended){
        cc_homotopic_frontal_for_register_to_orig
            .concat(cc_homotopic_occipital_for_register_to_orig)
            .concat(cc_homotopic_temporal_for_register_to_orig)
            .concat(cc_homotopic_parietal_for_register_to_orig)
            .concat(cc_homotopic_insular_for_register_to_orig)
            .concat(cc_homotopic_cingulum_for_register_to_orig)
            .concat(corticostriatal_for_register_to_orig)
            .concat(coronaradiata_for_register_to_orig)
            .concat(optical_radiation_for_register_to_orig)
            .concat(ushape_for_register_to_orig)
            .concat(cing_for_register_to_orig)
            .concat(slf_for_register_to_orig)
            .concat(af_for_register_to_orig)
            .concat(corticopontine_frontal_for_register_to_orig)
            .concat(corticopontine_POT_for_register_to_orig)
            .concat(cst_for_register_to_orig)
            .concat(fornix_for_register_to_orig)
            .concat(ifof_for_register_to_orig)
            .concat(uf_for_register_to_orig)
            .concat(ilf_for_register_to_orig)
            .concat(brainstem_for_register_to_orig)
            .concat(cerebellum_for_register_to_orig)
            .concat(accx_for_register_to_orig)
            .concat(plausible_commissural_for_register_to_orig)
            .combine(transformation_and_t1_for_transformation_to_orig_bundles, by: 0)
            .set{bundles_for_register}

        trk_plausible_for_register_plausible_to_orig
            .concat(trk_unplausible_for_register_to_orig)
            .combine(transformation_and_t1_for_transformation_to_orig, by: 0)
            .set{trks_for_register}
      }
    else{
        trk_plausible_for_register_plausible_to_orig
            .concat(trk_unplausible_for_register_to_orig)
            .combine(transformation_and_t1_for_transformation_to_orig, by: 0)
            .set{trks_for_register}
    }
}
else{
  trks_for_register = Channel.create()
  trks_for_register.close()
}

process Register_to_orig{
  publishDir = params.final_output_orig_space
  cpus 1

  input:
    set sid, file(trk), file(t1), file(transfo), file(inv_deformation), file(deformation) from trks_for_register

  output:
    set sid, "${sid}__*_${params.orig_space}.trk"

  when:
    params.orig

  script:
  """
    scil_apply_transform_to_tractogram.py ${trk} ${t1} ${transfo}  ${trk.getSimpleName().replaceAll("mni_space", "orig_space")}.trk --in_deformation ${deformation} --reverse_operation --keep_invalid
  """
}

process Register_bundles_to_orig{
  publishDir = params.final_output_bundles_orig_space
  cpus 1

  input:
    set sid, file(trk), file(t1), file(transfo), file(inv_deformation), file(deformation) from bundles_for_register

  output:
    set sid, "${sid}__*_${params.orig_space}.trk"

  when:
    params.orig

  script:
  """
    scil_apply_transform_to_tractogram.py ${trk} ${t1} ${transfo}  ${trk.getSimpleName().replaceAll("mni_space", "orig_space")}.trk --in_deformation ${deformation} --reverse_operation --keep_invalid
  """
}

process Copy_t1_to_orig{
  publishDir = params.final_output_orig_space
  cpus 1

  input:
    tuple sid, file(t1) from t1s_for_copy_to_orig

  output:
    file("${sid}__t1_orig_space.nii.gz")

  when:
    params.orig

  script:
  """
    cp ${t1} ${sid}__t1_orig_space.nii.gz
  """
}
