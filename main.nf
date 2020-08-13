#!/usr/bin/env nextflow

params.root = false
params.help = false
params.keep = true

if(params.help) {
    usage = file("$baseDir/USAGE")
    engine = new groovy.text.SimpleTemplateEngine()
    template = engine.createTemplate(usage.text).make()

    print template.toString()
    return
    }

if (params.root){
    log.info "Input: $params.root"
    root = file(params.root)
    in_data = Channel
        .fromFilePairs("$root/**/*{.trk}",
                       size:1,
                       maxDepth:1,
                       flat: true) {it.parent.name}
}
else {
    error "Error ~ Please use --root for the input data."
}

(inputTractogram) = in_data
    .map{sid, tractogram -> [tuple(sid, tractogram)]}
    .separate(1)

sides = params.sides?.tokenize(',')

process rm_out_JHU {
    cpus 1
    tag = "Remove out of JHU"

    input:
    set sid, file(tractogram) from inputTractogram

    output:
    set sid, "${sid}_wb_in_JHU.trk" into wb_for_rm_crossing_gyri
    file "${sid}_wb_in_JHU.txt" optional true
    set sid, "${sid}_wb_out_JHU.trk" optional true into wb_out_JHU
    file "${sid}_wb_out_JHU.txt" optional true

    script:
    atlas=params.roisFolder+params.atlas.JHU
    mode="all"
    criteria="include"
    out_extension='wb_in_JHU'
    remaining_extension='wb_out_JHU'
    basename="${sid}"

    template "filter_with_atlas.sh"
}

process remove_crossing_gyri {
  cpus 1
  tag = "Remove crossing the gyri limits of the JHU template"

  input:
    set sid, file(tractogram) from wb_for_rm_crossing_gyri

  output:
   set sid, "${sid}_wb_rm_crossing_gyri.trk" into wb_for_pruning
   file "${sid}_wb_rm_crossing_gyri.txt" optional true
   set sid, "${sid}_wb_crossing_gyri.trk" optional true into wb_crossing_gyri
   file "${sid}_wb_crossing_gyri.txt" optional true

   script:
   atlas=params.roisFolder+params.atlas.shell_limits
   mode="any"
   criteria="exclude"
   out_extension='wb_rm_crossing_gyri'
   remaining_extension='wb_crossing_gyri'
   basename="${sid}"

   template "filter_with_atlas.sh"
}

process pruning {
    cpus 1
    tag = "Pruning min 20mm"

    input:
    set sid, file(tractogram) from wb_for_pruning

    output:
    set sid, "${sid}_wb_min20.trk" into wb_for_rmloop
    file "${sid}_wb_min20.txt" optional true
    set sid, "${sid}_wb_max20.trk" optional true into wb_max20
    file "${sid}_wb_max20.txt" optional true

    script:

      """
      scil_filter_streamlines_by_length.py ${tractogram} --minL 20 \
                                           --maxL 100000 ${sid}_wb_min20.trk \
                                           -f --display_counts > ${sid}_wb_min20.txt
      if ${params.keep}
      then
        scil_streamlines_math.py difference ${tractogram} ${sid}_wb_min20.trk ${sid}_wb_max20.trk -f
      fi
      """
}

process removing_loops {
  cpus 1
  tag = "Remove loop"

  input:
  set sid, file(wb_min20) from wb_for_rmloop

  output:
  set sid, "${sid}_wb_min20_noloop.trk" into wb_for_rm_end_in_cc_dwm
  set sid, "${sid}_wb_loops.trk" optional true into wb_loops
  file "${sid}_wb_min20_noloop.txt" optional true
  file "${sid}_wb_loops.txt" optional true

  script:

    """
    scil_detect_streamlines_loops.py ${wb_min20} ${sid}_wb_min20_noloop.trk \
                                     -a ${params.angle}  \
                                     --looping_tractogram ${sid}_wb_loops.trk \
                                     -f
    if ${params.keep};
    then
      scil_count_streamlines.py ${sid}_wb_loops.trk > ${sid}_wb_loops.txt
      scil_count_streamlines.py ${sid}_wb_min20_noloop.trk > ${sid}_wb_min20_noloop.txt
    fi
    """
}

process removing_end_in_cc_DWM {
  cpus 1
  tag = "Remove end in CC DWM"

  input:
  set sid, file(wb_min20_noloop) from wb_for_rm_end_in_cc_dwm

  output:
  set sid, "${sid}_wb_clean01.trk" into wb_for_cerebellum_split
  file "${sid}_wb_clean01.txt" optional true

  script:
  """
  scil_filter_tractogram.py ${wb_min20_noloop} ${sid}_wb_clean01.trk \
                    --drawn_roi ${params.roisFolder}${params.atlas.cc} either_end exclude \
                    --drawn_roi ${params.roisFolder}${params.atlas.dwm} either_end exclude \
                    -f --display_count > ${sid}_wb_clean01.txt
  """
}

process split_either_end_in_cerebellum {

  cpus 1
  tag = "Split either end in cerebellum"

  input:
    set sid, file(tractogram) from wb_for_cerebellum_split

  output:
  set sid, "${sid}_wb_clean01_nocereb.trk" into notEndInCerebellum
  set sid, "${sid}_all_cereb.trk" into endInCerebellum
  file "${sid}_all_cereb.txt" optional true
  file "${sid}_wb_clean01_nocereb.txt" optional true

  script:
  atlas=params.roisFolder+params.atlas.cerebellum
  mode=params.mode.either_end
  criteria=params.criteria.exclude
  out_extension='wb_clean01_nocereb'
  remaining_extension='all_cereb'
  basename="${sid}"

  template "filter_with_atlas.sh"
}

process split_either_end_in_brainstem {

  cpus 1
  tag = "Split either end in brainstem"

  input:
    set sid, file(tractogram) from notEndInCerebellum

  output:
    set sid, "${sid}_wb_clean02.trk" into notEndInBrainstem
    set sid, "${sid}_all_brainstem.trk" into endInBrainstem
    file "${sid}_wb_clean02.txt" optional true
    file "${sid}_all_brainstem.txt" optional true

  script:
    atlas=params.roisFolder+params.atlas.brainstem
    mode=params.mode.either_end
    criteria=params.criteria.exclude
    out_extension='wb_clean02'
    remaining_extension='all_brainstem'
    basename="${sid}"

    template "filter_with_atlas.sh"
}

//Valide
process removing_GM {

  cpus 1
  tag = "Brain Remove GM"

  input:
    set sid, file(tractogram) from notEndInBrainstem

  output:
  set sid, "${sid}_wb_either_CGM_SWM.trk" into endInCGMSWI
  set sid, "${sid}_no_CGM_SWM.trk" into endNotInCGMSWI
  file "${sid}_wb_either_CGM_SWM.txt" optional true
  file "${sid}_no_CGM_SWM.txt" optional true

  script:
    atlas=params.roisFolder+params.atlas.CGMSWM
    mode=params.mode.either_end
    criteria=params.criteria.include
    out_extension='wb_either_CGM_SWM'
    remaining_extension='no_CGM_SWM'
    basename="${sid}"

    template "filter_with_atlas.sh"
}

//Valide
process split_CC {
  cpus 1
  tag = "split CC"

  input:
  set sid, file(tractogram) from endInCGMSWI

  output:
  set sid, "${sid}_tmp_CC.trk" into inCC4BG, inCC4Other
  set sid, "${sid}_wb_either_CGM_SWM_noCC.trk" into notInCC
  file "${sid}_wb_either_CGM_SWM_noCC.txt" optional true
  file "${sid}_tmp_CC.txt" optional true

  script:
  atlas=params.roisFolder+params.atlas.midline
  mode=params.mode.any
  criteria=params.criteria.exclude
  out_extension="noCC"
  remaining_extension='tmp_CC'
  basename=tractogram.getSimpleName()

  template "filter_with_atlas.sh"
}

//Valide
process split_CC_BG {
  cpus 1
  tag = "Split Basal Ganglia"

  input:
  set sid, file(tractogram) from inCC4BG
  each side from sides

  output:
  set sid, "${sid}_contra_BG_${side}.trk" into inCCBG
  file "${sid}_contra_BG_${side}.txt" optional true

  script:
  atlas=params.roisFolder+params.atlas.subcortical+"_${side}.nii.gz"
  mode=params.mode.either_end
  criteria=params.criteria.include
  out_extension="contra_BG_" + "${side}"
  remaining_extension=params.ext.notUsed
  basename="${sid}"

  template "filter_with_atlas.sh"
}

//Valide sauf -> lost
process removing_unplausible_CC {
  cpus 1
  tag = "Removing unplausible CC"

  input:
  set sid, file(tractogram) from inCC4Other

  output:
  set sid, "${sid}_CC_Cx.trk" into ccCleaned
  set sid, "${sid}_CC_lost.trk" optional true into CC_lost
  file "${sid}_CC_Cx.txt" optional true
  file "${sid}_CC_lost.txt" optional true


  script:
  """
  scil_filter_tractogram.py ${tractogram} ${sid}_CC_Cx.trk \
    --drawn_roi ${params.roisFolder}${params.atlas.allsubcortical} either_end exclude \
    --drawn_roi ${params.roisFolder}${params.atlas.brainstemINF} any exclude \
    --drawn_roi ${params.roisFolder}${params.atlas.ic} any exclude \
    --drawn_roi ${params.roisFolder}${params.atlas.allThal} any exclude \
    --drawn_roi ${params.roisFolder}${params.atlas.midline} either_end exclude \
    --drawn_roi ${params.roisFolder}${params.atlas.allL} both_ends exclude \
    --drawn_roi ${params.roisFolder}${params.atlas.allR} both_ends exclude -f;
  if ${params.keep}
  then
    scil_streamlines_math.py difference ${tractogram} ${sid}_CC_Cx.trk ${sid}_CC_lost.trk
    scil_count_streamlines.py ${sid}_CC_lost.trk > ${sid}_CC_lost.txt
    scil_count_streamlines.py ${sid}_CC_Cx.trk > ${sid}_CC_Cx.txt
  fi
  """
}

//Valide
process split_noCC_asso_BG {
  cpus 1
  tag = "Split not CC in asso BG and not BG"

  input:
  set sid, file(tractogram) from notInCC

  output:
  set sid, "${sid}_all_subcortical_from_CGM_SWM_noCC_f.trk" into assoBG
  file "${sid}_all_subcortical_from_CGM_SWM_noCC_f.txt" optional true
  set sid, "${sid}_asso_noBG.trk" into assoNoBG
  file "${sid}_asso_noBG.txt" optional true

  script:
  atlas=params.roisFolder+params.atlas.allsubcortical
  mode=params.mode.either_end
  criteria=params.criteria.include
  out_extension="all_subcortical_from_CGM_SWM_noCC_f"
  remaining_extension='asso_noBG'
  basename="${sid}"

  template "filter_with_atlas.sh"
}

//Valide
process split_asso_in_hemi {
  cpus 1
  tag = "Split asso in hemispheres"

  input:
  set sid, file(tractogram) from assoNoBG
  each side from sides

  output:
  set sid, val(side), "${sid}_asso_${side}.trk" into asso
  set sid, "${sid}_asso_${side}_lost.trk" optional true into assoLost
  file "${sid}_asso_${side}.txt" optional true
  file "${sid}_asso_${side}_lost.txt" optional true

  script:
  if (side=='L')
    atlas=params.roisFolder+params.atlas.all+"_R.nii.gz"
  else
    atlas=params.roisFolder+params.atlas.all+"_L.nii.gz"

  mode=params.mode.any
  criteria=params.criteria.exclude
  out_extension="asso_${side}"
  remaining_extension="asso_${side}_lost"
  basename="${sid}"

  template "filter_with_atlas.sh"
}

assoLost.groupTuple().map{it.flatten().toList()}.set{assoLostGroup}

//Valide
process extract_lost {
  cpus 1
  tag = "Extract asso lost"

  input:
  set sid, file(trk1), file(trk2) from assoLostGroup

  output:
  file "${sid}_asso_lost.txt"
  set sid, "${sid}_asso_lost.trk" into allassoLost


  script:
  """
  scil_streamlines_math.py intersection ${trk1} \
                                        ${trk2} \
                                        ${sid}_asso_lost.trk
  scil_count_streamlines.py ${sid}_asso_lost.trk > ${sid}_asso_lost.txt
  """
}

//Valide
process split_ushape_cgm_asso {
  cpus 1
  tag = "Extracting U-shaped and streamlines restricted to Cortical GM and removing them from asso"

  input:
  set sid, val(side), file(tractogram) from asso

  output:
  set sid, "${sid}_asso_only_in_CGM_${side}.trk" into assoCGM
  set sid, "${sid}_asso_Ushape_${side}.trk" into assoUShape
  set sid, val(side), "${sid}_asso_f_${side}.trk" into assof
  file "${sid}_asso_only_in_CGM_${side}.txt" optional true
  file "${sid}_asso_Ushape_${side}.txt" optional true
  file "${sid}_asso_f_${side}.txt" optional true

  script:
  """
    scil_filter_tractogram.py ${tractogram} ${sid}_tmp1_${side}.trk \
      --drawn_roi ${params.roisFolder}${params.atlas.CGM}_${side}.nii.gz ${params.mode.both_ends} include -f

    scil_streamlines_math.py difference ${tractogram} ${sid}_tmp1_${side}.trk \
                             ${sid}_asso_SWM_${side}.trk -f

    scil_filter_tractogram.py ${sid}_tmp1_${side}.trk ${sid}_asso_only_in_CGM_${side}.trk \
      --drawn_roi ${params.roisFolder}${params.atlas.SWM}_${side}.nii.gz ${params.mode.any} exclude

    scil_streamlines_math.py difference ${sid}_tmp1_${side}.trk ${sid}_asso_only_in_CGM_${side}.trk \
                                 ${sid}_tmp2_${side}.trk -f

    scil_filter_tractogram.py ${sid}_tmp2_${side}.trk ${sid}_asso_Ushape_${side}.trk \
      --drawn_roi ${params.roisFolder}${params.atlas.DWM}_${side}.nii.gz ${params.mode.any} exclude

    scil_streamlines_math.py difference ${sid}_tmp2_${side}.trk ${sid}_asso_Ushape_${side}.trk \
                               ${sid}_asso_DWM_${side}.trk -f

    scil_streamlines_math.py concatenate ${sid}_asso_DWM_${side}.trk ${sid}_asso_SWM_${side}.trk ${sid}_asso_f_${side}.trk -f

    if ${params.keep}
    then
      scil_count_streamlines.py ${sid}_asso_only_in_CGM_${side}.trk > ${sid}_asso_only_in_CGM_${side}.txt
      scil_count_streamlines.py ${sid}_asso_Ushape_${side}.trk > ${sid}_asso_Ushape_${side}.txt
      scil_count_streamlines.py ${sid}_asso_f_${side}.trk > ${sid}_asso_f_${side}.txt
    fi
  """
}

//Valide
process removing_unplausible_long_range_asso {
  cpus 1
  tag = "Extracting unplausible long-range association streamlines passing through subcortical structures (Cd, Put, GP, Thal, Amyg)"

  input:
  set sid, val(side), file(tractogram) from assof

  output:
  set sid, val(side), "${sid}_asso_all_intra_inter_${side}.trk" into assoAllIntraInter
  set sid, "${sid}_asso_lost2_${side}.trk" into assoLost2
  file "${sid}_asso_all_intra_inter_${side}.txt" optional true
  file "${sid}_asso_lost2_${side}.txt" optional true

  script:
  atlas=params.roisFolder+params.atlas.allsubcortical
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

process cereb_removing_GM {

  cpus 1
  tag = "Cereb - Remove GM"

  input:
    set sid, file(tractogram) from inCerebellum

  output:
  set sid, "${sid}_all_cereb_nocx.trk" into ccEndNotInCGMSWI
  set sid, "${sid}_cereb_bin_01.trk" optional true into ccEndInCGMSWI
  file "${sid}_all_cereb_nocx.txt" optional true
  file "${sid}_cereb_bin_01.txt" optional true

  script:
    atlas=params.roisFolder+params.atlas.CGMSWM
    mode=params.mode.either_end
    criteria=params.criteria.exclude
    out_extension='all_cereb_nocx'
    remaining_extension='cereb_bin_01'
    basename="${sid}"

    template "filter_with_atlas.sh"
}

process cereb_removing_out_cerebellar_hemisphere {
  cpus 1
  tag = "Cereb - Remove GM"

  input:
    set sid, file(tractogram) from ccEndNotInCGMSWI

  output:
  set sid, "${sid}_all_cereb_nocx_f.trk" into noCXf
  set sid, "${sid}_cereb_bin_02.trk" optional true into bin02
  file "${sid}_all_cereb_nocx_f.txt" optional true
  file "${sid}_cereb_bin_02.txt" optional true

  script:
    atlas=params.roisFolder+params.atlas.cerebellumGWM
    mode=params.mode.either_end
    criteria=params.criteria.include
    out_extension='all_cereb_nocx_f'
    remaining_extension='cereb_bin_02'
    basename="${sid}"

    template "filter_with_atlas.sh"
}

process cereb_split {
  cpus 1
  tag = "Cereb - Split"

  input:
    set sid, file(tractogram) from noCXf

  output:
  set sid, "${sid}_all_cereb_nocx_in_cereb.trk" into noCX_inCereb
  set sid, "${sid}_all_cereb_nocx_out_cereb.trk" into noCX_outCereb
  file "${sid}_all_cereb_nocx_in_cereb.txt" optional true
  file "${sid}_all_cereb_nocx_out_cereb.txt" optional true

  script:
  atlas=params.roisFolder+params.atlas.cerebellum
  mode=params.mode.both_ends
  criteria=params.criteria.include
  out_extension='all_cereb_nocx_in_cereb'
  remaining_extension='all_cereb_nocx_out_cereb'
  basename="${sid}"

  template "filter_with_atlas.sh"
}

process cereb_no_loop{
  cpus 1
  tag = "Cereb - no loop"

  input:
    set sid, file(tractogram) from noCX_inCereb

  output:
  set sid, "${sid}_all_cereb_nocx_in_cereb_f.trk" into noC
  set sid, "${sid}_cereb_bin_03.trk" into bin03
  file "${sid}_all_cereb_nocx_in_cereb_f.txt" optional true

  script:
  """
    scil_filter_tractogram.py ${tractogram} ${sid}_all_cereb_nocx_in_cereb_f.trk \
      --drawn_roi ${params.roisFolder}${params.atlas.medulla} any exclude \
      --drawn_roi ${params.roisFolder}${params.atlas.midbrainNoSCP} any exclude \
      -f \
      --display_count > ${sid}_all_cereb_nocx_in_cereb_f.txt

    scil_streamlines_math.py difference ${tractogram} ${sid}_all_cereb_nocx_in_cereb_f.trk ${sid}_cereb_bin_03.trk
  """
}

process cereb_inMedulla {
  cpus 1
  tag = "Cereb - In Medulla"

  input:
    set sid, file(tractogram) from noCX_outCereb

  output:
    set sid, "${sid}_all_cereb_nocx_out_cereb_in_medulla.trk" into inMedulla
    set sid, "${sid}_all_cereb_nocx_out_cereb_no_medulla.trk" into noMedulla
    file "${sid}_all_cereb_nocx_out_cereb_in_medulla.txt" optional true
    file "${sid}_all_cereb_nocx_out_cereb_no_medulla.txt" optional true

  script:
  """
  scil_filter_tractogram.py ${tractogram} ${sid}_all_cereb_nocx_out_cereb_in_medulla.trk \
    --drawn_roi ${params.roisFolder}${params.atlas.medulla} either_end include \
    --drawn_roi ${params.roisFolder}${params.atlas.thalCPMidBrain} any exclude \
    --drawn_roi ${params.roisFolder}${params.atlas.MCPant} any exclude \
    --drawn_roi ${params.roisFolder}${params.atlas.CGMSWMDWM} any exclude -f

  scil_streamlines_math.py difference ${tractogram} ${sid}_all_cereb_nocx_out_cereb_in_medulla.trk ${sid}_all_cereb_nocx_out_cereb_no_medulla.trk

  if ${params.keep}
  then
    scil_count_streamlines.py ${sid}_all_cereb_nocx_out_cereb_in_medulla.trk > ${sid}_all_cereb_nocx_out_cereb_in_medulla.txt
    scil_count_streamlines.py ${sid}_all_cereb_nocx_out_cereb_no_medulla.trk > ${sid}_all_cereb_nocx_out_cereb_no_medulla.txt
  fi
  """
}

process cereb_inPons{
  cpus 1
  tag = "Cereb - In Pons"

  input:
    set sid, file(tractogram) from noMedulla

  output:
  set sid, "${sid}_all_cereb_nocx_out_cereb_in_pons.trk" into inPons
  set sid, "${sid}_all_cereb_nocx_out_cereb_no_pons.trk" into noPons
  file "${sid}_all_cereb_nocx_out_cereb_in_pons.txt" optional true
  file "${sid}_all_cereb_nocx_out_cereb_no_pons.txt" optional true


  script:
  """
  scil_filter_tractogram.py ${tractogram} ${sid}_all_cereb_nocx_out_cereb_in_pons.trk \
    --drawn_roi ${params.roisFolder}${params.atlas.pons} either_end include \
    --drawn_roi ${params.roisFolder}${params.atlas.medulla} any exclude \
    --drawn_roi ${params.roisFolder}${params.atlas.thalCPMidBrain} any exclude \
    --drawn_roi ${params.roisFolder}${params.atlas.CGMSWMDWM} any exclude \
    --drawn_roi ${params.roisFolder}${params.atlas.MCPant} any include -f

  scil_streamlines_math.py difference ${tractogram} ${sid}_all_cereb_nocx_out_cereb_in_pons.trk ${sid}_all_cereb_nocx_out_cereb_no_pons.trk -f

  if ${params.keep}
  then
    scil_count_streamlines.py ${sid}_all_cereb_nocx_out_cereb_in_pons.trk > ${sid}_all_cereb_nocx_out_cereb_in_pons.txt
    scil_count_streamlines.py ${sid}_all_cereb_nocx_out_cereb_no_pons.trk > ${sid}_all_cereb_nocx_out_cereb_no_pons.txt
  fi
  """
}

process cereb_inMidBrain{
  cpus 1
  tag = "Cereb - In the midbrain"

  input:
    set sid, file(tractogram) from noPons

  output:
  set sid, "${sid}_all_cereb_nocx_out_cereb_in_midbrain.trk" into inMidBrain
  set sid, "${sid}_all_cereb_nocx_out_cereb_no_midbrain.trk" into noMidBrain
  file "${sid}_all_cereb_nocx_out_cereb_in_midbrain.txt" optional true
  file "${sid}_all_cereb_nocx_out_cereb_no_midbrain.txt" optional true


  script:
  """
    scil_filter_tractogram.py ${tractogram} ${sid}_all_cereb_nocx_out_cereb_in_midbrain.trk \
      --drawn_roi ${params.roisFolder}${params.atlas.midbrain} either_end include \
      --drawn_roi ${params.roisFolder}${params.atlas.thalCPMedulla} any exclude \
      --drawn_roi ${params.roisFolder}${params.atlas.MCPant} any exclude \
      --drawn_roi ${params.roisFolder}${params.atlas.CGMSWMDWM} any exclude -f

    scil_streamlines_math.py difference ${tractogram} ${sid}_all_cereb_nocx_out_cereb_in_midbrain.trk ${sid}_all_cereb_nocx_out_cereb_no_midbrain.trk -f

    if ${params.keep}
    then
      scil_count_streamlines.py ${sid}_all_cereb_nocx_out_cereb_in_midbrain.trk > ${sid}_all_cereb_nocx_out_cereb_in_midbrain.txt
      scil_count_streamlines.py ${sid}_all_cereb_nocx_out_cereb_no_midbrain.trk > ${sid}_all_cereb_nocx_out_cereb_no_midbrain.txt
    fi
  """
}

process cereb_inRedNucleusThalamus{
  cpus 1
  tag = "Cereb - One termination in the red nucleus and/or the thalamus"

  input:
    set sid, file(tractogram) from noMidBrain

  output:
  set sid, "${sid}_all_cereb_nocx_out_cereb_in_redN_and_thalamus.trk" into inRedNThal
  set sid, "${sid}_cereb_bin_04.trk" into cereb_bin_04
  file "${sid}_all_cereb_nocx_out_cereb_in_redN_and_thalamus.txt" optional true
  file "${sid}_cereb_bin_04.txt" optional true

  script:
  """
  scil_filter_tractogram.py ${tractogram} ${sid}_all_cereb_nocx_out_cereb_in_redN_and_thalamus.trk \
    --drawn_roi ${params.roisFolder}${params.atlas.BGICCPMCPantMedulla} any exclude \
    --drawn_roi ${params.roisFolder}${params.atlas.CGMSWMDWM} any exclude -f

  scil_streamlines_math.py difference ${tractogram} ${sid}_all_cereb_nocx_out_cereb_in_redN_and_thalamus.trk ${sid}_cereb_bin_04.trk -f

  if ${params.keep}
  then
    scil_count_streamlines.py ${sid}_all_cereb_nocx_out_cereb_in_redN_and_thalamus.trk > ${sid}_all_cereb_nocx_out_cereb_in_redN_and_thalamus.txt
    scil_count_streamlines.py ${sid}_cereb_bin_04.trk > ${sid}_cereb_bin_04.txt
  fi
  """
}

/*
  Brainstem
*/

process brainstem_endInBrainstem {

  cpus 1
  tag = "Brainstem - Split both_ends and either_end"

  input:
    set sid, file(tractogram) from inBrainstem

  output:
  set sid, "${sid}_all_brainstem_both_ends.trk" into bothEndInBrainstem
  set sid, "${sid}_all_brainstem_either_end.trk" into eitherEndInBrainstem
  file "${sid}_all_brainstem_both_ends.txt" optional true
  file "${sid}_all_brainstem_either_end.txt" optional true

  script:
    atlas=params.roisFolder+params.atlas.brainstem
    mode=params.mode.both_ends
    criteria=params.criteria.include
    out_extension='all_brainstem_both_ends'
    remaining_extension='all_brainstem_either_end'
    basename="${sid}"

    template "filter_with_atlas.sh"
}

process brainstem_fastMedulla{

  cpus 1
  tag = "Brainstem - Fast segmentation Medulla"

  input:
    set sid, file(tractogram) from bothEndInBrainstem

  output:
    set sid, "${sid}_all_brainstem_both_ends_Medulla.trk" into EndInMedulla
    set sid, "${sid}_all_brainstem_both_ends_noMedulla.trk" into NotEndInMedulla
    file "${sid}_all_brainstem_both_ends_Medulla.txt" optional true
    file "${sid}_all_brainstem_both_ends_noMedulla.txt" optional true

  script:
  """
  scil_filter_tractogram.py ${tractogram} ${sid}_all_brainstem_both_ends_Medulla.trk \
    --drawn_roi ${params.roisFolder}${params.atlas.CGMSWMsubcortical} both_ends exclude \
    --drawn_roi ${params.roisFolder}${params.atlas.midbrain} both_ends exclude \
    --drawn_roi ${params.roisFolder}${params.atlas.CP} both_ends exclude \
    --drawn_roi ${params.roisFolder}${params.atlas.pons} both_ends exclude \
    --drawn_roi ${params.roisFolder}${params.atlas.medulla} both_ends include \
    --drawn_roi ${params.roisFolder}${params.atlas.pons} any exclude -f

  scil_filter_tractogram.py ${tractogram} ${sid}_all_brainstem_both_ends_noMedulla.trk \
    --drawn_roi ${params.roisFolder}${params.atlas.CGMSWMsubcortical} both_ends exclude \
    --drawn_roi ${params.roisFolder}${params.atlas.midbrain} both_ends exclude \
    --drawn_roi ${params.roisFolder}${params.atlas.CP} both_ends exclude \
    --drawn_roi ${params.roisFolder}${params.atlas.pons} both_ends exclude \
    --drawn_roi ${params.roisFolder}${params.atlas.medulla} both_ends exclude -f
  """
}

process brainstem_noCP{

  cpus 1
  tag = "Brainstem - No  ends in cereberal peduncles"

  input:
    set sid, file(tractogram) from NotEndInMedulla

  output:
    set sid, "${sid}_all_brainstem_both_ends_noMedulla_noCP.trk" into brainstemNoCP
    set sid, "${sid}_brainstem_bin_02.trk" into brainstemBin02
    file "${sid}_all_brainstem_both_ends_noMedulla_noCP.txt" optional true
    file "${sid}_brainstem_bin_02.txt" optional true

  when:
    params.keep

  script:
    atlas=params.roisFolder+params.atlas.CP
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

process cc_homotopic {
  cpus 1
  tag = "CC Homo AGWM"

  input:
  set sid, file(tractogram) from CC_for_homotopic
  each pair from cc_homotopic_pairs

  output:
  set sid, "${sid}_cc_homotopic_${pair}.trk" into CC_Homotopic

  script:
    filtering_list=params.filteringListsFolder+"CC_homo_${pair}_filtering_list_f.txt"
    out_extension="cc_homotopic_${pair}"
    remaining_extension="garbage_${pair}"
    basename="${sid}"

    template "filter_with_list.sh"
}

/*
  ASSO
*/

assoLists=params.assoLists?.tokenize(',')

process asso_part1 {
  cpus 1
  tag = "Asso filtering to get asso_SLS->SLF+AF and asso_ILS->IFOF+UF"

  input:
  set sid, val(side), file(tractogram) from assoAllIntraInter_for_filtering
  each assoList from assoLists

  output:
  set sid, val(side), "${sid}_asso_${assoList}_${side}.trk" into assoAllIntraInter_filtered
  set sid, "${sid}_asso_lost_${assoList}_${side}.trk" into assoLost3
  file "${sid}_asso_${assoList}_${side}.txt" optional true
  file "${sid}_asso_lost_${assoList}_${side}.txt" optional true

  script:
  filtering_list=params.filteringListsFolder+"ASSO_${assoList}_${side}_filtering_list.txt"
  out_extension="asso_${assoList}_${side}"
  remaining_extension="asso_lost_${assoList}_${side}"
  basename="${sid}"

  template "filter_with_list.sh"
}
