The `src` directory includes local and external input data used in the analyses.

|File/Directory Name|Description|Source|Comments|
|---|---|---|---|
|`parcellations/*h.schaefer-*.annot`|Schaefer-N parcellation in fsaverage5 space|[ThomasYeoLab / CBIG](https://github.com/ThomasYeoLab/CBIG/tree/master/stable_projects/brain_parcellation/Schaefer2018_LocalGlobal/Parcellations/FreeSurfer5.3/fsaverage5/label)|Originally named `*h.Schaefer2018_*Parcels_7Networks_order.annot`|
|`parcellations/*h.schaefer-*.fsa.annot`|Schaefer-N parcellation in fsaverage space|[ThomasYeoLab / CBIG](https://github.com/ThomasYeoLab/CBIG/tree/master/stable_projects/brain_parcellation/Schaefer2018_LocalGlobal/Parcellations/FreeSurfer5.3/fsaverage/label)|Originally named `*h.Schaefer2018_*Parcels_7Networks_order.annot`|
|`parcellations/lut_*_mics.csv`|Parcel infos from micapipe|[MICA-MNI / micapipe](https://github.com/MICA-MNI/micapipe/)|Cortical parcels are renamed to match FC parcel names|
|`parcellations/tpl-MNI152_desc-schaefer-*_parcellation_1mm.nii.gz`|Schaefer-N parcellation in MNI152 space|[ThomasYeoLab / CBIG](https://github.com/ThomasYeoLab/CBIG/tree/master/stable_projects/brain_parcellation/Schaefer2018_LocalGlobal/Parcellations/MNI)|Originally named `Schaefer2018_*Parcels_7Networks_order_FSLMNI152_1mm.nii.gz`|
|`parcellations/Schaefer2018_*Parcels_7Networks_order.dlabel.nii`|Schaefer-N parcellation in fsLR space|[ThomasYeoLab / CBIG](https://github.com/ThomasYeoLab/CBIG/tree/master/stable_projects/brain_parcellation/Schaefer2018_LocalGlobal/Parcellations/HCP/fslr32k/cifti)||
|`parcellations/Schaefer2018_*Parcels_7Networks_order.txt`|Schaefer-N parcellation labels|[ThomasYeoLab / CBIG](https://github.com/ThomasYeoLab/CBIG/tree/master/stable_projects/brain_parcellation/Schaefer2018_LocalGlobal/Parcellations/MNI/freeview_lut)||
|`fsaverage_*h.pial_semi_inflated`|Semi-inflated fsaverage mesh|Freesurfer 7.4||
|`PET_nifti_images`|PET images of NMDA and GABA_A/BZ receptors|[netneurolab / hansen_receptors](https://github.com/netneurolab/hansen_receptors/tree/main/data)||
|`tpl-fs_LR_hemi-L_den-32k_desc-LTC_G1.shape.gii`|Laminar thickness covariance G1 in fsLR space|[amnsbr / laminar_organization](https://github.com/amnsbr/laminar_organization)||
|`ahba_parc-schaefer100_hemi-L_ibf_threshold-0.5_missing-centroids.npz`|AHBA gene expression data of Schaefer-100 parcels in left hemisphere|[amnsbr / laminar_organization](https://github.com/amnsbr/laminar_organization)||