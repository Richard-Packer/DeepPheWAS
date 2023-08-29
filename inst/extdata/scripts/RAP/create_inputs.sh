#!/bin/bash

DATASET=/app43027_20230422174636.dataset
DEST=/deep_phewas/inputs

dx run table-exporter \
    -idataset_or_cohort_or_dashboard=$DATASET \
    -ientity="participant" \
    -ifield_names_file_txt="/deep_phewas/inputs/fields_use.txt" \
    -icoding_option=RAW \
    -iheader_style=UKB-FORMAT \
    -ioutput="ukb_test_tab" \
    --destination $DEST \
    --name extract_fields \
    --brief -y

for ENTITY in death death_cause gp_clinical gp_scripts hesin hesin_diag hesin_oper; do
    dx run table-exporter \
        -idataset_or_cohort_or_dashboard=$DATASET \
        -ientity="$ENTITY" \
        -ioutput_format=TSV \
        -icoding_option=RAW \
        -iheader_style=UKB-FORMAT \
        -ioutput="$ENTITY" \
        --destination $DEST \
        --name extract_$ENTITY \
        --brief -y
done
