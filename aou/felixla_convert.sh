#!/bin/bash
# Convert phased VCF + RFMix .msp local-ancestry calls into the FELIXla format
# consumed by step 2. Positional args: phased VCF, .msp, number of ancestries,
# number of samples, output prefix.
set -o pipefail
set -o errexit

rfmix_msp_to_tractor_hybrid \
    "${PHASED_VCF}" \
    "${MSP_FILE}" \
    5 \
    7669 \
    "${prefix}"
