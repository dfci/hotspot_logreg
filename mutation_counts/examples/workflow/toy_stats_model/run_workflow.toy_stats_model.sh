#!/bin/bash -eux

/usr/bin/java \
    -jar /data/pipeline/v1/tools/cromwell/v0.21/cromwell-0.21.jar run \
    workflow.wdl \
    workflow_inputs.toy_stats_model.json
