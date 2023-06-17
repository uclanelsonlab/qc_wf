version 1.0

import "tasks/multiqc.wdl" as multiqc_task

task find_files {
    meta {
        author: "George Carvalho"
        email: "gcarvalhoneto@mednet.ucla.edu"
        description: "## Search for the files to run MultiQC for multiple samples"
    }
    input {
        String dir_path
        String? project
    }
    String DX_PROJECT_CONTEXT_ID=select_first([project, "project-F292qGQ02k8Y48QfFzJqX2j0"])
    command <<<
        set -uexo pipefail
        dx find data --path "~{DX_PROJECT_CONTEXT_ID}:~{dir_path}" --class file --delimiter "," | cut -d, -f4 | grep qc | cut -d, -f2 > files.list
        mkdir stats_files/
        for file in $(cat files.list); do 
            dx download "~{DX_PROJECT_CONTEXT_ID}:${file}" -f -o stats_files/; 
        done
        ls stats_files/
    >>>
    runtime {
        dx_instance_type: "mem1_ssd1_v2_x16"
        dx_ignore_reuse: true
        dx_restart: object {
            default: 1,
            max: 1,
            errors: object {
                UnresponsiveWorker: 2,
                ExecutionError: 2,
            }
        }
        dx_timeout: "5H30M"
        dx_access: object {
            network: ["*"],
            developer: true
        }
    }
    output {
        Array[File] stats_files = glob("stats_files/*")
    }
}

workflow find_data_multiqc_wf {
    parameter_meta {
        dir_path: {
            description: "Directory path to search for stats files"
        }
        prefix: {
            description: "Sample prefix to be used in output creation"
        }
        multiqc_image: {
            description: "TAR zip docker image from MultiQC stored on DNAnexus (file-GX16qZ802k8xKbkF797xKqP4)"
        }
    }
    input {
        String dir_path
        String? project
        # MultiQC
        String prefix
        String? multiqc_image
    }
    call find_files {
        input:
            dir_path=dir_path,
            project=project
    }
    call multiqc_task.multiqc_array as multiqc {
        input:
            stats_files=find_files.stats_files,
            prefix=prefix,
            multiqc_image=multiqc_image
    }
    output {
        File multiqc_html = multiqc.multiqc_html
        File multiqc_data = multiqc.multiqc_data
    }
}