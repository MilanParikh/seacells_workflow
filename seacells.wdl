version 1.0
workflow seacells {
    input {
    	String output_directory
        File anndata_file
        String sample_col
        Int n_cells_per_seacell = 75
        Int n_waypoint_eigenvalues = 10
        Int min_iter = 10
        Int max_iter = 100
        #general parameters
        Int cpu = 16
        String memory = "64G"
        Int extra_disk_space = 32
        String docker = "mparikhbroad/seacells_workflow:latest"
        Int preemptible = 2
    }

    String output_directory_stripped = sub(output_directory, "/+$", "")

    call split_anndata {
        input:
            anndata_file = anndata_file,
            sample_col = sample_col,
            cpu=cpu,
            memory=memory,
            extra_disk_space = extra_disk_space,
            docker=docker,
            preemptible=preemptible
    }
    scatter(sample_anndata_file in split_anndata.sample_anndata_files) {
        call run_seacells {
            input:
                sample_anndata_file = sample_anndata_file,
                n_cells_per_seacell = n_cells_per_seacell,
                n_waypoint_eigenvalues = n_waypoint_eigenvalues,
                min_iter = min_iter,
                max_iter = max_iter,
                cpu=cpu,
                memory=memory,
                extra_disk_space = extra_disk_space,
                docker=docker,
                preemptible=preemptible
        }
    }

    call combine_anndatas {
        input:
            output_dir = output_directory_stripped,
            seacell_anndata_files = run_seacells.seacell_anndata_file,
            cpu=cpu,
            memory=memory,
            extra_disk_space = extra_disk_space,
            docker=docker,
            preemptible=preemptible
    }

    output {
        File combined_anndata_file = combine_anndatas.combined_anndata_file
    }
}

task split_anndata {

    input {
        File anndata_file
        String sample_col
        String memory
        Int extra_disk_space
        Int cpu
        String docker
        Int preemptible
    }

    command <<<
        set -e

        mkdir -p outputs

        python <<CODE
        import os
        import scanpy as sc

        adata = sc.read_h5ad("~{anndata_file}")

        for sample in adata.obs['~sample_col'].unique():
            temp = adata[adata.obs['~{sample_col}'] == sample].copy()
            temp.write(f"outputs/{sample}.h5ad")

        CODE
    >>>

    output {
        Array[File] sample_anndata_files = glob('outputs/*.h5ad')
    }

    runtime {
        docker: docker
        memory: memory
        bootDiskSizeGb: 12
        disks: "local-disk " + (ceil(size(anndata_file, "GB")*2) + extra_disk_space) + " HDD"
        cpu: cpu
        preemptible: preemptible
    }

}

task run_seacells {

    input {
        File sample_anndata_file
        Int n_cells_per_seacell
        Int n_waypoint_eigenvalues
        Int min_iter
        Int max_iter
        String memory
        Int extra_disk_space
        Int cpu
        String docker
        Int preemptible
    }

    command <<<
        set -e

        mkdir -p outputs

        python <<CODE
        import os
        import scanpy as sc
        import matplotlib
        import SEACells

        adata = sc.read_h5ad("~{sample_anndata_file}")

        n_SEACells = np.ceil(adata.shape[0]/~{n_cells_per_seacell})
        print(f"***Number of SEACells = {n_SEACells}***")
        n_waypoint_eigs = ~{n_waypoint_eigenvalues}

        model = SEACells.core.SEACells(adata, 
                                build_kernel_on='X_pca', 
                                n_SEACells=n_SEACells, 
                                n_waypoint_eigs=n_waypoint_eigs,
                                convergence_epsilon = 1e-5)
        model.construct_kernel_matrix()
        M = model.kernel_matrix

        model.initialize_archetypes()

        model.fit(min_iter=~{min_iter}, max_iter=~{max_iter})

        adata.write("outputs/sample_anndata.h5ad")

        CODE
    >>>

    output {
        File seacell_anndata_file = 'outputs/sample_anndata.h5ad'
    }

    runtime {
        docker: docker
        memory: memory
        bootDiskSizeGb: 12
        disks: "local-disk " + (ceil(size(sample_anndata_file, "GB")*2) + extra_disk_space) + " HDD"
        cpu: cpu
        preemptible: preemptible
    }

}

task combine_anndatas {

    input {
        String output_dir
        Array[File] seacell_anndata_files
        String memory
        Int extra_disk_space
        Int cpu
        String docker
        Int preemptible
    }

    command <<<
        set -e

        mkdir -p outputs

        python <<CODE
        import os
        import scanpy as sc
        import anndata

        list_of_seacell_anndata_files = ["~{sep='","' seacell_anndata_files}"]
        adata_list = []
        for seacell_anndata_file in list_of_seacell_anndata_files:
            temp = sc.read_h5ad(seacell_anndata_file)
            adata_list.append(temp)

        adata = anndata.concat(adata_list, join='outer', merge='same', uns_merge='same', index_unique=None)
        adata.write('outputs/combined_anndata.h5ad')

        CODE

        gsutil -m cp outputs/combined_anndata.h5ad ~{output_dir}/
    >>>

    output {
        File combined_anndata_file = 'outputs/combined_anndata.h5ad'
    }

    runtime {
        docker: docker
        memory: memory
        bootDiskSizeGb: 12
        disks: "local-disk " + (ceil(size(seacell_anndata_files, "GB")*2) + extra_disk_space) + " HDD"
        cpu: cpu
        preemptible: preemptible
    }

}