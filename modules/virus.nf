params.vs2_groups = 'dsDNAphage,NCLDV,RNA,ssDNA,lavidaviridae'
params.vs2_minlen = 300
params.vb_sens = false

process VirSorter2 {
    conda '/shared/homes/120274/miniconda3/envs/virsorter'
    publishDir params.outdir, mode: 'copy', saveAs: {fn -> "$name/${asm_run}spades/virus/$fn"}
    scratch '/scratch/u120274'
    cpus 40
    memory '64 GB'
    
    input:
    tuple val(name), val(asm_run), path(contigs)

    output:
    tuple val(name), val(asm_run), path('virsorter2_out')

    """
    virsorter run \
           --min-length ${params.vs2_minlen} \
           --include-groups ${params.vs2_groups} \
           -j ${task.cpus-2} \
           -i $contigs \
           -w virsorter2_out \
           --tmpdir ./localtmp
    """
}

process VirBot {
    conda '/shared/homes/120274/miniconda3/envs/virbot'
    publishDir params.outdir, mode: 'copy', saveAs: {fn -> "$name/${asm_run}spades/virus/$fn"}
    cpus 8
    memory '64 GB'

    input:
    tuple val(name), val(asm_run), path(contigs)

    output:
    tuple val(name), val(asm_run), path('virbot_*out')

    script:
    if (params.vb_sens) {
        """
        virbot --sen \
               --threads ${task.cpus} \
               --input $contigs \
               --output virbot_sensitive_out
        """    
    }
    else {
        """
        virbot --threads ${task.cpus} \
               --input $contigs \
               --output virbot_default_out
        """
    }
}

workflow find_rnavirus {
    take:
    contig_sets
    
    main:
    VirBot(contig_sets)
    VirSorter2(contig_sets)

    emit:
    virbot = VirBot.out
    virsorter = VirSorter2.out
}

workflow find_rnavirus2 {
    take:
    contig_sets

    main:
    VirBot(contig_sets)

    emit:
    VirBot.out
}

workflow {

    contig_sets = Channel.fromPath(params.contig_sets).splitCsv()

    VirBot(contig_sets)
    VirSorter2(contig_sets)
}
