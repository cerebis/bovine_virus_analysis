params.vs2_groups = 'dsDNAphage,NCLDV,RNA,ssDNA,lavidaviridae'
params.vs2_minlen = 300

process VirSorter2 {
    conda '/shared/homes/120274/miniconda3/envs/vs2'
    publishDir params.outdir, mode: 'copy', saveAs: {fn -> "$name/virus/$asm_run/$fn"}
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
           --tmpdir localtmp
    """
}

process VirBot {
    conda '/shared/homes/120274/miniconda3/envs/virbot'
    publishDir params.outdir, mode: 'copy', saveAs: {fn -> "$name/virus/$asm_run/$fn"}
    cpus 16
    memory '64 GB'

    input:
    tuple val(name), val(asm_run), path(contigs)

    output:
    tuple val(name), val(asm_run), path('virbot_out')

    """
    virbot --sen \
           --threads ${task.cpus} \
           --input $contigs \
           --output virbot_out
    """
}
 
workflow find_rnavirus {
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
