process GzipValidator {
    publishDir 'validate_out'
    cpus 1
    memory '8 GB'

    input:
    tuple val(acc), path(r1), path(r2)    

    output:
    path("${acc}.out"), optional: true

    """
    gzip -t $r1 || echo "$r1"  > ${acc}.out 
    gzip -t $r2 || echo "$r2" >> ${acc}.out
    """
}

workflow {

    read_sets = Channel.fromPath(params.input).splitCsv()

    GzipValidator(read_sets)
}
