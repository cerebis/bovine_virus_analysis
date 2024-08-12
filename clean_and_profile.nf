include { reads_profile; contig_assign } from './modules/tax_profile.nf'
include { trim; clean_up               } from './modules/clean_reads.nf'
include { kmer_spades; meta_spades     } from './modules/assemble.nf'
include { find_rnavirus                } from './modules/virus.nf'

params.trim_only = false
params.rank = "genus"
params.min_len = 5000

/*
 * global args:
 * --outdir [PATH] --read_sets [CSV]
 *
 * clean_up args:
 * --human_index [PATH] --host_index [PATH]
 *
 * profile args:
 * --sm_ref [PATH] --sm_tax [PATH] --cat_ref [PATH] --cat_tax [PATH] --min_len [INT] --rank [RANK]
 *
 * assembly args:
 * --asm_mode [STR] Optional: ( --extra_kmers [INT,INT,INT,...] | --kmer_series [INT,INT,INT,...] )
 */

process FasterqDump {
    publishDir params.outdir, mode: 'copy', saveAs: {fn -> "$name/raw_reads/$fn"}
    cpus 8
    memory '32 GB'

    input:
    tuple val(name), path(run) from runs

    output:
    tuple name, path('*_1.fastq.gz'), path('*_2.fastq.gz')

    """
    /shared/homes/120274/sratoolkit.3.1.1-ubuntu64/bin/sratools.3.1.1/fasterq-dump \
        -O . -e ${task.cpus} -S -v
    pigz -p ${task.cpus} *fastq
    """
}

workflow extract_fastq {

    take:
    runs

    main:
    FasterqDump(runs)

    emit:
    FasterqDump.out
}

workflow {

    read_sets = extract_fastq(Channel.fromPath(params.runs).splitCsv())
//    read_sets = Channel.fromPath(params.read_sets).splitCsv()
    
    if (params.trim_only == true) {
        read_sets = trim(read_sets)
    }
    else {
        read_sets = clean_up(read_sets)
    }

    reads_profile(read_sets,
                  Channel.value(params.sm_ref),
                  Channel.value(params.sm_tax),
                  Channel.value(params.min_len),
                  Channel.value(params.rank))

//    asm = kmer_spades(read_sets,
//                      Channel.value(params.kmer_series))
    asm = meta_spades(read_sets)

//    if (params.asm_mode == "rna") {
//        asm = rna_spades(read_sets)
//    }
//    else if (params.asm_mode == "meta") {
//        asm = meta_spades(read_sets)
//    }

//    contig_assign(asm.asm_seqs,
//                  Channel.fromPath(params.cat_db),
//                  Channel.fromPath(params.cat_tax))

    find_rnavirus(asm.asm_seqs)
}
