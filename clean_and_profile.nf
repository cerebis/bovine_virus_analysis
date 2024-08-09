include { reads_profile; contig_assign } from './modules/tax_profile.nf'
include { trim; clean_up               } from './modules/clean_reads.nf'
include { rna_spades; meta_spades      } from './modules/assemble.nf'
include { find_rnavirus                } from './modules/virus.nf'

/*
 * global args:
 * --outdir [PATH] --read_sets [CSV]
 *
 * clean_up args:
 * --human_index [PATH] --host_index [PATH]
 *
 * profile args:
 * --reference_db [PATH] --minlen [INT] --rank [RANK] --taxonomy_db [PATH]
 *
 * assembly args:
 * --asm_mode [STR] --extra_kmers [INT,INT,INT,...]
 */

workflow {

    read_sets = Channel.fromPath(params.read_sets).splitCsv()
    
    if (params.trim_only == true) {
        read_sets = trim(read_sets)
    }
    else {
        read_sets = clean_up(read_sets)
    } 

    reads_profile(read_sets)
    
    if (params.asm_mode == "rna") {
        asm = rna_spades(read_sets)
    }
    else if (params.asm_mode == "meta") {
        asm = meta_spades(read_sets)
    }

    contig_assign(asm.asm_seqs, Channel.fromPath(params.cat_db), Channel.fromPath(params.cat_tax))
    find_rnavirus(asm.asm_seqs)
}
