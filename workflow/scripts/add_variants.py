import os
import sys
import datetime
import pandas as pd
import numpy as np
from itertools import chain, combinations
import concurrent.futures
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from intervaltree import IntervalTree  
import copy
import re
from pathlib import Path
import time
from concurrent.futures import as_completed, TimeoutError

# ------------------------------------------------------------
# 1) Configuration du log : redirige stdout & stderr vers un fichier
# ------------------------------------------------------------

# Utilisation de l'output Snakemake pour générer le log 
output_path = Path(snakemake.output.mutated_output_fasta)
log_path = output_path.with_suffix(".log")  # ex: 1_to_9_mutated_transcripts.log
log_fh = open(log_path, "w", buffering=1)

def log(msg: str):
    ts = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    log_fh.write(f"{ts} - {msg}\n")
    log_fh.flush()

# Redirige stdout & stderr vers le log
sys.stdout = log_fh
sys.stderr = log_fh

# ------------------------------------------------------------
# 2) Fonctions du pipeline
# ------------------------------------------------------------

def read_vcf(vcf_file: str) -> pd.DataFrame:
    """
    Lit un VCF simple (pas bgzippé) et retourne un DataFrame avec
    colonnes ['chromosome','position','ref_nucleotide','alt_nucleotide']
    """
    mutations = []
    with open(vcf_file, 'r') as v:
        for line in v:
            if line.startswith('#'):
                continue
            chrom, pos, _, ref, alt, *rest = line.strip().split('\t')
            alts = alt.split(',')
            for a in alts:
                mutations.append((chrom, int(pos), ref, a))

    df = pd.DataFrame(mutations,
                      columns=['chromosome','position','ref_nucleotide','alt_nucleotide'])
    df['position'] = df['position'].astype(np.uint32)
    log(f"1) VCF chargé avec succès → {len(df)} mutations")
    return df

# Génère toutes les combinaisons non vides
def all_subsets(iterable, max_size=None):
    s = list(iterable)
    max_len = len(s) if max_size is None else min(max_size, len(s))
    return chain.from_iterable(combinations(s, r) for r in range(1, max_len + 1))


# Variables globales partagées par les workers
global_genome = None
global_gtf_exons = None
global_trees = None

def init_worker(genome_pkl: str, tree_pkl: str, exon_parquet: str):
    global global_genome, global_gtf_exons, global_trees
    import pickle
    with open(genome_pkl, "rb") as f:
        genome_dict = pickle.load(f)
    # Wrap into SeqRecord objects
    global_genome = {
        tid: SeqRecord(Seq(seq), id=tid, description="")
        for tid, seq in genome_dict.items()
    }
    with open(tree_pkl, "rb") as f:
        global_trees = pickle.load(f)
    global_gtf_exons = pd.read_parquet(exon_parquet)

    for exon in global_gtf_exons.itertuples():
        chrom = exon.chromosome
        start = exon.start
        end = exon.end
        if start > end:
            start, end = end, start  # inverse si start > end
        end += 1  # intervaltree half-open
        tree = global_trees.setdefault(chrom, IntervalTree())
        tree.addi(start, end, exon.transcript_id)
    log(f"Clés global_trees : {list(global_trees.keys())}")
    chrom = "1"
    df1 = global_gtf_exons[global_gtf_exons["chromosome"] == chrom]
    log(f"Exons chr{chrom}: start min={df1['start'].min()}, start max={df1['start'].max()}, end min={df1['end'].min()}, end max={df1['end'].max()}")

# Récupère et assemble les exons d'un transcript
def _print_transcript_exons(transcript_id):
    global global_genome, global_gtf_exons
    exons = global_gtf_exons[global_gtf_exons['transcript_id'] == transcript_id]
    if exons.empty or transcript_id not in global_genome:
        return None, None
    strand = exons.iloc[0]['strand']
    sorted_exons = exons.sort_values('start', ascending=(strand == '+'))
    seq = str(global_genome[transcript_id].seq)
    mapping, recon = {}, []
    offset = 0
    for _, exon in sorted_exons.iterrows():
        length = exon['end'] - exon['start'] + 1
        coords = (range(exon['start'], exon['end']+1)
                  if strand=='+' else
                  range(exon['end'], exon['start']-1, -1))
        for i, gpos in enumerate(coords):
            mapping[gpos] = offset + i
        recon.append(seq[offset:offset+length])
        offset += length
    return "".join(recon), mapping

# Applique les mutations sur la séquence du transcript
def _apply_mutations_to_transcript(transcript_id, seq, mapping, muts):
    muts_in = [copy.deepcopy(m) for m in muts if m['position'] in mapping]

    applied_mutations = []
    for m in muts_in:
        m['transcript_pos'] = mapping[m['position']]
    muts_in.sort(key=lambda x: x['transcript_pos'])
    new_seq, shift = seq, 0
    for m in muts_in:
        p   = m['transcript_pos'] + shift
        ref = m['ref_nucleotide']
        alt = m['alt_nucleotide']
        if new_seq[p:p+len(ref)] != ref:
            log(f"[WARN] Mutation ignorée pour {transcript_id} à pos transcript {p+1}: attendu {ref}, trouvé {new_seq[p:p+len(ref)]}")
            continue
        log(f"[DEBggggUG] Applying: pos={p}, ref={ref}, alt={alt}")

        new_seq = new_seq[:p] + alt + new_seq[p+len(ref):]
        
        shift += len(alt) - len(ref)
        applied_mutations.append(m)
    return new_seq, applied_mutations

# Traite un transcript dans un worker
def _process_transcript_worker(transcript_id, mutations):
    global global_genome, global_gtf_exons
    log(f"Début du traitement pour {transcript_id} ({len(mutations)} mutations)")
    exons = global_gtf_exons[global_gtf_exons['transcript_id'] == transcript_id]
    # skip lncRNA
#    if exons.empty or exons.iloc[0]['gene_biotype'] == 'lncRNA':
#        return []
    seq, mapping = _print_transcript_exons(transcript_id)
    if seq is None:
        return []
    out = []

    wt_rec = SeqRecord(
        Seq(seq),
        id=transcript_id,
        description=global_genome[transcript_id].description
    )
    out.append(wt_rec)

    for combo in all_subsets(mutations, max_size=4):

        if any(m['position'] not in mapping for m in combo):
            continue
        ms, applied = _apply_mutations_to_transcript(transcript_id, seq, mapping, combo)
        if len(applied) != len(combo):
            log(f"[INFO] Combo ignoré pour {transcript_id} : demandé {len(combo)} mutations, appliquées {len(applied)}")
            continue

        signature = "_".join(f"{m['ref_nucleotide']}({m['alt_nucleotide']})" for m in applied)
        positions = ";".join(str(mapping[m['position']] + 1) for m in applied)

        new_id = f"{transcript_id}_mut_{signature}_pos={positions}"
        new_description = global_genome[transcript_id].description
        log(f"[DEBUG] type(ms)={type(ms)} value(ms)={ms}")
        rec = SeqRecord(ms if isinstance(ms, Seq) else Seq(ms), id=new_id, description=new_description)
        out.append(rec)
    return out 

# Génère les transcrits mutés en parallèle en utilisant l'IntervalTree
def _generate_mutated_transcripts_parallel(genome_pkl, tree_pkl, exon_parquet, vcf_df, max_workers=12):

    log("2) Chargement des exons pour le main…")
    # Reconstruit les arbres pour le thread principal
    init_worker(genome_pkl, tree_pkl, exon_parquet)
    for chrom, tree in global_trees.items():
        intervals = list(tree)
        if intervals:
            min_start = min(iv.begin for iv in intervals)
            max_end = max(iv.end for iv in intervals)
            log(f"IntervalTree {chrom}: {len(intervals)} intervalles, start min={min_start}, end max={max_end}")
        else:
            log(f"IntervalTree {chrom}: AUCUN intervalle.")

    log(f"   → {len(global_gtf_exons)} exons chargés.")
    # 3) Mapping transcript → mutations via IntervalTree
    tx_muts = {}
    for mut in vcf_df.to_dict('records'):
        chrom = mut['chromosome']
        pos   = mut['position']
        if chrom not in global_trees:
            continue
        for iv in global_trees[chrom][pos]:
            tx = iv.data
            tx_muts.setdefault(tx, []).append(mut)
    # Filtre pour ne garder que les transcrits qui ont <= 5 mutations
    tx_muts = {tx: muts for tx, muts in tx_muts.items() if len(muts) <= 5}
    for tx, muts in tx_muts.items():
        log(f"→ {tx} a {len(muts)} mutations (positions uniques : {len({m['position'] for m in muts})})")

    total = len(tx_muts)
    log(f"3) Mapping construit → {total} transcripts à traiter (≤5 mutations chacun).")
    # 4) Exécution parallèle
    log("4) Lancement du ProcessPoolExecutor…")

    with concurrent.futures.ProcessPoolExecutor(
        max_workers=max_workers,
        initializer=init_worker,
        initargs=(genome_pkl, tree_pkl, exon_parquet)
    ) as exe:
        futures = {
            exe.submit(_process_transcript_worker, tx, muts): tx
            for tx, muts in tx_muts.items()
        }

        for fut in as_completed(futures):
            tx_name = futures[fut]
            try:
                start = time.time()
                results = fut.result(timeout=300)  # timeout : 300 sec (5 min)
                duration = time.time() - start
                log(f"[OK] {tx_name} traité en {duration:.2f} s")
                for rec in results:
                    yield rec
            except TimeoutError:
                log(f"[TIMEOUT] {tx_name} a été ignoré après 5 minutes.")
            except Exception as e:
                log(f"[ERROR] {tx_name} a échoué : {e}")

# Ecrit les transcrits mutés par batch pour économiser la mémoire
def save_mutated_transcripts_in_batches_parallel(
    genome_pkl, tree_pkl, exon_parquet, vcf_df, out_fasta,
    batch_size=5, max_workers=12
):

    batch = []
    with open(out_fasta, "w") as out:
        for rec in _generate_mutated_transcripts_parallel(
            genome_pkl, tree_pkl, exon_parquet, vcf_df, max_workers
        ):
            batch.append(rec)
            if len(batch) >= batch_size:
                SeqIO.write(batch, out, "fasta")
                out.flush(); os.fsync(out.fileno())
                batch.clear()
        if batch:
            SeqIO.write(batch, out, "fasta")
            out.flush(); os.fsync(out.fileno())

def filter_transcripts(input_fasta, output_fasta):
    """
    Parcourt le FASTA d'entrée et regroupe les séquences par identifiant de transcript.
    Ne conserve que les groupes qui contiennent à la fois :
      - Un transcript WT (pas de "_mut_" dans l'ID)
      - Au moins un transcript muté (ID contenant "_mut_")
    """
    groups = {}
    total_records = 0

    for record in SeqIO.parse(input_fasta, "fasta"):
        total_records += 1
        m = re.search(r"(ENST\d+)", record.id)
        if not m:
            continue
        tid = m.group(1)
        groups.setdefault(tid, []).append(record)

    filtered = []
    kept = 0
    for tid, recs in groups.items():
        has_wt  = any("_mut_" not in r.id for r in recs)
        has_mut = any("_mut_" in r.id     for r in recs)
        if has_wt and has_mut:
            filtered.extend(recs)
            kept += 1

    SeqIO.write(filtered, output_fasta, "fasta")
    print(f" Total lus : {total_records}")
    print(f" Groupes gardés (WT+mut) : {kept}")
    print(f" Total écrits : {len(filtered)}")

# ------------------------------------------------------------
# 3) Fonction main()
# ------------------------------------------------------------
def main():
    log("=== DÉMARRAGE DU PIPELINE ===")
    vcf_file          = snakemake.input.vcf
    transcripts_fasta = snakemake.input.transcripts_fasta
    mutated_output_fa = snakemake.output.mutated_output_fasta
    exon_parquet      = snakemake.input.exon_parquet
    genome_pkl = snakemake.input.genome_pkl
    tree_pkl = snakemake.input.tree_pkl

#    vcf_file = "20QC_variant_tronque.vcf"
#    transcripts_fasta= "breast_cancer/pickle/transcriptome.fa"
#    transcripts_fasta = "breast_cancer/pickle/gencode.v47.transcripts.fa"
#    mutated_output_fa  = "mutated_transcripts.fa"
#    exon_parquet = "exon_data.parquet"  
    final_fasta = "final_transcriptome.fa"
    print("bob")
    try:
        df = read_vcf(vcf_file)
        save_mutated_transcripts_in_batches_parallel(
            genome_pkl, tree_pkl, exon_parquet, df,
            mutated_output_fa,
            batch_size=5,
            max_workers=12
        )
        filter_transcripts(mutated_output_fa, final_fasta)
        log(f"=== FASTA filtré écrit dans {final_fasta} ===")
        log("=== PIPELINE TERMINÉ AVEC SUCCÈS ===")
    except Exception as e:
        log(f"!!! ERREUR CRITIQUE : {e}")
    finally:
        log_fh.close()

if __name__ == "__main__":
    main()
