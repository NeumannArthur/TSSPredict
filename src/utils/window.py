import numpy as np
import pandas as pd
from sklearn import preprocessing as pp
import os
from typing import List, Dict 

def load_wiggle(wiggle_path):
    cov = {}
    with open(wiggle_path, "r") as f:
        for line in f:
            if not line or line.startswith("track") or line.startswith("variableStep"):
                continue
            p, v = line.strip().split()[:2]
            cov[int(p)] = float(v)
    return cov

def build_TSS_window(
    tss_df: pd.DataFrame, 
    window_size: int=1000,
    conditions=None,
    strand:str="all",
    preview:bool=False,
    preview_n:int=1001, 
    preview_path:str="/ceph/ibmi/it/thesis_data/neumannarthur/results/previews",
    output_dir:str="/ceph/ibmi/it/thesis_data/neumannarthur/results", 
    wiggle_forward:str="/ceph/ibmi/it/thesis_data/neumannarthur/data/ecoli_data/wiggles/rnaseq_files/SRR21871309_VCE_Non_enriched_ctrl_div_by_8609484.0_multi_by_1000000.0_forward.wig",
    wiggle_backward:str="/ceph/ibmi/it/thesis_data/neumannarthur/data/ecoli_data/wiggles/rnaseq_files/SRR21871309_VCE_Non_enriched_ctrl_div_by_8609484.0_multi_by_1000000.0_reverse.wig"
):
    """ 
    Buid TSS windows per condition and generate strand/condition specific coverage tables
    """

    conditions = list(conditions) if conditions else ["Control"]

    os.makedirs(output_dir, exist_ok=True)
    
    print(f"[Window] Checking for preview: {preview}")
    preview_dir = None
    if preview:
        preview_dir = preview_path
        os.makedirs(preview_dir, exist_ok=True)

    create_master = (strand in ("master", "all"))
    create_plus = (strand in ("+", "all"))
    create_minus = (strand in ("-", "all"))

    coverage_plus = load_wiggle(wiggle_forward) if create_plus else None
    coverage_minus = load_wiggle(wiggle_backward) if create_minus else None

    offsets = np.arange(-window_size, window_size + 1, dtype=int)
    window_total_size = offsets.size 

    sequence = {}
    name, chunks = None, []
    with open("/ceph/ibmi/it/thesis_data/neumannarthur/data/ecoli_data/genome_annotation/GCF_000005845.2_ASM584v2_genomic.fna") as fh:
        for line in fh:
            if line.startswith(">"):
                if name is not None:
                    sequence[name] = "".join(chunks).upper()
                name = line[1:].split()[0]
                chunks = []
            else:
                chunks.append(line.strip())
    if name is not None:
        sequence[name] = "".join(chunks).upper()

    reverse_complement_trans_table = str.maketrans("ACGTNacgtn", "TGCANtgcan")
    def reverse_complement(sequence):
        return sequence.translate(reverse_complement_trans_table)[::-1]

    for cond in conditions:  
        print(f"[Window.{cond}] Starting window creation with size {window_size}") 

        sub = tss_df[tss_df["Genome"] == cond].copy()
        if sub.empty:
            print(f"[Window.{cond}] No rows found! Skipping...")
            continue

        if create_master:
            print(f"[Window.{cond}] Creating master windows...")
            master_record = []
            for _, row in sub.iterrows():
                contig = row["contigID"]
                pos = int(row["contigPos"])
                contig_key = str(contig.split("|")[0].split()[0]) # ? War nötig, weil in tss_df der key manchmal z.B. "NC_000913.3|NC_000913.3" is anstatt nur "NC_000913.3"
                seq = sequence.get(contig_key)
                L = len(seq)
                abs_pos = ((pos + offsets - 1) % L) + 1

                window_characters = [seq[(pos + off - 1) % L] for off in offsets]
                window = "".join(window_characters)
                base_series = window if row["SuperStrand"] == "+" else reverse_complement(window)
                master_record.append({
                        "contigPos": pos,
                        "contigID": contig,
                        "windowBegin": int(abs_pos[0]),
                        "windowEnd": int(abs_pos[-1]),
                        "Genome": row["Genome"],
                        "SuperStrand": row["SuperStrand"],
                        "Bases": base_series
                    })
            master_df = pd.DataFrame(master_record, columns=["contigPos", "contigID", "windowBegin", "windowEnd", "Genome", "SuperStrand", "Bases"])
            master_path = os.path.join(output_dir + f"/{cond}/master/", f"[Window.{cond}]_master.tsv")
            master_df.to_csv(master_path, sep="\t", index=False)
            print(f"[Window.{cond}] Saved master table to {master_path}")
            
            print(f"[Window.{cond}] Creating negative TSS windows...")
            neg_master_records = []
            for (contig, sstrand), grp in sub.groupby(["contigID", "SuperStrand"]):
                contig_key = str(contig.split("|")[0].split()[0])
                seq = sequence.get(contig_key)
                if seq == None:
                    print("WARNING")
                    continue
                L = len(seq)
                tss_positions = sorted(set(int(p) for p in grp["contigPos"].tolist()))
                if len(tss_positions) < 2:
                    continue
                pairs = [(tss_positions[i], tss_positions[(i+1) % len(tss_positions)]) for i in range(len(tss_positions))]
                for first, second in pairs:
                    d_forward = (second - first) if second >= first else (L - first + second)
                    mid = ((first - 1 + d_forward // 2) % L) + 1

                    abs_pos = ((mid + offsets - 1) % L) + 1
                    window_chars = [seq[(mid + off - 1) % L] for off in offsets]
                    window = "".join(window_chars)
                    base_series = window if sstrand == "+" else reverse_complement(window)

                    neg_master_records.append({
                        "contigPos": int(mid),
                        "contigID": contig,
                        "windowBegin": int(abs_pos[0]),
                        "windowEnd": int(abs_pos[-1]),
                        "Genome": cond,
                        "SuperStrand": sstrand,
                        "Bases": base_series
                    })
            if neg_master_records:
                neg_master_df = pd.DataFrame(
                    neg_master_records,
                    columns=["contigPos", "contigID", "windowBegin", "windowEnd", "Genome", "SuperStrand", "Bases"]
                )
                neg_master_path = os.path.join(output_dir + f"/{cond}/master/", f"[Window.{cond}]_master_negative.tsv")
                neg_master_df.to_csv(neg_master_path, sep="\t", index=False)
                print(f"[Window.{cond}] Saved negative master table to {neg_master_path}")

                if preview:
                    preview_path = os.path.join(preview_dir + f"/{cond}/master/", f"[Window.{cond}]_negative_master_preview.tsv")
                    neg_master_df.head(window_total_size).to_csv(preview_path, sep="\t", index=False)
                    print(f"[Window.{cond}] Saved negative master preview to {preview_path}")

            if preview:
                preview_path = os.path.join(preview_dir + f"/{cond}/master/", f"[Window.{cond}]_master_preview.tsv")
                master_df.head(window_total_size).to_csv(preview_path, sep="\t", index=False)
                print(f"[Window.{cond}] Saved master preview to {preview_path}")

        if create_minus:
            print(f"[Window.{cond}] Creating minus strand windows")
            minus_record: List[dict] = []
            minus_sub = sub[sub["SuperStrand"] == "-"]
            for _, row in minus_sub.iterrows():
                contig = row["contigID"]
                pos = int(row["contigPos"])
                contig_key = str(contig.split("|")[0].split()[0])
                seq = sequence.get(contig_key)
                L = len(seq)
                abs_pos = ((pos + offsets - 1) % L) + 1

                window_characters = [seq[(pos + off - 1) % L] for off in offsets]
                window = "".join(window_characters)
                base_series = reverse_complement(window)
                cov_list = [coverage_minus.get(int(p), 0.0) for p in abs_pos]
                minus_record.append({
                        "contigPos": pos,
                        "contigID": contig,
                        "windowBegin": int(abs_pos[0]),
                        "windowEnd": int(abs_pos[-1]),
                        "Genome": row["Genome"],
                        "SuperStrand": row["SuperStrand"],
                        "Bases": base_series,
                        "coverage": pp.normalize(np.array(cov_list).reshape(1, -1), norm="max").tolist()
                    })
            minus_df = pd.DataFrame(minus_record, columns=["contigPos", "contigID", "windowBegin", "windowEnd", "Genome", "SuperStrand", "Bases", "coverage"])
            minus_path = os.path.join(output_dir + f"/{cond}/minus/", f"[Window.{cond}]_minusStrand.tsv")
            minus_df.to_csv(minus_path, sep="\t", index=False)
            print(f"[Window.{cond}] Saved plus strand table to {minus_path}")

            print(f"[Window.{cond}] Creating negative TSS windows...")
            neg_minus_records = []
            for (contig, sstrand), grp in minus_sub.groupby(["contigID", "SuperStrand"]):
                contig_key = str(contig.split("|")[0].split()[0])
                seq = sequence.get(contig_key)
                if seq == None:
                    print("WARNING")
                    continue
                L = len(seq)
                tss_positions = sorted(set(int(p) for p in grp["contigPos"].tolist()))
                if len(tss_positions) < 2:
                    continue
                pairs = [(tss_positions[i], tss_positions[(i+1) % len(tss_positions)]) for i in range(len(tss_positions))]
                for first, second in pairs:
                    d_forward = (second - first) if second >= first else (L - first + second)
                    mid = ((first - 1 + d_forward // 2) % L) + 1

                    abs_pos = ((mid + offsets - 1) % L) + 1
                    window_chars = [seq[(mid + off - 1) % L] for off in offsets]
                    window = "".join(window_chars)
                    base_series = reverse_complement(window)

                    coverage = [coverage_minus.get(int(p), 0.0) for p in abs_pos]

                    neg_minus_records.append({
                        "contigPos": int(mid),
                        "contigID": contig,
                        "windowBegin": int(abs_pos[0]),
                        "windowEnd": int(abs_pos[-1]),
                        "Genome": cond,
                        "SuperStrand": sstrand,
                        "Bases": base_series,
                        "coverage": pp.normalize(np.array(coverage).reshape(1, -1), norm="max").tolist()
                    })
            if neg_minus_records:
                neg_minus_df = pd.DataFrame(
                    neg_minus_records,
                    columns=["contigPos", "contigID", "windowBegin", "windowEnd", "Genome", "SuperStrand", "Bases", "coverage"]
                )
                neg_minus_path = os.path.join(output_dir + f"/{cond}/minus/", f"[Window.{cond}]_minusStrand_negative.tsv")
                neg_minus_df.to_csv(neg_minus_path, sep="\t", index=False)
                print(f"[Window.{cond}] Saved negative minus table to {neg_minus_path}")

                if preview:
                    preview_path = os.path.join(preview_dir + f"/{cond}/minus/", f"[Window.{cond}]_negative_minus_preview.tsv")
                    neg_minus_df.head(window_total_size).to_csv(preview_path, sep="\t", index=False)
                    print(f"[Window.{cond}] Saved negative minus preview to {preview_path}")

            if preview:
                preview_path = os.path.join(preview_dir + f"/{cond}/minus/", f"[Window.{cond}]_minusStrand_preview.tsv")
                minus_df.head(window_total_size).to_csv(preview_path, sep="\t", index=False)
                print(f"[Window.{cond}] Saved minus strand preview to {preview_path}")
        
        if create_plus:
            print(f"[Window.{cond}] Creating plus strand windows")
            plus_record: List[dict] = []
            plus_sub = sub[sub["SuperStrand"] == "+"]
            for _, row in plus_sub.iterrows():
                contig = row["contigID"]
                pos = int(row["contigPos"])
                contig_key = str(contig.split("|")[0].split()[0]) # ? War nötig, weil in tss_df der key manchmal z.B. "NC_000913.3|NC_000913.3" war aber im FASTA nur "NC_000913.3"
                seq = sequence.get(contig_key)
                L = len(seq)
                abs_pos = ((pos + offsets - 1) % L) + 1

                window_characters = [seq[(pos + off - 1) % L] for off in offsets]
                window = "".join(window_characters)
                base_series = window
                cov_list = [coverage_plus.get(int(p), 0.0) for p in abs_pos]
                plus_record.append({
                        "contigPos": pos,
                        "contigID": contig,
                        "windowBegin": int(abs_pos[0]),
                        "windowEnd": int(abs_pos[-1]),
                        "Genome": row["Genome"],
                        "SuperStrand": row["SuperStrand"],
                        "Bases": base_series,
                        "coverage": pp.normalize(np.array(cov_list).reshape(1, -1), norm="max").tolist()
                    })
            plus_df = pd.DataFrame(plus_record, columns=["contigPos", "contigID", "windowBegin", "windowEnd", "Genome", "SuperStrand", "Bases", "coverage"])
            plus_path = os.path.join(output_dir + f"/{cond}/plus/", f"[Window.{cond}]_plusStrand.tsv")
            plus_df.to_csv(plus_path, sep="\t", index=False)
            print(f"[Window.{cond}] Saved plus strand table to {plus_path}")



            print(f"[Window.{cond}] Creating negative TSS windows...")
            neg_plus_records = []
            for (contig, sstrand), grp in plus_sub.groupby(["contigID", "SuperStrand"]):
                contig_key = str(contig.split("|")[0].split()[0])
                seq = sequence.get(contig_key)
                if seq == None:
                    print("WARNING")
                    continue
                L = len(seq)
                tss_positions = sorted(set(int(p) for p in grp["contigPos"].tolist()))
                if len(tss_positions) < 2:
                    continue
                pairs = [(tss_positions[i], tss_positions[(i+1) % len(tss_positions)]) for i in range(len(tss_positions))]
                for first, second in pairs:
                    d_forward = (second - first) if second >= first else (L - first + second)
                    mid = ((first - 1 + d_forward // 2) % L) + 1

                    abs_pos = ((mid + offsets - 1) % L) + 1
                    window_chars = [seq[(mid + off - 1) % L] for off in offsets]
                    window = "".join(window_chars)
                    base_series = reverse_complement(window)

                    coverage = [coverage_plus.get(int(p), 0.0) for p in abs_pos]

                    neg_plus_records.append({
                        "contigPos": int(mid),
                        "contigID": contig,
                        "windowBegin": int(abs_pos[0]),
                        "windowEnd": int(abs_pos[-1]),
                        "Genome": cond,
                        "SuperStrand": sstrand,
                        "Bases": base_series,
                        "coverage": pp.normalize(np.array(coverage).reshape(1, -1), norm="max").tolist()
                    })
            if neg_plus_records:
                neg_plus_df = pd.DataFrame(
                    neg_plus_records,
                    columns=["contigPos", "contigID", "windowBegin", "windowEnd", "Genome", "SuperStrand", "Bases", "coverage"]
                )
                neg_plus_path = os.path.join(output_dir + f"/{cond}/plus/", f"[Window.{cond}]_plusStrand_negative.tsv")
                neg_plus_df.to_csv(neg_plus_path, sep="\t", index=False)
                print(f"[Window.{cond}] Saved negative plus table to {neg_plus_path}")

                if preview:
                    preview_path = os.path.join(preview_dir + f"/{cond}/plus/", f"[Window.{cond}]_negative_plus_preview.tsv")
                    neg_plus_df.head(window_total_size).to_csv(preview_path, sep="\t", index=False, )
                    print(f"[Window.{cond}] Saved negative plus preview to {preview_path}")



            if preview:
                preview_path = os.path.join(preview_dir + f"/{cond}/plus/", f"[Window.{cond}]_plusStrand_preview.tsv")
                plus_df.head(window_total_size).to_csv(preview_path, sep="\t", index=False)
                print(f"[Window.{cond}] Saved plus strand preview to {preview_path}")
        

    return True

