#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import csv
import os
import subprocess
import time
from pathlib import Path
from Bio import Align, SeqIO
from Bio.Seq import Seq
import pandas as pd
import sys

ref_seq_corr = os.environ["ref_seq_corr"]
reference_set = os.path.expandvars(os.path.expanduser(ref_seq_corr))
REFERENCE = { rec.id: str(rec.seq).upper() for rec in SeqIO.parse(reference_set, "fasta") }

# ----------------------------
# Helpers
# ----------------------------
def _vsearch_best_refs(
    qry_file: str,
    reference_set: str,
    out_tsv: str,
    vsearch_exe: str = "vsearch",
    min_id: float = 0.7,  # newer script used 0.7 to allow more distant matches
    maxaccepts: int = 1,
    threads: int = 1,
    maxhits: int = 1,
) -> dict:
    """
    Runs vsearch --usearch_global and returns {query_id: top_hit_id}.
    """
    cmd = [
        vsearch_exe,
        "--usearch_global", str(qry_file),
        "--db", str(reference_set),
        "--blast6out", str(out_tsv),
        "--id", str(min_id),
        "--strand", "both",
        "--maxaccepts", str(maxaccepts),
        "--maxhits", str(maxhits),
        "--threads", str(threads),
    ]

    proc = subprocess.run(cmd, check=False, capture_output=True, text=True)
    if proc.returncode != 0:
        raise RuntimeError(
            "vsearch failed.\n"
            f"CMD: {' '.join(cmd)}\n"
            f"STDOUT:\n{proc.stdout}\n"
            f"STDERR:\n{proc.stderr}\n"
        )

    vsearch_dict = {}
    out_path = Path(out_tsv)
    if not out_path.exists() or out_path.stat().st_size == 0:
        return vsearch_dict

    with out_path.open("r", newline="") as f:
        reader = csv.reader(f, delimiter="\t")
        for row in reader:
            if not row:
                continue
            qid, sid = row[0], row[1]
            vsearch_dict[qid] = sid

    return vsearch_dict

def collect_hit_seqs(reference_fasta: str, hit_ids: set) -> dict:
    """
    Returns {hit_id: UPPERCASE_SEQUENCE} for hit_ids present in reference_fasta.
    """
    hit_seqs = {}
    for rec in SeqIO.parse(str(reference_fasta), "fasta"):
        if rec.id in hit_ids:
            hit_seqs[rec.id] = str(rec.seq).upper()
    return hit_seqs

def end_correction(ref_aln: str, qry_aln: str):
    """
    - If reference starts/ends with gaps: trim those positions from BOTH.
    - If query starts/ends with gaps: replace those gaps with Ns.
    Returns:
        ref_fixed, qry_fixed, endfixed_bool, startgapcounter, endgapcounter
    """
    ref_fixed = ref_aln
    qry_fixed = qry_aln
    endfixed = False
    startgapcounter = 0
    endgapcounter = 0

    # Trim leading gaps in reference (query is longer at the start)
    if ref_fixed.startswith("-"):
        i = 0
        while i < len(ref_fixed) and ref_fixed[i] == "-":
            i += 1
        ref_fixed = ref_fixed[i:]
        qry_fixed = qry_fixed[i:]
        startgapcounter = i
        endfixed = True

    # Trim trailing gaps in reference (query is longer at the end)
    if ref_fixed.endswith("-"):
        j = len(ref_fixed) - 1
        while j >= 0 and ref_fixed[j] == "-":
            j -= 1
            endgapcounter += 1
        ref_fixed = ref_fixed[: j + 1]
        qry_fixed = qry_fixed[: j + 1]
        endfixed = True

    # Replace leading gaps in query with Ns (query is missing bases at start)
    if qry_fixed.startswith("-"):
        i = 0
        while i < len(qry_fixed) and qry_fixed[i] == "-":
            i += 1
        qry_fixed = ("n" * i) + qry_fixed[i:]
        endfixed = True

    # Replace trailing gaps in query with Ns (query missing bases at end)
    if qry_fixed.endswith("-"):
        j = len(qry_fixed) - 1
        while j >= 0 and qry_fixed[j] == "-":
            j -= 1
        qry_fixed = qry_fixed[: j + 1] + ("n" * (len(qry_fixed) - (j + 1)))
        endfixed = True

    return ref_fixed, qry_fixed, endfixed, startgapcounter, endgapcounter

def homopolymer_length(seq: str, index: int):

    if index < 0 or index >= len(seq):
        return "", "", "", "", ""

    base = seq[index]

    # nearest non-gap on left
    left_index = index - 1
    while left_index >= 0 and seq[left_index] == "-":
        left_index -= 1
    base_left = seq[left_index] if left_index >= 0 else ""

    # nearest non-gap on right
    last_gap_position = index
    right_index = index + 1
    while right_index < len(seq) and seq[right_index] == "-":
        last_gap_position = right_index
        right_index += 1
    base_right = seq[right_index] if right_index < len(seq) else ""

    # left run
    l = left_index
    left_start_index = index
    while l >= 0 and (seq[l] == base_left or seq[l] == "-" or seq[l] == "N"):
        left_start_index = l
        l -= 1
    homopolymer_left = seq[(l + 1): index].replace("-", "")

    # right run
    r = right_index
    right_end_index = index
    while r < len(seq) and (seq[r] == base_right or seq[r] == "-" or seq[r] == "N"):
        right_end_index = r
        r += 1
    homopolymer_right = seq[right_index:r].replace("-", "")

    if base == "-":
        if base_left == base_right:
            homopolymer = homopolymer_left + homopolymer_right
            placement = "part"
        else:
            if len(homopolymer_left) >= len(homopolymer_right):
                homopolymer = homopolymer_left
                right_end_index = index - 1
            else:
                homopolymer = homopolymer_right
                left_start_index = last_gap_position + 1
            placement = "part"

    elif base == "N":
        if (base_left == base_right) or (base_left == "N") or (base_right == "N"):
            homopolymer = homopolymer_left + base + homopolymer_right
            placement = "part"
        else:
            if len(homopolymer_left) >= len(homopolymer_right):
                homopolymer = homopolymer_left + base
                right_end_index = index
            else:
                homopolymer = base + homopolymer_right
                left_start_index = index
            placement = "part"

    else:
        if base == base_left == base_right:
            homopolymer = homopolymer_left + base + homopolymer_right
            placement = "part"
        elif base == base_left != base_right:
            homopolymer = homopolymer_left + base
            placement = "part"
            right_end_index = index - 1
        elif base == base_right != base_left:
            homopolymer = base + homopolymer_right
            placement = "part"
            left_start_index = last_gap_position + 1
        else:
            if len(homopolymer_left) >= len(homopolymer_right):
                homopolymer = homopolymer_left
                placement = "adjacent-L"
                right_end_index = index - 1
            else:
                homopolymer = homopolymer_right
                placement = "adjacent-R"
                left_start_index = last_gap_position + 1

    return base, placement, homopolymer, left_start_index, right_end_index

def second_tier_check(seq_name: str, qry_aln_edited: str, indels_not_corrected_list: list):

    for run in indels_not_corrected_list:
        query_or_reference = run[0]
        positions = run[1:]

        if len(positions) <= 2:
            if query_or_reference == "r":
                qry_aln_edited = (
                    qry_aln_edited[:positions[0]]
                    + ("-" * len(positions))
                    + qry_aln_edited[(positions[-1] + 1):]
                )
                print(f"Deleted bases from sequence {seq_name} at position(s) {positions} (non-homopolymer)")

            elif query_or_reference == "q":
                qry_aln_edited = (
                    qry_aln_edited[:positions[0]]
                    + ("n" * len(positions))
                    + qry_aln_edited[(positions[-1] + 1):]
                )
                print(f"Added Ns to sequence {seq_name} at position(s) {positions} (non-homopolymer)")

    return qry_aln_edited

def hmm_check(input_sequence: str, trans_table: int):
    """
    Translate frame starting at index 1 (frame 2),
    pad with Ns so frame length is divisible by 3.
    """
    seq_nogaps = input_sequence.replace("-", "").upper()
    frame = seq_nogaps[1:]
    length_issue = False

    while len(frame) % 3 != 0:
        length_issue = True
        frame += "N"

    protein_seq = str(Seq(frame).translate(table=trans_table))
    stop_codon = ("*" in protein_seq)
    return stop_codon, length_issue

def scan_gap_runs(aln_string: str, tag: str):
    """
    Returns list like [['q', i, i+1], ['q', j], ...]
    """
    runs = []
    current = [tag]
    for i, ch in enumerate(aln_string):
        if ch == "-":
            current.append(i)
        else:
            if len(current) > 1:
                runs.append(current)
            current = [tag]
    if len(current) > 1:
        runs.append(current)
    return runs

def safe_align_filename(seq_name: str) -> str:
    """
    Avoid crashing if seq_name lacks pipes.
    """
    parts = seq_name.split("|")
    if len(parts) >= 2:
        return f"{parts[0]}_{parts[1]}_align.fasta"
    return f"{parts[0]}_align.fasta"

# ----------------------------
# Main pipeline
# ----------------------------
def run_autocorrect(sampleid, wd):  
    qry_file = os.path.join(wd, f"{sampleid}_consensus.fasta")
    problem_file=   os.path.join(wd, f"{sampleid}_problem_seqs.fasta")
    output_file_complete = os.path.join(wd, f"{sampleid}_corrected_otus.fas")
    output_file = os.path.join(wd, f"{sampleid}_edited_seqs.fasta")
    uncertain_file = os.path.join(wd, f"{sampleid}_uncertain_edits.fasta")
    hmm_out = os.path.join(wd, f"{sampleid}_hmm_issues.csv")                         
    alignment_directory=os.path.expanduser("sequence_correction_alignment_dir")
    vsearch_ex = os.environ["vsearch_exe"]
    vsearch_exe = os.path.expandvars(os.path.expanduser(vsearch_ex))
    threads = int(os.environ.get("THREADS", "1"))
    start_time = time.perf_counter()
    print(f"Using {threads} threads") 

    aligner = Align.PairwiseAligner()
    aligner.mode = "global"

    aligner.match_score = 2
    aligner.mismatch_score = -4
    aligner.open_internal_gap_score = -20
    aligner.extend_internal_gap_score = -5
    aligner.open_end_deletion_score = 0
    aligner.extend_end_deletion_score = 0
    aligner.open_end_insertion_score = -10
    aligner.extend_end_insertion_score = -2

    qry_path = Path(qry_file)
    #ref_path = Path(reference_set)
    ref_path = reference_set
    
    qry_seq_list = list(SeqIO.parse(str(qry_path), "fasta"))
    output_dict = {}
    uncertain_dict = {}
    problem_dict = {}
    output_dict_complete = {}

    align_dir = Path(alignment_directory)
    align_dir.mkdir(parents=True, exist_ok=True)

    # ==========================================================
    # Run vsearch to get the best batches of reference sequences
    # ==========================================================
    vsearch_out =  os.path.join(wd, f"{sampleid}_vsearch_hit.tsv")
    vsearch_dict = _vsearch_best_refs(
        qry_path, ref_path, vsearch_out, vsearch_exe=vsearch_exe, threads=threads
    )
    top_hits = set(vsearch_dict.values())
    hit_seqs_dict = collect_hit_seqs(ref_path, top_hits)

    output_no_change_count = 0
    excluded_for_length_count = 0
    excluded_for_nonmatch_count = 0
    endfixed_count = 0
    no_changes_except_endfix_count = 0
    autocorrect_count = 0
    non_homopolymer_autocorrect_count = 0
    indel_not_corrected_count = 0
    stop_count = 0
    no_match = 0

    row_list = []

    for qry_record in qry_seq_list:
        seq_name = qry_record.id
        qry_seq = str(qry_record.seq).upper()

        autocorrect_bool_list = []
        indels_not_corrected_list = []
        no_issues = None
        hmm_check_needed = True
        temp_name = None
        qry_aln_second_tier = None
        non_hp_edit = False

        hit_id = vsearch_dict[seq_name] if seq_name in vsearch_dict else None
        ref_seq = hit_seqs_dict[hit_id] if hit_id in hit_seqs_dict else None

        try:
            trans_table = int(hit_id.rsplit("|", 1)[-1])
        except Exception:
            trans_table = 5
            print(f"Warning: could not parse translation table from {hit_id}; using {trans_table}")

        edited_name = seq_name

        # Must have a vsearch match
        if seq_name not in vsearch_dict:
            no_match += 1
            print(f"Sequence {seq_name} - no match found!")
            excluded_for_nonmatch_count += 1

            qry_aln_edited_nogaps = qry_seq
            stop_codon, length_issue = hmm_check(qry_aln_edited_nogaps, trans_table)
            edited_name = f"{seq_name}|no_match"

            problem_dict[edited_name] = qry_aln_edited_nogaps.upper()
            output_dict_complete[seq_name] = qry_aln_edited_nogaps.upper()
            row_list.append([seq_name, stop_codon, length_issue, qry_aln_edited_nogaps.upper()])
            continue

        if hit_id not in hit_seqs_dict:
            print(f"Sequence {seq_name} - hit {hit_id} not found in reference FASTA headers.")
            excluded_for_nonmatch_count += 1

            qry_aln_edited_nogaps = qry_seq
            stop_codon, length_issue = hmm_check(qry_aln_edited_nogaps, trans_table)
            edited_name = f"{seq_name}|no_ref_hit"

            problem_dict[edited_name] = qry_aln_edited_nogaps.upper()
            output_dict_complete[seq_name] = qry_aln_edited_nogaps.upper()
            row_list.append([seq_name, stop_codon, length_issue, qry_aln_edited_nogaps.upper()])
            continue

        alignment = aligner.align(ref_seq, qry_seq)[0]
        ref_aln, qry_aln = alignment

        # Case 1: exact match (no gaps)
        if "-" not in ref_aln and "-" not in qry_aln:
            print(f"Sequence {seq_name} - no gaps in alignment, skipping endfix/homopolymer checks.")
            qry_aln_edited = qry_seq
            qry_aln_edited_nogaps = qry_seq
            qry_aln_edited_for_alignment = qry_seq
            no_issues = True
            output_no_change_count += 1
            startgapcounter = 0
            endgapcounter = 0

        # Case 2: gaps exist -> attempt endfix + simple homopolymer edits
        else:
            print(f"Sequence {seq_name} - gaps detected in alignment, applying endfix and homopolymer checks.")
            ref_fix, qry_fix, endfixed, startgapcounter, endgapcounter = end_correction(ref_aln, qry_aln)
            qry_aln_edited = qry_fix

            if endfixed:
                edited_name = f"{seq_name}|endfixed"
                endfixed_count += 1

            # if any gaps remain, try your homopolymer logic
            if ("-" in ref_fix) or ("-" in qry_fix):
                no_issues = False

                gap_runs_ref = scan_gap_runs(ref_fix, "r")
                gap_runs_qry = scan_gap_runs(qry_fix, "q")
                gap_runs = gap_runs_qry + gap_runs_ref

                # Process right-to-left so index edits don't shift upcoming positions
                gap_runs.sort(key=lambda x: x[1], reverse=True)

                for run in gap_runs:
                    query_or_reference = run[0]
                    positions = run[1:]

                    if len(positions) <= 2:
                        N_index_list = []
                        deletion_index_list = []
                        Ns_deleted = 0

                        # If gap is in the reference and the query has N at that aligned position,
                        # delete N first because it is already ambiguous and safest to remove.
                        if query_or_reference == "r":
                            for gap_index in positions:
                                if qry_aln_edited[gap_index] == "N":
                                    qry_aln_edited = (
                                        qry_aln_edited[:gap_index]
                                        + "-"
                                        + qry_aln_edited[gap_index + 1:]
                                    )
                                    Ns_deleted += 1
                                    N_index_list.append(gap_index)

                            if Ns_deleted > 0:
                                print(f"Deleted ambiguous nucleotides from sequence {seq_name} at position(s) {N_index_list}")
                                autocorrect_bool_list.append(True)

                            if Ns_deleted == len(positions):
                                continue

                        (
                            letter,
                            placement,
                            homopolymer,
                            left_start_index,
                            right_end_index,
                        ) = homopolymer_length(qry_aln_edited, positions[0])

                        if query_or_reference == "r":
                            same_letter = True

                            if len(homopolymer) >= 4:
                                # Newer logic:
                                # If the homopolymer contains Ns, try deleting those before deleting confident bases.
                                if "N" in homopolymer:
                                    updated_hp = homopolymer

                                    for idx, char in enumerate(homopolymer):
                                        if char == "N" and Ns_deleted < len(positions):
                                            updated_hp = updated_hp[:idx] + "-" + updated_hp[idx + 1:]
                                            N_index_list.append(idx + left_start_index)
                                            Ns_deleted += 1

                                    if placement == "part":
                                        if left_start_index < positions[0]:
                                            qry_aln_edited = (
                                                qry_aln_edited[:left_start_index]
                                                + updated_hp
                                                + qry_aln_edited[positions[-1]:]
                                            )
                                        else:
                                            qry_aln_edited = (
                                                qry_aln_edited[:positions[0]]
                                                + updated_hp
                                                + qry_aln_edited[right_end_index + 1:]
                                            )

                                    elif placement == "adjacent-L":
                                        qry_aln_edited = (
                                            qry_aln_edited[:left_start_index]
                                            + updated_hp
                                            + qry_aln_edited[positions[0]:]
                                        )

                                    elif placement == "adjacent-R":
                                        qry_aln_edited = (
                                            qry_aln_edited[:positions[-1] + 1]
                                            + updated_hp
                                            + qry_aln_edited[right_end_index + 1:]
                                        )

                                    print(
                                        f"Deleted ambiguous nucleotides from sequence {seq_name} at position(s) {N_index_list} (homo) {placement}"
                                    )
                                    autocorrect_bool_list.append(True)

                                if Ns_deleted == len(positions):
                                    continue

                                if placement == "part":
                                    for b in positions:
                                        if qry_aln_edited[b] != letter:
                                            same_letter = False

                                    if same_letter:
                                        qry_aln_edited = (
                                            qry_aln_edited[:positions[0]]
                                            + ("-" * len(positions))
                                            + qry_aln_edited[(positions[-1] + 1):]
                                        )
                                        print(f"Deleted bases from sequence {seq_name} at position(s) {positions} {placement}")
                                        autocorrect_bool_list.append(True)
                                    else:
                                        # Newer logic:
                                        # If the gap run letters are not identical, try deleting from the nearby
                                        # homopolymer flank instead of directly at the gap positions.
                                        if left_start_index < positions[0]:
                                            qry_aln_edited = (
                                                qry_aln_edited[:positions[0] - len(positions)]
                                                + ("-" * len(positions))
                                                + qry_aln_edited[positions[0]:]
                                            )
                                            deletion_index_list = [x - len(positions) for x in positions]
                                        else:
                                            qry_aln_edited = (
                                                qry_aln_edited[:positions[-1] + 1]
                                                + ("-" * len(positions))
                                                + qry_aln_edited[positions[-1] + len(positions):]
                                            )
                                            deletion_index_list = [x + len(positions) for x in positions]

                                        print(
                                            f"Deleted bases from sequence {seq_name} at position(s) {deletion_index_list} {placement}"
                                        )
                                        autocorrect_bool_list.append(True)

                                elif placement == "adjacent-L":
                                    for d in positions:
                                        deletion_index_list.append(d - len(positions))

                                    qry_aln_edited = (
                                        qry_aln_edited[:deletion_index_list[0]]
                                        + ("-" * len(positions))
                                        + qry_aln_edited[(deletion_index_list[-1] + 1):]
                                    )
                                    print(
                                        f"Deleted bases from sequence {seq_name} at position(s) {deletion_index_list} {placement}"
                                    )
                                    autocorrect_bool_list.append(True)

                                elif placement == "adjacent-R":
                                    for d in positions:
                                        deletion_index_list.append(d + len(positions))

                                    qry_aln_edited = (
                                        qry_aln_edited[:deletion_index_list[0]]
                                        + ("-" * len(positions))
                                        + qry_aln_edited[(deletion_index_list[-1] + 1):]
                                    )
                                    print(
                                        f"Deleted bases from sequence {seq_name} at position(s) {deletion_index_list} {placement}"
                                    )
                                    autocorrect_bool_list.append(True)

                            else:
                                # Match original behavior: only flag unresolved non-codon indels.
                                # Codon-sized gaps (3, 6, ...) do not shift the reading frame.
                                if len(positions) % 3 != 0:
                                    indels_not_corrected_list.append([query_or_reference, *positions])
                                    autocorrect_bool_list.append(False)

                        elif query_or_reference == "q":
                            # Missing bases in query:
                            # if part of a long homopolymer, fill with Ns.
                            if len(homopolymer) >= 4:
                                qry_aln_edited = (
                                    qry_aln_edited[:positions[0]]
                                    + ("n" * len(positions))
                                    + qry_aln_edited[(positions[-1] + 1):]
                                )
                                print(f"Added Ns to sequence {seq_name} at position(s) {positions}")
                                autocorrect_bool_list.append(True)
                            else:
                                # Match original behavior: only flag unresolved non-codon indels.
                                # Codon-sized gaps (3, 6, ...) do not shift the reading frame.
                                if len(positions) % 3 != 0:
                                    indels_not_corrected_list.append([query_or_reference, *positions])
                                    autocorrect_bool_list.append(False)

                    else:
                        # More than 2 consecutive gaps are not auto-corrected.
                        # Match original behavior: only flag unresolved non-codon indels.
                        if len(positions) % 3 != 0:
                            indels_not_corrected_list.append([query_or_reference, *positions])
                            autocorrect_bool_list.append(False)

            else:
                no_issues = True
                if endfixed:
                    no_changes_except_endfix_count += 1
                else:
                    output_no_change_count += 1

            qry_aln_edited_for_alignment = ("-" * startgapcounter) + qry_aln_edited + ("-" * endgapcounter)
            qry_aln_edited_nogaps = qry_aln_edited.replace("-", "")

        if True in autocorrect_bool_list:
            edited_name = f"{edited_name}|autoedit"
            autocorrect_count += 1

        # if there are unresolved short indels, try a second-tier correction and keep it
        # only if it resolves the coding/length issue safely.
        if False in autocorrect_bool_list:
            # Match original behavior:
            # 1) First test whether unresolved non-codon indels actually cause a coding/length issue.
            # 2) If not, keep the current sequence but label it as containing an unresolved indel.
            # 3) If second-tier alignment-only edits fix the issue, DO NOT put them in the main
            #    corrected FASTA. Store them in uncertain_edits.fasta and mark the main output |check.
            stop_codon, length_issue = hmm_check(qry_aln_edited, trans_table)

            if stop_codon is False and length_issue is False:
                indel_not_corrected_count += 1
                hmm_check_needed = False
                edited_name = f"{edited_name}|indel"

            else:
                qry_aln_second_tier = second_tier_check(
                    seq_name,
                    qry_aln_edited,
                    indels_not_corrected_list,
                )
                stop_codon_second, length_issue_second = hmm_check(qry_aln_second_tier, trans_table)

                if stop_codon_second is True or length_issue_second is True:
                    # Second-tier edit did not solve the issue, so discard it.
                    indel_not_corrected_count += 1
                    edited_name = f"{edited_name}|indel"

                else:
                    # Keep main output conservative; put suggested non-HP edit into uncertain file.
                    non_homopolymer_autocorrect_count += 1
                    non_hp_edit = True
                    temp_name = f"{edited_name}|autoedit_nonHP"
                    edited_name = f"{edited_name}|check"
                    uncertain_dict[temp_name] = qry_aln_second_tier.replace("-", "").upper()

                hmm_check_needed = False

        if hmm_check_needed:
            stop_codon, length_issue = hmm_check(qry_aln_edited, trans_table)

        if stop_codon:
            edited_name = f"{edited_name}|STOP"
            stop_count += 1

        if length_issue:
            edited_name = f"{edited_name}|INDEL"

        if no_issues:
            edited_name = f"{edited_name}|OK"

        # do final length filtering AFTER all correction attempts.
        final_seq_length = len(qry_aln_edited_nogaps)
        if not (640 <= final_seq_length <= 670):
            print(f"Sequence {seq_name} excluded due to length of {final_seq_length} bp.")
            problem_dict[edited_name] = qry_aln_edited_nogaps
            excluded_for_length_count += 1
            output_dict_complete[seq_name] = qry_aln_edited_nogaps.upper()
            row_list.append([seq_name, stop_codon, length_issue, qry_aln_edited_nogaps.upper()])
            continue

        if (no_issues is False) or stop_codon or length_issue:
            align_file_name = safe_align_filename(seq_name)
            out_path = os.path.join(align_dir, align_file_name)
            with open(out_path, "w") as f:
                f.write(f">{hit_id}|REF\n{ref_aln}\n")

                if qry_aln_edited_nogaps == qry_seq and non_hp_edit is False:
                    f.write(f">{edited_name}\n{qry_aln_edited_for_alignment}\n")
                else:
                    f.write(f">{seq_name}|orig\n{qry_aln}\n")
                    f.write(f">{edited_name}\n{qry_aln_edited_for_alignment}\n")

                    if qry_aln_second_tier and temp_name:
                    # If there is a sequence with suggested edits based on the second tier check, add this as a final item
                        qry_aln_second_tier = ('-' * startgapcounter) + qry_aln_second_tier + ('-' * endgapcounter)
                        f.write(f">{temp_name}\n{qry_aln_second_tier}\n")

        output_dict[edited_name] = qry_aln_edited_nogaps
        output_dict_complete[seq_name] = qry_aln_edited_nogaps.upper()
        row_list.append([seq_name, stop_codon, length_issue, qry_aln_edited_nogaps.upper()])

    if no_match > len(qry_seq_list) / 2:
        print(
            "Warning!: a large proportion of sequences had no vsearch match. "
            "Check your reference database and sequence orientation of All.fasta or vsearch parameters."
        )

    # Problem file now writes the actual problematic sequences with edited headers.
    with open(problem_file, "w") as f:
        for header, seq in problem_dict.items():
            f.write(f">{header}\n{seq}\n")
    
    with open(uncertain_file, "w") as f:
        for header, seq in uncertain_dict.items():
            f.write(f">{header}\n{seq}\n")

    with open(output_file, "w") as f:
        for header, seq in output_dict.items():
            f.write(f">{header}\n{seq}\n")

    print(f"\nFinal corrected sequences written to {output_file_complete}")
    print(f"output_dict_complete keys: {list(output_dict_complete.keys())[:10]}")
    print(f"output_dict_complete values: {[seq[:10] for seq in output_dict_complete.values()][:10]}")
    
    with open(output_file_complete, "w") as f:
        for header, seq in output_dict_complete.items():
            f.write(f">{header}\n{seq}\n")

    df = pd.DataFrame(row_list, columns=["SeqName", "Stop_Codon", "Indel_Issue", "Seq"])
    df.to_csv(hmm_out, index=False)

    print("\nDone.")
    print(f"Excluded (length): {excluded_for_length_count}")
    print(f"Excluded (no match / missing ref): {excluded_for_nonmatch_count}")
    print(f"No change: {output_no_change_count}")
    print(f"Endfixed: {endfixed_count}")
    print(f"No changes except endfix: {no_changes_except_endfix_count}")
    print(f"Autocorrected by homopolymer logic: {autocorrect_count}")
    print(f"have corrections suggested based on alignment only (should be checked): {non_homopolymer_autocorrect_count}")
    print(f"Indel not corrected (flagged): {indel_not_corrected_count}")
    print(f"Stop codon flagged: {stop_count}")

    end_time = time.perf_counter()
    print(f"Elapsed time: {end_time - start_time:.4f} seconds")


if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Running sequence correction <sampleid> <working_directory>")
        sys.exit(1)

    sampleid = sys.argv[1]
    wd = sys.argv[2]

    run_autocorrect(sampleid, wd)
