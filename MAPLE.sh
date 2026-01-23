#!/usr/bin/env bash
clear
echo -e '\n\n\n########## STARTING MAPLE ANALYSIS ##########'
# v 0.9

# REQUIREMENTS:
# 1) FASTQ files from sequencer
# 2) Parameters file

#######################################################################
################# STEP 1: Set user-specific arguments #################
#######################################################################

componentreads=1                        # 0 = do not generate ASV component reads; 1 = generate and keep ASV component reads (takes longer)
working_dir="$HOME/Metabarcoding"       # Set the location of your working directory here
reference_lib_dir="$HOME/REFS/MAPLE"    # Set the location of your sintax-foramtted reference library here
scripts_dir="$HOME/SCRIPTS/MAPLE"       # Set the location of your MAPLE scripts (bash + R) here
sintax_cutoff=0.6                       # Set vsearch --sintax cutoff here (0.6 = retain identifications with >= 60% confidence)

#######################################################################
################# STEP 2: Collect run parametrs #######################
#######################################################################
# Detect number of cores
cores=$(getconf _NPROCESSORS_ONLN 2>/dev/null || sysctl -n hw.ncpu)

# Collect and reformat parameters information
cd $working_dir
Rscript "$scripts_dir/1-MAPLE_getparameters.R" "$PWD"

# Remove any Windows carriage returns from mapping file
sed -i.bak 's/\r$//' mapping*.txt && rm -f mapping*.txt.bak

# Extract info from parameters file
read -r runid < <(sed -n '1p' runinfo.txt)
read -r pe_reads < <(sed -n '2p' runinfo.txt)
read -r numreps  < <(sed -n '3p' runinfo.txt)
numreps=${numreps%.*}   # convert to integer
read -r minreads < <(sed -n '4p' runinfo.txt)
minreads=${minreads%.*}   # convert to integer
read -r asv_dist1  < <(sed -n '5p' runinfo.txt)
read -r asv_dist2  < <(sed -n '6p' runinfo.txt)

# Extract forward and reverse UMIs from mapping file
awk -F'\t' -v plate="$runid" '
    BEGIN {OFS="\t"}
    NR > 1 {
      if ($6 != "NA" && $6 != "") umis[$6] = 1   # forward UMI
    }
    END {
      for (u in umis) {
        print ">" u "\n" u
      }
    }' mapping_"$runid".txt > fwd_umis.fasta

awk -F'\t' -v plate="$runid" '
    BEGIN {OFS="\t"}
    NR > 1 {
      if ($7 != "NA" && $7 != "") umis[$7] = 1   # reverse UMI
    }
    END {
      for (u in umis) {
        print ">" u "\n" u
      }
    }' mapping_"$runid".txt > rev_umis.fasta

seqtk seq -r -c rev_umis.fasta > rev_umis_rc.fasta
rm rev_umis.fasta

# Calculate minimum length of fwd and rev UMIs
min_umi_len_fwd=$(awk 'NR>1 {print length($0)}' fwd_umis.fasta | sort -n | head -1)
min_umi_len_rev=$(awk 'NR>1 {print length($0)}' rev_umis_rc.fasta | sort -n | head -1)

#######################################################################
############ STEP 3a: Merge PE reads (Illumina only) ##################
#######################################################################

# Merge paired-end reads (Illumina only)
if [ "$pe_reads" = "Yes" ]; then
    
    # Get read1 and read2 from fastq files
    read1=$(ls *.gz | sed -n 1p)
    read2=$(ls *.gz | sed -n 2p)

    # Merge paired end reads
    echo -e "******** Merging paried-end reads..."
    pear -j $((cores - 1)) -f $read1 -r $read2 -o $runid > log.txt

    # Delete read1 and read2 FASTQ files to save space
    rm $read1 $read2

    # Delete discarded and unassembled paired-end reads
    rm $runid".discarded.fastq" $runid".unassembled.forward.fastq" $runid".unassembled.reverse.fastq"

    mv $runid.assembled.fastq all.fastq
fi

#######################################################################
####### STEP 3b: Concatenate raw fastq file (if not pair-ends) ########
#######################################################################

# Merge FASTQ files into single file
if [ "$pe_reads" = "No" ]; then
    echo -e "******** Merging FASTQ files..."

    # Decompress FASTQ files
    ulimit -n 65536
    find . -name "*.gz" | parallel -j $((cores-1)) gzip -d

    # Concatenate all raw FASTQ files into a single file and delete originals
    find . -name "*.fastq" -print0 | xargs -0 cat > all.fastq2
    find . -name "*.fastq" -print0 | xargs -0 rm
    
    # Remove any extra text from sequence headers
    awk '{if(NR%4==1) sub(/\t.*/, "", $0); print}' all.fastq2 > all.fastq
    rm all.fastq2
fi

#######################################################################
################## STEP 4: Filter and demultiplex #####################
#######################################################################

# Count raw reads
seqkit stats -T all.fastq | awk 'NR==2 {print "Raw Reads\t"$4}' > "${runid}_readcounts.txt"

# Compress raw FASTQ
echo -e "******** Compressing raw reads..."
pigz -p $((cores-1)) all.fastq

# Count reads after removing primer dimers
cutadapt -j $((cores - 1)) -m 100 -o all_filtered.fastq all.fastq.gz
seqkit stats -T all_filtered.fastq | awk 'NR==2 {print "After removing primer dimers\t"$4}' >> "${runid}_readcounts.txt"

# Demultiplex forward UMIs at 5' end 
echo -e "******** Demultiplexing reads..."

if [ "$min_umi_len_fwd" -ge 12 ]; then
    cutadapt -j $((cores-1)) \
        -e 0.125 \
        -O  $(( (min_umi_len_fwd * 75) / 100 )) \
        -g file:fwd_umis.fasta \
        -o "{name}.fwd5" \
        -m 100 \
        all_filtered.fastq
else
    cutadapt -j $((cores-1)) \
        -e 0.0 \
        -O  $min_umi_len_fwd \
        -g file:fwd_umis.fasta \
        -o "{name}.fwd5" \
        -m 100 \
        all_filtered.fastq
fi

rm all_filtered.fastq

# Demultiplex fwd UMIs at 3' end
if [ "$pe_reads" = "No" ]; then
    
    # Rev-comp reads that don't have 5'UMI
    seqtk seq -r unknown.fwd5 > unknown.fwd5.rc
    rm unknown.fwd5

    # Search for all fwd UMIs at new 5' end
    if [ "$min_umi_len_fwd" -ge 12 ]; then
        cutadapt -j $((cores-1)) \
            -e 0.25 \
            -O $(( (min_umi_len_fwd * 75) / 100 )) \
            -g file:fwd_umis.fasta \
            -o "{name}.fwd3" \
            -m 100 \
            --discard-untrimmed \
            unknown.fwd5.rc
    else
        cutadapt -j $((cores-1)) \
            -e 0.0 \
            -O $min_umi_len_fwd \
            -g file:fwd_umis.fasta \
            -o "{name}.fwd3" \
            -m 100 \
            --discard-untrimmed \
            unknown.fwd5.rc
    fi
    rm unknown.fwd5.rc

    # Merge fwd5 and fwd3 files into a single fasta file
    for fwd in *.fwd5; do
        sample="${fwd%.fwd5}"
        rev="${sample}.fwd3"
        
        if [[ -f "$rev" ]]; then
            cat "$fwd" "$rev" > "${sample}.fwd53"
        else
            echo "Warning: $rev not found, skipping $sample"
        fi
    done

    rm *.fwd5 *.fwd3
    rename 's/.fwd53/.fwd5/g' *.fwd53
fi

# Get count of reads with forward UMI
sum=0
for k in *.fwd5; do
    sum=$((sum + $(echo $(cat "$k" | wc -l)/4 | bc)))
done
echo $sum | awk '{print "With forward UMI" "\t" $1}' >> "${runid}_readcounts.txt"

# Demultiplex reverse UMIs at 3' end
for f in *.fwd5; do
    base=${f%.fwd5}

    # Extract the relevant 3′ UMIs for this 5′ UMI
    awk -F'\t' -v umi="$base" 'BEGIN{FS="\t"} NR>1 && $7!="" && $6==umi {print ">"$7"\n"$7}' mapping_"$runid".txt > rev_umis_temp.fasta

    # Convert to reverse complement
    seqkit seq -r -p rev_umis_temp.fasta > rev_umis_rc_temp.fasta

    if [ "$min_umi_len_rev" -ge 12 ]; then
        cutadapt -j $((cores-1)) \
            -e 0.25 \
            -O $(( (min_umi_len_rev * 75) / 100 )) \
            -a file:rev_umis_rc_temp.fasta \
            -o "${base}-{name}.fastq" \
            -m 100 \
            --discard-untrimmed \
            "$f"
    else
        cutadapt -j $((cores-1)) \
            -e 0.0 \
            -O $min_umi_len_rev \
            -a file:rev_umis_rc_temp.fasta \
            -o "${base}-{name}.fastq" \
            -m 100 \
            --discard-untrimmed \
            "$f"
    fi
done

rm fwd_umis.fasta rev_umis_rc.fasta *.fwd5 rev_umis_rc_temp.fasta rev_umis_temp.fasta

# Remove any FASTQ files with 0 reads
find . -type f -name "*.fastq" -size 0 -delete

# Process the mapping file and rename FASTQ files
sanitized_name_index="name_map.tsv" > "$sanitized_name_index"
tail -n +2 "mapping_${runid}.txt" | while IFS=$'\t' read -r plt well samp fwdp revp fwdu revu marker minlen maxlen amplen _rest; do
    # Original filename based on UMI mapping
    old="$(echo "${fwdu}-${revu}.fastq" | tr -d '[:space:]')"

    # New filename: append marker and amplicon length to make unique
    new="${samp}_${marker}_${amplen}bp.fastq.tmp"
    new_samp_only="${samp}"

    if [[ -f "$old" ]]; then
        mv -- "$old" "$new"
    fi

    # Create sanitized filename for mapping: replace # with -
    sanitized="${new_samp_only//#/-}"

    # Append to name_map.tsv: first column is samp, second is sanitized filename
    echo -e "${samp}\t${sanitized}" >> "$sanitized_name_index"
done


# Remove any leftover FASTA files that were not renamed
[ "$(find . -maxdepth 1 -name '*.fastq' | wc -l)" -gt 0 ] && rm -- *.fastq

# Remove the temporary ".tmp" suffix (will keep marker_amp in filename for uniqueness)
rename 's/\.tmp$//' *.tmp

# Get read count of demultiplexed reads
sum=0
for k in *.fastq; do
    sum=$((sum + $(echo $(cat "$k" | wc -l)/4 | bc)))
done
echo $sum | awk '{print "With forward and reverse UMI" "\t" $1}' >> "${runid}_readcounts.txt"

# Archive demultiplexed FASTQ files
mkdir -m 777 Individual_Raw_Fastq_Files
cp *.fastq Individual_Raw_Fastq_Files
tar -czf Individual_Raw_Fastq_Files.tar.gz Individual_Raw_Fastq_Files && rm -rf Individual_Raw_Fastq_Files

# Convert to FASTA
for f in *.fastq;do paste - - - - < $f | cut -f 1,2 | sed 's/^@/>/' | tr "\t" "\n" > $f.fasta;done&
wait
rm *.fastq
rename 's/.fastq//g' *.fastq.fasta

# Rename read headers to shorter versions (if spaces in raw read headers)
for f in *bp.fasta;do cut -d" " -f1 $f > $f.fas;done&
wait
rm *bp.fasta
rename 's/.fasta.fas/.fasta/g' *.fasta.fas

# Generate heat map of demultiplexed reads by well per plate
echo -e "******** Generating demultiplexing report..."
Rscript "$scripts_dir/2-MAPLE_heatmap.R" "$PWD" "mapping_${runid}.txt" "$runid"

# Organize FASTA files by marker/amplicon length
tail -n +2 "mapping_${runid}.txt" |
while IFS=$'\t' read -r plt well samp fwdp revp fwdu revu marker min max amplen _rest; do
    fasta_file="${samp}_${marker}_${amplen}bp.fasta"
    dest="${marker}_${amplen}bp"

    if [[ -f "$fasta_file" ]]; then
        mkdir -p "$dest"
        mv "$fasta_file" "$dest"/
    else
        echo "Warning: $fasta_file does not exist and likely failed to yield any reads"
    fi
done

# Replace '#' in filenames with '-' to avoid issues, and store original and sanitized names in an index file
# Rename files: replace # with - in filenames
find . -depth -name '*#*' | while read -r f; do
    newname="${f//#/-}"   # replace all # with -
    mv "$f" "$newname"
done


#######################################################################
################## STEP 5: Process each sample ########################
#######################################################################

# Function to process each sample
task1() (
    # Arguments
    local fasta="$1"
    local runid="$2"
    local cores="$3"
    local pe_reads="$4"
    local minreads="$5"
    local user="$6"
    local marker="$7"
    local ampsize="$8"
    local sampleid="$(basename "$fasta" .fasta)"
    local samplemap=$(awk -F'\t' -v f="$sampleid" '$2==f {print $1; exit}' "../$sanitized_name_index")

    # Get primer information
    local tmp_fwd="temp_fwd_primers.fa"
    local tmp_rev="temp_rev_primers.fa"
    rm -f "$tmp_fwd" "$tmp_rev" "$tmp_fwd".rc "$tmp_rev".rc

    # Create a temporary reduced mapping file for this marker/ampsize
    local tmp_mapping="temp_mapping_${runid}_${marker}_${ampsize}.txt"
    awk -F'\t' -v m="$marker" -v a="$ampsize" '$8 == m && $11 == a' "../mapping_${runid}.txt" > "$tmp_mapping"

    # Forward primer info
    local fprimers
    fprimers=$(awk -F'\t' -v s="$samplemap" '$3 == s {print $4}' "$tmp_mapping")
    IFS=',' read -ra primers_array <<< "$fprimers"
    for idx in "${!primers_array[@]}"; do
        local primer="${primers_array[$idx]//[[:space:]]/}"
        printf ">fprimer%02d\n%s\n" $((idx+1)) "$primer" >> "$tmp_fwd"
    done
    seqtk seq -r "$tmp_fwd" > "${tmp_fwd}.rc"

    local minlen
    minlen=$(awk 'BEGIN{min=1e9} /^>/ {next} {if(length($0)<min) min=length($0)} END{print min}' "$tmp_fwd")
    local O_arg_f=$(( minlen * 75 / 100 ))

    # Reverse primers info
    local rprimers
    rprimers=$(awk -F'\t' -v s="$samplemap" '$3==s {print $5}' "$tmp_mapping")
    IFS=',' read -ra rprimers_array <<< "$rprimers"
    for idx in "${!rprimers_array[@]}"; do
        local primer="${rprimers_array[$idx]//[[:space:]]/}"
        printf ">rprimer%02d\n%s\n" $((idx+1)) "$primer" >> "$tmp_rev"
    done
    seqtk seq -r "$tmp_rev" > "${tmp_rev}.rc"

    minlen=$(awk 'BEGIN{min=1e9} /^>/ {next} {if(length($0)<min) min=length($0)} END{print min}' "$tmp_rev")
    local O_arg_r=$(( minlen * 75 / 100 ))

    # Size & marker info
    local minreadlen maxreadlen ampsize marker
    minreadlen=$(awk -F'\t' -v s="$samplemap" '$3==s {print $9}' "$tmp_mapping")
    maxreadlen=$(awk -F'\t' -v s="$samplemap" '$3==s {print $10}' "$tmp_mapping")
    ampsize=$(awk -F'\t' -v s="$samplemap" '$3==s {print $11}' "$tmp_mapping")
    marker=$(awk -F'\t' -v s="$samplemap" '$3==s {print $8}' "$tmp_mapping")

    # Primer trimming
    cutadapt -j $((cores-1)) -e 0.2 -O $O_arg_f -g file:"$tmp_fwd" \
        --untrimmed-output "${sampleid}_fwd5_noprimer.fas" \
        -o "${sampleid}_fwd5_withprimer.fas" "$fasta"

    cutadapt -j $((cores-1)) -e 0.2 -O $O_arg_r -a file:"${tmp_rev}.rc" \
        --discard-untrimmed \
        -o "${sampleid}_rev3_withprimer.fas" "${sampleid}_fwd5_withprimer.fas"

    cutadapt -j $((cores-1)) -e 0.2 -O $O_arg_r -g file:"$tmp_rev" \
        --discard-untrimmed \
        -o "${sampleid}_rev5_withprimer.fas" "${sampleid}_fwd5_noprimer.fas"

    cutadapt -j $((cores-1)) -e 0.2 -O $O_arg_f -a file:"${tmp_fwd}.rc" \
        --discard-untrimmed \
        -o "${sampleid}_fwd3_withprimer.fas" "${sampleid}_rev5_withprimer.fas"

    awk 'NR%4==2 && /^$/ {$0="N"}
        NR%4==0 && /^$/ {$0="S"}
        {print}' "${sampleid}_fwd3_withprimer.fas" > "${sampleid}_fwd3_withprimer_noempty.fas"

    seqtk seq -r "${sampleid}_fwd3_withprimer_noempty.fas" > "${sampleid}_fwd3_withprimer.rc.fas"

    cat "${sampleid}_rev3_withprimer.fas" "${sampleid}_fwd3_withprimer.rc.fas" > "${sampleid}.fast"

    rm "$fasta" "${sampleid}_fwd5_noprimer.fas" "${sampleid}_fwd5_withprimer.fas" \
        "${sampleid}_rev5_withprimer.fas" "${sampleid}_fwd3_withprimer.fas" \
        "${sampleid}_fwd3_withprimer.rc.fas" "${sampleid}_rev3_withprimer.fas" \
        "${sampleid}_fwd3_withprimer_noempty.fas"

    # Size filter again
    cutadapt -j $((cores-1)) -m $minreadlen -M $maxreadlen -o "$sampleid".fasta "$sampleid".fast
    rm "$sampleid".fast

    # Cluster into ASVs
    echo -e "******** ASV clustering $sampleid..."
    if [ "$componentreads" -eq 0 ]; then
        vsearch --cluster_fast "$sampleid".fasta \
        --id "$(awk -v d="$asv_dist1" 'BEGIN { printf "%.6f", (100 - d)/100 }')" \
        --consout "${sampleid}_consensus.fasta" \
        --iddef 3 \
        --threads $((cores-1))
    else
        vsearch --cluster_fast "$sampleid".fasta \
        --id "$(awk -v d="$asv_dist1" 'BEGIN { printf "%.6f", (100 - d)/100 }')" \
        --consout "${sampleid}_consensus.fasta" \
        --clusters "${sampleid}|ASV" \
        --iddef 3 \
        --threads $((cores-1))
    fi

    # Delete pre-clustered file
    rm "$sampleid".fasta

    # Rename ASV consensus sequence headers to correct ASV name
    consensus_fasta="${sampleid}_consensus.fasta"
    renamed_fasta="${sampleid}_consensus_renamed.fasta"

    awk -v sampleid="$sampleid" '
        BEGIN { asv=0 }
        /^>/ {
            reads=1
            # extract number after ";seqs=" if present
            if ($0 ~ /;seqs=[0-9]+/) {
                # split on ";seqs=" and take the number part
                split($0, parts, ";seqs=")
                reads=parts[2] + 0   # convert to number
            }
            print ">" sampleid "|ASV" asv "|reads-" reads
            asv++
            next
        }
        { print }
        ' "$consensus_fasta" > "$renamed_fasta"
    mv "$renamed_fasta" "$consensus_fasta"
    
    # Filter ASV sequences with fewer than minreads (in both consensus and component read files if necessary)
    echo "******** Filtering ASVs below $minreads reads..."
    consensus_fasta="${sampleid}_consensus.fasta"
    filtered_fasta="${sampleid}_consensus_filtered.fasta"

    # Export variables for parallel jobs
    export minreads cores

    # Function to process one chunk
    process_chunk() {
        local chunk="$1"
        local outfile="$2"

        local tmp_out="${outfile}_$(basename "$chunk").tmp"

        echo "  [START] Processing chunk $(basename "$chunk")..."

        awk -v min="$minreads" -v out="$tmp_out" '
            BEGIN { keep=0; written=0 }
            /^>/ {
                reads = 0
                pos = index($0, "|reads-")
                if (pos > 0) {
                    reads = substr($0, pos + 7) + 0
                }
                keep = (reads >= min)
            }
            { if (keep) { print > out; written++ } }
            END { printf "[DONE] %s: wrote %d sequences\n", FILENAME, written > "/dev/stderr" }
        ' "$chunk"
    }



    export -f process_chunk

    # Function to process the full FASTA in parallel chunks
    process_fasta() {
        local fasta="$1"
        local filtered="$2"

        echo "Processing $fasta..."

        # Count total sequences
        local num_headers
        num_headers=$(grep -c '^>' "$fasta")
        echo "  $num_headers sequences"

        # Determine number of chunks
        local n_chunks=$((cores - 2))
        n_chunks=$(( n_chunks>0 ? n_chunks : 1 ))

        # Sequences per chunk
        local seqs_per_chunk=$(( (num_headers + n_chunks - 1) / n_chunks ))

        # Temporary directory for chunks
        local tmp_dir
        tmp_dir=$(mktemp -d)

        # Split FASTA into chunks
        awk -v seqs="$seqs_per_chunk" -v out="$tmp_dir/chunk_" '
            BEGIN {c=0; n=0; file=sprintf("%s%03d.fasta", out, c)}
            /^>/ {
                if(n++ % seqs == 0 && n>1){c++; file=sprintf("%s%03d.fasta", out, c)}
            }
            {print > file}
        ' "$fasta"

        # Run chunks in parallel
        local chunk_files=( "$tmp_dir"/chunk_*.fasta )
        printf "%s\n" "${chunk_files[@]}" | parallel --bar -j "$n_chunks" process_chunk {} "$filtered"

        # Combine filtered sequences
        cat "${filtered}"_chunk_*.tmp > "$filtered"

        # Cleanup
        rm -rf "$tmp_dir" "${filtered}"_chunk_*.tmp
    }

    export -f process_fasta

    # Run the parallel filtering
    process_fasta "$consensus_fasta" "$filtered_fasta"

    # Filter ASV component reads if necessary
    if [ "$componentreads" -eq 1 ]; then

        export minreads sampleid

        process_asv() {
            local f="$1"

            # skip if file doesn't exist (race conditions)
            [ -e "$f" ] || return

            # count sequences
            local seq_count
            seq_count=$(grep -c "^>" "$f")

            if (( seq_count < minreads )); then
                echo "Skipping $f (only $seq_count sequences < minreads=$minreads)"
                rm -f -- "$f"
            else
                local newname="${f}_ComponentReads.fasta"
                echo "Keeping $f → $newname"
                cp -- "$f" "$newname"
                rm -f -- "$f"
            fi
        }

        export -f process_asv

        echo "Finding ASV component files for sample '$sampleid'..."

        # Stream into GNU Parallel safely (handles unlimited files)
        find . -maxdepth 1 -type f -name "${sampleid}|ASV*" -print0 \
            | parallel -0 --bar -j $((cores - 2)) process_asv {}
    fi

    # Replace consensus file with filtered one
    mv "${sampleid}_consensus_filtered.fasta" "${sampleid}_consensus.fasta"
    
    # Auto-trim if necessary
    if [[ "$marker" == "COI-5P" && "$ampsize" -eq 500 ]]; then
        # Auto-trim sequences
        echo -e "******** Auto-trimming $sampleid ASV consensus sequences..."
        Rscript "$scripts_dir/3-MAPLE_autotrim.R" "$sampleid" "$PWD"
        rm "${sampleid}_consensus.fasta"
    else
        # Skip auto-trimming but rename file
        mv "${sampleid}_consensus.fasta" "${sampleid}_finalasvs.fas"
    fi

    # Remove sequences that fall outside of expected size range
    echo -e "******** Size filtering $sampleid ASV consensus sequences..."
    cutadapt \
      -j $((cores-1)) \
      -m $((ampsize-30)) \
      -M $((ampsize+30)) \
      -o "${sampleid}_FinalASVs.fast" \
      --too-short-output "${sampleid}_finalasvs_tooshort.fast" \
      --too-long-output "${sampleid}_finalasvs_toolong.fast" \
      "${sampleid}_finalasvs.fas"

    rm "$sampleid"_finalasvs.fas

     if [ "$componentreads" -eq 1 ]; then
        for f in "${sampleid}"_finalasvs_tooshort.fast "${sampleid}"_finalasvs_toolong.fast; do
            if [[ -s "$f" ]]; then
                echo "Processing $f ..."

                # Extract sampleid|ASV# from headers and deduplicate
                grep "^>" "$f" | sed -E 's/^>//; s/\|reads-[0-9]+$//' | sort -u |
                while read -r prefix; do
                    comp="${prefix}_ComponentReads.fasta"
                    if [[ -e "$comp" ]]; then
                        echo "Deleting $comp (discarded $prefix)"
                        rm -f "$comp"
                    else
                        echo "Warning: $comp not found"
                    fi
                done
            else
                echo "Skipping $f (missing or empty)"
            fi
        done
    fi

    rm "${sampleid}_finalasvs_tooshort.fast" "${sampleid}_finalasvs_toolong.fast"

    # Remove ASVs that are obvious NUMTs and try to correct polymerase indels (only if >= 500 bp COI-5P)
    if [[ "$marker" == "COI-5P" && "$ampsize" -eq 500 ]]; then
        echo -e "******** Removing NUMTs for $sampleid..."
        Rscript "$scripts_dir/4-MAPLE_numt_screen_indel_correction.R" "$sampleid" $cores $ampsize "$PWD"
        rm "$sampleid"_FinalASVs.fast

         if [ "$componentreads" -eq 1 ]; then
            # Extract the remaining ASVs (without the reads count)
            grep "^>" "${sampleid}_FinalASVs.fasta" | sed -E 's/^>//; s/\|reads-[0-9]+$//' | sort -u > kept_asvs.txt

            # Loop over all ComponentReads.fasta files
            for comp in "${sampleid}"*ComponentReads.fasta; do
                # Extract prefix
                prefix=$(basename "$comp" | sed -E 's/_ComponentReads\.fasta$//')

                # Check if that prefix is in the kept ASV list
                if grep -Fxq "$prefix" kept_asvs.txt; then
                    echo "Keeping $comp (found $prefix in kept_asvs.txt)"
                else
                    echo "Deleting $comp (not found in kept_asvs.txt)"
                    rm -f "$comp"
                fi
            done

            # Cleanup
            rm -f kept_asvs.txt

        fi
    else
        mv "$sampleid"_FinalASVs.fast "$sampleid"_FinalASVs.fasta
    fi

    # Convert final ASV FASTA to single-line
    awk '{if(NR==1) {print $0} else {if($0 ~ /^>/) {print "\n"$0} else {printf $0}}}' "${sampleid}_FinalASVs.fasta" > "${sampleid}_FinalASVs.fasta2"
    rm "${sampleid}_FinalASVs.fasta"
    mv "${sampleid}_FinalASVs.fasta2" "${sampleid}_FinalASVs.fasta"

    # Delete temp primer files
    rm -f "$tmp_fwd" "$tmp_rev" "$tmp_fwd".rc "$tmp_rev".rc

    if [ "$componentreads" -eq 1 ]; then
        # Move ASV component read files into their own folder
        mv ./*_ComponentReads.fasta ./"$runid"_"$marker"_"$ampsize"_ASV_Component_Reads/
    fi
)

# Main loop: process per marker
for marker_dir in */; do
     # Skip unwanted directory
    [[ "$marker_dir" == "Individual_Raw_Fasta_Files/" ]] && continue

    echo "=== Processing marker folder: ${marker_dir%/} ==="
    
    # Move into the marker folder, initialize some variables, and make ASV component read folder
    cd "$marker_dir" || continue
    dirname=$(basename "$PWD")
    marker="${dirname%_*}"
    ampsize="${dirname##*_}"
    ampsize="${ampsize%bp}" 
    
    # Strip the marker_amp tags from filenames, then lookup the refernece library for this marker_ampsize comination
    for f in *.fasta; do
        newname="$(echo "$f" | sed -E "s/_${marker}_${ampsize}bp\.fasta$/.fasta/")"
        mv "$f" "$newname"
    done

    reflib="$(awk -F'\t' -v m="$marker" -v a="$ampsize" 'NR>1 && $8==m && $11==a {print $12; exit}' "../mapping_${runid}.txt")"
    db_path=$(ls -d $reference_lib_dir/"$reflib"*.fasta)

    # Create directory to store ASV component reads if necessary
    if [ "$componentreads" -eq 1 ]; then
        mkdir -m 777 "$runid"_"$marker"_"$ampsize"_ASV_Component_Reads
    fi

    # Process all FASTA files in this folder
    for fasta in *.fasta; do
        [ -e "$fasta" ] || continue
        task1 "$(realpath "$fasta")" "$runid" "$cores" "$pe_reads" "$minreads" "$user" "$marker" "$ampsize"
    done

    # Delete files with 0 bytes
    find . -type f -name "*.fasta" -size 0 -exec rm -f {} +

    if [ "$componentreads" -eq 1 ]; then
        # Rename ASV component read files
        find ./"$runid"_"$marker"_"$ampsize"_ASV_Component_Reads/ -type f -name '*_Reads.fasta' -exec rename 's/[=]/-/g' {} +
        find ./"$runid"_"$marker"_"$ampsize"_ASV_Component_Reads/ -type f -name '*_Reads.fasta' -exec rename 's/[|]/_/g' {} +
        find ./"$runid"_"$marker"_"$ampsize"_ASV_Component_Reads/ -type f -name '*_Reads.fasta' -exec rename 's/_Reads.fasta/.fasta/g' {} +
    fi

    #######################################################################
    ################## STEP 6: Generate run-wide ASVs #####################
    #######################################################################
    # Merge all ASVs into a single master file
    awk 'NF{if(/^>/){if(seen++)printf "\n"; print $0} else print $0}' *FinalASVs.fasta > all_asvs.fasta
    rm *FinalASVs.fasta
    
    # Change ASV sequence names from |reads-n to ;size=n for subsequent size sorting
    sed '/^>/ s/|reads-/;size=/g' all_asvs.fasta > all_asvs.fasta2
    rm all_asvs.fasta

    # Cluster into run-wide ASVs
    vsearch --cluster_size all_asvs.fasta2 \
        --id "$(awk -v d="$asv_dist2" 'BEGIN { printf "%.6f", (100 - d)/100 }')" \
        --consout all_asvs_consensus.fasta \
        --iddef 3 \
        --threads $((cores-1)) \
        --uc all_asvs_cluster_info.uc

    # Reformat run-wide ASV info
    Rscript "$scripts_dir/5-MAPLE_runwideasvformat.R" "$PWD"
    rm all_asvs_cluster_info.uc all_asvs.fasta2

    # Convert run-wide ASVs to single-line FASTA
    awk '{if(NR==1) {print $0} else {if($0 ~ /^>/) {print "\n"$0} else {printf $0}}}' all_asvs_consensus.fasta > all_asvs_consensus.fasta2
    rm all_asvs_consensus.fasta
    mv all_asvs_consensus.fasta2 all_asvs_consensus.fasta

    #######################################################################
    ################## STEP 7: Identify run-wide ASVs #####################
    #######################################################################

    # Indetify ASVs using sintax
    vsearch --sintax all_asvs_consensus.fasta \
        -db $db_path \
        -tabbedout temp_sintax_output.txt \
        -strand plus \
        -sintax_cutoff $sintax_cutoff \
        -threads $((cores-2))

    #######################################################################
    ######################## STEP 8: Denoise ASVs #########################
    #######################################################################

    # De-noise data
    echo -e "******** Denoising data..."
    Rscript "$scripts_dir/6-MAPLE_denoise.R" $minreads $numreps "$PWD" "$runid" "$marker" $ampsize
       
    # Chimera screen using UCHIME (only if ampsize is >= 200 bp)
    if [[ "$ampsize" -ge 200 ]]; then
        echo -e "******** Performing chimera screen..."
        vsearch --uchime_denovo temp.fasta \
            --abskew 10 \
            --chimeras chimeras.fas \
            --fasta_width 0 \
            --mindiv 0.0005
    else
        > chimeras.fas
    fi

    echo -e "******** Finalizing data..."
    Rscript "$scripts_dir/7-MAPLE_final_results.R" "$runid" "$marker" $minreads $numreps $ampsize "$PWD" "$reflib" $componentreads
    rm all_asvs_consensus.fasta all_asvs_run_wide_ASV_info.txt chimeras.fas temp.fasta temp_sintax_output.txt df.txt

    if [ "$componentreads" -eq 1 ]; then
        # Compress ASV component read files
        zip -r -0 "$runid"_"$marker"_"$ampsize"_ASV_Component_Reads.zip "$runid"_"$marker"_"$ampsize"_ASV_Component_Reads && rm -rf "$runid"_"$marker"_"$ampsize"_ASV_Component_Reads
    fi

    #######################################################################
    ############ STEP 9: Perform BIN analysis (if applicable) #############
    #######################################################################

    if [[ "$marker" == "COI-5P" && "$ampsize" -ge 300 ]]; then
        
        vsearch_db=$(ls -d $reference_lib_dir/"$reflib"*vsearch)
        
        if [ ! -f "$vsearch_db" ]; then
          # make vsearch reference library from sintax
          echo -e "******** Creating vsearch reference library..."
          cp $db_path $reference_lib_dir/temp_sintax_no_spaces.fasta
          sed -i '' '/^>/ s/ /_/g' "$reference_lib_dir/temp_sintax_no_spaces.fasta"
          vsearch --makeudb_usearch $reference_lib_dir/temp_sintax_no_spaces.fasta --output $reference_lib_dir/"$reflib".vsearch
          vsearch_db=$(ls -d $reference_lib_dir/"$reflib"*vsearch)
          rm $reference_lib_dir/temp_sintax_no_spaces.fasta
        fi

        echo -e "******** Performing BIN match..."
        # Convert FASTA to single-line
        awk '{if(NR==1) {print $0} else {if($0 ~ /^>/) {print "\n"$0} else {printf $0}}}' temp_toBIN.fasta > single_line.fasta
        rm temp_toBIN.fasta
        mv single_line.fasta temp_toBIN.fasta

        # Identify sequences using VSEARCH
        vsearch --usearch_global temp_toBIN.fasta \
            --db  $vsearch_db \
            --blast6out temp_vsearch_output.txt \
            --id 0.75 \
            --maxhits 3 \
            --maxaccepts 3 \
            --threads $((cores-1))

        # Remove low-level hits from VSEARCH results
        min_overlap=$(printf "%.0f" "$(echo "$ampsize * 0.75" | bc -l)")
        awk -v min_overlap="$min_overlap" '
            BEGIN { OFS="\t" }
            $4 >= min_overlap { print }
            ' temp_vsearch_output.txt | sort -k1,1 -k3nr -k4nr | awk '!seen[$1]++' > filtered_vsearch_output.txt
        
        # Parse hit column into ProcessID and BIN;tax=
        awk -F'\t' 'BEGIN { OFS="\t" } $2 ~ /\|/ { split($2, hit_parts, "|"); print $1, hit_parts[1], hit_parts[2], $3, $4 }' filtered_vsearch_output.txt > parsed_vsearch_output.txt
       
        # Parse BIN;tax= column into BIN and tax=...
        awk -F'\t' '
        {
            # Split column 3 on semicolon
            split($3, parts, ";")

            bold = parts[1]              # Keep BIN name
            taxpart = ""

            # Find the part beginning with tax=
            for (i in parts) {
                if (parts[i] ~ /^tax=/) taxpart = parts[i]
            }

            gsub(/^tax=/, "", taxpart)   # Remove "tax="

            # Split into individual ranks
            n = split(taxpart, ranks, ",")

            # Extract only the taxon names from each rank
            for (i = 1; i <= n; i++) {
                sub(/^[a-z]:/, "", ranks[i])   # Remove sintax rank prefixes (e.g., k:, p:, etc)
            }

            # Print: original col1/col2, BOLD code, then taxonomy fields, then col4/col5
            printf "%s\t%s\t%s\t", $1, $2, bold

            for (i = 1; i <= n; i++) {
                printf "%s\t", ranks[i]
            }

            printf "%s\t%s\n", $4, $5
        }
        ' parsed_vsearch_output.txt > final_output.txt

        # Add BIN_MATCH column
        awk -F'\t' 'BEGIN { OFS="\t" }
        {
            # Check if column 11 (% ID) >= 97.7
            if ($11 >= 97.7) {
                bin_match = "BIN_MATCH"
            } else {
                bin_match = "NO_MATCH"
            }

            # Add new column at the end
            $(NF+1) = bin_match

            print $0
        }' final_output.txt > final_output_with_bin_match.txt

        # Add BIN match results to "by sample" results tab in Excel output
        Rscript "$scripts_dir/8-MAPLE_binmatch.R" "$runid" "$PWD" "$marker" $ampsize
        
        #######################################################################
        ################### STEP 10: Generate html report #####################
        #######################################################################
        
        # Generate interactive html file using Quarto
       cat > params.yaml <<EOF
runid: "$runid"
working_dir: "$PWD"
marker: "$marker"
ampsize: $ampsize
numreps: $numreps
EOF

        quarto render "$scripts_dir/9b-MAPLE_reporting_with_bins.qmd" --execute-params params.yaml
        mv "$scripts_dir/9b-MAPLE_reporting_with_bins.html" "$PWD/MAPLE Report - "${runid}"_"${marker}"_${ampsize}.html"

        # Tidy up directory
        rm params.yaml temp_vsearch_output.txt filtered_vsearch_output.txt parsed_vsearch_output.txt final_output.txt final_output_with_bin_match.txt temp_toBIN.fasta temp_mapping*.txt   
    else
        # Generate interactive html file using Quarto
       cat > params.yaml <<EOF
runid: "$runid"
working_dir: "$PWD"
marker: "$marker"
ampsize: $ampsize
numreps: $numreps
EOF

        quarto render "$scripts_dir/9a-MAPLE_reporting.qmd" --execute-params params.yaml
        mv "$scripts_dir/9a-MAPLE_reporting.html" "$PWD/MAPLE Report - "${runid}"_"${marker}"_${ampsize}.html"

        # tidy up directory
        rm params.yaml temp_mapping*.txt
    fi

    # Tidy up directory and move back to main working directory
    mkdir -m 777 "1-Results and Report"
    mkdir -m 777 "2-TSV Versions of Results"
    mkdir -m 777 "3-Negative Control ASVs"
    mv *.xlsx *.html "1-Results and Report"
    mv Metabarcoding_Results*.tsv "2-TSV Versions of Results"
    mv *NegativeControlASVs.tsv "3-Negative Control ASVs"
    cd $working_dir
done

#######################################################################
##################### STEP 11: Clean up directory #####################
#######################################################################

cd $working_dir
rm runinfo.txt "$runid"_readcounts.txt "mapping_$runid.txt" "metadata_$runid.txt" name_map.tsv
mkdir -m 777 "Input_Files"
mkdir -m 777 "Demultiplexing Results"
mkdir -m 777 "$runid"
mv parameters_*.xlsx "./Input_Files"
mv log.txt "./Input_Files"
mv all.fastq.gz "./Input_Files"
mv Demultiplexing_Results*"$runid".pdf Individual_Raw_Fastq_Files.tar.gz "./Demultiplexing Results"
mv * "$runid"
echo -e '\n\n\n########## MAPLE ANALYSIS COMPLETE ##########'