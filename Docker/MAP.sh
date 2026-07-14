#!/usr/bin/env bash
clear
START_TIME=$(date +%s)
start_human=$(date "+%Y-%m-%d %H:%M:%S")
echo -e '\n\n\n########## STARTING MAP ANALYSIS ##########'
# v 1.0


# REQUIREMENTS:
# 1) FASTQ files from sequencer
# 2) Parameters file

#######################################################################
############################## OUTLINE ################################

########################### FOR FASTQ FILE ############################
################# STEP 1: Set user-specific arguments #################
################# STEP 2: Collect run parameters ######################
############ STEP 3a: Merge PE reads (Illumina only) ##################
####### STEP 3b: Concatenate raw fastq file (if not pair-ends) ########
############# STEP 4: Filter, demultiplex, primer trim ################

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~# TASK 1 #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
####### STEP 5: Produce low sequence variant clusters by sample #######
#~~~~~~~~~~~~~~~~~~~~~~~~~~~# END TASK 1 #~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

################# STEP 6: Run-wide OTU clustering #####################
######### STEP 7: Chimera check (only for long read) ##################
################### STEP 8: Sequence correction #######################
################### STEP 9: Identify run-wide OTUs ####################
#################### STEP 10: Remove contaminants #####################
############ STEP 11: Perform BIN analysis (if applicable) ############
############# STEP 12: Generate output and html report ################
######################### END MARKER FOR LOOP #########################

############################ END OUTLINE ##############################
#######################################################################



#######################################################################
################# STEP 1: Set user-specific arguments #################
#######################################################################                     


fastq_file="${fastq_file:-/MAP/Metabarcoding/PHAUS_1K_RawReads.fastq.gz}" # You may use a wildcard here to refer to multiple files, but please assign a similar prefix! (e.g., PHAUS_Illumina_*.fastq.gz)
params_file="${params_file:-/MAP/Metabarcoding/parameters.xlsx}"
reference_lib_dir="${reference_lib_dir:-/MAP/REFS}"
working_dir="${working_dir:-/MAP/Metabarcoding}" 
scripts_dir="${scripts_dir:-/MAP/SCRIPTS}"     
sintax_cutoff=0.6 #0-1
componentreads=0 #1 for yes, 0 for no
cores_to_leave=2 #How many cores to leave free. MAP will use the rest.
ref_seq_corr="${ref_seq_corr:-/MAP/REFS/reference_seqs_327K.fasta}" # File used for sequence correction

#~#~#~#~#~#~#~#~#~#~#
#Advanced parameters
#~#~#~#~#~#~#~#~#~#~#

#### Ia. UMI and primer overlap minimum proportion. 
# This value is used as a multiplier with UMI lengths during demultiplexing. We recommend 0.75 for UMIs greather than or equal to 12 nucleotides, and 1.0 for UMIs with fewer than 12 nucleotides.
umi_overlap_min=0.75
# This value is used as a multiplier with primer lengths. We recommend 0.75.
primer_overlap_min=0.75
#### Ib. Maximum error rate (Cutadapt parameter) for UMIs (default is 0.125, could be lower for shorter UMI sequences [e.g., < 10bp])
error_umi1=0.125
error_umi2=0.125
#### Ic. Maximum error rate (Cutadapt parameter) for primers (default is 0.2)
error_primer1=0.2
error_primer2=0.2
#### Id. Minimum quality score permitted with Chopper. We recommend 10 for long-read sequences with lower overall predicted quality (Oxford Nanopore) for better retention.
# Quality score closer to 20 is likely more useful for sequences associated with higher overall predicted quality (e.g., HiFi PacBio).
minqual=10
# ... or (if illumina PE data) Maximum expected error value permitted.
maxee=1

#### Ie. minimum and maximum lengths of sequence+primers to be removed with Chopper. 
# These are mainly in place to remove primer dimers and other sequencing errors. We keep these loose to account for variable markers with variable polymorphic lengths. 
min_read_and_primer_length=100
max_read_and_primer_length=1000

# Chimera removal (Illumina)
Ill_abskew=10 # Only for paired-end (e.g., Illumina) data. Parameter to be used for the VSEARCH's uchime_denovo command (abskew)
Ill_mindiv=0.0005 #Only for paired-end (e.g., Illumina) data. Parameter to be used for the VSEARCH's uchime_denovo command (mindiv)

# Chimera removal (long read data)
LR_abskew=10 # Only for long read (e.g., Oxford Nanopore) data. Parameter to be used for the VSEARCH's uchime_denovo command (abskew)
LR_mindiv=0.0005 #Only for long read (e.g., Oxford Nanopore) data. Parameter to be used for the VSEARCH's uchime_denovo command (mindiv)

#### II. Low sequence variant clustering (by sample)
minsize_unoise=2 # Minimum cluster size (within samples) for paired-end reads using VSEARCH cluster_unoise

#### III. Remove Contaminants
min_rep_prop=0.65 #minimum proportion of replicate samples for OTU to be present in in order for OTUs to be retained
high_dens_prop=0.65 #minimum proportion of replicate samples (of the replicates where OTU is actually present) for OTU sequence count within high kernel density area for OTUs to be retained
alpha_default=0.001 #only if alpha_quantile not defined.
alpha_quantile=0.9 #which quantile of the calculated alpha (contaminated reads floor estimate) to use? 
k_multiplier=1 #factor multiplied by alpha to determine contaminant floor. 1 by default, less than 1 will be more permissive, over 1 will be stricter

#### III. BIN identification (optional)
BIN_percent_ID=0.85 # We recommend 0.85 for short-read data. This is the threshold to keep matches to BINs in final database, regardless of 'formal' BIN assignment. Thresholds below the universal threshold will be labeled "NO MATCH" in final data. 
BIN_maxaccepts=3 #parameter to feed VSEARCH's usearch_global command.
BIN_maxhits=3 #parameter to feed VSEARCH's usearch_global command.


while [[ $# -gt 0 ]]; do
 case "$1" in 
    --fastq)
        fastq_file="$2" 
        shift 2 
        ;;
    --params)
        params_file="$2"
        shift 2
        ;;
    --refs)
        reference_lib_dir="$2" 
        shift 2 
        ;;
    --wd)
        working_dir="$2" 
        shift 2 
        ;; 
    --scripts)
        scripts_dir="$2" 
        shift 2 
        ;;
    --sintax_cutoff)
        sintax_cutoff="$2"
        shift 2
        ;;   
    --componentreads)
        componentreads="$2"
        shift 2
        ;;
    --cores_to_leave)
        cores_to_leave="$2"
        shift 2
        ;;
    --ref_seq_corr)
        ref_seq_corr="$2" 
        shift 2 
        ;;
    --umi_overlap_min)
        umi_overlap_min="$2"
        shift 2
        ;;  
    --primer_overlap_min)
        primer_overlap_min="$2"
        shift 2
        ;;      
    --error_umi1)
        error_umi1="$2"
        shift 2
        ;;
    --error_umi2)
        error_umi2="$2"
        shift 2
        ;; 
    --error_primer1)
        error_primer1="$2"
        shift 2
        ;;
    --error_primer2)
        error_primer2="$2"
        shift 2
        ;;
    --maxee)
        maxee="$2"
        shift 2
        ;;
    --minqual)
        minqual="$2"
        shift 2
        ;;
    --min_read_and_primer_length)
        min_read_and_primer_length="$2"
        shift 2
        ;;
    --max_read_and_primer_length)
        max_read_and_primer_length="$2"
        shift 2
        ;;     
    --Ill_abskew)
        Ill_abskew="$2"
        shift 2
        ;;
    --Ill_mindiv)
        Ill_mindiv="$2"
        shift 2
        ;;      
    --LR_abskew)
        LR_abskew="$2"
        shift 2
        ;;
    --LR_mindiv)
        LR_mindiv="$2"
        shift 2
        ;;
    --minsize_unoise)
        minsize_unoise="$2"
        shift 2
        ;;     
    --min_rep_prop)
        min_rep_prop="$2"
        shift 2
        ;; 
    --high_dens_prop)
        high_dens_prop="$2"
        shift 2
        ;; 
    --alpha_default)
        alpha_default="$2"
        shift 2
        ;; 
    --alpha_quantile)
        alpha_quantile="$2"
        shift 2
        ;;
    --k_multiplier)
        k_multiplier="$2"
        shift 2
        ;;
    --BIN_percent_ID)
        BIN_percent_ID="$2"
        shift 2
        ;;
    --BIN_maxaccepts)
        BIN_maxaccepts="$2"
        shift 2
        ;;
    --BIN_maxhits)
        BIN_maxhits="$2"
        shift 2
        ;;
    *)
        echo "Unknown option: $1" 
        exit 1 
        ;; 
    esac
done

#### Resolve relative input paths against the launch directory (BEFORE any cd) (to avoid using the full $HOME path). 
make_abs() { case "$1" in /*) printf '%s' "$1" ;; "") printf '' ;; *) printf '%s/%s' "$PWD" "$1" ;; esac; }
fastq_file="$(make_abs "$fastq_file")"
params_file="$(make_abs "$params_file")"
reference_lib_dir="$(make_abs "$reference_lib_dir")"
working_dir="$(make_abs "$working_dir")"
scripts_dir="$(make_abs "$scripts_dir")"
ref_seq_corr="$(make_abs "$ref_seq_corr")"

echo "Using fastq file: $fastq_file"
echo "Using parameters file: $params_file"
echo "Using reference directory: $reference_lib_dir"
echo "Using working directory: $working_dir"
echo "Using scripts directory: $scripts_dir"
echo "Sintax cutoff: $sintax_cutoff"
echo "Using component reads (1 for yes, 0 for no): $componentreads"

#### Fetch the BOLDdistilled sintax reference library on first use.

if ! ls "$reference_lib_dir"/BOLDistilled*.fasta >/dev/null 2>&1; then
    echo -e "\n****** BOLDdistilled sintax reference library not found in $reference_lib_dir — downloading latest..."
    mkdir -p "$reference_lib_dir"
    reflib_tmp="$(mktemp -d)"
    curl -fSL https://us-sea-1.linodeobjects.com/boldistilled/sintax.zip -o "$reflib_tmp/sintax.zip"
    python -m zipfile -e "$reflib_tmp/sintax.zip" "$reflib_tmp"
    mv "$reflib_tmp"/sintax/* "$reference_lib_dir"/
    rm -rf "$reflib_tmp"
    echo "****** Reference library download complete."
fi

#######################################################################
############################# FUNCTIONS ###############################
#######################################################################

#### Function to process one chunk
process_chunk() {
    local chunk="$1"
    local outfile="$2"

    local tmp_out="${outfile}_$(basename "$chunk").tmp"

    echo "  [START] Processing chunk $(basename "$chunk")..."

    # ensure file exists
    : > "$tmp_out"
    
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

process_fasta() {
    local fasta="$1"
    local filtered="$2"
    echo "Processing $fasta..."
    awk -v min="$minreads" '
        /^>/ {
            reads = 0
            pos = index($0, "|reads-")
            if (pos > 0) reads = substr($0, pos + 7) + 0
            keep = (reads >= min)
        }
        keep { print }
    ' "$fasta" > "$filtered"
}


export -f process_fasta
############################ END FUNCTIONS ##############################
#########################################################################

#########################################################################
########################### PROCESS (task1) #############################

#~#~#~#~#~#~#~#~#~#~#~#~#~#~##~#~#~#
#~#~#~#~#~#~#~#~ TASK 1 ~#~#~#~#~#~#
#~#~#~#~#~~#~#~#~#~#~#~#~#~#~#~#~#~#

task1() (
#######################################################################
####### STEP 5: Produce low sequence variant clusters by sample #######
#######################################################################
    echo "Using $cores cores"
    echo "Using $minreads minreads" 
    local fasta="$1"
    local sampleid="$(basename "$fasta" .fasta)"
    

    #### FOR PAIRED END READS only - Chimera screen using UCHIME (only if ampsize is >= 200 bp)
    if [[ "$ampsize" -ge 200 && "$pe_reads" = "Yes" ]]; then        
        # Dereplicate for chimera screen
        vsearch --derep_fulllength "${sampleid}.fasta" \
        --output "${sampleid}.d.fasta" \
        --sizeout 
        wait
        echo -e "******** Performing chimera screen..."
        vsearch --uchime_denovo "${sampleid}.d.fasta" \
        --chimeras "${sampleid}.chimeras.fasta" \
        --nonchimeras "${sampleid}.nonchimeras.fasta" \
        --fasta_width 0 \
        --abskew $Ill_abskew \
        --mindiv $Ill_mindiv
    else
        mv "${sampleid}.fasta" "${sampleid}".nonchimeras.fasta
    fi

    #### Cluster 
    echo -e "******** Making low sequence variant clusters for $sampleid..."
    if [[ "$pe_reads" = "Yes" ]]; then
        if [ "$componentreads" -eq 0 ]; then
           vsearch --cluster_unoise "${sampleid}.nonchimeras.fasta" \
            --consout "${sampleid}_consensus.fasta" \
            --iddef 3 \
            --sizein \
            --sizeout \
            --minsize $minsize_unoise \
            --threads 1
        else
            vsearch --cluster_unoise "${sampleid}.nonchimeras.fasta" \
            --consout "${sampleid}_consensus.fasta" \
            --clusters "${sampleid}|OTU" \
            --iddef 3 \
            --sizein \
            --sizeout \
            --minsize $minsize_unoise \
            --threads 1
        fi
    else
        if [ "$componentreads" -eq 0 ]; then
           vsearch --cluster_fast "${sampleid}.nonchimeras.fasta" \
            --id "$(awk -v d="$otu_dist1" 'BEGIN { printf "%.6f", (100 - d)/100 }')" \
            --consout "${sampleid}_consensus.fasta" \
            --iddef 3 \
            --sizeout \
            --threads 1
        else
            vsearch --cluster_fast "${sampleid}.nonchimeras.fasta" \
            --id "$(awk -v d="$otu_dist1" 'BEGIN { printf "%.6f", (100 - d)/100 }')" \
            --consout "${sampleid}_consensus.fasta" \
            --clusters "${sampleid}|OTU" \
            --iddef 3 \
            --sizeout \
            --threads 1
        fi
    fi


    #### Delete pre-clustered file
    rm "$sampleid".fasta "$sampleid".d.fasta "$sampleid".nonchimeras.fasta

    #### Rename OTU consensus sequence headers to correct OTU name
    consensus_fasta="${sampleid}_consensus.fasta"
    renamed_fasta="${sampleid}_consensus_renamed.fasta"

    awk -v sampleid="$sampleid" '
        BEGIN { otu=0 }
        /^>/ {
            reads=1
            # extract number after ";size=" if present
            if ($0 ~ /;size=[0-9]+/) {
                # split on ";size=" and take the number part
                split($0, parts, ";size=")
                reads=parts[2] + 0   # convert to number
            }
            print ">" sampleid "|OTU" otu "|reads-" reads
            otu++
            next
        }
        { print }
        ' "$consensus_fasta" > "$renamed_fasta"
    mv "$renamed_fasta" "$consensus_fasta"

    #### Filter OTU sequences with fewer than minreads (in both consensus and component read files if necessary)
    echo "******** Filtering OTUs below $minreads reads..."
    consensus_fasta="${sampleid}_consensus.fasta"
    filtered_fasta="${sampleid}_consensus_filtered.fasta"

    #### Export variables for parallel jobs
    export minreads cores

    #### Run filtering
    process_fasta "$consensus_fasta" "$filtered_fasta"

    if [ "$componentreads" -eq 1 ]; then

        export minreads sampleid

        process_otu() {
            local f="$1"

            # skip if file doesn't exist (race conditions)
            [ -e "$f" ] || return

            # count sequences
            local seq_count
            seq_count=$(grep -c "^>" "$f")

            local newname="${f}_ComponentReads.fasta"
            echo "Keeping $f → $newname"
            cp -- "$f" "$runid"_"$marker"_"$ampsize"_OTU_Component_Reads/"$newname"
            rm -f -- "$f"
        }

        export -f process_otu

        echo "Finding OTU component files for sample '$sampleid'..."

        # Stream into GNU Parallel safely (handles unlimited files)
        find . -maxdepth 1 -type f -name "${sampleid}|OTU*" -print0 \
            | parallel -0 --bar -j 1 process_otu {}
    fi

    #### Replace consensus file with filtered one
    mv "${sampleid}_consensus_filtered.fasta" "${sampleid}_consensus.fasta"

    #### Convert final OTU FASTA to single-line
    awk '{if(NR==1) {print $0} else {if($0 ~ /^>/) {print "\n"$0} else {printf $0}}}' "${sampleid}_consensus.fasta" > "${sampleid}_consensus2.fasta"

    rm "${sampleid}_consensus.fasta"

    #### Prepare global cleanup step
    echo -e "Preparing global cleanup step..."
    grep -h ">" "${sampleid}_consensus2.fasta" \
        | sed 's/^>//; s/|reads-.*$//' >> keep.txt

)


############################ END PROCESSES ##############################
#########################################################################


#########################################################################
###################### EXECUTE PIPELINE: START ##########################

#######################################################################
################# STEP 2: Collect run parameters #######################
#######################################################################

#### Detect number of cores
cores=$(($(getconf _NPROCESSORS_ONLN 2>/dev/null || sysctl -n hw.ncpu) - $cores_to_leave))

#### Collect and reformat parameters information
cd $working_dir

python3.12 <<XL1
import pandas as pd
import re
import sys



PF = "${params_file}"
xls = pd.ExcelFile(PF)

def fmt(x):
    s = str(x)
    if re.fullmatch(r"[0-9.]+", s):
        return f"{round(float(x), 1):.1f}"
    return s

def cell(x):
    return "" if (x is None or (isinstance(x, float) and pd.isna(x))) else x

# 'UMIs and Primers': Run Name (top header) + per-sample table 
up_raw = pd.read_excel(PF, sheet_name="UMIs and Primers", header=None)
plate_row = None
for i in range(up_raw.shape[0]):
    if str(cell(up_raw.iloc[i, 0])).strip().lower() == "plate":
        plate_row = i
        break
if plate_row is None:
    sys.exit("ERROR: no 'Plate' header row found in 'UMIs and Primers'.")

run_id = None
for i in range(plate_row):
    vals = [str(cell(v)) for v in up_raw.iloc[i, :].tolist()]
    for j, v in enumerate(vals):
        if v.strip().rstrip(":").strip().lower() == "run name":
            for k in range(j + 1, len(vals)):
                if vals[k].strip() != "":
                    run_id = vals[k].strip()
                    break
            break
    if run_id:
        break
if not run_id:
    sys.exit("ERROR: 'Run Name' not found in the header block atop 'UMIs and Primers'.")

up = pd.read_excel(PF, sheet_name="UMIs and Primers", header=plate_row)
up.columns = [str(c).strip() for c in up.columns]
up = up[up["Sample"].notna() & (up["Sample"].astype(str).str.strip() != "")]

# 'Dictionary Update': primer lookup + run-level params 
dic = pd.read_excel(PF, sheet_name="Dictionary Update")
dic.columns = [str(c).strip() for c in dic.columns]
if dic.shape[0] == 0:
    sys.exit("ERROR: 'Dictionary Update' tab is empty.")

run_cols = ["Paired-End Reads", "Min Reads per OTU", "Replicates per Sample",
            "Intra-OTU Clustering Threshold", "Inter-OTU Clustering Threshold"]
miss_run = [c for c in run_cols if c not in dic.columns]
if miss_run:
    sys.exit(f"ERROR: 'Dictionary Update' missing run-level column(s): {miss_run}")

pe_reads = str(cell(dic.iloc[0]["Paired-End Reads"])).strip()
minreads = cell(dic.iloc[0]["Min Reads per OTU"])
numreps  = cell(dic.iloc[0]["Replicates per Sample"])
otu1     = cell(dic.iloc[0]["Intra-OTU Clustering Threshold"])
otu2     = cell(dic.iloc[0]["Inter-OTU Clustering Threshold"])

dict_need = ["Forward Primer Name", "Reverse Primer Name",
             "Forward Primer Sequence", "Reverse Primer Sequence",
             "Marker", "Reference Library",
             "Min Amplicon Length", "Max Amplicon Length", "Target Amplicon Length"]
miss_d = [c for c in dict_need if c not in dic.columns]
if miss_d:
    sys.exit(f"ERROR: 'Dictionary Update' missing column(s): {miss_d}")

lut = {}
for _, d in dic.iterrows():
    lut[(str(d["Forward Primer Name"]).strip(),
         str(d["Reverse Primer Name"]).strip())] = d

# runinfo.txt : runid, pe_reads, numreps, minreads, otu_dist1, otu_dist2 
with open("runinfo.txt", "w") as f:
    for v in [run_id, pe_reads, fmt(numreps), fmt(minreads), fmt(otu1), fmt(otu2)]:
        f.write(f"{v}\n")

#  mapping_<run>.txt : expand primer NAMES to sequences + dictionary fields 
MAP_COLS = ["Plate", "Well", "Sample", "Forward Primer", "Reverse Primer",
            "Forward UMI", "Reverse UMI", "Marker",
            "Min Amplicon Length", "Max Amplicon Length", "Target Amplicon Length",
            "Reference Library", "Negative Control"]
rows = []
for _, r in up.iterrows():
    fn = str(cell(r.get("Forward Primer Name", ""))).strip()
    rn = str(cell(r.get("Reverse Primer Name", ""))).strip()
    if (fn, rn) not in lut:
        sys.exit(f"ERROR: primer pair not in 'Dictionary Update': {fn} / {rn}")
    d = lut[(fn, rn)]
    rows.append([
        cell(r.get("Plate", "")), cell(r.get("Well", "")), cell(r.get("Sample", "")),
        d["Forward Primer Sequence"], d["Reverse Primer Sequence"],
        cell(r.get("Forward UMI", "")), cell(r.get("Reverse UMI", "")),
        d["Marker"], d["Min Amplicon Length"], d["Max Amplicon Length"],
        d["Target Amplicon Length"], d["Reference Library"],
        cell(r.get("Negative Control", "")),
    ])
out = pd.DataFrame(rows, columns=MAP_COLS)
with open(f"mapping_{run_id}.txt", "w") as fh:
    fh.write("\t".join(MAP_COLS) + "\n")
    out.to_csv(fh, sep="\t", header=False, index=False)

# metadata_<run>.txt : from 'Bulk Sample Metadata' 
meta_sheet = "Bulk Sample Metadata" if "Bulk Sample Metadata" in xls.sheet_names else "Sample Metadata"
metadata = pd.read_excel(PF, sheet_name=meta_sheet)
metadata.to_csv(f"metadata_{run_id}.txt", sep="\t", header=True, index=False)
XL1

#### Remove any Windows carriage returns from mapping file
sed -i.bak 's/\r$//' mapping*.txt && rm -f mapping*.txt.bak

#### Extract info from parameters file
read -r runid < <(sed -n '1p' runinfo.txt)
read -r pe_reads < <(sed -n '2p' runinfo.txt)
read -r numreps  < <(sed -n '3p' runinfo.txt)
numreps=${numreps%.*}   # convert to integer
read -r minreads < <(sed -n '4p' runinfo.txt)
minreads=${minreads%.*}   # convert to integer
read -r otu_dist1  < <(sed -n '5p' runinfo.txt)
read -r otu_dist2  < <(sed -n '6p' runinfo.txt)


#### Extract forward and reverse UMIs from mapping file

awk -F'\t' -v err1="$error_umi1" -v err2="$error_umi2" -v umi_ov=$umi_overlap_min '
NR>1 && $6!="" && $6!="NA" && $7!="" && $7!="NA" {
    key = $6 FS $7
    if (!(key in seen)) {
        seen[key]=1
        f[++n]=$6
        r[n]=$7
        fnum[n]= int(length($6) * umi_ov)
        rnum[n]= int(length($7) * umi_ov)
    }
}
END {
    for (i=1;i<=n;i++)
        print ">"f[i]"\n"f[i]";min_overlap="fnum[i]";max_error_rate=" err1 > "fwd_umis.fasta"
    close("fwd_umis.fasta")
    for (i=1;i<=n;i++)
        print ">"r[i]";min_overlap="rnum[i]";max_error_rate=" err2 "\n"r[i] > "rev_umis.fasta"
}
' mapping_${runid}.txt

seqtk seq -r -c rev_umis.fasta > rev_umis_rc_raw.fasta
rm rev_umis.fasta

awk '
NR % 2 == 1 {
    # header line
    split($0, a, ";")
    print a[1]
    suffix = (length(a) > 1 ? ";" a[2] ";" a[3] : "")
    next
}
{
    # sequence line
    print $0 suffix
}
' rev_umis_rc_raw.fasta > rev_umis_rc.fasta
rm rev_umis_rc_raw.fasta

#### Make linked primers file.
paste fwd_umis.fasta rev_umis_rc.fasta \
|   awk '
    {
    print $1 "..." $2
    }
' > linked_umis.fasta

#### Calculate minimum length of fwd and rev UMIs
min_umi_len_fwd=$(awk 'NR>1 {print length($0)}' fwd_umis.fasta | sort -n | head -1)
min_umi_len_rev=$(awk 'NR>1 {print length($0)}' rev_umis_rc.fasta | sort -n | head -1)

#######################################################################
############ STEP 3a: Merge PE reads (Illumina only) ##################
#######################################################################

#### Check if any file matching the pattern shares the same prefix as all.fastq.gz (to save time for repeat analysis)
match=0
if [ -f all.fastq.gz ]; then
    prefix_all=$(gzip -dc all.fastq.gz 2>/dev/null | head -c 5)
    for f in $fastq_file; do
        [ -f "$f" ] || continue
        prefix_f=$(gzip -dc "$f" 2>/dev/null | head -c 5)
        if [ "$prefix_all" = "$prefix_f" ]; then
            match=1
            break
        fi
    done
fi

#### Merge paired-end reads (PAIRED-END READS only)
if [ "$pe_reads" = "Yes" ]; then #if already a copy of all.fastq.gz present, skip this step. must be merged fastq of same input files!
    if [ "$match" = "1" ]; then
        #### Decompress, but keep all.fastq.gz 
        echo -e ****** "decompressing all.fastq.gz"
        pigz -d -k -p $cores all.fastq.gz
    else
        echo "WARNING: all.fastq.gz exists but does not match input files. Renaming to OLD_all.fastq.gz and re-merging."
        mv all.fastq.gz OLD_all.fastq.gz
        #### Get read1 and read2 from fastq files. must follow either _R1* / _R2* OR _1.fastq.gz /_2.fastq.gz convention.
        read1=( $(ls *.gz | grep -E '_R1[_.]|_1\.fastq') )
        read2=( $(ls *.gz | grep -E '_R2[_.]|_2\.fastq') )

        #### Merge paired end reads
        echo -e "******** Merging paired-end reads..."
        pear -j $cores -f $read1 -r $read2 -o $runid > log.txt

        #### Delete discarded and unassembled paired-end reads
        rm $runid".discarded.fastq" $runid".unassembled.forward.fastq" $runid".unassembled.reverse.fastq"
        mv $runid.assembled.fastq all.fastq
    fi
fi

#######################################################################
####### STEP 3b: Concatenate raw fastq file (if not pair-ends) ########
#######################################################################

#### Merge FASTQ files into single file
if [ "$pe_reads" = "No" ]; then
    echo -e "******** Merging FASTQ files..."

    # Decompress and concatenate all FASTQ files into single file.
    ulimit -n 65536
    echo $fastq_file | tr ' ' '\n' | xargs pigz -p $cores -dc > all.fastq2
    
    # Remove any extra text from sequence headers
    awk '{if(NR%4==1) sub(/\t.*/, "", $0); print}' all.fastq2 > all.fastq
    rm all.fastq2
fi

#######################################################################
############## STEP 4: Filter, demultiplex, primer trim ###############
#######################################################################

#### Count raw reads, write to readcounts file.
seqkit stats -T all.fastq | awk 'NR==2 {print "Raw Reads\t"$4}' > "${runid}_readcounts.txt"

#### Remove reads with low quality scores, and primer dimer forming reads that are outliers on the size distribution.
if [ "$pe_reads" = "Yes" ]; then 
    vsearch --fastq_filter all.fastq --fastq_maxee $maxee --fastqout all_filt.fastq
else
    chopper -q $minqual --minlength $min_read_and_primer_length --maxlength $max_read_and_primer_length < all.fastq > all_filt.fastq
fi

#### Save a copy of all.fastq.gz (for now)
if [ "$match" = "1" ]; then
    echo -e "******** all.fastq.gz found, skipping..."
else
    pigz -p $cores -k all.fastq
fi
rm all.fastq

#### Demultiplex forward UMIs at 5' end 
echo -e "******** Demultiplexing reads..."

##### Evaluate whether the UMIs provided are symmetrical or asymmetrical.
if [ $(diff <(awk 'NR>1 {print $6}' mapping_"$runid".txt) \
           <(awk 'NR>1 {print $7}' mapping_"$runid".txt) | wc -l) -gt 0 ]; then
    symm="asymmetrical UMIs"
    echo "$symm"
else
    symm="symmetrical UMIs"
    echo "$symm"
fi

#### Extract forward and reverse primers from mapping file

awk -F'\t' -v err1="$error_primer1" -v err2="$error_primer2" -v prim_ov=$primer_overlap_min '
NR>1 && $4!="" && $4!="NA" && $5!="" && $5!="NA" {
    key = $4 FS $5
    if (!(key in seen)) {
        seen[key]=1
        f[++n]=$4
        r[n]=$5
        fnum[n]= int(length($4) * prim_ov)
        rnum[n]= int(length($5) * prim_ov)
    }
}
END {
    for (i=1;i<=n;i++)
        print ">"f[i]"\n"f[i]";min_overlap="fnum[i]";max_error_rate=" err1 > "fwd_primers.fasta"

    for (i=1;i<=n;i++)
        print ">"r[i]";min_overlap="rnum[i]";max_error_rate=" err2 "\n"r[i] > "rev_primers.fasta"
}
' mapping_"${runid}.txt"

seqtk seq -r -c rev_primers.fasta > rev_primers_rc_raw.fasta
rm rev_primers.fasta

awk '
NR % 2 == 1 {
    # header line
    split($0, a, ";")
    print a[1]
    suffix = (length(a) > 1 ? ";" a[2] ";" a[3] : "")
    next
}
{
    # sequence line
    print $0 suffix
}
' rev_primers_rc_raw.fasta > rev_primers_rc.fasta
rm rev_primers_rc_raw.fasta

#### Make linked primers file.
paste fwd_primers.fasta rev_primers_rc.fasta \
|   awk '
    {
    print $1 "..." $2
    }
' > linked_primers.fasta

#### Orient sequences by searching the primers. This will make primers, umis, and other non-amplicon seqs lowercase. 
if [ "$pe_reads" = "Yes" ]; then # No reverse complementing for Illumina
    cutadapt -j $cores \
        -g file:linked_primers.fasta \
        --action=lowercase \
        -o "${runid}"_g_link_primer.fastq \
        all_filt.fastq 
    else #### Allow reverse complementing for non-Illumina
        cutadapt -j $cores \
        -g file:linked_primers.fasta \
        --action=lowercase \
        -o "${runid}"_g_link_primer.fastq \
        all_filt.fastq \
        --revcomp    
fi

rm all_filt.fastq

#### Add info to readcounts file.
seqkit stats -T "${runid}"_g_link_primer.fastq | awk 'NR==2 {print "After removing primer dimers\t"$4}' >> "${runid}_readcounts.txt"

#### Demultiplex using linked UMIs.
cutadapt -j $cores \
    -g file:linked_umis.fasta \
    --action=trim \
    -o {name}.fastq \
    "${runid}"_g_link_primer.fastq \
    --untrimmed-output "${runid}"_unt.fwd.fastq 

#### Get read count of demultiplexed reads with linked UMIs
sum=0
for k in *'>'*.fastq; do
    sum=$((sum + $(echo $(cat "$k" | wc -l)/4 | bc)))
done
echo $sum | awk '{print "With linked UMIs" "\t" $1}' >> "${runid}_readcounts.txt"


##################################
##### if symmetrical UMIS... #####
##################################

if [ "$symm" = "symmetrical UMIs" ]; then
    echo -e "******** Symmetrical UMIs detected ... "
#### Here, prepare the single, unlinked symmetrical UMIs (0-X or X-0). Start with the fwd UMIs.
    cutadapt -j $cores \
        -g file:fwd_umis.fasta \
        --action=none \
        -o mark_fwd_untrim.fastq \
        "${runid}"_unt.fwd.fastq \
        --rename {adapter_name}_fwd_sing \
        --untrimmed-output "${runid}"_fwd_leftover_untrim.fasta 

#### Repeat UMI prep for demultiplexing on reverse UMIs.
    cutadapt -j $cores \
        -a file:rev_umis_rc.fasta \
        --action=none \
        -o mark_rev_untrim.fastq \
        "${runid}"_unt.fwd.fastq \
        --rename {adapter_name}_rev_sing \
        --untrimmed-output "${runid}"_rev_leftover_untrim.fasta
    
#### Combine all seqs where UMIs were found. 
    cat mark_fwd_untrim.fastq mark_rev_untrim.fastq > mark_all_untrim.fastq
    rm mark_fwd_untrim.fastq mark_rev_untrim.fastq


#### Since all the linked symmetrical UMIs have been demultiplexed, all other "linked" UMIs are tag switched. Remove any duplicate (tag-switched) reads. 


#### Remove the first of the duplicates and create file showing which sequences these are.
    seqkit seq -s -w 0 mark_all_untrim.fastq | sort | uniq -d > dup_seqs.txt

#### Remove any read whose sequence matches a duplicate sequence
    awk '
    NR==FNR { dups[$1]=1; next }
    /^@/ { header=$0; getline seq; getline plus; getline qual }
    !(seq in dups) { print header "\n" seq "\n" plus "\n" qual }
    ' dup_seqs.txt mark_all_untrim.fastq |   seqkit replace -p '^.*_([^_]+_[^_]+_[^_]+)$' -r '$1'   > mark_uniq_ref.fastq

    seqkit split -i mark_uniq_ref.fastq 
    
    #### Get count of reads with single forward UMI
    (
    sum=0
    for k in mark_uniq_ref.fastq.split/*_fwd_sing.fastq; do
        sum=$((sum + $(echo $(cat "$k" | wc -l)/4 | bc)))
    done
    echo $sum | awk '{print "Forward-only symmetric UMIs" "\t" $1}' >> "${runid}_readcounts.txt"
    ) &

    #### Get count of reads with single reverse UMI
    (
    sum=0
    for k in mark_uniq_ref.fastq.split/*_rev_sing.fastq; do
        sum=$((sum + $(echo $(cat "$k" | wc -l)/4 | bc)))
    done
    echo $sum | awk '{print "Reverse-only symmetric UMIs" "\t" $1}' >> "${runid}_readcounts.txt"
    ) &
    wait

    sanitized_name_index="name_map.tsv"
    tail -n +2 "mapping_${runid}.txt" |
    while IFS=$'\t' read -r plt well samp fwdp revp fwdu revu marker minlen maxlen amplen _rest; do

        #### List all possible original filenames
        demux_f=(
            "${fwdu}...>${revu}.fastq"
            "mark_uniq_ref.fastq.split/mark_uniq_ref.part_${fwdu}_fwd_sing.fastq"
            "mark_uniq_ref.fastq.split/mark_uniq_ref.part_${revu}_rev_sing.fastq"
        )

        #### Assign new filename
        new="${samp}_${marker}_${amplen}bp.fastq.tmp"
        > "$new"  # truncate/create

        for old in "${demux_f[@]}"; do
            if [[ -f "$old" ]]; then
                cat "$old" >> "$new"
                rm "$old"
            fi
        done

        #### Sanitize sample name
        sanitized="${samp//#/-}"

        #### Append mapping
        printf "%s\t%s\n" "$samp" "$sanitized" >> "$sanitized_name_index"

    done
    rmdir mark_uniq_ref.fastq.split

    rm  *_sing_UMI_unt_marked.fasta *_all_untrim.fastq  \
    *_fwd_leftover_untrim.fasta *_rev_leftover_untrim.fasta mark*.fastq dup_seqs.txt


else

##################################
##### if asymmetrical UMIS... ####
##################################
    echo -e "******** Asymmetrical UMIs detected ... "

#### Change the fastq file names to match the samples
    sanitized_name_index="name_map.tsv"
    
    tail -n +2 "mapping_${runid}.txt" | while IFS=$'\t' read -r plt well samp fwdp revp fwdu revu marker minlen maxlen amplen _rest; do
    #### Original filename based on UMI mapping
        old="$(echo "${fwdu}...>${revu}.fastq" | tr -d '[:space:]')"

    #### New filename: append marker and amplicon length to make unique
        new="${samp}_${marker}_${amplen}bp.fastq.tmp"
        new_samp_only="${samp}"

        if [[ -f "$old" ]]; then
           mv -- "$old" "$new"
        fi

    #### Create sanitized filename for mapping: replace # with -
        sanitized="${new_samp_only//#/-}"

    #### Append to name_map.tsv: first column is sample, second is sanitized filename
        echo -e "${samp}\t${sanitized}" >> "$sanitized_name_index"
    done
fi

#### Remove the temporary ".tmp" suffix (will keep marker_amp in filename for uniqueness)
rename 's/\.tmp$//' *.tmp
rm *g_link_primer.fastq *.fasta *unt.fwd.fastq

#### Remove any FASTQ files with 0 reads
find . -type f -name "*.fastq" -size 0 -delete
tar -cf - *.fastq | pigz -p $cores > Individual_Raw_Fastq_Files.tar.gz

#### Convert to fasta and rename read headers to shorter versions (if spaces in raw read headers)
parallel -j $cores '
    seqkit fq2fa {} | cut -d" " -f1 > {.}.fasta
' ::: *.fastq
rm *.fastq

#### Generate heat map of demultiplexed reads by well per plate
echo -e "******** Generating demultiplexing report..."
Rscript "$scripts_dir/1-MAP_heatmap.R" "$PWD" "mapping_${runid}.txt" "$runid"


#### Build job table from mapping file
tail -n +2 "mapping_${runid}.txt" |
awk -F'\t' '{
    fasta = $3 "_" $8 "_" $11 "bp.fasta"
    print fasta, $9, $10
}' OFS='\t' > cutadapt_jobs.tsv

#### Strip lowercase nucleotide function:
strip_lower() {
    infile="$1"
    min="$2"
    max="$3"

    base="${infile%.fasta}"

    #Do UMI and primer trimming by removing all lowercase DNA seqs (these are umis and primers)
    awk 'NR % 2 == 0 { gsub(/[a-z]/, "") }1' "$infile" > "${base}.tr.fasta"

    cutadapt -j 1 \
        -m "$min" \
        -M "$max" \
        -o "${base}.fasta" \
        "${base}.tr.fasta"

}
export -f strip_lower
export min_read_and_primer_length max_read_and_primer_length

#### Strip lowercase adapters using the TSV file.
parallel -j "$cores" --colsep '\t' \
    strip_lower {1} {2} {3} :::: cutadapt_jobs.tsv

rm *.tr.fasta cutadapt_jobs.tsv

#### Organize FASTA files by marker/amplicon length
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

#### Replace '#' in filenames with '-' to avoid issues, and store original and sanitized names in an index file
find . -depth -name '*#*' | while read -r f; do
    newname="${f//#/-}"   # replace all # with -
    mv "$f" "$newname"
done

#### Make output dir
mkdir ./output

#~#~#~#~#~#~#~#~#~#~#~#~#~#~##~#~#~#
#~#~#~#~ MARKER FOR LOOP #~#~#~#~#~#
#~#~#~#~#~~#~#~#~#~#~#~#~#~#~#~#~#~#

for marker_dir in */; do
     #### Skip unwanted directory
    [[ "$marker_dir" == "Individual_Raw_Fasta_Files/" ]] || \
    [[ "$marker_dir" == "output/" ]] && continue

    echo "=== Processing marker folder: ${marker_dir%/} ==="
    
    #### Move into the marker folder, initialize some variables, and make OTU component read folder
    cd "$marker_dir" || continue
    dirname=$(basename "$PWD")
    marker="${dirname%_*}"
    ampsize="${dirname##*_}"
    ampsize="${ampsize%bp}" 

    #### Strip the marker_amp tags from filenames, then lookup the reference library for this marker_ampsize comination
    for f in *.fasta; do
        newname="$(echo "$f" | sed -E "s/_${marker}_${ampsize}bp\.fasta$/.fasta/")"
        mv "$f" "$newname"
    done

    reflib="$(awk -F'\t' -v m="$marker" -v a="$ampsize" 'NR>1 && $8==m && $11==a {print $12; exit}' "../mapping_${runid}.txt")"
    db_path=$(ls -d $reference_lib_dir/"$reflib"*.fasta)

    #### Create directory to store OTU component reads if necessary
    if [ "$componentreads" -eq 1 ]; then
        mkdir -m 777 "$runid"_"$marker"_"$ampsize"_OTU_Component_Reads
    fi
    
    # Print correct reference lib at the end of the run
    ref_file=$(ls -1 "$reference_lib_dir/$reflib"*.fasta 2>/dev/null | head -n1)

    # filename only, then strip extensions (e.g. .fasta, .SEQUENCES_sintax.fasta -> base name)
    ref_name=$(basename "$ref_file")          # e.g. BOLDistilled_COI_Jan2026_SEQUENCES_sintax.fasta
    ref_name=${ref_name%.fasta}               # drop trailing .fasta or .fa
    ref_name=${ref_name%.fa}                  

    # write a marker logfile
    marker_log="${marker}_${ampsize}bp.log"
    {
        echo "Run name:           ${runid}"
        echo "Marker:             ${marker}"
        echo "Amplicon length:    ${ampsize} bp"
        echo "Reference library:  ${ref_name:-unknown}"
    } > "$marker_log"

#####################
##### Run task1 #####
#####################


    #### Export everything task1 needs
    export -f task1 process_fasta process_otu
    export cores ampsize minreads componentreads otu_dist1 runid marker pe_reads Ill_abskew Ill_mindiv minsize_unoise

    fasta_files=( *.fasta )
    last_fasta="${fasta_files[-1]}"

    #### Run task1 in parallel
    parallel -j "$cores" task1 {} ::: *.fasta

    #### Run the global cleanup once
    echo "Running global cleanup step..."
    ls | grep -E '^[^\.]+$' | grep -v -f keep.txt | xargs rm
    rm keep.txt

    
    #### Delete files with 0 bytes
    find . -type f -name "*.fasta" -size 0 -exec rm -f {} +

    if [ "$componentreads" -eq 1 ]; then
        # Rename OTU component read files
        find ./"$runid"_"$marker"_"$ampsize"_OTU_Component_Reads/ -type f -name '*ComponentReads.fasta' -exec rename 's/[=]/-/g' {} +
        find ./"$runid"_"$marker"_"$ampsize"_OTU_Component_Reads/ -type f -name '*ComponentReads.fasta' -exec rename 's/[|]/_/g' {} +
    fi

    #### Merge all OTUs into a single master file
    seqkit seq *consensus2.fasta > all_otus_raw_consensus.fasta
    rm *consensus2.fasta *.chimeras.fasta

    #### Change OTU sequence names from |reads-n to ;size=n for subsequent size sorting
    sed '/^>/ s/|reads-/;size=/g' all_otus_raw_consensus.fasta > all_otus1_corrected_otus.fasta2

#######################################################################
################## STEP 6: Run-wide OTU clustering #####################
#######################################################################

    #### Cluster into run-wide OTUs
    vsearch --cluster_size all_otus1_corrected_otus.fasta2 \
        --id "$(awk -v d="$otu_dist2" 'BEGIN { printf "%.6f", (100 - d)/100 }')" \
        --consout all_otus_consensus.fasta \
        --iddef 3 \
        --threads $cores \
        --sizein \
        --sizeout \
        --uc all_otus_cluster_info.uc
  
    #### Reformat run-wide OTU info
    Rscript "$scripts_dir/2-MAP_runwideotuformat.R" "$PWD"

    rm all_otus_cluster_info.uc 

    #### Convert run-wide OTUs to single-line FASTA
    awk '{if(NR==1) {print $0} else {if($0 ~ /^>/) {print "\n"$0} else {printf $0}}}' all_otus_consensus.fasta > all_otus_consensus.fasta2
    rm all_otus_consensus.fasta

#######################################################################
############ STEP 7: Chimera check (only for long read) ###############
#######################################################################

    #### Chimera screen using UCHIME (only if ampsize is >= 200 bp)
    if [[ "$ampsize" -ge 200 && "$pe_reads" = "No" ]]; then
        echo -e "******** Performing chimera screen..."
        vsearch --uchime_denovo all_otus_consensus.fasta2 \
        --chimeras all_otus_consensus.chimeras.fasta \
        --nonchimeras all_otus_consensus.fasta \
        --fasta_width 0 \
        --abskew $LR_abskew \
        --mindiv $LR_mindiv
    else
        mv all_otus_consensus.fasta2 all_otus_consensus.fasta
    fi

    #### Remove from master table the chimeric OTUs
    if [ -s all_otus_consensus.chimeras.fasta ]; then
        echo -e "******** Removing chimeric OTUs from master table..."
        grep '^>' all_otus_consensus.chimeras.fasta | sed 's/^>//; s/;.*//' | sort -u > chimera_otu_ids.txt
        awk -F'\t' '
            NR==FNR { chim[$0]=1; next }
            FNR==1  { print; next }
            { key=$6; sub(/;.*/, "", key); if (!(key in chim)) print }
        ' chimera_otu_ids.txt all_otus_run_wide_OTU_info.txt > all_otus_run_wide_OTU_info.nochim.txt
        mv all_otus_run_wide_OTU_info.nochim.txt all_otus_run_wide_OTU_info.txt
        echo "   Removed $(wc -l < chimera_otu_ids.txt) chimeric run-wide OTUs from master table"
        rm -f chimera_otu_ids.txt
    else
        echo "******** No chimeras to remove from master table."
    fi

#######################################################################
################### STEP 8: Sequence correction #######################
#######################################################################

    threads=$cores
    export THREADS=$threads
    export ref_seq_corr
    vsearch_exe=$(which vsearch)
    export vsearch_exe

    if [[ "$marker" == "COI-5P" && "$ampsize" -ge 500 ]]; then
        echo "******** Sequence-correcting OTU consensus sequences..."
        python "$scripts_dir/3-MAP_sequence_correction.py" "all_otus" "$PWD"
    else
        mv "all_otus_consensus.fasta" "all_otus_corrected_otus.fas"
    fi 
#######################################################################
################## STEP 9: Identify run-wide OTUs #####################
#######################################################################

    db_path=$(ls -d $reference_lib_dir/"$reflib"*.fasta)

    #### Identify OTUs using sintax
    vsearch --sintax all_otus_corrected_otus.fas \
        -db $db_path \
        -tabbedout temp_sintax_output.txt \
        -strand plus \
        -sintax_cutoff $sintax_cutoff \
        -threads $cores

#######################################################################
######################## STEP 10: Remove contaminants #########################
#######################################################################

    #### De-noise data
    echo -e "******** Removing contaminants ..."
    Rscript "$scripts_dir/4-RemoveContaminants.R" "$PWD" "$runid" "$marker" $ampsize $min_rep_prop $high_dens_prop $alpha_default $alpha_quantile $k_multiplier "../mapping_${runid}.txt"
       
    #### Create a temporary reduced mapping file for this marker/ampsize
    tmp_mapping="temp_mapping_${runid}_${marker}_${ampsize}.txt"
    awk -F'\t' -v m="$marker" -v a="$ampsize" '$8 == m && $11 == a' "../mapping_${runid}.txt" > "$tmp_mapping"
    echo -e $marker marker $runid runid $ampsize ampsize $tmp_mapping tmp_mapping $minreads minreads $numreps numreps

    echo -e "******** Finalizing data..."
    Rscript "$scripts_dir/5-MAP_final_results.R" "$runid" "$marker" $minreads $numreps $ampsize "$PWD" "$reflib" $componentreads $params_file
    rm all_otus_raw_consensus.fasta all_otus1_corrected_otus.fasta2 all_otus_consensus.fasta2 all_otus_consensus.chimeras.fasta
    rm all_otus_consensus.fasta all_otus_run_wide_OTU_info.txt temp.fasta temp_sintax_output.txt df.txt

    if [ "$componentreads" -eq 1 ]; then
        #### Compress OTU component read files
        zip -r -0 "$runid"_"$marker"_"$ampsize"_OTU_Component_Reads.zip "$runid"_"$marker"_"$ampsize"_OTU_Component_Reads && rm -rf "$runid"_"$marker"_"$ampsize"_OTU_Component_Reads
    fi

    #######################################################################
    ############ STEP 11: Perform BIN analysis (if applicable) ############
    #######################################################################

    if [[ "$marker" == "COI-5P" && "$ampsize" -ge 300 ]]; then
        
        vsearch_db=$(ls -d $reference_lib_dir/"$reflib"*vsearch)
        
        if [ ! -f "$vsearch_db" ]; then
          #### make vsearch reference library from sintax
          echo -e "******** Creating vsearch reference library..."
          #### Name the vsearch DB after the actual dated reference fasta
          #### (e.g. BOLDistilled_COI_Apr2026.vsearch), not the bare $reflib prefix.
          vsearch_name=$(basename "$db_path")
          vsearch_name=${vsearch_name%.fasta}
          vsearch_name=${vsearch_name%.fa}
          vsearch_name=${vsearch_name%_SEQUENCES_sintax}
          cp $db_path $reference_lib_dir/temp_sintax_no_spaces.fasta
          sed '/^>/ s/ /_/g' "$reference_lib_dir/temp_sintax_no_spaces.fasta" > "$reference_lib_dir/temp_sintax_no_spaces.fasta.tmp" \
                && mv "$reference_lib_dir/temp_sintax_no_spaces.fasta.tmp" "$reference_lib_dir/temp_sintax_no_spaces.fasta" || rm -f "$reference_lib_dir/temp_sintax_no_spaces.fasta.tmp"
          vsearch --makeudb_usearch $reference_lib_dir/temp_sintax_no_spaces.fasta --output $reference_lib_dir/"$vsearch_name".vsearch
          vsearch_db=$(ls -d $reference_lib_dir/"$reflib"*vsearch)
          rm $reference_lib_dir/temp_sintax_no_spaces.fasta
        fi

        echo -e "******** Performing BIN match..."
        #### Convert FASTA to single-line
        awk '{if(NR==1) {print $0} else {if($0 ~ /^>/) {print "\n"$0} else {printf $0}}}' temp_toBIN.fasta > single_line.fasta
        rm temp_toBIN.fasta
        mv single_line.fasta temp_toBIN.fasta

        #### Identify sequences using VSEARCH
        vsearch --usearch_global temp_toBIN.fasta \
            --db  $vsearch_db \
            --blast6out temp_vsearch_output.txt \
            --id $BIN_percent_ID \
            --maxhits $BIN_maxhits \
            --maxaccepts $BIN_maxaccepts \
            --threads $cores

        #### Remove low-level hits from VSEARCH results
        min_overlap=$(printf "%.0f" "$(echo "$ampsize * 0.75" | bc -l)")
        awk -v min_overlap="$min_overlap" '
            BEGIN { OFS="\t" }
            $4 >= min_overlap { print }
            ' temp_vsearch_output.txt | sort -k1,1 -k3nr -k4nr | awk '!seen[$1]++' > filtered_vsearch_output.txt
        
        #### Parse hit column into ProcessID and BIN;tax=
        awk -F'\t' 'BEGIN { OFS="\t" } $2 ~ /\|/ { split($2, hit_parts, "|"); print $1, hit_parts[1], hit_parts[2], $3, $4 }' filtered_vsearch_output.txt > parsed_vsearch_output.txt
       
        #### Parse BIN;tax= column into BIN and tax=...
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

        #### Add BIN_MATCH column
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

        #### Add BIN match results to "by sample" results tab in Excel output
        Rscript "$scripts_dir/6-MAP_binmatch.R" "$runid" "$PWD" "$marker" $ampsize
        
#######################################################################
############# STEP 12: Generate output and html report ################
#######################################################################
        
        #### Generate interactive html file using Quarto
       cat > params.yaml <<EOF
runid: "$runid"
working_dir: "$PWD"
marker: "$marker"
ampsize: $ampsize
numreps: $numreps
param_file: "$params_file"
EOF

        quarto render "$scripts_dir/7b-MAP_reporting_with_bins.qmd" --execute-params params.yaml
        mv "$scripts_dir/7b-MAP_reporting_with_bins.html" "$PWD/MAP Report - "${runid}"_"${marker}"_${ampsize}.html"

        #### Tidy up directory
        rm params.yaml temp_vsearch_output.txt filtered_vsearch_output.txt parsed_vsearch_output.txt final_output.txt final_output_with_bin_match.txt temp_toBIN.fasta temp_mapping*.txt   
    else
        #### Generate interactive html file using Quarto
       cat > params.yaml <<EOF
runid: "$runid"
working_dir: "$PWD"
marker: "$marker"
ampsize: $ampsize
numreps: $numreps
param_file: "$params_file"
EOF

        quarto render "$scripts_dir/7a-MAP_reporting.qmd" --execute-params params.yaml
        mv "$scripts_dir/7a-MAP_reporting.html" "$PWD/MAP Report - "${runid}"_"${marker}"_${ampsize}.html"

        #### Tidy up directory
        rm params.yaml temp_mapping*.txt
    fi

    #### Tidy up directory and move back to main working directory
    mkdir -m 777 "1-Results and Report"
    mkdir -m 777 "2-TSV Versions of Results"
    mkdir -m 777 "3-Negative Control OTUs"
    mv *.xlsx *.html "1-Results and Report"
    mv Metabarcoding_Results*.tsv "2-TSV Versions of Results"
    mv *NegativeControlOTUs.tsv "3-Negative Control OTUs"
    cd $working_dir
    mv $marker_dir ./output/
done

#~#~#~#~#~#~#~#~#~#~#~#~#~#~##~#~#~#
#~#~#~#~ END MARKER FOR LOOP #~#~#~#
#~#~#~#~#~~#~#~#~#~#~#~#~#~#~#~#~#~#

##################### EXECUTE PIPELINE: FINISH ##########################
#########################################################################

END_TIME=$(date +%s)
end_human=$(date "+%Y-%m-%d %H:%M:%S")
ELAPSED=$(( END_TIME - START_TIME ))
# Format elapsed time as HH:MM:SS
HOURS=$(( ELAPSED / 3600 ))
MINUTES=$(( (ELAPSED % 3600) / 60 ))
SECONDS=$(( ELAPSED % 60 ))

{
echo "==============================="
echo "Pipeline Complete"
echo "==============================="
echo "Machine   : $(hostname)"
echo "CPU Cores : $cores"
echo "Run name:   ${runid}"
echo "Started:    ${start_human}"
echo "Finished:   ${end_human}"
echo "Total Time: $(printf '%02d:%02d:%02d' $HOURS $MINUTES $SECONDS)"
echo "==============================="
} >> ./output/pipeline.log

mv Demultiplexing_Results_* Individual_Raw_Fastq_Files* *_readcounts.txt metadata* mapping* runinfo.txt name_map.tsv ./output/




