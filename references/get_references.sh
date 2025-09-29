#!/bin/bash

# NGS Pipeline - Reference Genome Setup Script
# Downloads and builds indices for hs1 (T2T/CHM13) reference genome

set -euo pipefail

# Default parameters
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REF_DIR="${SCRIPT_DIR}/hs1"
BUILD_INDICES=""
DOWNLOAD_DBSNP=false
THREADS=8
MEMORY=50  # GB for STAR indexing

# NCBI/Ensembl URLs for hs1 (T2T/CHM13)
HS1_FASTA_URL="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/009/914/755/GCF_009914755.1_T2T-CHM13v2.0/GCF_009914755.1_T2T-CHM13v2.0_genomic.fna.gz"
HS1_GTF_URL="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/009/914/755/GCF_009914755.1_T2T-CHM13v2.0/GCF_009914755.1_T2T-CHM13v2.0_genomic.gtf.gz"
DBSNP_URL="https://ftp.ncbi.nih.gov/snp/organisms/human_9606_b151_GRCh38p7/VCF/All_20180423.vcf.gz"

# Color codes for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# Logging functions
log_info() {
    echo -e "${BLUE}[INFO]${NC} $1"
}

log_warn() {
    echo -e "${YELLOW}[WARN]${NC} $1"
}

log_error() {
    echo -e "${RED}[ERROR]${NC} $1"
}

log_success() {
    echo -e "${GREEN}[SUCCESS]${NC} $1"
}

# Progress indicator
show_progress() {
    local pid=$1
    local message=$2
    local spin='|/-\'
    local i=0
    while kill -0 $pid 2>/dev/null; do
        local char=${spin:$((i%4)):1}
        printf "\r${BLUE}[INFO]${NC} $message $char"
        sleep 0.1
        ((i++))
    done
    printf "\r${GREEN}[SUCCESS]${NC} $message ✓\n"
}

# Function to check if command exists
command_exists() {
    command -v "$1" >/dev/null 2>&1
}

# Function to check available memory
check_memory() {
    local required_gb=$1
    local available_kb
    
    if [[ "$OSTYPE" == "linux-gnu"* ]]; then
        available_kb=$(grep MemAvailable /proc/meminfo | awk '{print $2}')
        available_gb=$((available_kb / 1024 / 1024))
    elif [[ "$OSTYPE" == "darwin"* ]]; then
        available_bytes=$(vm_stat | perl -ne '/free.*?(\d+)/ and $free+=$1; /speculative.*?(\d+)/ and $spec+=$1; END { print +($free+$spec)*4096 }')
        available_gb=$((available_bytes / 1024 / 1024 / 1024))
    else
        log_warn "Cannot detect available memory on this system"
        return 0
    fi
    
    if [[ $available_gb -lt $required_gb ]]; then
        log_error "Insufficient memory: ${available_gb}GB available, ${required_gb}GB required"
        return 1
    fi
    
    log_info "Memory check passed: ${available_gb}GB available, ${required_gb}GB required"
    return 0
}

# Function to download file with progress
download_file() {
    local url=$1
    local output=$2
    local description=$3
    
    log_info "Downloading $description..."
    
    if command_exists wget; then
        wget --progress=bar:force -O "$output" "$url" 2>&1 | \
        grep -o '[0-9]*%' | \
        while read percent; do
            printf "\r${BLUE}[INFO]${NC} Downloading $description: $percent"
        done
        printf "\n"
    elif command_exists curl; then
        curl -# -L -o "$output" "$url"
    else
        log_error "Neither wget nor curl found. Please install one of them."
        exit 1
    fi
    
    if [[ ! -f "$output" ]]; then
        log_error "Download failed: $output"
        exit 1
    fi
    
    log_success "Downloaded $description"
}

# Function to check file integrity
check_integrity() {
    local file=$1
    local description=$2
    
    log_info "Checking integrity of $description..."
    
    # Basic checks
    if [[ ! -f "$file" ]]; then
        log_error "$description not found: $file"
        return 1
    fi
    
    if [[ ! -s "$file" ]]; then
        log_error "$description is empty: $file"
        return 1
    fi
    
    # Check if it's a valid gzipped file
    if [[ "$file" == *.gz ]]; then
        if ! gzip -t "$file" 2>/dev/null; then
            log_error "$description is corrupted (invalid gzip): $file"
            return 1
        fi
    fi
    
    log_success "$description integrity check passed"
    return 0
}

# Function to build STAR index
build_star_index() {
    local fasta_file=$1
    local gtf_file=$2
    local index_dir=$3
    
    log_info "Building STAR index (this may take 30-60 minutes)..."
    
    if ! check_memory 30; then
        log_error "STAR indexing requires at least 30GB RAM for hs1 genome"
        return 1
    fi
    
    mkdir -p "$index_dir"
    
    # Check if index already exists
    if [[ -f "${index_dir}/SA" ]]; then
        log_warn "STAR index already exists in $index_dir"
        return 0
    fi
    
    STAR \
        --runMode genomeGenerate \
        --genomeDir "$index_dir" \
        --genomeFastaFiles "$fasta_file" \
        --sjdbGTFfile "$gtf_file" \
        --sjdbOverhang 100 \
        --runThreadN "$THREADS" \
        --limitGenomeGenerateRAM $((MEMORY * 1024 * 1024 * 1024)) \
        --genomeSAindexNbases 14 > "${index_dir}/star_build.log" 2>&1 &
    
    show_progress $! "Building STAR index"
    wait
    
    if [[ ! -f "${index_dir}/SA" ]]; then
        log_error "STAR index build failed. Check ${index_dir}/star_build.log"
        return 1
    fi
    
    log_success "STAR index built successfully"
}

# Function to build Bowtie2 index
build_bowtie2_index() {
    local fasta_file=$1
    local index_prefix=$2
    
    log_info "Building Bowtie2 index..."
    
    # Check if index already exists
    if [[ -f "${index_prefix}.1.bt2" ]]; then
        log_warn "Bowtie2 index already exists: $index_prefix"
        return 0
    fi
    
    bowtie2-build \
        --threads "$THREADS" \
        "$fasta_file" \
        "$index_prefix" > "${index_prefix}_build.log" 2>&1 &
    
    show_progress $! "Building Bowtie2 index"
    wait
    
    if [[ ! -f "${index_prefix}.1.bt2" ]]; then
        log_error "Bowtie2 index build failed. Check ${index_prefix}_build.log"
        return 1
    fi
    
    log_success "Bowtie2 index built successfully"
}

# Function to build Bismark index
build_bismark_index() {
    local ref_dir=$1
    
    log_info "Building Bismark index..."
    
    # Check if index already exists
    if [[ -d "${ref_dir}/Bisulfite_Genome" ]]; then
        log_warn "Bismark index already exists in $ref_dir"
        return 0
    fi
    
    bismark_genome_preparation \
        --parallel "$((THREADS / 4))" \
        --verbose \
        "$ref_dir" > "${ref_dir}/bismark_build.log" 2>&1 &
    
    show_progress $! "Building Bismark index"
    wait
    
    if [[ ! -d "${ref_dir}/Bisulfite_Genome" ]]; then
        log_error "Bismark index build failed. Check ${ref_dir}/bismark_build.log"
        return 1
    fi
    
    log_success "Bismark index built successfully"
}

# Function to process FASTA file
process_fasta() {
    local fasta_gz=$1
    local fasta_file=$2
    
    log_info "Processing FASTA file..."
    
    if [[ ! -f "$fasta_file" ]]; then
        log_info "Decompressing FASTA..."
        gunzip -c "$fasta_gz" > "$fasta_file"
    fi
    
    # Create FASTA index
    if [[ ! -f "${fasta_file}.fai" ]]; then
        log_info "Creating FASTA index..."
        samtools faidx "$fasta_file"
    fi
    
    # Create sequence dictionary for GATK
    if [[ ! -f "${fasta_file%.*}.dict" ]]; then
        log_info "Creating sequence dictionary..."
        if command_exists picard; then
            picard CreateSequenceDictionary \
                R="$fasta_file" \
                O="${fasta_file%.*}.dict"
        elif command_exists gatk; then
            gatk CreateSequenceDictionary \
                -R "$fasta_file" \
                -O "${fasta_file%.*}.dict"
        else
            log_warn "Neither picard nor gatk found. Skipping dictionary creation."
        fi
    fi
    
    log_success "FASTA processing complete"
}

# Function to process GTF file
process_gtf() {
    local gtf_gz=$1
    local gtf_file=$2
    
    log_info "Processing GTF file..."
    
    if [[ ! -f "$gtf_file" ]]; then
        log_info "Decompressing GTF..."
        gunzip -c "$gtf_gz" > "$gtf_file"
    fi
    
    # Create BED file for RSeQC
    local bed_file="${gtf_file%.*}.bed"
    if [[ ! -f "$bed_file" ]]; then
        log_info "Converting GTF to BED format..."
        if command_exists gtfToGenePred && command_exists genePredToBed; then
            gtfToGenePred "$gtf_file" "${gtf_file%.*}.genePred"
            genePredToBed "${gtf_file%.*}.genePred" "$bed_file"
            rm -f "${gtf_file%.*}.genePred"
        else
            log_warn "UCSC tools not found. Skipping GTF to BED conversion."
        fi
    fi
    
    log_success "GTF processing complete"
}

# Function to show usage
show_usage() {
    cat << EOF
NGS Pipeline - Reference Setup

USAGE:
    $0 [OPTIONS]

OPTIONS:
    --hs1                   Download hs1 (T2T/CHM13) reference genome
    --build-index TOOLS     Build indices for specified tools (comma-separated)
                           Available: star, bowtie2, bismark, all
    --dbsnp                 Download dbSNP VCF file
    --ref-dir DIR          Reference directory (default: ./hs1)
    --threads N            Number of threads for indexing (default: 8)
    --memory N             Memory limit in GB for STAR indexing (default: 50)
    --help                 Show this help message

EXAMPLES:
    # Download hs1 genome only
    $0 --hs1

    # Download and build all indices
    $0 --hs1 --build-index all

    # Download with specific indices
    $0 --hs1 --build-index star,bowtie2 --dbsnp

    # Custom reference directory
    $0 --hs1 --build-index all --ref-dir /opt/references/hs1

REQUIREMENTS:
    - wget or curl (for downloading)
    - STAR (for STAR index)
    - bowtie2 (for Bowtie2 index)
    - bismark (for Bismark index)
    - samtools (for FASTA indexing)
    - At least 30GB RAM for STAR indexing
    - At least 100GB free disk space

EOF
}

# Parse command line arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        --hs1)
            DOWNLOAD_HS1=true
            shift
            ;;
        --build-index)
            BUILD_INDICES="$2"
            shift 2
            ;;
        --dbsnp)
            DOWNLOAD_DBSNP=true
            shift
            ;;
        --ref-dir)
            REF_DIR="$2"
            shift 2
            ;;
        --threads)
            THREADS="$2"
            shift 2
            ;;
        --memory)
            MEMORY="$2"
            shift 2
            ;;
        --help|-h)
            show_usage
            exit 0
            ;;
        *)
            log_error "Unknown option: $1"
            show_usage
            exit 1
            ;;
    esac
done

# Check if no action specified
if [[ -z "${DOWNLOAD_HS1:-}" ]]; then
    log_error "No action specified. Use --hs1 to download reference genome."
    show_usage
    exit 1
fi

# Main execution
main() {
    log_info "NGS Pipeline - Reference Setup"
    log_info "=========================================="
    log_info "Reference directory: $REF_DIR"
    log_info "Threads: $THREADS"
    log_info "Memory limit: ${MEMORY}GB"
    echo
    
    # Create reference directory
    mkdir -p "$REF_DIR"
    cd "$REF_DIR"
    
    # Download hs1 reference
    if [[ "${DOWNLOAD_HS1:-}" == "true" ]]; then
        log_info "Setting up hs1 (T2T/CHM13) reference genome..."
        
        # Download FASTA
        if [[ ! -f "hs1.fasta.gz" ]]; then
            download_file "$HS1_FASTA_URL" "hs1.fasta.gz" "hs1 genome FASTA"
        fi
        check_integrity "hs1.fasta.gz" "hs1 genome FASTA"
        
        # Download GTF
        if [[ ! -f "hs1.gtf.gz" ]]; then
            download_file "$HS1_GTF_URL" "hs1.gtf.gz" "hs1 annotation GTF"
        fi
        check_integrity "hs1.gtf.gz" "hs1 annotation GTF"
        
        # Process files
        process_fasta "hs1.fasta.gz" "hs1.fasta"
        process_gtf "hs1.gtf.gz" "hs1.gtf"
    fi
    
    # Download dbSNP
    if [[ "$DOWNLOAD_DBSNP" == "true" ]]; then
        if [[ ! -f "dbsnp.vcf.gz" ]]; then
            download_file "$DBSNP_URL" "dbsnp.vcf.gz" "dbSNP VCF"
        fi
        check_integrity "dbsnp.vcf.gz" "dbSNP VCF"
        
        # Index VCF
        if [[ ! -f "dbsnp.vcf.gz.tbi" ]]; then
            log_info "Indexing dbSNP VCF..."
            tabix -p vcf dbsnp.vcf.gz
        fi
    fi
    
    # Build indices
    if [[ -n "$BUILD_INDICES" ]]; then
        log_info "Building indices: $BUILD_INDICES"
        
        # Create indices directory
        mkdir -p indices
        
        # Parse build indices
        IFS=',' read -ra INDICES <<< "$BUILD_INDICES"
        
        for index in "${INDICES[@]}"; do
            case $index in
                star|all)
                    if command_exists STAR; then
                        build_star_index "hs1.fasta" "hs1.gtf" "indices/star"
                    else
                        log_error "STAR not found. Please install STAR."
                    fi
                    ;;
                bowtie2|all)
                    if command_exists bowtie2-build; then
                        build_bowtie2_index "hs1.fasta" "indices/bowtie2/hs1"
                    else
                        log_error "Bowtie2 not found. Please install Bowtie2."
                    fi
                    ;;
                bismark|all)
                    if command_exists bismark_genome_preparation; then
                        # Copy FASTA to bismark directory for indexing
                        mkdir -p indices/bismark
                        cp hs1.fasta indices/bismark/
                        build_bismark_index "indices/bismark"
                    else
                        log_error "Bismark not found. Please install Bismark."
                    fi
                    ;;
                *)
                    log_warn "Unknown index type: $index"
                    ;;
            esac
        done
    fi
    
    # Generate summary
    log_info "Generating reference summary..."
    cat > reference_info.txt << EOF
NGS Pipeline - Reference Information
==============================================

Reference Genome: hs1 (T2T/CHM13v2.0)
Date Generated: $(date)
Directory: $(pwd)

Files Generated:
EOF
    
    for file in hs1.fasta hs1.gtf hs1.fasta.fai hs1.dict dbsnp.vcf.gz; do
        if [[ -f "$file" ]]; then
            size=$(du -h "$file" | cut -f1)
            echo "  ✓ $file ($size)" >> reference_info.txt
        fi
    done
    
    echo "" >> reference_info.txt
    echo "Indices Built:" >> reference_info.txt
    
    for index_dir in indices/*/; do
        if [[ -d "$index_dir" ]]; then
            index_name=$(basename "$index_dir")
            size=$(du -sh "$index_dir" | cut -f1)
            echo "  ✓ $index_name ($size)" >> reference_info.txt
        fi
    done
    
    echo
    log_success "Reference setup completed successfully!"
    log_info "Summary saved to: ${REF_DIR}/reference_info.txt"
    
    if [[ -f "reference_info.txt" ]]; then
        echo
        cat reference_info.txt
    fi
}

# Check prerequisites
check_prerequisites() {
    log_info "Checking prerequisites..."
    
    local missing_tools=()
    
    # Check for required tools
    if ! command_exists wget && ! command_exists curl; then
        missing_tools+=("wget or curl")
    fi
    
    if ! command_exists samtools; then
        missing_tools+=("samtools")
    fi
    
    if [[ "$BUILD_INDICES" == *"star"* ]] || [[ "$BUILD_INDICES" == *"all"* ]]; then
        if ! command_exists STAR; then
            missing_tools+=("STAR")
        fi
    fi
    
    if [[ "$BUILD_INDICES" == *"bowtie2"* ]] || [[ "$BUILD_INDICES" == *"all"* ]]; then
        if ! command_exists bowtie2-build; then
            missing_tools+=("bowtie2")
        fi
    fi
    
    if [[ "$BUILD_INDICES" == *"bismark"* ]] || [[ "$BUILD_INDICES" == *"all"* ]]; then
        if ! command_exists bismark_genome_preparation; then
            missing_tools+=("bismark")
        fi
    fi
    
    if [[ ${#missing_tools[@]} -gt 0 ]]; then
        log_error "Missing required tools:"
        for tool in "${missing_tools[@]}"; do
            echo "  - $tool"
        done
        echo
        log_info "Please install missing tools and try again."
        exit 1
    fi
    
    log_success "All prerequisites satisfied"
}

# Run prerequisite check and main function
check_prerequisites
main

log_success "Reference setup script completed!"