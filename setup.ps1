#!/usr/bin/env powershell
# NGS Pipeline Setup Script for Windows
# Equivalent to 'make setup' for Windows PowerShell

param(
    [string]$Target = "help"
)

function Show-Help {
    Write-Host "NGS Pipeline Setup" -ForegroundColor Green
    Write-Host "==============================" -ForegroundColor Green
    Write-Host ""
    Write-Host "Available commands:" -ForegroundColor Yellow
    Write-Host "  .\setup.ps1 setup        - Install dependencies and setup environment"
    Write-Host "  .\setup.ps1 install-deps - Install conda/mamba dependencies"
    Write-Host "  .\setup.ps1 build-refs   - Download and build reference files"
    Write-Host "  .\setup.ps1 test         - Run quick test with toy data"
    Write-Host "  .\setup.ps1 demo         - Run full demo pipeline"
    Write-Host "  .\setup.ps1 check-env    - Check environment and dependencies"
    Write-Host "  .\setup.ps1 clean        - Clean up temporary files"
    Write-Host ""
    Write-Host "Example usage:" -ForegroundColor Cyan
    Write-Host "  .\setup.ps1 setup"
    Write-Host ""
}

function Install-Dependencies {
    Write-Host "Installing conda environment for Windows..." -ForegroundColor Yellow
    Write-Host "Note: Using Windows-compatible environment (bioinformatics tools via Docker/WSL)" -ForegroundColor Cyan
    
    # Check if mamba or conda is available
    $mambaAvailable = Get-Command mamba -ErrorAction SilentlyContinue
    $condaAvailable = Get-Command conda -ErrorAction SilentlyContinue
    
    # Use Windows-specific environment file
    $envFile = "env/environment_minimal.yml"
    
    if ($mambaAvailable) {
        Write-Host "Using mamba for faster installation..." -ForegroundColor Green
        try {
            & mamba env create -f $envFile -n ngs_pipeline
        } catch {
            Write-Host "Environment already exists, updating..." -ForegroundColor Yellow
            & mamba env update -f $envFile -n ngs_pipeline
        }
    } elseif ($condaAvailable) {
        Write-Host "Using conda for installation..." -ForegroundColor Green
        try {
            & conda env create -f $envFile -n ngs_pipeline
        } catch {
            Write-Host "Environment already exists, updating..." -ForegroundColor Yellow
            & conda env update -f $envFile -n ngs_pipeline
        }
    } else {
        Write-Error "Neither mamba nor conda found. Please install Anaconda/Miniconda first."
        Write-Host "Download from: https://docs.conda.io/en/latest/miniconda.html" -ForegroundColor Cyan
        exit 1
    }
    
    if ($LASTEXITCODE -ne 0) {
        Write-Error "Failed to create conda environment"
        exit 1
    }
    
    Write-Host "Installing Python package in development mode..." -ForegroundColor Yellow
    
    # Activate environment and install package
    $condaBase = & conda info --base
    & cmd /c "call `"$condaBase\Scripts\activate.bat`" ngs_pipeline && pip install -e ."
    
    if ($LASTEXITCODE -eq 0) {
        Write-Host "Dependencies installed successfully!" -ForegroundColor Green
        Write-Host ""
        Write-Host "IMPORTANT: For full bioinformatics pipeline functionality on Windows:" -ForegroundColor Yellow
        Write-Host "1. Install Docker Desktop from: https://www.docker.com/products/docker-desktop" -ForegroundColor Cyan
        Write-Host "2. Or install WSL2 with Ubuntu from Microsoft Store" -ForegroundColor Cyan
        Write-Host "3. Install Nextflow: https://www.nextflow.io/docs/latest/getstarted.html#installation" -ForegroundColor Cyan
    } else {
        Write-Error "Failed to install dependencies"
        exit 1
    }
}

function Build-References {
    Write-Host "Setting up reference files..." -ForegroundColor Yellow
    
    # Create references directory if it doesn't exist
    if (-not (Test-Path "references")) {
        New-Item -ItemType Directory -Path "references" -Force
        Write-Host "Created references directory"
    }
    
    # Check if we're in WSL/Cygwin or native Windows
    if (Test-Path "references/get_references.sh") {
        Write-Host "Attempting to run reference download script..." -ForegroundColor Yellow
        
        # Try different approaches to run the bash script
        $bashAvailable = Get-Command bash -ErrorAction SilentlyContinue
        $wslAvailable = Get-Command wsl -ErrorAction SilentlyContinue
        
        if ($bashAvailable) {
            try {
                Write-Host "Using bash to download references..." -ForegroundColor Green
                & bash references/get_references.sh --hs1 --build-index star bowtie2 bismark --dbsnp
                if ($LASTEXITCODE -eq 0) {
                    Write-Host "Reference files downloaded successfully." -ForegroundColor Green
                    return
                }
            } catch {
                Write-Warning "Bash script execution failed."
            }
        }
        
        if ($wslAvailable) {
            try {
                Write-Host "Using WSL to download references..." -ForegroundColor Green
                & wsl bash references/get_references.sh --hs1 --build-index star bowtie2 bismark --dbsnp
                if ($LASTEXITCODE -eq 0) {
                    Write-Host "Reference files downloaded successfully." -ForegroundColor Green
                    return
                }
            } catch {
                Write-Warning "WSL script execution failed."
            }
        }
    }
    
    # If bash script fails, provide manual instructions
    Write-Host "Automatic reference download not available on Windows." -ForegroundColor Yellow
    Write-Host "Please download reference files manually or use Docker/WSL:" -ForegroundColor Cyan
    Write-Host ""
    Write-Host "Required reference files:" -ForegroundColor White
    Write-Host "1. Human T2T CHM13v2.0 genome:" -ForegroundColor Cyan
    Write-Host "   https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/009/914/755/GCF_009914755.1_T2T-CHM13v2.0/" -ForegroundColor Blue
    Write-Host ""
    Write-Host "2. Gene annotation (GTF):" -ForegroundColor Cyan
    Write-Host "   https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/" -ForegroundColor Blue
    Write-Host ""
    Write-Host "3. dbSNP variants:" -ForegroundColor Cyan
    Write-Host "   https://ftp.ncbi.nlm.nih.gov/snp/organisms/human_9606/VCF/" -ForegroundColor Blue
    Write-Host ""
    Write-Host "Alternative: Use Docker for full bioinformatics pipeline:" -ForegroundColor Yellow
    Write-Host "  docker run -it --rm -v `${PWD}:/workspace ngs_pipeline:latest" -ForegroundColor Green
}

function Run-Test {
    Write-Host "Running quick test with toy data..." -ForegroundColor Yellow
    
    # Check environment first
    Check-Environment
    
    $condaBase = & conda info --base
    & cmd /c "call `"$condaBase\Scripts\activate.bat`" ngs_pipeline && nextflow run main.nf -profile test -with-report test_report.html -with-timeline test_timeline.html"
    
    if ($LASTEXITCODE -eq 0) {
        Write-Host "Test completed successfully!" -ForegroundColor Green
        Write-Host "Check test_report.html for detailed results." -ForegroundColor Cyan
    } else {
        Write-Error "Test failed. Check the output above for errors."
    }
}

function Run-Demo {
    Write-Host "Running full demo pipeline..." -ForegroundColor Yellow
    
    Check-Environment
    
    $condaBase = & conda info --base
    & cmd /c "call `"$condaBase\Scripts\activate.bat`" ngs_pipeline && nextflow run main.nf -profile local --samples config/samplesheet.tsv --enable_tcga false --enable_string false -with-report demo_report.html -with-timeline demo_timeline.html"
    
    if ($LASTEXITCODE -eq 0) {
        Write-Host "Demo completed successfully!" -ForegroundColor Green
        Write-Host "Check demo_report.html for detailed results." -ForegroundColor Cyan
    } else {
        Write-Error "Demo failed. Check the output above for errors."
    }
}

function Check-Environment {
    Write-Host "Checking environment and dependencies..." -ForegroundColor Yellow
    
    # Check if conda environment exists
    $envExists = & conda env list | Select-String "ngs_pipeline"
    if (-not $envExists) {
        Write-Error "ngs_pipeline conda environment not found. Run 'install-deps' first."
        return
    }
    
    # Check Python package
    $condaBase = & conda info --base
    & cmd /c "call `"$condaBase\Scripts\activate.bat`" ngs_pipeline && python -c `"import ngs_pipeline; print('Python package: OK')`""
    
    # Check Nextflow
    $nextflowAvailable = Get-Command nextflow -ErrorAction SilentlyContinue
    if ($nextflowAvailable) {
        & nextflow -version
        Write-Host "Environment check: PASSED" -ForegroundColor Green
    } else {
        Write-Warning "Nextflow not found in PATH. Make sure it's installed and accessible."
        Write-Host "Install from: https://www.nextflow.io/docs/latest/getstarted.html#installation" -ForegroundColor Cyan
    }
}

function Clean-Files {
    Write-Host "Cleaning temporary files..." -ForegroundColor Yellow
    
    # Remove Nextflow work directory
    if (Test-Path "work") {
        Remove-Item -Path "work" -Recurse -Force
        Write-Host "Removed work directory"
    }
    
    # Remove Nextflow hidden files
    Get-ChildItem -Path "." -Name ".nextflow*" | ForEach-Object {
        Remove-Item -Path $_ -Force
        Write-Host "Removed $_"
    }
    
    # Remove HTML reports
    Get-ChildItem -Path "." -Name "*.html" | ForEach-Object {
        Remove-Item -Path $_ -Force
        Write-Host "Removed $_"
    }
    
    # Remove other temporary files
    $tempFiles = @("*.svg", "trace.txt", "timeline.html", "report.html")
    foreach ($pattern in $tempFiles) {
        Get-ChildItem -Path "." -Name $pattern | ForEach-Object {
            Remove-Item -Path $_ -Force
            Write-Host "Removed $_"
        }
    }
    
    # Remove Python cache
    Get-ChildItem -Path "." -Recurse -Name "*.pyc" | Remove-Item -Force
    Get-ChildItem -Path "." -Recurse -Directory -Name "__pycache__" | Remove-Item -Recurse -Force
    
    Write-Host "Cleanup complete." -ForegroundColor Green
}

function Setup-Complete {
    Write-Host "Setting up NGS Pipeline..." -ForegroundColor Green
    Install-Dependencies
    Build-References
    Write-Host ""
    Write-Host "Setup complete! Ready to run pipeline." -ForegroundColor Green
    Write-Host ""
    Write-Host "Next steps:" -ForegroundColor Yellow
    Write-Host "1. Run a quick test: .\setup.ps1 test"
    Write-Host "2. Or run the full demo: .\setup.ps1 demo"
    Write-Host "3. Check the Jupyter notebook: notebooks/ngs_pipeline.ipynb"
}

# Main execution logic
switch ($Target.ToLower()) {
    "help" { Show-Help }
    "setup" { Setup-Complete }
    "install-deps" { Install-Dependencies }
    "build-refs" { Build-References }
    "test" { Run-Test }
    "demo" { Run-Demo }
    "check-env" { Check-Environment }
    "clean" { Clean-Files }
    default {
        Write-Error "Unknown target: $Target"
        Show-Help
    }
}