# =============================================================================
#  00_arrange_files_by_cancer.ps1
#
#  Organises GDC STAR-Counts TSV files into cancer-type subfolders by
#  querying the GDC API to resolve each file UUID to a TCGA project ID.
#
#  Input  : GDC manifest (.txt) + a flat folder of UUID-named subfolders
#  Output : ByCancer/{CANCER_TYPE}/*.tsv
#
#  Usage  : Edit the three config paths below, then run in PowerShell:
#           .\00_arrange_files_by_cancer.ps1
# =============================================================================

# ── CONFIGURE THESE THREE PATHS ──────────────────────────────────────────────
$manifestPath = "path\to\gdc_manifest.txt"     # <-- edit
$sourceDir    = "path\to\RNA-Seq"              # <-- edit (flat UUID folders live here)
$destParent   = Join-Path $sourceDir "ByCancer" # output: RNA-Seq\ByCancer\

# ── 0. REMOVE INCORRECTLY CREATED FOLDERS (idempotency guard) ────────────────
$wrongFolders = @(
    Join-Path $destParent "RNA-Seq",
    Join-Path $destParent "ByCancer"
)
foreach ($folder in $wrongFolders) {
    if (Test-Path $folder) {
        Write-Host "Removing incorrect folder: $folder"
        Remove-Item $folder -Recurse -Force
    }
}

# ── 1. CREATE ByCancer DIRECTORY ─────────────────────────────────────────────
if (-not (Test-Path $destParent)) {
    New-Item -ItemType Directory -Path $destParent | Out-Null
}

# ── 2. READ MANIFEST — extract file UUIDs (column 0) ─────────────────────────
$manifestLines = Get-Content $manifestPath | Select-Object -Skip 1
$uuids = $manifestLines | ForEach-Object { ($_ -split "\t")[0] }
Write-Host "Total UUIDs in manifest: $($uuids.Count)"

# ── 3. QUERY GDC API IN BATCHES OF 100 ───────────────────────────────────────
$batchSize = 100
$uuidMap   = @{}   # UUID -> cancer type abbreviation (e.g. BRCA)

for ($i = 0; $i -lt $uuids.Count; $i += $batchSize) {
    $batch = $uuids[$i..([Math]::Min($i + $batchSize - 1, $uuids.Count - 1))]

    $body = @{
        filters = @{
            op      = "in"
            content = @{
                field  = "file_id"
                value  = $batch
            }
        }
        fields = "file_id,cases.project.project_id"
        size   = $batchSize
    } | ConvertTo-Json -Depth 5

    try {
        $response = Invoke-RestMethod `
            -Uri         "https://api.gdc.cancer.gov/files" `
            -Method      POST `
            -ContentType "application/json" `
            -Body        $body

        foreach ($hit in $response.data.hits) {
            $projectId  = $hit.cases[0].project.project_id   # e.g. TCGA-BRCA
            $cancerType = $projectId -replace "TCGA-", ""    # e.g. BRCA
            $uuidMap[$hit.file_id] = $cancerType
        }

        $batchNum   = [Math]::Ceiling(($i + 1) / $batchSize)
        $totalBatch = [Math]::Ceiling($uuids.Count / $batchSize)
        Write-Host "  Batch $batchNum / $totalBatch done"
    }
    catch {
        Write-Warning "Batch at index $i failed: $_"
    }
}

Write-Host "Cancer types resolved for $($uuidMap.Count) / $($uuids.Count) files."

# ── 4. MOVE TSV FILES INTO ByCancer/{CANCER_TYPE}/ ───────────────────────────
$moved   = 0
$skipped = 0

Get-ChildItem $sourceDir -Directory | ForEach-Object {
    $uuidFolder = $_.Name
    if ($uuidFolder -eq "ByCancer") { return }

    $cancerType = $uuidMap[$uuidFolder]
    if (-not $cancerType) {
        Write-Warning "No cancer type for UUID: $uuidFolder — skipping"
        $skipped++
        return
    }

    $cancerFolder = Join-Path $destParent $cancerType
    if (-not (Test-Path $cancerFolder)) {
        New-Item -ItemType Directory -Path $cancerFolder | Out-Null
    }

    $tsv = Get-ChildItem $_.FullName -Filter "*.tsv" | Select-Object -First 1
    if ($tsv) {
        Move-Item $tsv.FullName (Join-Path $cancerFolder $tsv.Name) -Force
        $moved++
    }
}

Write-Host ""
Write-Host "Done!"
Write-Host "  Files moved  : $moved"
Write-Host "  Files skipped: $skipped"
Write-Host "  Output folder: $destParent"
