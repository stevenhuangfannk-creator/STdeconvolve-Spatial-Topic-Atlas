# Windows 串行运行器：每个模块独立 R 进程，记录运行时间、父进程峰值内存和完整日志。
# 输入：在仓库根目录执行；Rscript 参数选择已经审查的 R，不安装、不升级环境。
param(
    [Parameter(Mandatory=$true)][string]$Rscript,
    [string[]]$Scripts = @('01_Official_Reproduction.R','02_Topic_Number_Selection.R',
        '03_Spatial_Visualization.R','04_Topic_Annotation.R','05_Functional_Enrichment.R',
        '06_Multi_Sample_Analysis.R','07_Spatial_Neighborhood.R')
)
$ErrorActionPreference = 'Stop'
if (!(Test-Path -LiteralPath '.\00_Start.R')) { throw '请从仓库根目录运行。' }
$taskRoot = (Get-Location).Path
$taskR = (Resolve-Path -LiteralPath $Rscript).Path
New-Item -ItemType Directory -Force -Path 'environment/logs' | Out-Null
foreach ($script in $Scripts) {
    if ($script -notmatch '^[0-9A-Za-z_]+\.R$') { throw '脚本必须是项目根目录的 R 文件名。' }
    if (!(Test-Path -LiteralPath $script)) { throw "缺少脚本 $script" }
    $taskStart = Get-Date
    $taskInfo = [System.Diagnostics.ProcessStartInfo]::new()
    $taskInfo.FileName = $taskR
    $taskInfo.Arguments = "--vanilla $script"
    $taskInfo.WorkingDirectory = $taskRoot
    $taskInfo.UseShellExecute = $false
    $taskInfo.CreateNoWindow = $true
    $taskInfo.RedirectStandardOutput = $true
    $taskInfo.RedirectStandardError = $true
    $taskInfo.Environment['OMP_NUM_THREADS'] = '1'
    $taskInfo.Environment['OPENBLAS_NUM_THREADS'] = '1'
    $taskProcess = [System.Diagnostics.Process]::Start($taskInfo)
    $taskOut = $taskProcess.StandardOutput.ReadToEndAsync()
    $taskErr = $taskProcess.StandardError.ReadToEndAsync()
    $taskPeak = 0L
    while (!$taskProcess.HasExited) {
        $taskProcess.Refresh()
        try { $taskPeak = [Math]::Max($taskPeak, $taskProcess.PeakWorkingSet64) } catch { }
        Start-Sleep -Milliseconds 500
    }
    $taskProcess.WaitForExit()
    $taskLog = $taskOut.Result + "`n--- STDERR / warnings ---`n" + $taskErr.Result
    # 发布日志去除机器私人目录，保留错误、参数和科学结果。
    $taskLog = $taskLog.Replace($taskRoot, '<PROJECT>').Replace($env:USERPROFILE, '<USER>')
    $taskLog | Set-Content -LiteralPath "environment/logs/$($script -replace '\.R$','').log" -Encoding utf8
    [pscustomobject]@{
        script = $script
        elapsed_seconds = [Math]::Round(((Get-Date) - $taskStart).TotalSeconds, 3)
        parent_peak_working_set_mb = [Math]::Round($taskPeak / 1MB, 2)
        exit_code = $taskProcess.ExitCode
        memory_scope = 'sampled R parent only; excludes child workers and other applications'
    } | Export-Csv -LiteralPath 'environment/process_metrics.csv' -NoTypeInformation -Append -Encoding utf8
    Write-Output "$script : exit $($taskProcess.ExitCode), peak parent $([Math]::Round($taskPeak / 1MB)) MB"
    if ($taskProcess.ExitCode -ne 0) { throw "模块失败；见 environment/logs/$($script -replace '\.R$','').log" }
}
